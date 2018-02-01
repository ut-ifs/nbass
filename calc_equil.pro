@lib/efit_psi2b.pro
@lib/read_geqdsk.pro

pro calc_equil,chord,field,plot=plot
    @equil_param.idl

    ;Construct equilibrium and B field
    case Bmodel of
	'millerold': begin
	;Old calculation of B(grid point) based on Miller model [Miller et al, PoP 5 4 (1998)]
	;It's not so easy to invert Miller coordinates, but we need to know the r,theta coordinates that map to chr, chz.
	;So we will evaluate the Miller coordinates on a high resolution grid of of r,theta and pick the values that get numerically closest to chr, chz
	n_theta_p = 1700 ;number of theta points in grid
	n_rho_s = 300 ;number of rho points in grid
	theta_p = maken(0d,2*!dpi,n_theta_p) ;dim [n_theta_p]
	rho_s = maken(0d,aminor,n_rho_s)     ;dim [n_rho_s]
	rovera = rho_s/aminor               ;dim [n_rho_s]
	;assume delta(r) = delta95/0.95^2 * (r/a)^2
	triang = dblarr(n_rho_s,n_theta_p)  ;dim [n_rho_s,n_theta_p]
	triang[*,where(theta_p le !dpi,count)] = cmreplicate(uppertri/0.95d^2*rovera^2,count)
	triang[*,where(theta_p gt !dpi,count)] = cmreplicate(lowertri/0.95d^2*rovera^2,count)
	x_sp = asin(triang)                  ;dim [n_rho_s,n_theta_p]
	;Assume that kappa(r) = kappa95  . constant elongation
	elong = dblarr(n_rho_s,n_theta_p)   ;dim [n_rho_s,n_theta_p]
	elong[*,where(theta_p le !dpi)] = upperelong
	elong[*,where(theta_p gt !dpi)] = lowerelong
	shafranov = cmreplicate(shafranov0*(1-(rho_s/aminor)^2),n_theta_p) ;dim [n_rho_s,n_theta_p] ;assume a parabolic Shafranov shift for simplicity.

	;grow some matrices so the calculation can be vectorized
	rho_sp = cmreplicate(rho_s,n_theta_p) ;minor radius in meters, not normalized
	theta_sp = extend(n_rho_s,theta_p)

	;eliminate some common subexpressions
	sin_theta_sp = sin(theta_sp)
	cos_theta_sp = cos(theta_sp)

	R_sp = Rmajor + shafranov + rho_sp*cos(theta_sp + x_sp*sin_theta_sp)
	Z_sp = elong*rho_sp*sin_theta_sp + z0
	dRdtheta = -rho_sp*(1 + x_sp*cos_theta_sp)*sin(theta_sp + x_sp*sin_theta_sp)
	dZdtheta = elong*rho_sp*cos_theta_sp

	;plot for verification
	if keyword_set(plot) then begin
	    for i=0,n_rho_s-1,20 do begin
		oploti,R_sp[i,*],Z_sp[i,*],/isotropic,icolor=1
	    endfor
	    stop
	endif

	;Search for the rho,theta coordinates closest to the R,Z coordinates of the chord
	;This will fail miserably for R,Z outside the outermost flux surface
	;Miller geometry only works for inner surfaces.
	npoint = chord.npoint
	chdist = dblarr(npoint)
	chrhoi = lonarr(npoint)
	chrmid = dblarr(npoint)
	chthetai = lonarr(npoint)
	print,'chord to RZ mapping'
	for i=0l,npoint-1 do begin
	    chdist0 = (R_sp - chord.r[i])^2 + (Z_sp - chord.z[i])^2    ;slow!
	    chdist[i] = min(chdist0,index)
	    index = array_indices(chdist0,index)
	    chrhoi[i] = index[0]
	    chthetai[i] = index[1]
	    chrmid[i] = R_sp[index[0],0]
	    if i mod 100 eq 0 then print,string(13b),i,npoint,format='(A,"iter ",I10," of ",I10,$)'
	endfor
	print,'done'

	;Restrict results to points in the plasma where Miller reconstruction is valid
	valid = where(chdist le 0.001) ;hardcoded distance.
	npoint = n_elements(valid)
	chx = chord.x[valid]
	chy = chord.y[valid]
	chz = chord.z[valid]
	chr = chord.r[valid]
	chphi = chord.phi[valid]
	chrmid = chrmid[valid]
	chrhoi = chrhoi[valid]
	chthetai = chthetai[valid]
	chvect = chord.vect[valid,*]
	chord = {npoint:npoint, x:chx, y:chy, z:chz, r: chr, phi: chphi, rmid:chrmid, rho:rovera[chrhoi], vect:chvect, ds:chord.ds, nray:chord.nray} ;update chord structure

	;calculate the B field at each mesh point
	;This part is confusing and needs to be looked over for mistakes
	B_tor = b_tor0*Rmajor/chr ;Basic vacuum toroidal field. (if you want to include poloidal current effects, better use EFIT)
	;Calculate B_poloidal from q profile using Miller eq 37
	s_k = 0 ;assume kappa(r) = kappa95
	s_delta = 2*triang / sqrt(1-triang^2)
	denom = cos(x_sp*sin_theta_sp) - shafranov0*2*rho_sp/aminor^2*cos_theta_sp + (s_k - s_delta*cos_theta_sp + (1+s_k)*x_sp*cos_theta_sp)*sin_theta_sp*sin(theta_sp + x_sp*sin_theta_sp)
	bpolcoef = 1 / elong * sqrt(sin(theta_sp + x_sp*sin_theta_sp)^2 * (1 + x_sp*cos_theta_sp)^2 + elong^2*cos_theta_sp^2)/denom    ;bpolcoef is Bpol*R/(dpsi/dr)
	;equation 16: q = f/(2*pi) * integral{dl/(R^2*Bpol)}
	; f ~= B0*R0
	; q = B0*R0/(2*pi) * total{1/(bpolcoef*R*dpsi/drho)}*2*pi/n_theta_p
	arclen = sqrt(dRdtheta^2 + dZdtheta^2)
	arclen[0,*] = 1d ; Calculation blows up at magnetic axis. Remove the divide by zero which causes problems later
	bpolint = total(arclen/(bpolcoef*R_sp^2),2)/n_theta_p * (B_tor0*Rmajor) ;q = bpolint / (dpsi/drho)  ; rho in meters
	;calculate (dpsi/drho) based on input q profile
	qinterpol = interpol(qprof,qprof_rho,rovera)
	dpsi = bpolint / qinterpol ;dim [n_rho_s]
	;calculate B_poloidal
	Bpol = bpolcoef/R_sp*cmreplicate(dpsi,n_theta_p)
	;calculate projections of B_poloidal along R and Z directions
	bR_sp = Bpol * dRdtheta/arclen
	bZ_sp = Bpol * dZdtheta/arclen
	Bz = bZ_sp[chrhoi,chthetai]
	Br = bR_sp[chrhoi,chthetai]
	end

	'miller': begin
	;New calculation of B(grid point) based on Miller model [Miller et al, PoP 5 4 (1998)]
	;This is quicker and more elegant than the old calculation
	elong0 = (chord.z ge z0)*(upperelong-lowerelong) + lowerelong
	tri0 = ((chord.z ge z0)*(uppertri-lowertri) + lowertri) / 0.95d^2
	z_map = (chord.z-z0)/elong0
	;rhoguess = sqrt((chord.r - Rmajor - 0.75d*shafranov0 + 0.25d*tri0) + z_map^2)
	rhoguess = sqrt((chord.r - Rmajor - 0.75d*shafranov0 + 0.25d*tri0)^2 + z_map^2)
	sinthetaguess = z_map/(chord.r - Rmajor - 0.75d*shafranov0 + 0.25d*tri0)
	if keyword_set(plot) then begin
	    figure,1
	    clf,/all
	    figure,2
	    clf,/all
	endif
	for i=0,10 do begin
	    shafranovguess = shafranov0*(1-(rhoguess/aminor)^2)
	    triangterm = asin(tri0*(rhoguess/aminor)^2)*sinthetaguess
	    ;triangshift = z_map*triangterm ;first order Taylor on rho*cos(theta+x*sin(theta)) about rho(cos(theta))
	    triangshift = z_map*triangterm + 0.5d*(chord.r-Rmajor-shafranovguess)*triangterm^2 ;second order
	    ;triangshift = z_map*triangterm + 0.5d*(chord.r-Rmajor-shafranovguess)*triangterm^2 - z_map*triangterm^3/6d
	    r_map = chord.r - Rmajor - shafranovguess + triangshift
	    rhoguess = sqrt(r_map^2 + z_map^2)
	    sinthetaguess = z_map/rhoguess
	    if keyword_set(plot) then begin
		figure,1
		oploti,rhoguess,icolor=i
		figure,2
		oploti,sinthetaguess,icolor=i
	    endif
	endfor

	;Restrict results to points in the plasma where Miller reconstruction is valid
	valid = where(rhoguess lt aminor)
	npoint = n_elements(valid)
	chx = chord.x[valid]
	chy = chord.y[valid]
	chz = chord.z[valid]
	chr = chord.r[valid]
	chphi = chord.phi[valid]
	chrho = rhoguess[valid]/aminor
	chtheta = asin(sinthetaguess[valid])
	chrmid = Rmajor + rhoguess[valid] + shafranovguess[valid]
	chvect = chord.vect[valid,*]
	chord = {npoint:npoint, x:chx, y:chy, z:chz, r: chr, phi: chphi, rmid:chrmid, rho:chrho, vect:chvect, ds:chord.ds, nray:chord.nray} ;update chord structure
	triangterm = triangterm[valid]

	;Need to calculate the B field on a grid to get dpsi
	n_theta_p = 400
	n_rho_s = 400
	theta_p = maken(0d,2*!dpi,n_theta_p) ;dim [n_theta_p]
	rho_s = maken(0d,aminor,n_rho_s)     ;dim [n_rho_s]
	rovera = rho_s/aminor               ;dim [n_rho_s]
	;assume delta(r) = delta95/0.95^2 * (r/a)^2
	triang = dblarr(n_rho_s,n_theta_p)  ;dim [n_rho_s,n_theta_p]
	triang[*,where(theta_p le !dpi,count)] = cmreplicate(uppertri/0.95d^2*rovera^2,count)
	triang[*,where(theta_p gt !dpi,count)] = cmreplicate(lowertri/0.95d^2*rovera^2,count)
	x_sp = asin(triang)                  ;dim [n_rho_s,n_theta_p]
	;Assume that kappa(r) = kappa95  . constant elongation
	elong = dblarr(n_rho_s,n_theta_p)   ;dim [n_rho_s,n_theta_p]
	elong[*,where(theta_p le !dpi)] = upperelong
	elong[*,where(theta_p gt !dpi)] = lowerelong
	shafranov = cmreplicate(shafranov0*(1-(rho_s/aminor)^2),n_theta_p) ;dim [n_rho_s,n_theta_p] ;assume a parabolic Shafranov shift for simplicity.
	;grow some matrices so the calculation can be vectorized
	rho_sp = cmreplicate(rho_s,n_theta_p) ;minor radius in meters, not normalized
	theta_sp = extend(n_rho_s,theta_p)
	;eliminate some common subexpressions
	sin_theta_sp = sin(theta_sp)
	cos_theta_sp = cos(theta_sp)
	R_sp = Rmajor + shafranov + rho_sp*cos(theta_sp + x_sp*sin_theta_sp)
	dRdtheta = -rho_sp*(1 + x_sp*cos_theta_sp)*sin(theta_sp + x_sp*sin_theta_sp)
	dZdtheta = elong*rho_sp*cos_theta_sp
	;Calculate B_poloidal from q profile using Miller eq 37
	s_k = 0 ;assume kappa(r) = kappa95
	s_delta = 2*triang / sqrt(1-triang^2)
	denom = cos(x_sp*sin_theta_sp) - shafranov0*2*rho_sp/aminor^2*cos_theta_sp + (s_k - s_delta*cos_theta_sp + (1+s_k)*x_sp*cos_theta_sp)*sin_theta_sp*sin(theta_sp + x_sp*sin_theta_sp)
	bpolcoef = 1 / elong * sqrt(sin(theta_sp + x_sp*sin_theta_sp)^2 * (1 + x_sp*cos_theta_sp)^2 + elong^2*cos_theta_sp^2)/denom    ;bpolcoef is Bpol*R/(dpsi/dr)
	;equation 16: q = f/(2*pi) * integral{dl/(R^2*Bpol)}
	; f ~= B0*R0
	; q = B0*R0/(2*pi) * total{1/(bpolcoef*R*dpsi/drho)}*2*pi/n_theta_p
	arclen = sqrt(dRdtheta^2 + dZdtheta^2)
	arclen[0,*] = 1d ; Calculation blows up at magnetic axis. Remove the divide by zero which causes problems later
	bpolint = total(arclen/(bpolcoef*R_sp^2),2)/n_theta_p * (B_tor0*Rmajor) ;q = bpolint / (dpsi/drho)  ; rho in meters
	;calculate (dpsi/drho) based on input q profile
	qinterpol = interpol(qprof,qprof_rho,rovera)
	dpsi = bpolint / qinterpol ;dim [n_rho_s]

	;now recalculate B field on the chord using dpsi above
	B_tor = b_tor0*Rmajor/chr ;Basic vacuum toroidal field. (if you want to include poloidal current effects, better use EFIT)
	;Calculate B_poloidal from q profile using Miller eq 37
	s_k = 0 ;assume kappa(r) = kappa95
	triang = tri0*chrho^2
	s_delta = 2*triang / sqrt(1-triang^2)
	costheta = cos(chtheta)
	sintheta = sinthetaguess[valid]
	denom = cos(triangterm) - shafranov0*2*chrho/aminor*cos(chtheta) + (s_k - s_delta*costheta + (1+s_k)*asin(triang)*costheta)*sintheta*sin(chtheta + triangterm)
	bpolcoef = 1 / elong0 * sqrt(sin(chtheta + triangterm)^2 * (1 + asin(triang)*costheta)^2 + elong0^2*costheta^2)/denom    ;bpolcoef is Bpol*R/(dpsi/dr)
	;equation 16: q = f/(2*pi) * integral{dl/(R^2*Bpol)}
	; f ~= B0*R0
	; q = B0*R0/(2*pi) * total{1/(bpolcoef*R*dpsi/drho)}*2*pi/n_theta_p
	dRdtheta = -chrho*(1 + asin(triang)*costheta)*sin(chtheta + triangterm)
	dZdtheta = elong0*chrho*costheta
	arclen = sqrt(dRdtheta^2 + dZdtheta^2)
	invalid = where(arclen eq 0)
	if invalid[0] ne -1 then arclen[invalid] = 1
	dpsi = interpol(dpsi,rho_s,chrho*aminor)
	Bpol = bpolcoef/chord.r*dpsi
	;calculate projections of B_poloidal along R and Z directions
	bR = Bpol * dRdtheta/arclen
	bZ = Bpol * dZdtheta/arclen
	end

	'efit_file': begin
	if file_search(efit_file) eq efit_file then begin
	    gdata = read_geqdsk(filename=efit_file,/quiet)
	endif else message,'Cannot find EFIT geqdsk file (efit_file): '+efit_file

	zind = mindex(gdata.zgrid,0)
	psimidplane = gdata.psirz[*,zind]
	chordr_ind = interpol(dindgen(n_elements(gdata.rgrid)),gdata.rgrid,chord.r,/spline)
	chordz_ind = interpol(dindgen(n_elements(gdata.zgrid)),gdata.zgrid,chord.z,/spline)
	psichord = interpolate(reform(gdata.psirz[*,*]),chordr_ind,chordz_ind,/cubic)
	chrmid = interpol(gdata.rgrid,psimidplane,psichord)

	sidif = gdata.ssibdry-gdata.ssimag
	psichordnorm = (psichord-gdata.ssimag)/sidif
	chrho = sqrt(psichordnorm)

	;Restrict results to points inside the LCFS
	valid = where(chrho le 1)
	if valid[0] eq -1 then begin
	    message,'Viewing chord misses plasma!'
	endif
	npoint = n_elements(valid)
	chx = chord.x[valid]
	chy = chord.y[valid]
	chz = chord.z[valid]
	chr = chord.r[valid]
	chphi = chord.phi[valid]
	chrmid = chrmid[valid]
	chrho = chrho[valid]
	chvect = chord.vect[valid,*]
	chord = {npoint:npoint, x:chx, y:chy, z:chz, r: chr, phi: chphi, rmid:chrmid, rho:chrho, vect:chvect, ds:chord.ds, nray:chord.nray} ;update chord structure

	efit_psi2b,chr,chz,br=br,bz=bz,bt=b_tor,bmod=b_tot,rgrid=gdata.rgrid,zgrid=gdata.zgrid,tgrid=[1],psirz=gdata.psirz,fpol=gdata.fpol,ssibry=gdata.ssibdry,ssimag=gdata.ssimag
	end

	'efit_mds': begin
	;if file_search('efit_east_61973.sav') eq 'efit_east_61973.sav' then begin ;mds access too slow for testing, so we saved it
	;    restore,'efit_east_61973.sav'
	;endif else begin
	    print,'getting efit from mdsplus'
	    connect_efit_mdsplus
	    psirz = mdsvalue('\psirz')
	    zdim = mdsvalue('\z')
	    rdim = mdsvalue('\r')
	    tdim = mdsvalue('\time')
	    fpol = mdsvalue('\fpol') ;array order is opposite on EAST and C-Mod for some reason
	    si0 = mdsvalue('\ssimag')  ;what is the difference between ssimag and simagx ?
	    sib = mdsvalue('\ssibry')  ;what is the difference between ssibry and psibdy ?
	;    save,psirz,zdim,rdim,tdim,fpol,si0,sib,filename='efit_east_61973.sav'
	;endelse
	tind = mindex(tdim,efit_time)
	zind = mindex(zdim,0)
	psimidplane = psirz[*,zind,tind]
	chordr_ind = interpol(dindgen(n_elements(rdim)),rdim,chord.r,/spline)
	chordz_ind = interpol(dindgen(n_elements(zdim)),zdim,chord.z,/spline)
	psichord = interpolate(reform(psirz[*,*,tind]),chordr_ind,chordz_ind,/cubic)
	chrmid = interpol(rdim,psimidplane,psichord)

	sidif = sib-si0
	psichordnorm = (psichord-si0[tind])/sidif[tind]
	chrho = sqrt(psichordnorm)

	;Restrict results to points inside the LCFS
	valid = where(chrho le 1)
	if valid[0] eq -1 then begin
	    message,'Viewing chord misses plasma!'
	endif
	npoint = n_elements(valid)
	chx = chord.x[valid]
	chy = chord.y[valid]
	chz = chord.z[valid]
	chr = chord.r[valid]
	chphi = chord.phi[valid]
	chrmid = chrmid[valid]
	chrho = chrho[valid]
	chvect = chord.vect[valid,*]
	chord = {npoint:npoint, x:chx, y:chy, z:chz, r: chr, phi: chphi, rmid:chrmid, rho:chrho, vect:chvect, ds:chord.ds, nray:chord.nray} ;update chord structure

	efit_psi2b,chr,chz,t1=efit_time,br=br,bz=bz,bt=b_tor,bmod=b_tot,rgrid=rdim,zgrid=zdim,tgrid=tdim,psirz=psirz,fpol=fpol,ssibry=sib,ssimag=si0
	end
    endcase

    Bx = Br*cos(chphi)-B_tor*sin(chphi)
    By = Br*sin(chphi)+B_tor*cos(chphi)
    B_tot = sqrt(b_tor^2+br^2+bz^2)
    B_hat = [[bx],[by],[bz]]/rebin(B_tot,npoint,3) ; dim [npoint,3]: unit vectors in B direction
;    B_vect_full = [[bx],[by],[bz]] ; dim [npoint,3]: B field vector

    ;unit vector in purely toroidal direction with sign given by B field
    ;b_tor_vect = [-sign(by[0])*sin(spot_pos[1]),sign(by[0])*cos(spot_pos[1]),0]

    ;Construct E field. Only consider radial E fields. Toroidal E field from loop voltage is negligible and greatly complicates analysis
    case Ermodel of
	'none': begin
	Ex = fltarr(npoint)
	Ey = fltarr(npoint)
	Ez = fltarr(npoint)
	Er = fltarr(npoint)
	E_tot = fltarr(npoint)
	E_hat = fltarr(npoint,3)
	end
	'tabulated': begin
	Ermid = interpol(Er_val,Er_r,chrmid) ;[npoint] Er_val should be negative for inward Er
	bsign = sign(by[(where(chrmid ge Rmajor+aminor/2))[0]])
	Er = Ermid*bz/b_tot*bsign ;build Er and Ez perpendicular to B_r, B_z
	Ez = Ermid*br/b_tot*bsign
	Ex = Er*cos(chphi)
	Ey = Er*sin(chphi)
	E_tot = sqrt(Er^2+Ez^2)
	E_hat = [[Ex],[Ey],[Ez]]/rebin(E_tot,npoint,3) ; dim [nstep,3]: unit vectors in E direction
	end
	'forcebalance': begin
	; given by 1/(n*Z*e) * dp_i/dr - v_{theta,i}*B_phi + v_{phi,i}*B_theta  (c.f. Wesson 2004 4.19.4)
	diamagnetic = interpol(Er_diamagnetic,Er_r,chrmid)
	v_pol = interpol(Er_vtheta,Er_r,chrmid)
	v_tor = interpol(Er_vphi,Er_r,chrmid)
	stop,'need to look at this code'
	Ermid = diamagnetic - v_pol*B_tor + v_tor*B_pol
	bsign = sign(by[(where(chrmid ge Rmajor+aminor/2))[0]])
	Er = Ermid*bz/b_tot*bsign ;build Er and Ez perpendicular to B_r, B_z
	Ez = Ermid*br/b_tot*bsign
	Ex = Er*cos(chphi)
	Ey = Er*sin(chphi)
	E_tot = sqrt(Er^2+Ez^2)
	E_hat = [[Ex],[Ey],[Ez]]/rebin(E_tot,npoint,3) ; dim [nstep,3]: unit vectors in E direction
	end
	'romannikov': begin
	    message,'TODO'
	end
    endcase

    field = {Bx:Bx, By:By, Bz:Bz, Br:Br, B_tor:B_tor, B_tot:B_tot, B_hat:B_hat, Ex:Ex, Ey:Ey, Ez:Ez, Er:Er, E_tot:E_tot, E_hat:E_hat}
end
