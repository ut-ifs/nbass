@starkzeeman.pro

;Read the nonstatistical line amplitude calculation results from an external file
;external file should be a text file with 13 or more columns
;each row is a point along the beam-line, with distance to grid given by column 0
function get_pop_ratio,filename,chdist
    ;statistical_amps is the line amplitudes for a pure Stark effect with statistically populated upper states
                       ;pi4        sigma6       pi8           sigma1     pi3        sigma5       -pi2        sigma0     pi2         -sigma5      -pi3       -sigma1    -pi8          -pi6         -pi4
    statistical_amps = [0.17826087,0.0019088017,0.00010604454,0.20530223,0.24432662,0.0016967126,0.077306469,0.58218452,0.077306469,0.0016967126,0.24432662,0.20530223,0.00010604454,0.0019088017,0.17826087]
    data = read_ascii(filename,comment_symbol='#')
    distance = reform(data.field01[0,*]) ;meters
    ;density = reform(data.field01[1,*])
    ;temperature = reform(data.field01[2,*])
    ;time = reform(data.field01[3,*])
    ;n0 = reform(data.field01[4,*])
    sigma0 = reform(data.field01[5,*])
    sigma1 = reform(data.field01[6,*])
    sigma5 = reform(data.field01[7,*])
    sigma6 = reform(data.field01[8,*])
    pi2 = reform(data.field01[9,*])
    pi3 = reform(data.field01[10,*])
    pi4 = reform(data.field01[11,*])
    pi8 = reform(data.field01[12,*])

    chdist2 = ((chdist > min(distance)) < max(distance)) ;cap results so interpol doesn't try to extrapolate
    s0 = interpol(sigma0,distance,chdist2)
    s1 = interpol(sigma1,distance,chdist2)
    s5 = interpol(sigma5,distance,chdist2)
    s6 = interpol(sigma6,distance,chdist2)
    p2 = interpol(pi2,distance,chdist2)
    p3 = interpol(pi3,distance,chdist2)
    p4 = interpol(pi4,distance,chdist2)
    p8 = interpol(pi8,distance,chdist2)

    ;sort them into the same order that the Stark-Zeeman function uses for the 15 lines
    spectrum = [[p4],[s6],[p8],[s1],[p3],[s5],[p2],[s0],[p2],[s5],[p3],[s1],[p8],[s6],[p4]]
    ratio = transpose(spectrum / extend(n_elements(chdist),statistical_amps))
    return,ratio
end

function synth_bes,chord,field,beam,angle,lorentzE_norm,dangle_fg,dangle_ap,dangle_sp,profiles,bes_data
    @param.idl
    @detector_param.idl
    @beam_param.idl
    @view_param.idl
    @mesh_param.idl
    n_comp = beam.n_comp

    case beam.atom of
	'H': H_alpha = 6562.7952d        ; A        wavelength of H-alpha line (Air)
	'D': H_alpha = 6561.010d         ; A        wavelength of D-alpha line (Air)
    endcase

    ;Only do the calculation on chord points where there is significant beam
    ;n_bpoints is the number of points along the chord with beam
    fb = where(beam.n[*,0]/max(beam.n[*,0],fb0) gt 1d-5,n_bpoints)  ;fb stands for found beam

    n_bin = n_elements(wve)
    if keyword_set(autocenter) then begin
	;wve += H_alpha*beam.gamma[0]*(1-beam.beta[0]*cos(mean(angle.alpha[fb]))) - mean(wve) ;rough recentering
	wve += H_alpha/beam.gamma[0]/(1+beam.beta[0]*cos(mean(angle.alpha[fb]))) - mean(wve) ;rough recentering
	print,'Moving detector wavelength center to ',num2str(mean(wve),3)
    endif
    pixboundl = [wve[0],(wve[0:n_bin-1]+wve[1:*])/2]
    pixboundr = [(wve[0:n_bin-1]+wve[1:*])/2,wve[n_bin-1]]
    pixwidth = abs(pixboundr-pixboundl)
    ds1 = chord.ds / double(chord.nray)

    ;choice of models for relative line amplitudes, and choice of models for total H-alpha emission amplitude
    print,'Calculate emission amplitude'
    ;Apply beam excited population to line amplitudes
    if keyword_set(beamingas) then begin
	message,'cannot do Miller calculation for beam into gas'
    endif else begin
	amplitude = extend(15,transpose(beam.n3[fb,*]*1d6))  ;[15,n_comp,n_bpoints]. 1/m^3
	if keyword_set(nonstatistical) then begin
	    chdist = sqrt((chord.x-beam.grid_pos[0])^2+(chord.y-beam.grid_pos[1])^2+(chord.z-beam.grid_pos[2])^2)
	    ns_ratio = get_pop_ratio(nonstat_file,chdist[fb])
	    ;We only use the nonstatistical model for the full energy component (component 0)
	    ;The other components are assumed to be statistical due to low energy
	    amplitude[*,0,*] *= ns_ratio
	endif
    endelse
    amplitude /= 4*!dpi ; [15,n_comp,n_bpoints] 1/m^3/sr

    if linemodel eq 'delta' then needwidth=0 else begin
	needwidth=1
	case chordmodel of
	'ray'         : begin dalpha_ap = dangle_ap.alpha[fb]
			      dalpha_sp = dangle_sp.alpha[fb]
			end
	'grid'        : begin
		case ap_shape of
		    'point': dalpha_ap = dangle_ap.alpha[fb]
		    'hex7':  dalpha_ap = dangle_ap.alpha[fb]/3
		    'hex19': dalpha_ap = dangle_ap.alpha[fb]/5
		    'grid':  dalpha_ap = dangle_ap.alpha[fb]
		endcase
		case sp_shape of
		    'point': dalpha_sp = dangle_sp.alpha[fb]
		    'hex7':  dalpha_sp = dangle_sp.alpha[fb]/3
		    'hex19': dalpha_sp = dangle_sp.alpha[fb]/5
		    'grid':  dalpha_sp = dangle_sp.alpha[fb]
		endcase
	    end
	endcase
	circularfix = sqrt(3d)/2 ;effective width of a circular distribution
	dalpha_ap *= circularfix
	dalpha_sp *= circularfix
	dalpha_fg = dangle_fg.alpha[fb] * circularfix
	if num_pini gt 1 then dalpha_fg = 0
    endelse
;    sigmas = [1,3,5,7,9,11,13]
;    pis =    [0,2,4,6,8,10,12,14]
;    print,'starkzeeman'
    shifts = dblarr(15,num_pini)
    stokes = dblarr(15,4,num_pini)
    spectrum = dblarr(n_bin,n_comp,4) ;last dim is stokes vector
    B_lor = rebin(field.b_tot[fb]*1d4,n_bpoints,n_comp)*sqrt(1+outer(sin(angle.theta[fb])^2,beam.gamma^2-1)) ;for simplicity, I used angle.theta here instead of calculating for each pini
    E_cgs = lorentzE_norm[fb,*] * 1d4 / 299792458d0 ;[npoint,n_comp]
    print,'loop size (plasma points, beam energies, source points):'+strtrim(n_bpoints,2)+'x'+strtrim(n_comp,2)+'x'+strtrim(num_pini,2)
    for i=0,n_bpoints-1 do begin
    ;for i=fb0-fb[0],fb0-fb[0] do begin ;special test
	;print,'.',format='(A,$)'
	if i mod 50 eq 0 then print,'i='+strtrim(i,2)+'/'+strtrim(n_bpoints,2)+' : '+strtrim(100.0*i/n_bpoints,2)+'%'
	if num_pini gt 1 then begin
	    pos = [chord.x[fb[i]],chord.y[fb[i]],chord.z[fb[i]]]  ;[3]
	    beam_vect = extend([num_pini],pos) - beam.pini_pos ;[num_pini, 3]
	    beam_vect_norm = sqrt(total(beam_vect^2,2))
	    beam_vect /= rebin(beam_vect_norm,num_pini,3) ;[num_pini,3]: unit vectors
	    vxB0 = -crossp_1d2d(field.b_tot[fb[i]]*field.b_hat[fb[i],*],beam_vect) * 1d4 ;vxb / relativisticbeta, cgs
	    piniEr = rebin([[field.Ex[fb[i]]],[field.Ey[fb[i]]],[field.Ez[fb[i]]]],num_pini,3) / 299792458d0 * 1d4 ;cgs
	    piniErperp = piniEr - beam_vect*rebin(total(piniEr*beam_vect,2),num_pini,3)
	    cos_alpha = total(rebin(chord.vect[fb[i],*],num_pini,3)*beam_vect,2) ;[num_pini]
	    for j=0,n_comp-1 do begin
		vxB = vxB0*beam.gamma[j]*beam.beta[j] ;[num_pini,3]
		piniE = vxB + piniEr + (beam.gamma[j]-1)*piniErperp ;[num_pini,3]
		piniE_norm = sqrt(total(piniE^2,2))
		piniE /= rebin(piniE_norm,num_pini,3)
		piniEx = total(piniE*rebin(angle.xvect[fb[i],*],num_pini,3),2) ;project into Isler coord
		piniEy = total(piniE*rebin(angle.yvect[fb[i],*],num_pini,3),2)
		piniphi = sign(piniEx)*acos(-piniEy)
		piniE_cgs = piniE_norm ;[npini]
		for k=0,num_pini-1 do begin
		    starkzeeman,B_lor[i,j],piniE_cgs[k],piniphi[k],angle.psi[fb[i]],shifts_one,amplitude_one,/quadratic,H_alpha=H_alpha,stokes=stokes_one,vectH=angle.vectH[fb[i],*],vectV=angle.vectV[fb[i],*]
		    shifts[*,k] = reform(shifts_one,15)
		    stokes[*,*,k] = reform(stokes_one,15,4) / 9D
		endfor
		;line_wve = (H_alpha+shifts)*extend([15],beam.gamma[j]*(1-beam.beta[j]*cos_alpha)) ;[15,num_pini]
		line_wve = (H_alpha+shifts)/extend([15],beam.gamma[j]*(1+beam.beta[j]*cos_alpha)) ;[15,num_pini]
		stokes *= outer(rebin(amplitude[*,j,i],15,4),reform(beam.pini_w[fb[i],*,j])) ;[15,n_comp,num_pini] ;pht/sec/m^3/sr

		if needwidth then begin
		    sin_alpha = sqrt(1-cos_alpha^2)
		    ;dline_wve_dalpha_old = (H_alpha+shifts)*extend([15],1/beam.gamma[j]*beam.beta[j]*sin_alpha) ;[15,num_pini]
		    dline_wve_dalpha = line_wve*extend([15],(beam.beta[j]*sin_alpha)/(1+beam.beta[j]*cos_alpha))
		    ;dline_wve_dbeta_old = (H_alpha+shifts)*extend([15],1/beam.gamma[j]*cos_alpha) ;[15,num_pini]
		    dline_wve_dbeta = line_wve*extend([15],(cos_alpha)/(1+beam.beta[j]*cos_alpha))
		    width2 = dline_wve_dalpha^2*(dalpha_ap[i]^2 + dalpha_sp[i]^2) + (dline_wve_dbeta*beam.beta[j]*beam_ripple/2)^2 + l_instr^2
		    width = sqrt(width2)
		endif

		case linemodel of
		'gaussian': begin
		    for bin=0,n_bin-1 do begin
			exp_fact = exp(-(wve[bin] - line_wve)^2/width2) / (width*sqrt(!dpi))
			for si=0,3 do begin
			    spectrum[bin,j,si] += total(stokes[*,si,*]*exp_fact*ds1)
			endfor
		    endfor
		end
		'erf': begin
		    ;save some time by assuming pixboundl[bin] = pixboundr[bin-1]
		    last = erf((pixboundl[0] - line_wve)/width)
		    for bin=0,n_bin-1 do begin
			;cur = erf((pixboundr[bin] - line_wve)/width) ;This is the computation bottleneck
			 erfin = ((pixboundr[bin] - line_wve)/width)               ; these four lines do the same thing
			 cur = double(sign(erfin))                                 ; as the line above, but tiny
			 valid = where(erfin gt -5.89d and erfin lt 5.89d,nvalid)  ; bit faster
			 if nvalid gt 0 then cur[valid] = erf(erfin[valid])        ;
			erf_diff = abs(0.5d*(last - cur) / pixwidth[bin]) ; [15,n_comp,n_bpoints] divide by pixel width to change pht/m^2/sec/sr/pix to pht/m^2/sec/sr/A to be consistent with linemodel='erf'.
			last = cur
			for si=0,3 do begin
			    spectrum[bin,j,si] += total(stokes[*,si,*]*erf_diff*ds1)
			endfor
		    endfor
;if total(finite(spectrum)) lt 8192 then stop
		end
		'delta': begin
		    line_idx = value_locate(pixboundl,line_wve) > 0 ;[15,n_comp,n_bpoints]
		    for bi=0,14 do begin
			for bj=0,num_pini-1 do begin
			    spectrum[line_idx[bi,bj],j,*] += stokes[bi,*,bj] / pixwidth[line_idx[bi,bj]] * ds1 ;pht/s/m^2/st
			endfor
		    endfor
		end
		endcase
	    endfor
	endif else begin ;num_pini = 1, so just use the beam centerline
	    for j=0,n_comp-1 do begin
		starkzeeman,B_lor[i,j],E_cgs[i,j],angle.phiL[fb[i]],angle.psi[fb[i]],shifts_one,amplitude_one,/quadratic,H_alpha=H_alpha,stokes=stokes_one,vectH=angle.vectH[fb[i],*],vectV=angle.vectV[fb[i],*]
		shifts[*,0] = reform(shifts_one,15)
		stokes[*,*,0] = reform(stokes_one,15,4) / 9D ; rescale to accept total population of n=3 level instead of population of sublevel
		;line_wve = (H_alpha+shifts)*extend([15],beam.gamma[j]*(1-beam.beta[j]*cos(angle.alpha[fb[i]]))) ;[15,num_pini]
		line_wve = (H_alpha+shifts)/extend([15],beam.gamma[j]*(1+beam.beta[j]*cos(angle.alpha[fb[i]]))) ;[15,num_pini]
		stokes *= rebin(amplitude[*,j,i],15,4) ;[15,n_comp,num_pini] ;pht/sec/m^3/sr

		if needwidth then begin
		    sin_alpha = sin(angle.alpha[fb[i]])
		    cos_alpha = cos(angle.alpha[fb[i]])
;		    dline_wve_dalpha_old = (H_alpha+shifts)*extend([15],1/beam.gamma[j]*beam.beta[j]*sin_alpha) ;[15,num_pini]
		    dline_wve_dalpha = line_wve*extend([15],(beam.beta[j]*sin_alpha)/(1+beam.beta[j]*cos_alpha))
;		    dline_wve_dbeta_old = (H_alpha+shifts)*extend([15],1/beam.gamma[j]*cos_alpha) ;[15,num_pini]
		    dline_wve_dbeta = line_wve*extend([15],(cos_alpha)/(1+beam.beta[j]*cos_alpha))
		    width2 = dline_wve_dalpha^2*(dalpha_ap[i]^2 + dalpha_sp[i]^2 + dalpha_fg[i]^2) + (dline_wve_dbeta*beam.beta[j]*beam_ripple/2)^2 + l_instr^2
		    width = sqrt(width2)
		endif

		case linemodel of
		'gaussian': begin
		    for bin=0,n_bin-1 do begin
			exp_fact = exp(-(wve[bin] - line_wve)^2/width2) / (width*sqrt(!dpi))
			for si=0,3 do begin
			    spectrum[bin,j,si] += total(stokes[*,si,*]*exp_fact*ds1)
			endfor
		    endfor
		end
		'erf': begin
		    for bin=0,n_bin-1 do begin
			erf_diff = abs(0.5d*(erf((pixboundl[bin] - line_wve)/width) - erf((pixboundr[bin] - line_wve)/width)) / pixwidth[bin]) ; [15,n_comp,n_bpoints] divide by pixel width to change pht/m^2/sec/sr/pix to pht/m^2/sec/sr/A to be consistent with linemodel='erf'.
			for si=0,3 do begin
			    spectrum[bin,j,si] += total(stokes[*,si,*]*erf_diff*ds1)
			endfor
		    endfor
		end
		'delta': begin
		    line_idx = value_locate(pixboundl,line_wve) > 0 ;[15,n_comp,n_bpoints]
		    for bi=0,14 do begin
			for bj=0,num_pini-1 do begin
			    spectrum[line_idx[bi,bj],j,*] += stokes[bi,*,bj] / pixwidth[line_idx[bi,bj]] * ds1 ;pht/s/m^2/st
			endfor
		    endfor
		end
		endcase
	    endfor
	endelse
    endfor

    if keyword_set(plot) then begin
	figure,1
	clf,/all
	ploti,wve,spectrum[*,*,0],charsize=1.5,xtitle='wavelength (A)',ytitle='emissivity ph/m^2/s/sr'
    endif

;    if n_comp ne 1 then begin
;	stokes_s0 = total(spectrum[*,*,0],2)
;	stokes_s1 = total(spectrum[*,*,1],2)
;	stokes_s2 = total(spectrum[*,*,2],2)
;	stokes_s3 = total(spectrum[*,*,3],2)
;    endif else begin
	stokes_s0 = spectrum[*,*,0]
	stokes_s1 = spectrum[*,*,1]
	stokes_s2 = spectrum[*,*,2]
	stokes_s3 = spectrum[*,*,3]
;    endelse

    poldegree = sqrt(stokes_s1^2 + stokes_s2^2 + stokes_s3^2)/stokes_s0
    spsi = atan(stokes_s2,stokes_s1)/2
    schi = atan(stokes_s3,sqrt(stokes_s1^2+stokes_s2^2))/2
    ;bes_data = {spectrum:spectrum,stokes_s0:stokes_s0,stokes_s1:stokes_s1,stokes_s2:stokes_s2,stokes_s3:stokes_s3,poldegree:poldegree,spsi:spsi,schi:schi}
    bes_data = {stokes_s0:stokes_s0,stokes_s1:stokes_s1,stokes_s2:stokes_s2,stokes_s3:stokes_s3,poldegree:poldegree,spsi:spsi,schi:schi}

    if n_comp ne 1 then return,total(stokes_s0,2) else return,stokes_s0 ;pht/m^2/sR/sec/A
end
