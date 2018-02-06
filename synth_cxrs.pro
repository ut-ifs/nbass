@zeeman_helium.pro
@zeeman_boron.pro
@blom.pro

;returns the wavelength and amplitude of component lines
;input: species, B, psi
;output:
;   wavelength
;   amp
pro normalize_zeeman,species,wvezeeman,ampzeeman,B,psi
    switch species of
	'H_4_3':
	'D_4_3': begin
	    wvezeeman = dblarr(15,n_elements(B))
	    ampzeeman = dblarr(15,n_elements(B))
	    message,'incomplete'
	    break
	end
	'He3_4_3':
	'He4_4_3': begin
	    wvezeeman = dblarr(146,n_elements(B))
	    ampzeeman = dblarr(146,n_elements(B))
	    for i=0,n_elements(B)-1 do begin
		zeeman_helium,B[i],invwve,polarization,ampz,ampp,ampm,species=(species eq 'He3_4_3' ? 'He3' : 'He4')
		wavelength = 1d10 / invwve
		wavelength = wavelength / (1d0 + 2.735182d-4 + 131.4182d/wavelength^2+2.76249d8/wavelength^4) ;Vacuum to air. Morton (1991, ApJS, 77, 119)
		amp = sin(psi[i])^2*ampz+(1d0 + cos(psi[i])^2)*(ampp+ampm)/2D
		amp /= total(amp) ;normalize
		wvezeeman[*,i] = wavelength
		ampzeeman[*,i] = amp
	    endfor
	    break
	end
	'B_7_6': begin
	    wvezeeman = dblarr(560,n_elements(B))
	    ampzeeman = dblarr(560,n_elements(B))
	    for i=0,n_elements(B)-1 do begin
		zeeman_boron,B[i],wavelength,amplitude,pi_comp_arr
		amp = sin(psi[i])^2*amplitude*pi_comp_arr + (1d0 + cos(psi[i])^2)*amplitude*(1-pi_comp_arr)
		amp /= total(amp)
		wvezeeman[*,i] = wavelength
		ampzeeman[*,i] = amp
	    endfor
	    break
	end
    end
end

;returns the wavelength and amplitude of component lines
;input: species, B, psi
;vectH and vectV are vectors perpendicular to view direction, transformed into the local magnetic coordinate system (z = b_vect, view in xz plane)
;output:
;   wavelength
;   amp
pro normalize_zeeman_stokes,species,wvezeeman,s0,s1,s2,s3,B,vectH,vectV
    ;calculate the Stokes vectors for the zeeman pi and sigma components
    ;pi lines are linearly polarized in the z direction
    pi_s0 = vectH[*,2]^2 + vectV[*,2]^2 ;should be exactly equal to pi_s0
    pi_s1 = vectH[*,2]^2 - vectV[*,2]^2
    pi_s2 = 2*vectH[*,2]*vectV[*,2] ;not sure about sign
    ;pi_s0_1 = sqrt(pi_s1^2+pi_s2^2) ;should be exactly equal to pi_s0
    ;pi_s0_2 = sin(psi)^2 ;should be exactly equal to pi_s0

    ;sigma lines are circularly polarized in the xy (magnetic coordinates) plane
    ;The y component (magnetic coordinate) is phase delayed by +- pi/2 from the x component
    ;So, the y component should be understood as imaginary; the complex conjugate will flip the sign
    ;Stokes s2 = 2*Re(Ex Ey*), where Ex and Ey are in the (viewing coordinates)
    ;Stokes s3 = -2*Im(Ex Ey*), where Ex and Ey are in the (viewing coordinates)

    sigma_s0 = 0.5d*((vectH[*,0]^2+vectH[*,1]^2) + (vectV[*,0]^2+vectV[*,1]^2))
    sigma_s1 = 0.5d*((vectH[*,0]^2+vectH[*,1]^2) - (vectV[*,0]^2+vectV[*,1]^2))
    sigma_s2 = (vectH[*,0]*vectV[*,0] + vectH[*,1]*vectV[*,1])
    sigma_s3p = -(-vectH[*,0]*vectV[*,1] + vectH[*,1]*vectV[*,0])
    sigma_s3m = (-vectH[*,0]*vectV[*,1] + vectH[*,1]*vectV[*,0])
    ;sigma_s0_1 = sqrt(sigma_s1^2+sigma_s2^2+sigma_s3p^2) ;should be exactly equal to sigma_s0
    ;sigma_s0_2 = 0.5 + cos(psi)^2/2                      ;should be exactly equal to sigma_s0
    ;sigma_s3check = cos(psi) ;should be exactly equal to sigma_s3p

    switch species of
	'H_4_3':
	'D_4_3': begin
	    wvezeeman = dblarr(15,n_elements(B))
	    ampzeeman = dblarr(15,n_elements(B))
	    message,'incomplete'
	    break
	end
	'He3_4_3':
	'He4_4_3': begin
	    wvezeeman = dblarr(146,n_elements(B))
	    s0 = dblarr(146,n_elements(B))
	    s1 = dblarr(146,n_elements(B))
	    s2 = dblarr(146,n_elements(B))
	    s3 = dblarr(146,n_elements(B))
	    for i=0,n_elements(B)-1 do begin
		zeeman_helium,B[i],invwve,polarization,ampz,ampp,ampm,species=(species eq 'He3_4_3' ? 'He3' : 'He4')
		wavelength = 1d10 / invwve
		wavelength = wavelength / (1d0 + 2.735182d-4 + 131.4182d/wavelength^2+2.76249d8/wavelength^4) ;Vacuum to air. Morton (1991, ApJS, 77, 119)
		wvezeeman[*,i] = wavelength
		norma = total(ampp+ampm+ampz)/1.5d ;scale so that total(s0,1) = 1
		ampp /= norma
		ampm /= norma
		ampz /= norma
		s0[*,i] = ampz*pi_s0 + ampp*sigma_s0 + ampm*sigma_s0
		s1[*,i] = ampz*pi_s1 + ampp*sigma_s1 + ampm*sigma_s1
		s2[*,i] = ampz*pi_s2 + ampp*sigma_s2 + ampm*sigma_s2
		s3[*,i] = ampp*sigma_s3p + ampm*sigma_s3m
	    endfor
	    break
	end
	'B_7_6': begin
	    wvezeeman = dblarr(560,n_elements(B))
	    s0 = dblarr(560,n_elements(B))
	    s1 = dblarr(560,n_elements(B))
	    s2 = dblarr(560,n_elements(B))
	    s3 = dblarr(560,n_elements(B))
	    for i=0,n_elements(B)-1 do begin
		zeeman_boron,B[i],wavelength,amplitude,pi_comp_arr,pol
		;note, total amplitude of each sigma component seems half that of the pi component
		amp = sin(psi[i])^2*amplitude*pi_comp_arr + (1d0 + cos(psi[i])^2)*amplitude*(1-pi_comp_arr)
		amp /= total(amp)

		wvezeeman[*,i] = wavelength
		ind = where(pol eq 0)
		s0[ind,i] = amplitude[ind]*pi_s0
		s1[ind,i] = amplitude[ind]*pi_s1
		s2[ind,i] = amplitude[ind]*pi_s2
		s3[ind,i] = 0
		ind = where(pol eq 1)
		s0[ind,i] = amplitude[ind]*sigma_s0
		s1[ind,i] = amplitude[ind]*sigma_s1
		s2[ind,i] = amplitude[ind]*sigma_s2
		s3[ind,i] = amplitude[ind]*sigma_s3p
		ind = where(pol eq -1)
		s0[ind,i] = amplitude[ind]*sigma_s0
		s1[ind,i] = amplitude[ind]*sigma_s1
		s2[ind,i] = amplitude[ind]*sigma_s2
		s3[ind,i] = amplitude[ind]*sigma_s3m
	    endfor
	    break
	end
    end
end
;pro zeemanstokestest
;    angle = 2*!dpi*randomu(seed,10,/double)
;    zz = 2*randomu(seed,10,/double)-1
;    vectH = dblarr(10,3)
;    vectH[*,0] = sqrt(1-zz^2)*cos(angle)
;    vectH[*,1] = sqrt(1-zz^2)*sin(angle)
;    vectH[*,2] = zz
;    angle2 = 2*!dpi*randomu(seed,10,/double)
;    zz2 = 2*randomu(seed,10,/double)-1
;    vect2 = dblarr(10,3)
;    vect2[*,0] = sqrt(1-zz2^2)*cos(angle2)
;    vect2[*,1] = sqrt(1-zz2^2)*sin(angle2)
;    vect2[*,2] = zz2
;    vectV = crossp_2d2d(vect2,vectH)
;    vectV /= rebin(sqrt(total(vectV^2,2)),10,3)
;    vectS = crossp_2d2d(vectH,vectV)
;    psi = reform(acos(vectS[*,2]))
;    normalize_zeeman_stokes,'B_7_6',wve,amp,5.6,vectH,vectV,psi
;end

;result in ph/(m^2*s*sr)
;cs_effects = cross section effect corrections
;input:
;  wve (A)
;  ds  (cm)
;  n_e (1/cm^3)
;  t_e (eV)
;  alcbeam_results
;  l_instr (A)
;mass = atomic mass of beam
;todo, add H to D ratio and split line
function synth_cxrs,chord,field,beam,angle,profiles,cxrs_data,species
    @param.idl
    @view_param.idl
    @detector_param.idl
    @equil_param.idl

    ;Data from von Hellermann PPCF 1995, with some tweaked wavelengths
    ; coef = [wavelength(A), Q_0(10^-15 m^3 s^-1), E_m(keV/amu), p', q']
    case species of
	'H_3_2': coef = [6562.7952d, 3.54d, 33.3d, 1.39d, 3.75d]
	'D_3_2': coef = [6561.010d, 3.54d, 33.3d, 1.39d, 3.75d]
	'H_4_2': coef = [4861.3298d, 0.94d, 33.4d, 1.63d, 3.93d]
	'D_4_2': coef = [4860.007d, 0.94d, 33.4d, 1.63d, 3.93d]
	'He3_4_3': coef = [4685.9150d, 9.24d, 36.3d, 2.34d, 4.65d]
	'He4_4_3': coef = [4685.7048d, 9.24d, 36.3d, 2.34d, 4.65d]
	'B_7_6': coef = [4944.6d, 21.60d, 41.5d, 3.0d, 4.61d]
	'C_8_7': coef = [5290.5d, 23.70d, 49.2d, 3.05d, 4.76d]
	else: message,species+' cxrs cross section not implemented'
    endcase
    wved = coef[0]

    case species of
	'H_3_2' : mc2kev = 938262.0813d ;mass of proton * c^2 in keV (2014 CODATA recommended values)
	'D_3_2' : mc2kev = 1865612.928d ;mass of deuteron * c^2 in keV (2014 CODATA recommended values)
	'H_4_2'  : mc2kev = 938262.0813d ;mass of proton * c^2 in keV (2014 CODATA recommended values)
	'D_4_2'  : mc2kev = 1865612.928d ;mass of deuteron * c^2 in keV (2014 CODATA recommended values)
	'He3_4_3': mc2kev = 2808920.906d ;mass of He3 nucleus * c^2 in keV
	'He4_4_3': mc2kev = 3727379.24082d ;mass of He4 nucleus * c^2 in keV
	'B_7_6'  : mc2kev = 10.81d * 931494.0954d
		;boron 10 = 10.012937u, boron 11 = 11.0093055u
		;5 electron binding energy ~ 0.67098 keV. 5 electron mass = 2554.9947305 keV
	'C_8_7'  : mc2kev = 12.01d * 931494.0954d
	else: message,species+' ion mass not implemented'
    endcase
    c = 29979245800d            ; cm/s     speed of light in vacuum

    case species of
	'H_3_2': indion = (where(profiles.ions eq 'H'))[0]
	'D_3_2': indion = (where(profiles.ions eq 'D'))[0]
	'He3_4_3': indion = (where(profiles.ions eq 'He'))[0]
	'He4_4_3': indion = (where(profiles.ions eq 'He'))[0]
	'B_7_6': indion = (where(profiles.ions eq 'B'))[0]
	'C_8_7': indion = (where(profiles.ions eq 'C'))[0]
	else: message,species+' ion not implemented'
    endcase
    if indion eq -1 then message,species+' ion density missing from profiles data'

    case beam.atom of
	'H': mass = 1d
	'D': mass = 1.9990075d
    endcase

    ;Only do the calculation on chord points where there is significant beam
    ;n_bpoints is the number of points along the chord with beam
    fb = where(beam.n[*,0]/max(beam.n[*,0],fb0) gt 1d-5,n_bpoints)  ;fb stands for found beam
    print,'# chordbeam points',n_bpoints
    n_comp = beam.n_comp

    Ecol = beam.E / mass  ;keV/amu [4]
    xparm = Ecol/coef[2] ; Hellermann energy parameter x = energy/Em, Em ~= energy of max collision ~ 33.3keV

    ;need to divide Ecol by A, but A=1 for proton (from beam)
    qecol = coef[1]*(xparm^coef[3])/(1+xparm^coef[4])*1.0D-15 ; the rate coefficient (m^3/s)
    lineemis = profiles.n_ion[fb,indion]*1D12*(beam.n[fb,*]#qecol)/(4D*!dpi) ;emissivity in 1/(m^3*s*sR)

    ;Doppler broadening 1/e half width is lambda0*sqrt(t_e_rbm/m_d)/c
    ;width2 = 2*extend(n_zeeman,profiles.t_i[fb]/1000)*wvezeeman^2/mc2kev
    width2 = 2*profiles.t_i[fb]/1000*wved^2/mc2kev
    ds1 = chord.ds / double(chord.nray)

    ;Doppler shift
    ;calculate the ion velocity based on v_pol and v_tor and the magnetic field direction
    v_pol = interpol(Er_vtheta,Er_r,chord.rmid[fb])*100 ;convert to cm/s
    v_tor = interpol(Er_vphi,Er_r,chord.rmid[fb])*100   ;convert to cm/s
    sintheta = field.br[fb]/field.b_tot[fb]
    costheta = sqrt(1-sintheta^2)
    cosphi = cos(chord.phi[fb])
    sinphi = sin(chord.phi[fb])
    v_ion_x = -sinphi*v_tor-sintheta*cosphi*v_pol
    v_ion_y = cosphi*v_tor-sintheta*sinphi*v_pol
    v_ion_z = costheta*v_pol
    dopplerv = v_ion_x*chord.vect[fb,0] + v_ion_y*chord.vect[fb,1] + v_ion_z*chord.vect[fb,2]
    dopplerf = sqrt((1 + dopplerv/c)/(1 - dopplerv/c))
    wved = cmreplicate(wved*dopplerf,n_comp)

    ;Calculate zeeman spectrum
    if keyword_set(calczeeman) then begin
	case calczeeman of
	    1: begin
	    print,'Zeeman calculation'
	    ;normalize_zeeman,species,wvezeeman,ampzeeman,field.b_tot[fb],angle.psi[fb]
	    normalize_zeeman_stokes,species,wvezeeman,s0,s1,s2,s3,field.b_tot[fb],angle.vectH[fb,*],angle.vectV[fb,*]
	    n_zeeman = (size(wvezeeman,/dimensions))[0]
	    wvezeeman *= extend(n_zeeman,dopplerf)
	    end
	    'Blom': begin
	    print,'Blom Zeeman parameterization'
	    blom,species,wvezeeman,s0,s1,s2,s3,field.b_tot[fb],angle.vectH[fb,*],angle.vectV[fb,*],profiles.t_i[fb],widthfactor
	    width2 *= sqrt(1 + widthfactor^2)
	    n_zeeman = 3
	    wvezeeman *= extend(n_zeeman,dopplerf)
	    end
	endcase
    endif

    if keyword_set(cs_effects) then begin ;based on Von Hellerman PPCF 37 (1995) 71-94
	;Note, the second order Taylor expansion around vbeam-v is poor for the E/18 component.
	vbeam = beam.beta*c ;cm/s
	xp2 = sqrt(xparm) ; v_col/v_m
	Q0 = coef[1] ;3.54D
	p = coef[3]*2 ;1.39D*2
	q = coef[4]*2 ;3.75D*2
	deno = 1+xp2^q
	dqdv = Q0*(p*xp2^(p-1)/deno-q*xp2^(p+q-1)/deno^2)/vbeam
	d2qdv2 = Q0*(p*(p-1)*xp2^(p-2)/deno-q*(2*p+q-1)*xp2^(p+q-2)/deno^2+2*q^2*xp2^(p+2*q-2)/deno^3)/vbeam^2
	alpha1 = dqdv/qecol*1.0D-15
	beta = 1d/2d*(d2qdv2/qecol*1.0D-15 - alpha1^2)
	eps = alpha1/2/vbeam
	delta = !dpi/2-angle.alpha[fb] ;nchord
	Dparm = (sin(delta)^2)#eps+(cos(delta)^2)#beta
	vt2 = 2*profiles.t_e[fb]/1000/mc2kev*c^2 ;thermal vel^2, (cm/s)^2
	Gparm = 1 - vt2#(eps+beta)+(vt2^2)#(eps*beta) ; [n_fb, 4]
	Aparm = 1 - cmreplicate(vt2,n_comp)*Dparm  ;[n_fb, 4]
	Bparm = 1 - vt2#eps  ;[n_fb, 4]
	Cparm = (sin(delta)^2*cos(delta)^2*vt2)#(eps-beta)
	invflag = Bparm le 0.8
	invalid = where(invflag) ;Taylor expansion breaks down
	if invalid[0] ge 0 then begin
	    Dparm[invalid] = 0
	    Gparm[invalid] = 1
	    Bparm[invalid] = 1
	    Aparm[invalid] = 1
	endif
	width2 = cmreplicate(width2,n_comp)*(1-cmreplicate(vt2,n_comp)*Dparm)/Gparm ;from Hellerman eq 13   [n_fb, 4]
	vshift = -Bparm*((sin(delta)*vt2)#alpha1)/2/Gparm;-(Cparm*ry-Aparm*rz)/Gparm ;commented out plasma rotation sensitivity
	cswveshift = vshift*wved/c
	amp2 = extend(n_elements(fb),qecol)*exp(invflag*vt2#(alpha1^2)*(cmreplicate(sin(delta)^2,n_comp)*Bparm+cmreplicate(cos(delta)^2,n_comp)*Gparm)/(4*Aparm*Gparm))/sqrt(Gparm*Bparm)
	lineemis = cmreplicate(profiles.n_e[fb],n_comp)*1D12*beam.n[fb,*]*amp2/(4D*!dpi) ;emissivity in 1/(m^3*s*sR) [n_fb]
    endif
    width2 += l_instr^2
    width = sqrt(width2)

    case linemodel of
	'gaussian': begin
	    message,'fixme'
	    npixel = n_elements(wve)
	    A = lineemis/sqrt(!dpi*width2)*ds1
	    spectra = dblarr(npixel,n_comp)
	    for i=0,npixel-1 do begin
		spectra[i,*] = total(A*exp(-(wve[i]-wved)^2/width2),1)
	    endfor
	end
	'erf': begin
	    npixel = n_elements(wve)
	    pixboundl = [wve[0],(wve[0:npixel-1]+wve[1:*])/2]
	    pixboundr = [(wve[0:npixel-1]+wve[1:*])/2,wve[npixel-1]]
	    ;spectra = dblarr(npixel,n_comp)
	    stokes_s0 = dblarr(npixel,n_comp)
	    stokes_s1 = dblarr(npixel,n_comp)
	    stokes_s2 = dblarr(npixel,n_comp)
	    stokes_s3 = dblarr(npixel,n_comp)
	    if keyword_set(calczeeman) then begin
;figure,8
;clf,/all
;pol = (s3[*,0] gt 0) + 2*(s3[*,0] lt 0)
		for j=0,n_zeeman-1 do begin
		    ;print,'.',format='(A,$)'
		    print,'j='+strtrim(j,2)+'/'+strtrim(n_zeeman,2)+' : '+strtrim(100.0*j/n_zeeman,2)+'%'
		    wvenew = rebin(reform(wvezeeman[j,*]),n_bpoints,n_comp) + cswveshift
		    ;ampnew = lineemis*rebin(reform(ampzeeman[j,*]),n_bpoints,n_comp)
		    s0new = lineemis*rebin(reform(s0[j,*]),n_bpoints,n_comp)
		    s1new = lineemis*rebin(reform(s1[j,*]),n_bpoints,n_comp)
		    s2new = lineemis*rebin(reform(s2[j,*]),n_bpoints,n_comp)
		    s3new = lineemis*rebin(reform(s3[j,*]),n_bpoints,n_comp)
		    for i=0,n_elements(wve)-1 do begin
			erfvalue = 0.5d * (erf((pixboundl[i] - wvenew)/width) - erf((pixboundr[i] - wvenew)/width))
			stokes_s0[i,*] += total(s0new * erfvalue,1)
			stokes_s1[i,*] += total(s1new * erfvalue,1)
			stokes_s2[i,*] += total(s2new * erfvalue,1)
			stokes_s3[i,*] += total(s3new * erfvalue,1)
		    endfor
;wve1 = mean(wvezeeman[j,*])
;amp1 = s0new[j,0]
;oploti,[wve1,wve1],[0,amp1],icolor=pol[j]
		endfor
	    endif else begin
		wvenew = wved + cswveshift
		for i=0,n_elements(wve)-1 do begin
		    stokes_s0[i,*] = total(lineemis * 0.5d * (erf((pixboundl[i] - wvenew)/width) - erf((pixboundr[i] - wvenew)/width)),1)
		endfor
	    endelse
	    ;spectra *= rebin(ds1/(pixboundl-pixboundr),npixel,n_comp)
	    stokes_s0 *= rebin(ds1/(pixboundl-pixboundr),npixel,n_comp)
	    stokes_s1 *= rebin(ds1/(pixboundl-pixboundr),npixel,n_comp)
	    stokes_s2 *= rebin(ds1/(pixboundl-pixboundr),npixel,n_comp)
	    stokes_s3 *= rebin(ds1/(pixboundl-pixboundr),npixel,n_comp)
	end
	'delta': begin
	    message,'fixme'
	end
    endcase

    cxrs_data = {stokes_s0:stokes_s0,stokes_s1:stokes_s1,stokes_s2:stokes_s2,stokes_s3:stokes_s3}
    if n_comp eq 1 then return,stokes_s0 else return,total(stokes_s0,2)
end
