;Calculate the H Balmer-alpha line shifts and amplitudes with Stark+Zeeman effects
;Note, CGS units used
;Use coordinate system of R.C. Isler, Phys. Rev. A, vol 14.3, 1976, p1015
;   z is the magnetic field direction
;   view is in the xz plane
;   y is perpendicular to x and z
;Input:
;   B_tot: magnitude of B field, in Gauss
;   E_tot: magnitude of E field, in statvolt/cm = Gauss
;   gamma: relativistic Lorentz factor of the beam
;   psi:   the angle between the magnetic field and the viewing chord (no longer used! vectH and vectV have all the information!)
;   phi:   the angle between projection of viewing chord on to plane perpendicular to B,
;          and projection of beam on to plane perpendicular to B.
;          In other words, the angle between the view and beam in the xy plane
;   theta: the angle between the magnetic field and beam
;   jones: jones matrix to apply on result, to represent polarization effects of optics
;   H_alpha: wavelength of the unperturbed line. Made into a parameter to support D beams
;   vectH: vector perpendicular to the viewing direction in the horizontal plane, written in local Isler coordinate system
;   vectV: vector perpendicular to the viewing direction and vectH
;Output:
;   shifts:     array[15] of line shifts in angstroms
;   amplitude:  array[15] of line amplitudes in 1/s
;               Emission is amplitude * population of upper state
;               For statistical population, P(upper state) ~= P(n=3)/9
;Keywords:
;   stop:       stop the calculation for debugging
;   quadratic:  include quadratic stark shift as a perturbation to energy levels.
;               neglect any perturbation to the eigenvectors
;   onlystark:  set Zeeman splitting to zero to compare SZ with pure Stark effect
pro starkzeeman,B_tot,E_tot,phi,psi,shifts,amplitude,jones,stop=stop,quadratic=quadratic,onlystark=onlystark,stokes=stokes,H_alpha=H_alpha,vectH=vectH,vectV=vectV
    c = 29979245800d            ; cm/s     speed of light in vacuum
    ;H_alpha = 6562.7952d        ; A        wavelength of H-alpha line
    hbar = 1.054571726d-27      ; erg s
    ec = 4.8032068d-10          ; esu
    a0 = 0.52917721092d-8       ; cm       Bohr radius
    bohrmag = 9.27400968d-21    ; erg/G    Bohr magneton = e*hbar/(2*me*c)
    m_e = 9.10938291d-28        ; g        electron mass
    alpha = 7.2973525698d-3     ;          Fine structure constant

    gam = bohrmag*B_tot         ;gam = Bohr magneton times B_t (erg)
    eps = 3*ec*a0*E_tot         ;eps = Stark energy splitting coefficient (erg) = 2*E_L
    if keyword_set(onlystark) then gam = 0d

;testing
;vectv=transpose([-0.221525,-0.974146,0.0443431])
;vecth=transpose([-0.955197,0.225919,0.191204])  
;gam = 0.00018605221d*1.6022d-12
;eps=0.00104028d*1.6022d-12   
;phi=2.9825398
;stop=1

    tau=dcomplex(0,1)*exp(-dcomplex(0,1)*phi)
    tauc=conj(tau)

    q0=sqrt(gam^2+eps^2)        ;n=2 energy splitting parameter (erg)
    q1=sqrt(4*gam^2+9*eps^2)    ;n=3 energy splitting parameter (erg)

    ;electric dipole transition integral <psi_nlm|d|psi_n'l'm'> can be split into angular and radial parts
    ;<Y_lm|T|Y_l'm'><R_nl|r|R_n'l'>, where T is a spherical tensor
    ;Rnln'l' = <R_nl|r|R_n'l'> given by integral[R_nl[r]*R_n'l'[r]*r^3]
    R3221 = 165888D/15625D/sqrt(5D) ;results in units of Bohr radius
    R3120 = 27648D/15625D*sqrt(3D)
    R3021 = 10368D/15625D*sqrt(2D)

    ;numbering of states in H0 basis, nlm, n=3
    ; 1. 300
    ; 2. 311
    ; 3. 31-1
    ; 4. 320
    ; 5. 322
    ; 6. 32-2
    ; 7. 310
    ; 8. 321
    ; 9. 32-1
    ;numbering of states in H0 basis, nlm, n=2
    ; 1. 200
    ; 2. 211
    ; 3. 21-1
    ; 4. 210

    ;numbering of states in H' basis, n=3 (see Isler) (approximate matching nk|m| state, degenerate to sign of m)
    ; 1. q1      (6 mix)                                 3,+2,0
    ; 2. -q1     (6 mix)                                 3,-2,0
    ; 3. q1/2    (6 mix)                                 3,+1,+-1
    ; 4. -q1/2   (6 mix)                                 3,-1,+-1
    ; 5. 0       (6 mix)                                 3,0,blend  Isler eigenvectors don't seem to be eigenvectors of m^2
    ; 6. 0       (2 mix)                                 3,0,blend
    ; 7. 0       (3 mix)                                 3,0,blend
    ; 8. q1/2    (3 mix)                                 3,+1,+-1
    ; 9. -q1/2   (3 mix)                                 3,-1,+-1
    ;numbering of states in H' basis, n=2 (see Isler) (approximate matching nk|m| state, degenerate to sign of m)
    ; 1. 0       (3 mix)                                 2,0,+-1
    ; 2. q0                                              2,+1,0
    ; 3. -q0                                             2,-1,0
    ; 4. 0       (same as 210)                           2,0,+-1
    ;eigenvalues are given by [Isler, R.C., Phys. Rev. A, vol 14.3, 1976, p1015]

    ;transition matrix from H0 basis to H' basis, n=3, rows are H', columns H0
    t3 = [ $
	[3d*sqrt(3d)*eps^2, 3d/2d*tau*eps*(2*gam+q1), 3d/2d*tauc*eps*(2*gam-q1), -3d*sqrt(3d)*eps^2/2d/sqrt(2d), 1d/4d*tau^2*(8*gam^2+9*eps^2+4*gam*q1), 1d/4d*tauc^2*(8*gam^2+9*eps^2-4*gam*q1), 0, 0, 0]/q1^2, $ ;q1
	[3d*sqrt(3d)*eps^2, 3d/2d*tau*eps*(2*gam-q1), 3d/2d*tauc*eps*(2*gam+q1), -3d*sqrt(3d)*eps^2/2d/sqrt(2d), 1d/4d*tau^2*(8*gam^2+9*eps^2-4*gam*q1), 1d/4d*tauc^2*(8*gam^2+9*eps^2+4*gam*q1), 0, 0, 0]/q1^2, $ ;-q1
	[4d*sqrt(3d)*eps*gam, 1d/2d*tau*(4*gam^2-9*eps^2+2*gam*q1), 1d/2d*tauc*(4*gam^2-9*eps^2-2*gam*q1), -sqrt(6d)*eps*gam, -3d/2d*tau^2*eps*(2*gam+q1), -3d/2d*tauc^2*eps*(2*gam-q1), 0, 0, 0]/q1^2, $ ;q1/2
	[4d*sqrt(3d)*eps*gam, 1d/2d*tau*(4*gam^2-9*eps^2-2*gam*q1), 1d/2d*tauc*(4*gam^2-9*eps^2+2*gam*q1), -sqrt(6d)*eps*gam, -3d/2d*tau^2*eps*(2*gam-q1), -3d/2d*tauc^2*eps*(2*gam+q1), 0, 0, 0]/q1^2, $ ;-q1/2
	[1d/3d*sqrt(2d)*(8*gam^2-9*eps^2), -3d*sqrt(6d)*tau*gam*eps, -3d*sqrt(6d)*tauc*gam*eps, -1d/6d*(8*gam^2-9*eps^2), 9d*sqrt(3d)*tau^2*eps^2/2d/sqrt(2d), 9d*sqrt(3d)*tauc^2*eps^2/2d/sqrt(2d), 0, 0, 0]/q1^2, $ ;0
	[1d/3d*q1^2, 0, 0, 2d/3d*sqrt(2d)*q1^2, 0, 0, 0, 0, 0]/q1^2, $ ;0
	[0, 0, 0, 0, 0, 0, 2d*gam, -3d*tau*eps/sqrt(2d), -3d*tauc*eps/sqrt(2d)]/q1, $ ;0
	[0, 0, 0, 0, 0, 0, 3d*eps/sqrt(2d), 1d/2d*tau*(2*gam+q1), 1d/2d*tauc*(2*gam-q1)]/q1, $ ;q1/2
	[0, 0, 0, 0, 0, 0, 3d*eps/sqrt(2d), 1d/2d*tau*(2*gam-q1), 1d/2d*tauc*(2*gam+q1)]/q1] ;-q1/2
    t3 = reform(t3,9,9)

    t2 = [ $
	[gam, -tau*eps/sqrt(2d), -tauc*eps/sqrt(2d), 0], $ ;0
	[eps/sqrt(2d), 1d/2d*tau*(gam+q0), 1d/2d*tauc*(gam-q0), 0], $ ;q0
	[eps/sqrt(2d), 1d/2d*tau*(gam-q0), 1d/2d*tauc*(gam+q0), 0], $ ;-q0
	[0, 0, 0, q0]]/q0 ;0

    ;eigenenergies of the H' basis vectors
    ;e3 = [q1, -q1, q1/2, -q1/2, 0, 0, 0, q1/2, -q1/2]
    ;e2 = [0, q0, -q0, 0]

    ;transition amplitudes in H0 basis, (rows n=3, columns n=2)
    ;c.f. Howard Yuh's thesis
    ;row states: 200, 211, 21-1, 210
    dz = [ $ ;z polarized
	[0, 0, 0, sqrt(1d/3d)*R3021], $  ;300
	[0, 0, 0, 0], $                  ;311
	[0, 0, 0, 0], $                  ;31-1
	[0, 0, 0, sqrt(4d/15d)*R3221], $ ;320
	[0, 0, 0, 0], $                  ;322
	[0, 0, 0, 0], $                  ;32-2
	[sqrt(1d/3d)*R3120, 0, 0, 0], $  ;310
	[0, sqrt(1d/5d)*R3221, 0, 0], $  ;321
	[0, 0, sqrt(1d/5d)*R3221, 0]]    ;32-1
    dp = [ $ ;x+iy circularly polarized
	[0, sqrt(2d/3d)*R3021, 0, 0], $  ;300
	[0, 0, 0, 0], $                  ;311
	[-sqrt(2d/3d)*R3120, 0, 0, 0], $ ;31-1
	[0, -sqrt(2d/15d)*R3221, 0, 0], $ ;320
	[0, 0, 0, 0], $                  ;322
	[0, 0, -sqrt(4d/5d)*R3221, 0], $ ;32-2
	[0, 0, 0, 0], $                  ;310
	[0, 0, 0, 0], $                  ;321
	[0, 0, 0, -sqrt(2d/5d)*R3221]]   ;32-1
    dm = [ $ ;x-iy circularly polarized
	[0, 0, -sqrt(2d/3d)*R3021, 0], $ ;300
	[sqrt(2d/3d)*R3120, 0, 0, 0], $  ;311
	[0, 0, 0, 0], $                  ;31-1
	[0, 0, sqrt(2d/15d)*R3221, 0], $ ;320
	[0, sqrt(4d/5d)*R3221, 0, 0], $  ;322
	[0, 0, 0, 0], $                  ;32-2
	[0, 0, 0, 0], $                  ;310
	[0, 0, 0, sqrt(2d/5d)*R3221], $  ;321
	[0, 0, 0, 0]]                    ;32-1

    ;calculate transition amplitudes in H' basis
    ;result has (rows n=3 H' basis), (columns n=2 H' basis)
    dz2 = (t3##dz)##transpose(conj(t2))
    dp2 = (t3##dp)##transpose(conj(t2))
    dm2 = (t3##dm)##transpose(conj(t2))

    ;isotropy check
    ;total(dz2^2) == total(((dp2+dm2)/2)^2) == total(((dp2-dm2)/2)^2)

    ;calculate amplitude along two viewing polarizations (Jones vectors)
    ;multiply by the Jones matrix of the polarizer
    ;then take abs^2 to get intensities
    dx2 = (dp2 + dm2)/2
    dy2 = (dp2 - dm2)/2 * dcomplex(0,-1)
    asx = vectH[0]*dx2 + vectH[1]*dy2 + vectH[2]*dz2 ;amplitude of dipole aligned with the H direction
    asy = vectV[0]*dx2 + vectV[1]*dy2 + vectV[2]*dz2 ;amplitude of dipole aligned with the V direction
;    asx = sin(psi)*dz2 + cos(psi)*(dp2 + dm2)/2
;    asy = (dp2 - dm2)/2
    if n_elements(jones) gt 0 then begin
	a2 = abs(jones[0,0]*asx + jones[1,0]*asy)^2 + abs(jones[0,1]*asx + jones[1,1]*asy)^2
    endif else begin
	a2 = abs(asx)^2 + abs(asy)^2
    endelse
    if arg_present(stokes) then begin
	s0 = abs(asx)^2 + abs(asy)^2
	s1 = abs(asx)^2 - abs(asy)^2
	s2 = 2*real_part(asx*conj(asy)) ;same as (abs(asx + asy)^2 - abs(asx - asy)^2)/2
	;s3 = 2*imaginary(asx*conj(asy)) ;same as(abs(asx - dcomplex(0,1)*asy)^2 - abs(asx + dcomplex(0,1)*asy)^2)/2
	s3 = -2*imaginary(asx*conj(asy)) ;Is sign flipped?
    endif

    ;combine degenerate rows and columns
    ;numbering of degenerate states: n=3
    ; 1. -q1
    ; 2. -q1/2
    ; 3. 0
    ; 4. +q1/2
    ; 5. +q1
    ;numbering of degenerate states: n=2
    ; 1. -q0
    ; 2. 0
    ; 3. +q0

    ;matrices for combining degenerate states (could save a little bit of CPU by combining earlier but meh)
    ;n=3 (columns = H', rows = degenerate)
    dg3 = [ $
	[0,1,0,0,0,0,0,0,0], $    ; -q1     3,-2,0
	[0,0,0,1,0,0,0,0,1], $    ; -q1/2   3,-1,1
	[0,0,0,0,1,1,1,0,0], $    ; 0       3, 0,blend{0,2}
	[0,0,1,0,0,0,0,1,0], $    ; +q1/2   3,+1,1
	[1,0,0,0,0,0,0,0,0]]      ; +q1     3,+2,0

    ;n=2 (columns = degenerate, rows = H')
    dg2 = transpose([ $
	[0,0,1,0], $                  ; -q0     2,-1,1
	[1,0,0,1], $                  ;   0     2, 0,1
	[0,1,0,0]])                   ; +q0     2,+1,1
    amplitude = dg3##(a2##dg2)

    ;scale the dipole matrix elements into transition probabilities
    omega = 2*!dpi*c/(H_alpha*1d-8) ; 1/s
    spont = 4*alpha*omega^3/3/c^2*a0^2 ; 1/s/cm^2
    amplitude *= spont
    if arg_present(stokes) then begin
	s0 = dg3##(s0##dg2) * spont
	s1 = dg3##(s1##dg2) * spont
	s2 = dg3##(s2##dg2) * spont
	s3 = dg3##(s3##dg2) * spont
	stokes = [[[s0]],[[s1]],[[s2]],[[s3]]]
    endif

    ;eigenenergies of the degenerate vectors
    if keyword_set(quadratic) then begin
	;in nkm parabolic basis, the quadratic Stark shift is -1/16*E_L^2*n^4*(17n^2-3k^2-9m^2+19)
	;The Isler eigenvectors are _almost_ eigenstates of m^2, with the exception of the 320 and 322 states
	;Quadratic Stark breaks the degeneracy of the 320, 322 states, but they are close together so
	;the states are blended together here and a weighted average line center is used
	;ratio is calculated in  stark_limit_test.pro
	qc = -(eps/3)^2/(m_e*c^2*alpha^2)/16d ; = -1/16*E_L^2 'quadratic Stark constant'
	k3 = [-2, -1, 0, 1, 2]
	ratio322_320 = 24d/103d*(16d +5d*sqrt(2d)) ; ratio of line intensities of the sigma0+ transition to the sigma0- transition
	m3sq = [0, 1, 0, 1, 0]
	e3 = [-q1,-q1/2,0,q1/2,q1] + qc*3^4*(17*3^2-3*k3^2-9*m3sq+19)
	e3 = [1d,1d,1d]#e3
	e3[1,2] -= qc*3^4*9*4*ratio322_320 ;only apply the m^2=4 shift to the sigma_0 line. pi_2 lines are unaffected
	k2 = [-1, 0, 1]
	m2sq = [0, 1, 0]
	e2 = [-q0,0,q0] + qc*2^4*(17*2^2-3*k2^2-9*m2sq+19)
	e2 = e2#(make_array(n_elements(e3),value=1d))
    endif else begin
	e3 = [1d,1d,1d]#[-q1,-q1/2,0,q1/2,q1]
	e2 = [-q0,0,q0]#[1d,1d,1d,1d,1d]
    endelse
    shifts = -(e3-e2)*h_alpha^2/(hbar*2*!dpi*c)/1d8 ; Angstroms
;    test = dblarr(9,4)
;    test[0,*] = stokes[2,4,*]
;    test[1,*] = stokes[1,3,*]
;    test[2,*] = stokes[0,2,*]
;    test[3,*] = stokes[2,3,*]
;    test[4,*] = stokes[1,2,*]
;    test[5,*] = stokes[0,1,*]
;    test[6,*] = stokes[2,2,*]
;    test[7,*] = stokes[1,1,*]
;    test[8,*] = stokes[0,0,*]
;    figure,1
;    ploti,test
;    poltest = sqrt(total(test[*,1:3]^2,2))/test[*,0]
;    figure,2
;    ploti,poltest
    if keyword_set(stop) then stop
end

pro test_stark_zeeman,stop=stop
;reproduce the plots in Breton 1980
    physicalconstants
    !p.background='ffffff'x
    !p.color='000000'x

    Bt = 5d4
    vr = [1.5d7, 6d7, 2.4d8, 1d7, 3d7, 2.4d8]
    phi = [0,0,0,!dpi/2,!dpi/2,!dpi/2]
    psi = !dpi/2
    c = 29979245800d
    EL = vr*Bt/c

    vectH = [1,0,0]
    vectV = [0,0,1]

    for j=0,5 do begin
	window,j,xsize=500,ysize=350
	starkzeeman,Bt,EL[j],phi[j],psi,shifts,amplitudes,stop=stop,h_alpha=6562.7952d,vectH=vectH,vectV=vectV
	starkzeeman,Bt,EL[j],phi[j],psi,shifts,amplitudes1,[[1,0],[0,0]],h_alpha=6562.7952d,vectH=vectH,vectV=vectV
	starkzeeman,Bt,EL[j],phi[j],psi,shifts,amplitudes2,[[0,0],[0,1]],h_alpha=6562.7952d,vectH=vectH,vectV=vectV
	shifts = (shifts - vr[j]*cos(phi[j])/!const.c.cm_s*!const.Halpha)/10d
	amplitudes /= 9 * 1d6
	amplitudes1 /= -9 * 1d6
	amplitudes2 /= 9 * 1d6
	plot,[0],[0],/nodata,xrange=[1.05*min(shifts)-0.05*max(shifts),1.05*max(shifts)-0.05*min(shifts)],yrange=[min(amplitudes1),max(amplitudes2)],charsize=1.8,ytitle='A (10!U6!N s!U-1!N sr!U-1!N)',xtitle='shift (nm)'
	for i=0,14 do begin
	    oplot,[shifts[i],shifts[i]],[0,amplitudes1[i]],color='ff0000'x
	    oplot,[shifts[i],shifts[i]],[0,amplitudes2[i]],color='0000ff'x
	endfor
	oplot,[-100,100],[0,0],color='000000'x
    endfor
end

pro test_stokes
;    Bt = 5d4
;    vr = 3d8
;    phi = 0.4
;    psi = 1.5
;    c = 29979245800d
;    EL = vr*Bt/c
    Bt = 3.21414d4
    EL = 218.52389657d
    phi = 1.4117427d
    psi = 0.1d

    starkzeeman,Bt,EL,phi,psi,shifts,amplitudes,stokes=stokes,h_alpha=6562.7952d,/stop
    I = stokes[*,*,0]
    p = sqrt(total(stokes[*,*,1:3]^2,3))/stokes[*,*,0]
    spsi = atan(stokes[*,*,2],stokes[*,*,1])/2
    schi = atan(stokes[*,*,3],sqrt(stokes[*,*,1]^2+stokes[*,*,2]^2))/2
    order = [7,5,6,4,2,3,1,0] ;sigma0, sigma1, pi2, pi3, pi4, sigma5, sigma6, pi8

    stokes = reform(stokes,15,4)
    figure,1
    clf,/all
    for i=0,14 do begin
	oploti,[shifts[i],shifts[i]],[0,stokes[i,0]]
;	oploti,[shifts[i],shifts[i]],[0,stokes[i,1]],icolor=1
;	oploti,[shifts[i],shifts[i]],[0,stokes[i,2]],icolor=2
;	oploti,[shifts[i],shifts[i]],[0,stokes[i,3]],icolor=3
    endfor

    figure,2
    clf,/all
    for i=0,14 do begin
	oploti,[shifts[i],shifts[i]],[0,p[i]]
    endfor
;
;    figure,3
;    clf,/all
;    for i=0,14 do begin
;	oploti,[shifts[i],shifts[i]],[0,spsi[i]]
;    endfor
;
;    figure,4
;    clf,/all
;    for i=0,14 do begin
;	oploti,[shifts[i],shifts[i]],[0,schi[i]]
;    endfor
end

;test quadratic stark feature
pro test_quadratic
    Bt = 54000d
    m_p = 1.672721637e-24       ; g        proton mass
    v_beam = sqrt(2*50000d*1.60217646d-12/m_p)
    starkzeeman,54000,v_beam,1.2,0.93,1.4,sh1,amp1,/onlystark
    starkzeeman,54000,v_beam,1.2,0.93,1.4,sh2,amp2,/onlystark,/quadratic
    starkzeeman,54000,v_beam,1.2,0.93,1.4,sh3,amp3
    starkzeeman,54000,v_beam,1.2,0.93,1.4,sh4,amp4,/quadratic
    figure,1
    clf,/all
    ploti,/nodata,xtitle='\Delta\lambda (A)',ytitle='transition probability (10^8/s)',charsize=2
    for i=0,n_elements(sh1)-1 do begin
	oploti,[sh1[i],sh1[i]],[0,amp1[i]/1d8],icolor=0
	oploti,[sh3[i],sh3[i]],[0,amp3[i]/1d8],icolor=2
    endfor
    for i=0,n_elements(sh2)-1 do begin
	oploti,[sh2[i],sh2[i]],[0,amp2[i]/1d8],icolor=1
	oploti,[sh4[i],sh4[i]],[0,amp4[i]/1d8],icolor=3
    endfor
    legendi,label=['Linear Stark','Quadratic Stark','Linear Stark+Zeeman','Quadratic Stark+Zeeman'],icolor=[0,1,2,3]
end

pro statistical_amps,bt,eb
    defvalue,Bt,2.2d4
    m_p = 1.672721637e-24       ; g        proton mass
    defvalue,eb,30000d
    v_beam = sqrt(2*eb*1.60217646d-12/m_p)
    phi = !dpi/2
    psi = !dpi/2
    c = 29979245800d
    EL = v_beam*Bt/c

    vectH = [1,0,0]
    vectV = [0,0,1]

    starkzeeman,Bt,EL,phi,psi,shifts,amplitudes,/onlystark,h_alpha=6562.7952d,vectH=vectH,vectV=vectV
    stat_amps = flatten(amplitudes/total(amplitudes)*2)
    print,stat_amps
end
