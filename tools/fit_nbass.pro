;@/home/liao/bes_spectrum/bes_spectrum_stark_zeeman_miller.pro
@outer_sum.pro

;fit the stark zeeman pattern
;A[0] wavelength shift (due to miscalibration or misunderstanding of h-alpha center)
;A[1] |B|
;A[2] amplitude scaling on total E/1 component
;A[2+n] amplitude scaling on total E/(n+1) component
;A[3+n] width of line
;A[4+n] amplitude of CXRS gaussian
;A[5+n] center of CXRS gaussian
;A[6+n] width of CXRS gaussian
function fit_ls_erf,xfit,a,ampl=ampl,vsintheta=vsintheta,voverc=voverc,cosalpha=cosalpha,sens=sens,pixbound=pixbound,disp=disp
    c = 29979245800d            ; cm/s     speed of light in vacuum
    H_alpha = 6562.819d         ; A (air)
    hbar = 1.054571726d-27      ; erg s
    ec = 4.8032068d-10          ; esu
    a0 = 0.52917721092d-8       ; cm       Bohr radius
    bohrmag = 9.27400968d-21    ; erg/G    Bohr magneton = e*hbar/(2*me*c)
    m_e = 9.10938291d-28        ; g        electron mass
    alpha = 7.2973525698d-3     ;          Fine structure constant

    B_tot = A[1]
    gam = bohrmag*B_tot         ;gamma = Bohr magneton times B_t (erg)

    k3 = [-2, -1, 0, 1, 2]
    m3sq = [0, 1, 0, 1, 0]
    k2 = [-1, 0, 1]
    m2sq = [0, 1, 0]

    yout = dblarr(n_elements(xfit))
    ncomp = n_elements(a)-7 ;1 less than the number of beam components
    for bcomp=0,ncomp do begin
	eps = 3*ec*a0*vsintheta[bcomp]/c*B_tot ;epsilon = Stark energy splitting coefficient (erg) = 2*E_L
	qc = -(ec*a0*vsintheta[bcomp]/c*B_tot)^2/(m_e*c^2*alpha^2)/16d ; = -1/16*E_L^2 'quadratic Stark constant'

	q0=sqrt(gam^2+eps^2)        ;n=2 energy shift parameter (erg)
	q1=sqrt(4*gam^2+9*eps^2)    ;n=3 energy shift parameter (erg)

	e3 = [-q1,-q1/2,0,q1/2,q1] + qc*3^4*(17*3^2-3*k3^2-9*m3sq+19)

	;e3_adjust = -qc*3^4*9*mtwo/(mone+mzero)*4
	e3_adjust = -qc*3^4*9*0.58*4 ;from bes_spectrum_stark_zeeman simulation

	e2 = [-q0,0,q0] + qc*2^4*(17*2^2-3*k2^2-9*m2sq+19)

	;shifts = -outer_sum(e3,-e2)*h_alpha^2/(hbar*2*!dpi*c)/1d8 ; Angstroms
	shifts = -outer_sum(-e2,e3)
	shifts[1,2] -= e3_adjust
	shifts *= h_alpha^2/(hbar*2*!dpi*c)/1d8 ; Angstroms

	lines = (h_alpha + shifts + a[0])/(1+voverc[bcomp]*cosalpha)/(sqrt(1-voverc[bcomp]^2)) ;apply Doppler shift
	for linei=0,14 do begin
	    erf0 = erf((pixbound-lines[linei])/(a[3+ncomp]*sqrt(2)))
	    ;oneline = sqrt(!dpi/2)*ampl[linei,bcomp]*a[2+bcomp] * ([1d,erf0] - [erf0,-1d])
	    oneline = sqrt(!dpi/2)*ampl[linei,bcomp]*a[2+bcomp] * (erf0[0:n_elements(erf0)-2] - erf0[1:n_elements(erf0)-1])
	    yout += oneline
	endfor
	;try a more vectorized version. not really faster. actually a tiny bit slower. Maybe erf is internally not vectorized
;	erf0 = erf(outer_sum(pixbound,-reform(lines,15))/(a[3]*sqrt(2)))
;        oneline = extend(n_elements(xfit),sqrt(!dpi/2)*ampl[*,bcomp]*a[2+bcomp]) * ([ones(1,15),erf0] - [erf0,-ones(1,15)])
;	yout += total(oneline,2)
    endfor
    yout *= sens / (-disp)
    yout += A[4+ncomp]*exp(-(xfit-A[5+ncomp])^2/A[6+ncomp]) ;CXRS can be done better.. just a quick fix
    return,yout
end

pro MSE_LS_synth,run=run,synthfile=synthfile,fixwidth=fixwidth,fixcenter=fixcenter,fitgamma=fitgamma,domainstart=domainstart,domainend=domainend,width0=width0,polfilter=polfilter,plot=plot,addimpurity=addimpurity,fitimpurity=fitimpurity
    @../beam_param.idl
    @../mesh_param.idl

    filename='../runs/'+run+'/'+synthfile
    defvalue,frame_scl,1d0
    restore,filename
    sens=result.sens * parameters.detector_param.int_time
    ;scalefactor=result.scalefactor ;conversion from ph/m2/s/sr/A to counts. units:counts/(ph/m2/s/sr/A)
    size0 = (size(result.noisy,/dim))
    if n_elements(size0) lt 2 then n_gen=1 else n_gen=(size(result.noisy,/dim))[1]

    openw,lun,filename+'.dat',/get_lun,width=1000
    printf,lun,'generation,shift[A],|B|[G],A_E0,linewidth,CXRSamp,CXRScenter,CXRSwidth,V_beam,alpha'
    resultB = dblarr(n_gen)   ; measured magnetic field vs generation
    result2a = dblarr(n_gen) ; measured magnetic field vs generation (fit_gamma)
    result2b = dblarr(n_gen) ; measured cos(gamma) vs generation (fit_gamma)
    ;radii = dblarr(n_elements(chan_sel))

    x = result.wve

    alcbeam_file = parameters.beam_param.alcbeam_file
    get_beam_data,0,data.chord,beam

    n_beam = beam.n
    n_beam_n1 = beam.n1
    n_beam_n2 = beam.n2
    n_beam_n3 = beam.n3
    E_beam = beam.E_amu
    v_beam = beam.beta * 299792458l ;m/s
    print,'Beam energy:',E_beam
    print,'Beam velocity:',v_beam
    ;approx width of beam (m)
    ds = parameters.mesh_param.ds * total(n_beam[*,0])/max(n_beam[*,0],fb)
    ds /= data.chord.nray
    rbm = data.chord.r[fb]
    theta = data.angle.theta[fb]
    alpha = data.angle.alpha[fb]
    b_tot=data.field.b_tot

    shifts = dblarr(15,n_elements(v_beam))    ;A
    amplitude = dblarr(15,n_elements(v_beam)) ;1/s [15,4,num_beam]
    marchuk_ratios = dblarr(15,n_elements(v_beam))
    for i=0,n_elements(v_beam)-1 do begin
	starkzeeman,data.field.b_tot[fb]*1d4,v_beam*100,data.angle.phi[fb],data.angle.psi[fb],data.angle.theta[fb],shifts_one,amplitude_one,jones,/quadratic
	shifts[*,i] = reform(shifts_one,15)
	amplitude[*,i] = reform(amplitude_one,15) / 9D ; rescale to accept total population of n=3 level instead of population of sublevel

;	n_e = data.profiles.n_e[fb]
;	marchuk_n_e = n_e
;	marchuk_e_b = e_beam[i]/1000
;	marchuk_b_t = data.field.b_tot[fb]
;	marchuk_z = data.profiles.z_eff[fb]
;	marchuk_ratio = get_marchuk_ratios(marchuk_n_e,marchuk_e_b,marchuk_B_t,marchuk_z)
;	marchuk_ratios[*,i] = reform(marchuk_ratio,15,/overwrite)
    endfor

    ampl = amplitude; * marchuk_ratios ;in m^3 / s

    sigmas = [0,1,0,1,0,1,0,1,0,1,0,1,0,1,0]
    pis =    [1,0,1,0,1,0,1,0,1,0,1,0,1,0,1]
    if keyword_set(polfilter) then begin
	if polfilter eq 'cut_pi' then ampl *= rebin(sigmas,15,n_elements(v_beam))
	if polfilter eq 'cut_sigma' then ampl *= rebin(pis,15,n_elements(v_beam))
	if polfilter eq 'PEMs' then ampl *= rebin(0.6*sigmas + 1.3*pis,15,n_elements(v_beam)) ;roughly simulate the PEMs. Very coarse model
    endif
;    condense = [[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],$
;                [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],$
;                [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],$
;                [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],$
;                [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],$
;                [0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],$
;                [0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],$
;                [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],$
;                [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],$
;                [0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],$
;                [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],$
;                [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],$
;                [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],$
;                [0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],$
;                [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],$
;                [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],$
;                [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],$
;                [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]]
;    ampl = amplitude ## condense
    ;ampl *= rebin(n_beam_n3[fb,*],15,4) ;will be fitted

    ;estimate line width
    c = 29979245800d            ; cm/s     speed of light in vacuum
    H_alpha = 6562.819d         ; A (air)
    l_instr = parameters.detector_param.l_instr
    dline_wve_dalpha = H_alpha*outer(beam.beta[0]/c,sin(data.angle.alpha[fb]))  ;[4] A
    dopp_broad_fg = data.dangle_fg.alpha[fb]*dline_wve_dalpha ;[4]
    dopp_broad_ap = data.dangle_ap.alpha[fb]*dline_wve_dalpha
    dopp_broad_sp = data.dangle_sp.alpha[fb]*dline_wve_dalpha
    ripple_broad = [0.24498053]
    rad_w_ratio = 0.750846d ;empirical fit
    ;rad_w_ratio = 2d ;hack
    ;defvalue,width0, sqrt(l_instr^2 + (dopp_broad_ap*rad_w_ratio)^2 + (dopp_broad_fg*rad_w_ratio)^2 + $
    ;                 (dopp_broad_sp*rad_w_ratio)^2 + ripple_broad^2)
    if run eq 'CMOD_SMSE' then defvalue,width0,1.6d else defvalue,width0,1.0d
    ;defvalue,domainstart,H_alpha*(1+v_beam[0]/c*cos(data.angle.alpha[fb]))-45
    ;defvalue,domainend,H_alpha*(1+v_beam[0]/c*cos(data.angle.alpha[fb]))+45
    defvalue,domainstart,H_alpha*(1-beam.beta[0]*cos(data.angle.alpha[fb]))-45
    defvalue,domainend,H_alpha*(1-beam.beta[0]*cos(data.angle.alpha[fb]))+45

    wveimpurity = maken(6560,6610,n_gen+1)
    for gen=-1,n_gen-1 do begin
	if gen eq -1 then begin
	    yadd = result.pure
	    ysub = result.pure_bk
	endif else begin
	    yadd = result.noisy[*,gen]
	    ysub = result.noisy_bk[*,gen]
	endelse
	yadderr = sqrt(yadd)+1
	ysuberr = sqrt(ysub)

	y = yadd - ysub/frame_scl
	yerr2 = yadderr^2 + (ysuberr/frame_scl)^2

	if keyword_set(addimpurity) then begin ;quick and dirty hack to test impurity sensitivity
	    if isnum(addimpurity) then y += gen*exp(-(x-addimpurity)^2 / 2.1) else y += 500*exp(-(x-wveimpurity[gen+1])^2 / 2.1)
	endif

	nx = n_elements(x)-1
	scaling = 1d
	domain = where(x gt domainstart and x lt domainend) ;gets a bit of D alpha CXRS, perhaps introducing some error, but it still reduces spread.
	yerrfix = sqrt(yerr2[domain])
	disp = diffc(x[domain])
	pixbound = (x[0:nx-1] + x[1:nx])/2
	pixbound = [2*pixbound[0]-pixbound[1], pixbound, 2*pixbound[nx-1]-pixbound[nx-2]]
	pixbound = pixbound[domain[0]:domain[n_elements(domain)-1]+1]
;		loop1:
;		read,scaling,prompt='scaling['+num2str(scaling)+']:'
	randomize = randomn(seed,1)*300
	case run of
	    'ITER': a = [0, data.field.b_tot[fb]*10000, scaling*ds*n_beam_n3[fb], width0, 3200*scaling, 6562.8d, 1000d] ;note cgs units in fit_ls, need to fix width
	    'EAST': a = [0, data.field.b_tot[fb]*10000 + randomize, ds*n_beam_n3[fb,0], ds*n_beam_n3[fb,1], ds*n_beam_n3[fb,2], width0, 23000*scaling, 6562.8d, 80d] ;note cgs units in fit_ls, need to fix width
	    'CMOD_SMSE': a = [0, data.field.b_tot[fb]*10000, ds*n_beam_n3[fb,0], ds*n_beam_n3[fb,1], ds*n_beam_n3[fb,2], ds*n_beam_n3[fb,3], width0, 10000*scaling, 6562.8d, 80d]
	endcase
	if keyword_set(plot) then ploti,x[domain],y[domain]
	res = fit_ls_erf(x[domain],a,ampl=ampl,vsintheta=v_beam*sin(theta)*100,voverc=v_beam/299792458d,cosalpha=cos(alpha),sens=sens[domain],pixbound=pixbound,disp=disp)
;		oploti,x[domain],res,icolor=1
;		stop
;		goto,loop1
	ncomp = n_elements(v_beam)-1
	parinfo = replicate({fixed:0, limited:[1,1], limits:[0.D,0.D], step:0d},7+ncomp)
	parinfo[0].limits = [-4d,4d] ;halpha shift
	;parinfo[1].limits = [data.field.b_tot[fb]*9500,data.field.b_tot[fb]*11000] ;b
	parinfo[1].limits = [data.field.b_tot[fb]*8500,data.field.b_tot[fb]*11000] ; needed expanded error range for impurity sensitivity test
	for i=0,ncomp do begin
	    parinfo[2+i].limits = [a[2+i]/2,a[2+i]*2] ;n_beam
	endfor
	parinfo[3+ncomp].limits = [width0/2,width0*2] ;width
	parinfo[4+ncomp].limits = [0d,5d4] ;amplitude of cxrs
	parinfo[5+ncomp].limits = [6560d, 6565d] ;center of cxrs
	parinfo[6+ncomp].limits = [50d, 1600d] ;width of cxrs
	parinfo[0].step = 0.1
	parinfo[1].step = 500
	for i=0,ncomp do begin
	    parinfo[2+i].step = 0.2*a[2+i]
	endfor
	parinfo[3+ncomp].step = 0.03
	parinfo[4+ncomp].step = 1d2
	parinfo[5+ncomp].step = 0.05
	parinfo[6+ncomp].step = 100d
	if keyword_set(fixwidth) then parinfo[6].fixed=1
	if keyword_set(fixcenter) then parinfo[0].fixed=1
	b = mpfitfun('fit_ls_erf',x[domain],y[domain],yerrfix,a,MAXITER=300,YFIT=yout,perror=perror,FUNCTARGS={ampl:ampl,vsintheta:v_beam*sin(theta)*100,voverc:v_beam/299792458d,cosalpha:cos(alpha),sens:sens[domain],pixbound:pixbound,disp:disp},parinfo=parinfo,NPRINT=0,/quiet)
	;printf,lun,gen,b[0],b[1],b[2],b[3],b[4],b[5],b[6],e_beam[0],alpha
	printf,lun,gen,b,e_beam[0],alpha
	if keyword_set(fitgamma) then begin
	    sigmas = [1,3,5,7,9,11,13]
	    pis =    [0,2,4,6,8,10,12,14]
	    amplnorm = ampl
	    cosgamma2 = cos(data.angle.gamma[fb])^2
	    amplnorm[sigmas,*] /= (1d0 + cosgamma2)
	    amplnorm[pis,*] /= (1d0 - cosgamma2)
	    a2 = [b,cosgamma2]
	    parinfo2 = [parinfo,parinfo[0]]
	    parinfo2[7+ncomp].limits = [-0.1D,1.1D]
	    parinfo2[7+ncomp].step = 0.1

;	    res2 = fit_gamma(x[domain],a2,ampl=ampl,vsintheta=v_beam*sin(theta)*100,voverc=v_beam/299792458d,cosalpha=cos(alpha),sens=sens[domain],pixbound=pixbound,disp=disp)
	    b2 = mpfitfun('fit_gamma',x[domain],y[domain],yerrfix,a2,MAXITER=300,YFIT=yout2,perror=perror,FUNCTARGS={ampl:amplnorm,vsintheta:v_beam*sin(theta)*100,voverc:v_beam/299792458d,cosalpha:cos(alpha),sens:sens[domain],pixbound:pixbound,disp:disp},parinfo=parinfo2,NPRINT=0,/quiet,status=status)
	    if status le 0 then b2[*] = !values.d_nan
	    if gen ge 0 then result2a[gen] = b2[1]/10000 else result2a0 = b2[1]/10000
	    if gen ge 0 then result2b[gen] = b2[7+ncomp] else result2b0 = b2[7+ncomp]
	    printf,lun,'&',b2[0],b2[1],b2[2],b2[3],b2[4],b2[5],b2[6],b2[7],e_beam[0],alpha
	endif
	if keyword_set(fitimpurity) then begin
	    sigmas = [1,3,5,7,9,11,13]
	    pis =    [0,2,4,6,8,10,12,14]
	    amplnorm = ampl
	    cosgamma2 = cos(data.angle.gamma[fb])^2
	    amplnorm[sigmas,*] /= (1d0 + cosgamma2)
	    amplnorm[pis,*] /= (1d0 - cosgamma2)
	    a2 = [b,cosgamma2,200]
	    parinfo2 = [parinfo,parinfo[0],parinfo[0]]
	    parinfo2[7+ncomp].limits = [-0.1D,1.1D]
	    parinfo2[7+ncomp].step = 0.1
	    parinfo2[8+ncomp].limits = [0.0D,10000D]
	    parinfo2[8+ncomp].step = 100

;	    res2 = fit_gamma(x[domain],a2,ampl=ampl,vsintheta=v_beam*sin(theta)*100,voverc=v_beam/299792458d,cosalpha=cos(alpha),sens=sens[domain],pixbound=pixbound,disp=disp)
	    b2 = mpfitfun('fit_gamma_imp',x[domain],y[domain],yerrfix,a2,MAXITER=300,YFIT=yout2,perror=perror,FUNCTARGS={ampl:amplnorm,vsintheta:v_beam*sin(theta)*100,voverc:v_beam/299792458d,cosalpha:cos(alpha),sens:sens[domain],pixbound:pixbound,disp:disp,impwavelength:wveimpurity[gen+1]},parinfo=parinfo2,NPRINT=0,/quiet,status=status)
	    if status le 0 then b2[*] = !values.d_nan
	    if gen ge 0 then result2a[gen] = b2[1]/10000 else result2a0 = b2[1]/10000
	    if gen ge 0 then result2b[gen] = b2[7+ncomp] else result2b0 = b2[7+ncomp]
	    printf,lun,'&',b2[0],b2[1],b2[2],b2[3],b2[4],b2[5],b2[6],b2[7],e_beam[0],alpha
	endif
	if keyword_set(plot) then begin
	    oploti,x[domain],yout,icolor=2
stop
	endif
;	oplot,x[domain],yout
;	radii[chani] = su_rbm[su_ind]
	if gen ge 0 then resultB[gen] = b[1]/10000 else result0 = b[1]/10000
	if gen mod 100 eq 0 then begin
	    print,gen
	    figure,1
	    clf,/all
	    ploti,x[domain],y[domain],icolor=0,ptitle=ptitle,xtitle='Angstroms',ytitle='arb'
	    oploti,x[domain],res,icolor=1
	    oploti,x[domain],yout,icolor=2
	    if keyword_set(fitgamma) then begin
;		oploti,x[domain],res2,icolor=3
		oploti,x[domain],yout2,icolor=4
		legendi,label=['data','prefit','fit','fit2']
	    endif else legendi,label=['data','prefit','fit']
	endif
    endfor
    free_lun,lun
    filenameout='runs/'+run+'/'+synthfile+'.analyzed'
    gamma0 = data.angle.gamma[fb]
    save,result2a,result2b,b_tot,fb,n_gen,resultB,result0,e_beam,alpha,theta,gamma0,filename=filenameout
    saveallmain
;    plot_mse_ls_synth
end

pro plot_mse_ls_synth,special=special
    restoreallmain
    figure,3
    clf,/all
    iteration = dindgen(n_gen)+1

    ;calculate centroid b_tot
    b_tot0 = total(b_tot*n_beam_n3[*,0])/total(n_beam_n3[*,0])

    wveimpurity = maken(6560,6610,n_gen+1)
    if keyword_set(special) then oploti,ptitle='MSE-LS |B|',xtitle='Impurity wavelength (A)', ytitle='|B|, Tesla',wveimpurity[1:*],resultB,charsize=1.75,icolor=0 else $
    oploti,ptitle='MSE-LS |B|',xtitle='iteration', ytitle='|B|, Tesla',iteration,resultB,charsize=1.75,icolor=0
    plothl,result0,icolor=1
    plothl,b_tot[fb],icolor=2
    plothl,b_tot0,icolor=3
    legendi,icolor=[0,1,2,3],label=['Fitted synthetic','Fitted no-noise','Field at max beam','Field w average']

    B_err = stdev(resultB)
    B_ave = mean(resultB)
    print,'fitted |B|:'+strtrim(B_ave,2)+'  stddev:'+strtrim(B_err,2)
    print,'fitted |B| no-noise:'+strtrim(result0,2)
    print,'assumed |B| at beam center:'+strtrim(b_tot[fb],2)
    print,'beam weighted chord average |B|:'+strtrim(b_tot0,2)

    print,'Pitch angle used in generation (beam center):'+strtrim(data.angle.pitch[fb],2)+' = '+strtrim(data.angle.pitch[fb]*180/!dpi,2)+' deg'
    calcpitch = acos((B_tot[fb]*cos(data.angle.pitch[fb])/B_ave)<1)
    print,'Calculated pitch angle:'+strtrim(calcpitch,2)+' = '+strtrim(calcpitch*180/!dpi,2)+' deg'
    errpitch = B_err/B_ave/tan(data.angle.pitch[fb])
    print,'Propagated error in pitch:'+strtrim(errpitch,2)+' = '+strtrim(errpitch*180/!dpi,2)+' deg'
end

pro plot_mse_gamma,special=special
    restoreallmain
    figure,3
    clf,/all
    iteration = dindgen(n_gen)+1
    if keyword_set(special) then iteration = (maken(6560,6610,n_gen+1))[1:*]

    ;calculate centroid b_tot
    b_tot0 = total(b_tot*n_beam_n3[*,0])/total(n_beam_n3[*,0])

    oploti,ptitle='MSE-LS |B|',xtitle='iteration', ytitle='|B|, Tesla',iteration,result2a,charsize=1.75,icolor=0
    plothl,result2a0,icolor=1
    plothl,b_tot[fb],icolor=2
    plothl,b_tot0,icolor=3
    legendi,icolor=[0,1,2,3],label=['Fitted synthetic (fit_gamma)','Fitted no-noise','Field at max beam','Field w average']

    B_err = stdev(result2a)
    B_ave = mean(result2a)
    print,'fitted |B|:'+strtrim(B_ave,2)+'  stddev:'+strtrim(B_err,2)
    print,'fitted |B| no-noise:'+strtrim(result2a0,2)
    print,'assumed |B| at beam center:'+strtrim(b_tot[fb],2)
    print,'beam weighted chord average |B|:'+strtrim(b_tot0,2)

    figure,4
    clf,/all
    result2b[where(result2b eq cos(data.angle.gamma[fb])^2)] = !values.d_nan
    ploti,ptitle='MSE-LR cos(gamma)',xtitle='iteration', ytitle='cos(gamma)',iteration,result2b,charsize=1.75,icolor=0
    plothl,cos(data.angle.gamma[fb])^2
    legendi,icolor=[0,1],label=['Fitted cos^2(gamma)','Expected']
    print,'mean fitted cos^2(gamma):',mean(result2b)
    print,'actual cos^2(gamma) at beam center:'+strtrim(cos(data.angle.gamma[fb])^2,2)
    print,'Pitch angle used in generation (beam center):'+strtrim(data.angle.pitch[fb],2)+' = '+strtrim(data.angle.pitch[fb]*180/!dpi,2)+' deg'
    print,'mean LR fitted pitch angle (deg)',mean(acos(sqrt(result2b)),/nan)*180/!dpi
    paerror = stdev(acos(sqrt(result2b)))
    print,'stddev in LR fitted pitch angle (deg)',paerror*180/!dpi

;    figure,5
;    clf,/all
;    ploti,iteration,acos(sqrt(result2b))
end

pro radial_dependence
    radii=[6.2,6.5,6.8,7.1,7.4,7.7,8.0]
    b_err = dblarr(n_elements(radii))
    b_bias = dblarr(n_elements(radii))
    pa_error = dblarr(n_elements(radii))
    pa_bias = dblarr(n_elements(radii))
    for i=0,n_elements(radii)-1 do begin
;	nbass,'ITER',periscope='E-2',spot_radius=radii[i],profiles='kappatou'
	synthfile='synthetic_ITER_kappatou_E-2_R'+num2str(radii[i],3)+'.sav'
	MSE_LS_synth,run='ITER',synthfile=synthfile,fitgamma=1
	restoreallmain
	B_err[i] = stdev(result2a)
	B_bias[i] = mean(result2a) - b_tot[fb]
	pa_error[i] = stdev(acos(sqrt(abs(result2b)))) ;FIXME, absolute value
	pa_bias[i]  = acos(sqrt(abs(mean(result2b)))) - data.angle.gamma[fb] ;FIXME
    endfor
    figure,5
    clf,/all
    ploti,radii,b_err,charsize=1.8,xtitle='radius (m)',ytitle='|B| error (T)',ptitle='|B| Measurement Uncertainty',/ylog
    oploti,radii,pa_bias,icolor=1,psym=-4
    legendi,label=['std-dev(|B_{meas}|)','|mean(|B_{meas}| - |B_{true}|)|']
    figure,6
    clf,/all
    ploti,radii,pa_error,charsize=1.8,xtitle='radius (m)',ytitle='Pitch angle error (radians)',ptitle='Pitch Angle Measurement Uncertainty',/ylog,psym=-2
    oploti,radii,pa_bias,icolor=1,psym=-4
    legendi,label=['std-dev(\theta_{meas})','|mean(\theta_{meas} - \theta_{true})|']
    saveallmain
;    for i=0,n_elements(radii)-1 do begin
;	nbass,'ITER',periscope='E-2',spot_radius=radii[i],profiles='casper'
;    endfor
end

pro test_int_dependence,skip=skip,use_erf=use_erf
    if keyword_set(skip) then goto,skipped
    int_counts = makenlog(50d,30000d,10)
    sens_adj = int_counts/1500d
    result0s = dblarr(10)
    resultbias = dblarr(10)
    resultstd = dblarr(10)
    figure,3
    clf,/all
    for i=0,n_elements(sens_adj)-1 do begin
	MSE_LS_synth,synthfile='synth_intdep'+strtrim(i,2)+'.sav',use_erf=use_erf
	restoreallmain,except='I'
	result0s[i] = result0
	resultbias[i] = mean(resultB)
	resultstd[i] = stdev(resultB)
	figure,3
	oploti,x,yadd,color=gettrue(i+1)
    endfor
    saveallmain
    skipped:
    restoreallmain
    figure,3
    legendi,label=num2str(int_counts)
    figure,4
    clf,/all
    ploti,int_counts,resultbias-b_tot[fb],xtitle='spectrum peak intensity (counts)',ytitle='bias |B|',/xlog
    oploti,int_counts,result0s-b_tot[fb],icolor=1
    figure,5
    clf,/all
    ploti,int_counts,resultstd,xtitle='spectrum peak intensity (counts)',ytitle='stddev |B|',/xlog
    a = mpfitexpr('P[0]/sqrt(X)',int_counts,resultstd,replicate(1d0,10),[1d0],nprint=0)
    oploti,int_counts,a[0]/sqrt(int_counts),linestyle=2,icolor=1
end

pro test_broadening_dependence,skip=skip,use_erf=use_erf
    if keyword_set(skip) then goto,skipped
    broadening = dindgen(10)^2/10+0.1
    result0s = dblarr(10)
    resultbias = dblarr(10)
    resultstd = dblarr(10)
    figure,3
    clf,/all
    for i=0,n_elements(broadening)-1 do begin
	MSE_LS_synth,synthfile='synth_broaddep'+strtrim(i,2)+'.sav',width0=broadening[i],use_erf=use_erf
	restoreallmain,except='I'
	result0s[i] = result0
	resultbias[i] = mean(resultB)
	resultstd[i] = stdev(resultB)
	figure,3
	oploti,x,yadd,color=gettrue(i+1)
    endfor
    saveallmain
    skipped:
    restoreallmain
    figure,3
    legendi,label=num2str(broadening)
    figure,4
    clf,/all
    ploti,broadening,resultbias-b_tot[fb],xtitle='spectrum broadening (A)',ytitle='bias |B|',/xlog
    oploti,broadening,result0s-b_tot[fb],icolor=1
    figure,5
    clf,/all
    ploti,broadening,resultstd,xtitle='spectrum broadening (A)',ytitle='stddev |B|',/xlog
    a = mpfitexpr('P[0]*X^P[1]+P[2]',broadening,resultstd,replicate(1d0,10),[1d0,2d0,1d0],nprint=0)
    oploti,broadening,a[0]*broadening^a[1]+a[2],linestyle=2,icolor=1
end

pro test_broadening_dependence_7T,skip=skip,use_erf=use_erf
    if keyword_set(skip) then goto,skipped
    broadening = dindgen(10)^2/10+0.1
    result0s = dblarr(10)
    resultbias = dblarr(10)
    resultstd = dblarr(10)
    figure,3
    clf,/all
    for i=0,n_elements(broadening)-1 do begin
	MSE_LS_synth,synthfile='synth_broaddep_7T'+strtrim(i,2)+'.sav',width0=broadening[i],use_erf=use_erf,b_tor0=7d
	restoreallmain,except=['I','broadening']
	result0s[i] = result0
	resultbias[i] = mean(resultB)
	resultstd[i] = stdev(resultB)
	figure,3
	oploti,x,yadd,color=gettrue(i+1)
    endfor
    saveallmain
    skipped:
    restoreallmain
    figure,3
    legendi,label=num2str(broadening)
    figure,4
    clf,/all
    ploti,broadening,resultbias-b_tot[fb],xtitle='spectrum broadening (A)',ytitle='bias |B|',/xlog
    oploti,broadening,result0s-b_tot[fb],icolor=1
    figure,5
    clf,/all
    ploti,broadening,resultstd,xtitle='spectrum broadening (A)',ytitle='stddev |B|',/xlog
    a = mpfitexpr('P[0]*X^P[1]+P[2]',broadening,resultstd,replicate(1d0,10),[1d0,2d0,1d0],nprint=0)
    oploti,broadening,a[0]*broadening^a[1]+a[2],linestyle=2,icolor=1
end

pro test_field_dependence,skip=skip,use_erf=use_erf
    if keyword_set(skip) then goto,skipped
    field = maken(4,10,10)
    result0s = dblarr(10)
    resultbias = dblarr(10)
    resultstd = dblarr(10)
    figure,3
    clf,/all
    for i=0,n_elements(field)-1 do begin
	MSE_LS_synth,synthfile='synth_fielddep'+strtrim(i,2)+'.sav',b_tor0=field[i],use_erf=use_erf
	restoreallmain,except='I'
	result0s[i] = result0
	resultbias[i] = mean(resultB)
	resultstd[i] = stdev(resultB)
	figure,3
	oploti,x,yadd,color=gettrue(i+1)
    endfor
    saveallmain
    skipped:
    restoreallmain
    btotsav = dblarr(10)
    for i=0,9 do begin
	Calc_Geom_Angles_chord,shot=shot,row=row,col=col,t1=t1,t2=t2,angle=angle,dangle_fg=dangle_fg,dangle_ap=dangle_ap,dangle_sp=dangle_sp,chord=chord,b_tot=b_tot,aper_pos_shift=aper_pos_shift,ap_rad=ap_rad,b_tor0=field[i],r0=r0,pitch=pitch
	btotsav[i]=b_tot[fb]
    endfor
    figure,3
    legendi,label=num2str(field)
    figure,4
    clf,/all
    ploti,field,resultbias-btotsav,xtitle='B_{Tor} (T)',ytitle='bias |B|'
    oploti,field,result0s-btotsav,icolor=1
    figure,5
    clf,/all
    ploti,field,resultstd,xtitle='B_{Tor} (T)',ytitle='stddev |B|'
    a = mpfitexpr('P[0]/sqrt(X)',field,resultstd,replicate(1d0,10),[1d0],nprint=0)
    oploti,field,a[0]/sqrt(field),linestyle=2,icolor=1
end

pro test_wve_dependence,skip=skip,use_erf=use_erf
    if keyword_set(skip) then goto,skipped
    dwve = dindgen(10)^2/20+0.1
    result0s = dblarr(10)
    resultbias = dblarr(10)
    resultstd = dblarr(10)
    figure,3
    clf,/all
    for i=0,n_elements(dwve)-1 do begin
	MSE_LS_synth,synthfile='synth_wvedep_b1'+strtrim(i,2)+'.sav',use_erf=use_erf
	restoreallmain,except='I'
	result0s[i] = result0
	resultbias[i] = mean(resultB)
	resultstd[i] = stdev(resultB)
	figure,3
	oploti,x,yadd,color=gettrue(i+1)
    endfor
    saveallmain
    skipped:
    restoreallmain
    figure,3
    legendi,label=num2str(dwve)
    figure,4
    clf,/all
    ploti,dwve,resultbias-b_tot[fb],xtitle='disp (A)',ytitle='mean |B|'
    oploti,dwve,result0s-b_tot[fb],icolor=1
    figure,5
    clf,/all
    ploti,dwve,resultstd,xtitle='disp (A)',ytitle='stddev |B|'
    a = mpfitexpr('P[0]/sqrt(X)',dwve,resultstd,replicate(1d0,10),[1d0],nprint=0)
    oploti,dwve,a[0]/sqrt(dwve),linestyle=2,icolor=1
end

pro test_bgmismatch_dependence,skip=skip,use_erf=use_erf
    if keyword_set(skip) then goto,skipped
    bg_mismatch=maken(-0.5,0.5,10)
    result0s = dblarr(10)
    resultbias = dblarr(10)
    resultstd = dblarr(10)
    figure,3
    clf,/all
    for i=0,n_elements(bg_mismatch)-1 do begin
	MSE_LS_synth,synthfile='synth_bgmismatch'+strtrim(i,2)+'.sav',use_erf=use_erf
	restoreallmain,except=['I','bg_mismatch']
	result0s[i] = result0
	resultbias[i] = mean(resultB)
	resultstd[i] = stdev(resultB)
	figure,3
	oploti,x,yadd,color=gettrue(i+1)
    endfor
    saveallmain
    skipped:
    restoreallmain
    figure,3
    legendi,label=num2str(bg_mismatch)
    figure,4
    clf,/all
    ploti,bg_mismatch,resultbias-b_tot[fb],xtitle='background mismatch',ytitle='bias |B|',/xlog
    oploti,bg_mismatch,result0s-b_tot[fb],icolor=1
    figure,5
    clf,/all
    ploti,bg_mismatch,resultstd,xtitle='background mismatch',ytitle='stddev |B|',/xlog
    ;a = mpfitexpr('P[0]/sqrt(X)',bg_mismatch,resultstd,replicate(1d0,10),[1d0],nprint=0)
    ;oploti,bg_mismatch,a[0]/sqrt(bg_mismatch),linestyle=2,icolor=1
end

pro test_bglevel_dependence,skip=skip,use_erf=use_erf
    if keyword_set(skip) then goto,skipped
    bglevel = dindgen(10)^2/5
    result0s = dblarr(10)
    resultbias = dblarr(10)
    resultstd = dblarr(10)
    figure,3
    clf,/all
    for i=0,n_elements(bglevel)-1 do begin
	MSE_LS_synth,synthfile='synth_bglevel'+strtrim(i,2)+'.sav',use_erf=use_erf
	restoreallmain,except=['I','bglevel']
	result0s[i] = result0
	resultbias[i] = mean(resultB)
	resultstd[i] = stdev(resultB)
	figure,3
	oploti,x,yadd,color=gettrue(i+1)
    endfor
    saveallmain
    skipped:
    restoreallmain
    figure,3
    legendi,label=num2str(bglevel)
    figure,4
    clf,/all
    ploti,bglevel,resultbias-b_tot[fb],xtitle='bremsstrahlung level',ytitle='bias |B|',/xlog
    oploti,bglevel,result0s-b_tot[fb],icolor=1
    figure,5
    clf,/all
    ploti,bglevel,resultstd,xtitle='bremsstrahlung level',ytitle='stddev |B|',/xlog
    ;a = mpfitexpr('P[0]/sqrt(X)',bg_mismatch,resultstd,replicate(1d0,10),[1d0],nprint=0)
    ;oploti,bg_mismatch,a[0]/sqrt(bg_mismatch),linestyle=2,icolor=1
end

;fit the stark zeeman pattern with a pi to sigma ratio
;input ampl should already be normalized such that pi to sigma ratio is such that cos(gamma) = 1
;A[0] wavelength shift (due to miscalibration or misunderstanding of h-alpha center)
;A[1] |B|
;A[2] amplitude scaling on total E/1 component
;A[3] width of line
;A[4] amplitude of CXRS gaussian
;A[5] center of CXRS gaussian
;A[6] width of CXRS gaussian
;A[7] cos^2(gamma) == (1-r)/(1+r), where r is pi to sigma ratio
function fit_gamma,xfit,a,ampl=ampl,vsintheta=vsintheta,voverc=voverc,cosalpha=cosalpha,sens=sens,pixbound=pixbound,disp=disp
    c = 29979245800d            ; cm/s     speed of light in vacuum
    H_alpha = 6562.819d         ; A (air)
    hbar = 1.054571726d-27      ; erg s
    ec = 4.8032068d-10          ; esu
    a0 = 0.52917721092d-8       ; cm       Bohr radius
    bohrmag = 9.27400968d-21    ; erg/G    Bohr magneton = e*hbar/(2*me*c)
    m_e = 9.10938291d-28        ; g        electron mass
    alpha = 7.2973525698d-3     ;          Fine structure constant

    sigmas = [0,1,0,1,0,1,0,1,0,1,0,1,0,1,0]
    pis =    [1,0,1,0,1,0,1,0,1,0,1,0,1,0,1]
    ncomp = n_elements(a)-8 ;1 less than the number of beam components
    scaleratio = sigmas*(1d0+a[7+ncomp]) + pis*(1d0-a[7+ncomp])

    B_tot = A[1]
    gam = bohrmag*B_tot         ;gamma = Bohr magneton times B_t (erg)

    k3 = [-2, -1, 0, 1, 2]
    m3sq = [0, 1, 0, 1, 0]
    k2 = [-1, 0, 1]
    m2sq = [0, 1, 0]

    yout = dblarr(n_elements(xfit))
    for bcomp=0,ncomp do begin
	eps = 3*ec*a0*vsintheta[bcomp]/c*B_tot ;epsilon = Stark energy splitting coefficient (erg) = 2*E_L
	qc = -(ec*a0*vsintheta[bcomp]/c*B_tot)^2/(m_e*c^2*alpha^2)/16d ; = -1/16*E_L^2 'quadratic Stark constant'

	q0=sqrt(gam^2+eps^2)        ;n=2 energy shift parameter (erg)
	q1=sqrt(4*gam^2+9*eps^2)    ;n=3 energy shift parameter (erg)

	e3 = [-q1,-q1/2,0,q1/2,q1] + qc*3^4*(17*3^2-3*k3^2-9*m3sq+19)

	;e3_adjust = -qc*3^4*9*mtwo/(mone+mzero)*4
	e3_adjust = -qc*3^4*9*0.58*4 ;from bes_spectrum_stark_zeeman simulation

	e2 = [-q0,0,q0] + qc*2^4*(17*2^2-3*k2^2-9*m2sq+19)

	;shifts = -outer_sum(e3,-e2)*h_alpha^2/(hbar*2*!dpi*c)/1d8 ; Angstroms
	shifts = -outer_sum(-e2,e3)
	shifts[1,2] -= e3_adjust
	shifts *= h_alpha^2/(hbar*2*!dpi*c)/1d8 ; Angstroms

	lines = (h_alpha + shifts + a[0])/(1+voverc[bcomp]*cosalpha)/(sqrt(1-voverc[bcomp]^2)) ;apply Doppler shift
	for linei=0,14 do begin
	    erf0 = erf((pixbound-lines[linei])/(a[3+ncomp]*sqrt(2)))
	    ;oneline = sqrt(!dpi/2)*ampl[linei,bcomp]*a[2+bcomp] * ([1d,erf0] - [erf0,-1d])*scaleratio[linei]
	    oneline = sqrt(!dpi/2)*ampl[linei,bcomp]*a[2+bcomp] * (erf0[0:n_elements(erf0)-2] - erf0[1:n_elements(erf0)-1])*scaleratio[linei]
	    yout += oneline
	endfor
	;try a more vectorized version. not really faster. actually a tiny bit slower. Maybe erf is internally not vectorized
;	erf0 = erf(outer_sum(pixbound,-reform(lines,15))/(a[3]*sqrt(2)))
;        oneline = extend(n_elements(xfit),sqrt(!dpi/2)*ampl[*,bcomp]*a[2+bcomp]) * ([ones(1,15),erf0] - [erf0,-ones(1,15)])
;	yout += total(oneline,2)
    endfor
    yout *= sens / (-disp)
    yout += A[4+ncomp]*exp(-(xfit-A[5+ncomp])^2/A[6+ncomp]) ;CXRS can be done better.. just a quick fix
    return,yout
end

pro test_domain_dependence,skip=skip,use_erf=use_erf
    domainends=maken(6620d,6650d,10)
    if keyword_set(skip) then goto,skipped
    result0s = dblarr(10)
    resultbias = dblarr(10)
    resultstd = dblarr(10)
    figure,3
    clf,/all
    for i=0,n_elements(domainends)-1 do begin
	MSE_LS_synth,synthfile='synthetic1120621025_22,23_21,24_aprad1.sav',use_erf=use_erf,domainend=domainends[i]
	restoreallmain,except=['I','domainends']
	result0s[i] = result0
	resultbias[i] = mean(resultB)
	resultstd[i] = stdev(resultB)
	figure,3
	oploti,x,yadd,color=gettrue(i+1)
    endfor
    saveallmain
    skipped:
    restoreallmain
    figure,3
    legendi,label=num2str(domainends)
    figure,4
    clf,/all
    ploti,domainends,resultbias-b_tot[fb],xtitle='domain end',ytitle='bias |B|',/xlog
    oploti,domainends,result0s-b_tot[fb],icolor=1
    figure,5
    clf,/all
    ploti,domainends,resultstd,xtitle='domain end',ytitle='stddev |B|',/xlog
    ;a = mpfitexpr('P[0]/sqrt(X)',bg_mismatch,resultstd,replicate(1d0,10),[1d0],nprint=0)
    ;oploti,bg_mismatch,a[0]/sqrt(bg_mismatch),linestyle=2,icolor=1
end

pro test_domain_dependence2,skip=skip,use_erf=use_erf
    domainstarts=maken(6520d,6590d,10)
    if keyword_set(skip) then goto,skipped
    result0s = dblarr(10)
    resultbias = dblarr(10)
    resultstd = dblarr(10)
    figure,3
    clf,/all
    for i=0,n_elements(domainstarts)-1 do begin
	MSE_LS_synth,synthfile='synthetic1120621025_22,23_21,24_aprad1.sav',use_erf=use_erf,domainstart=domainstarts[i]
	restoreallmain,except=['I','domainstarts']
	result0s[i] = result0
	resultbias[i] = mean(resultB)
	resultstd[i] = stdev(resultB)
	figure,3
	oploti,x,yadd,color=gettrue(i+1)
    endfor
    saveallmain
    skipped:
    restoreallmain
    figure,3
    legendi,label=num2str(domainstarts)
    figure,4
    clf,/all
    ploti,domainstarts,resultbias-b_tot[fb],xtitle='domain end',ytitle='bias |B|',/xlog
    oploti,domainstarts,result0s-b_tot[fb],icolor=1
    figure,5
    clf,/all
    ploti,domainstarts,resultstd,xtitle='domain end',ytitle='stddev |B|',/xlog
    ;a = mpfitexpr('P[0]/sqrt(X)',bg_mismatch,resultstd,replicate(1d0,10),[1d0],nprint=0)
    ;oploti,bg_mismatch,a[0]/sqrt(bg_mismatch),linestyle=2,icolor=1
end

pro aps_2016_plot
    radii=[0.65,0.7,0.75,0.8,0.85]
    b_err = dblarr(n_elements(radii))
    b_bias = dblarr(n_elements(radii))
    pa_error = dblarr(n_elements(radii))
    pa_bias = dblarr(n_elements(radii))
    for i=0,n_elements(radii)-1 do begin
	synthfile='CMOD_'+num2str(radii[i],3)+'.sav'
	synthfile2='CMOD_'+num2str(radii[i],3)+'.sav.analyzed'
	if file_test('runs/CMOD_SMSE/'+synthfile2) then begin
	    print,'restoring '+synthfile2
	    restore,'runs/CMOD_SMSE/'+synthfile2
	endif else begin
	    print,'did not find '+synthfile2
	    MSE_LS_synth,run='CMOD_SMSE',synthfile=synthfile,polfilter='PEMs',/fitgamma
	    restore,'runs/CMOD_SMSE/'+synthfile2
;	    restoreallmain,except=['I','B_ERR','B_BIAS','PA_ERROR','PA_BIAS','RADII']
	endelse
	B_err[i] = stdev(result2a)
	B_bias[i] = mean(result2a) - b_tot[fb]
	pa_error[i] = stdev(acos(sqrt(abs(result2b)))) ;FIXME, absolute value
	;pa_bias[i]  = acos(sqrt(abs(mean(result2b)))) - data.angle.gamma[fb] ;FIXME
	pa_bias[i]  = acos(sqrt(abs(mean(result2b)))) - gamma0 ;FIXME
    endfor
    figure,5
    clf,/all
    ploti,radii,b_err,charsize=1.8,xtitle='radius (m)',ytitle='|B| error (T)',ptitle='|B| Measurement Uncertainty',/ylog
    oploti,radii,b_bias,icolor=1,psym=-4
    legendi,label=['std-dev(|B_{meas}|)','|mean(|B_{meas}| - |B_{true}|)|']
    figure,6
    clf,/all
    ploti,radii,pa_error*180/!dpi,charsize=1.8,xtitle='radius (m)',ytitle='Pitch angle error (deg.)',ptitle='Pitch Angle Measurement Uncertainty',/ylog,psym=-2
    oploti,radii,pa_bias*180/!dpi,icolor=1,psym=-4
    legendi,label=['std-dev(\theta_{meas})','|mean(\theta_{meas} - \theta_{true})|']
    saveallmain
end

;fit the stark zeeman pattern with a pi to sigma ratio
;input ampl should already be normalized such that pi to sigma ratio is such that cos(gamma) = 1
;A[0] wavelength shift (due to miscalibration or misunderstanding of h-alpha center)
;A[1] |B|
;A[2] amplitude scaling on total E/1 component
;A[3+nc] width of line
;A[4+nc] amplitude of CXRS gaussian
;A[5+nc] center of CXRS gaussian
;A[6+nc] width of CXRS gaussian
;A[7+nc] cos^2(gamma) == (1-r)/(1+r), where r is pi to sigma ratio
;A[8+nc] amplitude of impurity. We assume that the wavelength and width are known
function fit_gamma_imp,xfit,a,ampl=ampl,vsintheta=vsintheta,voverc=voverc,cosalpha=cosalpha,sens=sens,pixbound=pixbound,disp=disp,impwavelength=impwavelength
    c = 29979245800d            ; cm/s     speed of light in vacuum
    H_alpha = 6562.819d         ; A (air)
    hbar = 1.054571726d-27      ; erg s
    ec = 4.8032068d-10          ; esu
    a0 = 0.52917721092d-8       ; cm       Bohr radius
    bohrmag = 9.27400968d-21    ; erg/G    Bohr magneton = e*hbar/(2*me*c)
    m_e = 9.10938291d-28        ; g        electron mass
    alpha = 7.2973525698d-3     ;          Fine structure constant

    sigmas = [0,1,0,1,0,1,0,1,0,1,0,1,0,1,0]
    pis =    [1,0,1,0,1,0,1,0,1,0,1,0,1,0,1]
    ncomp = n_elements(a)-9 ;1 less than the number of beam components
    scaleratio = sigmas*(1d0+a[7+ncomp]) + pis*(1d0-a[7+ncomp])

    B_tot = A[1]
    gam = bohrmag*B_tot         ;gamma = Bohr magneton times B_t (erg)

    k3 = [-2, -1, 0, 1, 2]
    m3sq = [0, 1, 0, 1, 0]
    k2 = [-1, 0, 1]
    m2sq = [0, 1, 0]

    yout = dblarr(n_elements(xfit))
    for bcomp=0,ncomp do begin
	eps = 3*ec*a0*vsintheta[bcomp]/c*B_tot ;epsilon = Stark energy splitting coefficient (erg) = 2*E_L
	qc = -(ec*a0*vsintheta[bcomp]/c*B_tot)^2/(m_e*c^2*alpha^2)/16d ; = -1/16*E_L^2 'quadratic Stark constant'

	q0=sqrt(gam^2+eps^2)        ;n=2 energy shift parameter (erg)
	q1=sqrt(4*gam^2+9*eps^2)    ;n=3 energy shift parameter (erg)

	e3 = [-q1,-q1/2,0,q1/2,q1] + qc*3^4*(17*3^2-3*k3^2-9*m3sq+19)

	;e3_adjust = -qc*3^4*9*mtwo/(mone+mzero)*4
	e3_adjust = -qc*3^4*9*0.58*4 ;from bes_spectrum_stark_zeeman simulation

	e2 = [-q0,0,q0] + qc*2^4*(17*2^2-3*k2^2-9*m2sq+19)

	;shifts = -outer_sum(e3,-e2)*h_alpha^2/(hbar*2*!dpi*c)/1d8 ; Angstroms
	shifts = -outer_sum(-e2,e3)
	shifts[1,2] -= e3_adjust
	shifts *= h_alpha^2/(hbar*2*!dpi*c)/1d8 ; Angstroms

	lines = (h_alpha + shifts + a[0])/(1+voverc[bcomp]*cosalpha)/(sqrt(1-voverc[bcomp]^2)) ;apply Doppler shift
	for linei=0,14 do begin
	    erf0 = erf((pixbound-lines[linei])/(a[3+ncomp]*sqrt(2)))
	    ;oneline = sqrt(!dpi/2)*ampl[linei,bcomp]*a[2+bcomp] * ([1d,erf0] - [erf0,-1d])*scaleratio[linei]
	    oneline = sqrt(!dpi/2)*ampl[linei,bcomp]*a[2+bcomp] * (erf0[0:n_elements(erf0)-2] - erf0[1:n_elements(erf0)-1])*scaleratio[linei]
	    yout += oneline
	endfor
	;try a more vectorized version. not really faster. actually a tiny bit slower. Maybe erf is internally not vectorized
;	erf0 = erf(outer_sum(pixbound,-reform(lines,15))/(a[3]*sqrt(2)))
;        oneline = extend(n_elements(xfit),sqrt(!dpi/2)*ampl[*,bcomp]*a[2+bcomp]) * ([ones(1,15),erf0] - [erf0,-ones(1,15)])
;	yout += total(oneline,2)
    endfor
    yout *= sens / (-disp)
    yout += A[4+ncomp]*exp(-(xfit-A[5+ncomp])^2/A[6+ncomp]) ;CXRS can be done better.. just a quick fix
    yout += A[8+ncomp]*exp(-(xfit-impwavelength)^2/(A[3+ncomp]*sqrt(2))) ;CXRS can be done better.. just a quick fix
    return,yout
end
