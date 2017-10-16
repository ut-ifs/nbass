;Try to calculate the optimal bandwidth for an MSE filter
;
;This tool uses the tnmin optimization package,
;available from cow.physics.wisc.edu/~craigm/idl/idl.html

function bestfilter_sub,a,wve=wve,s0=s0,s1=s1,s2=s2,s3=s3
    msefilter = 0.8*exp(-(wve-a[0])^2/(a[1]^2/4/alog(2)))
    int_s0 = int_tabulated(wve,s0*msefilter)
    int_s1 = int_tabulated(wve,s1*msefilter)
    int_s2 = int_tabulated(wve,s2*msefilter)
    int_s3 = int_tabulated(wve,s3*msefilter)
    sig = sqrt(int_s1^2 + int_s2^2 + int_s3^2)/sqrt(int_s0)
    return,sig
end

pro bestfilter,run,filename,a0=a0
    if n_elements(filename) eq 0 then filename=run+'.sav'
    restore,getenv('NBASS_RUNS_PATH')+'/'+run+'/'+filename

    s0 = total(result.bes_data.stokes_s0,2) + total(result.cxrs_data.stokes_s0,2) + result.edge_spec + result.brem_spec
    s1 = total(result.bes_data.stokes_s1,2) + total(result.cxrs_data.stokes_s1,2)
    s2 = total(result.bes_data.stokes_s2,2) + total(result.cxrs_data.stokes_s2,2)
    s3 = total(result.bes_data.stokes_s3,2) + total(result.cxrs_data.stokes_s3,2)

    wve = result.wve
    if n_elements(a0) eq 0 then a0 = [6600.6d,2.0d]
    a = tnmin('bestfilter_sub', a0, functargs={wve:result.wve,s0:s0,s1:s1,s2:s2,s3:s3},/maximize,/autoderivative,/quiet)
    print,'New a',a

    bandwidth = a[1] ;FWHM
    filtercenter = a[0]
    filtertrans = 0.8
    ;filter shape. assume a gaussian
    msefilter = filtertrans*exp(-(wve-filtercenter)^2/(bandwidth^2/4/alog(2)))

    figure,3
    clf,/all
    ploti,result.wve,s0/max(s0),charsize=1.6,ptitle='Polarization',icolor=0
    poldegree = sqrt(s1^2 + s2^2 + s3^2)/s0
    linpoldegree = sqrt(s1^2 + s2^2)/s0
    oploti,result.wve,poldegree,icolor=1
    oploti,result.wve,linpoldegree,icolor=2
    oploti,result.wve,msefilter,icolor=4
    legendi,label=['I/Imax','polarization','linear pol.','Filter'],icolor=[0,1,2,4],charsize=1.6

    figure,4
    clf,/all
    ploti,result.wve,s0/max(s0)*msefilter,charsize=1.6,ptitle='Seen on spectrometer'

    int_s0 = int_tabulated(result.wve,s0*msefilter)
    int_s1 = int_tabulated(result.wve,s1*msefilter)
    int_s2 = int_tabulated(result.wve,s2*msefilter)
    int_s3 = int_tabulated(result.wve,s3*msefilter)
    print,'integrated s0,s1,s2,s3:',int_s0,int_s1,int_s2,int_s3

    signal = sqrt(int_s1^2 + int_s2^2 + int_s3^2)
    int_pol = signal/int_s0
    int_psi = atan(int_s2,int_s1)/2
    int_chi = atan(int_s3,sqrt(int_s1^2+int_s2^2))/2
    print,'pol:',int_pol
    print,'psi:',int_psi
    print,'chi:',int_chi
    print,'signal:',signal
    print,'SNR',signal/sqrt(int_s0)
end

pro filter_result,run,filename
    if n_elements(filename) eq 0 then filename=run+'.sav'
    restore,getenv('NBASS_RUNS_PATH')+'/'+run+'/'+filename

    s0 = result.bes_data.stokes_s0
    s1 = result.bes_data.stokes_s1
    s2 = result.bes_data.stokes_s2
    s3 = result.bes_data.stokes_s3

    figure,1
    clf,/all
    angles = maken(0,!dpi,500)
    pi2sigma = dblarr(500)
    for i=0,n_elements(angles)-1 do begin
	filterangle=angles(i)
	new_s0 = 0.5d*(s0 + cos(2*filterangle)*s1 + sin(2*filterangle)*s2)
	signal = new_s0 * result.scalefactor + result.brem_counts*0.5 + result.cxrs_counts*0.5
;	oploti,result.wve,signal,color=gradcolor(i,30)
	pi2sigma[i] = signal[mindex(result.wve,6519.3)]/signal[mindex(result.wve,6522.9)]
    endfor
    figure,2
    ploti,angles,pi2sigma
    a=max(pi2sigma,maxind)
    cut_sigma = angles[maxind]
    cut_pi = cut_sigma - !dpi/2
    ;max pi2sigma at 1.627

    sz = size(result.noisy)
    n_gen = sz[0] eq 1 ? 1 : sz[2]
    ;cut sigma
    new_s0 = 0.5d*(s0 + cos(2*cut_sigma)*s1 + sin(2*cut_sigma)*s2)
    pure = new_s0 * result.scalefactor + result.brem_counts*0.5d + result.cxrs_counts*0.5d
    pure_bk = result.brem_counts*0.5d
    noisy = dblarr(n_elements(result.wve),n_gen)
    noisy_bk = dblarr(n_elements(result.wve),n_gen)
    for i=0,n_elements(result.wve)-1 do begin
	noisy[i,*] = randomu(seed,n_gen,poisson=pure[i],/double)
	noisy_bk[i,*] = randomu(seed,n_gen,poisson=pure_bk[i],/double)
    endfor
    result.pure = pure
    result.pure_bk = pure_bk
    result.noisy = noisy
    result.noisy_bk = noisy_bk
    save,parameters,data,result,filename='runs/'+run+'/polfilter_cutsigma.sav',description='cut sigma components with filter'

    ;cut pi
    new_s0 = 0.5d*(s0 + cos(2*cut_pi)*s1 + sin(2*cut_pi)*s2)
    pure = new_s0 * result.scalefactor + result.brem_counts*0.5d + result.cxrs_counts*0.5d
    pure_bk = result.brem_counts*0.5d
    noisy = dblarr(n_elements(result.wve),n_gen)
    noisy_bk = dblarr(n_elements(result.wve),n_gen)
    for i=0,n_elements(result.wve)-1 do begin
	noisy[i,*] = randomu(seed,n_gen,poisson=pure[i],/double)
	noisy_bk[i,*] = randomu(seed,n_gen,poisson=pure_bk[i],/double)
    endfor
    result.pure = pure
    result.pure_bk = pure_bk
    result.noisy = noisy
    result.noisy_bk = noisy_bk
    save,parameters,data,result,filename='runs/'+run+'/polfilter_cutpi.sav',description='cut pi components with filter'

;    ;no filter. just high brem
;    new_s0 = s0
;    pure = new_s0 * result.scalefactor + result.brem_counts*0.5d + result.cxrs_counts*0.5d
;    pure_bk = result.brem_counts*0.5d
;    noisy = dblarr(n_elements(result.wve),n_gen)
;    noisy_bk = dblarr(n_elements(result.wve),n_gen)
;    for i=0,n_elements(result.wve)-1 do begin
;	noisy[i,*] = randomu(seed,n_gen,poisson=pure[i],/double)
;	noisy_bk[i,*] = randomu(seed,n_gen,poisson=pure_bk[i],/double)
;    endfor
;    result.pure = pure
;    result.pure_bk = pure_bk
;    result.noisy = noisy
;    result.noisy_bk = noisy_bk
;    save,parameters,data,result,filename='runs/'+run+'/highbrem.sav',description='cut pi components with filter'
end

