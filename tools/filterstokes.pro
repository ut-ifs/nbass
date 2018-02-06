;Apply a Gaussian bandpass filter to the results of a run
pro filterstokes,run,filename,bandwidth=bandwidth,center=center
    if n_elements(filename) eq 0 then filename=run+'.sav'
    restore,getenv('NBASS_RUNS_PATH')+'/'+run+'/'+filename

    defvalue,bandwidth,2.0 ;A, FWHM
    defvalue,center,6600.0

    wve = result.wve
    ;Assume that the bandpass filter has a Gaussian shaped transmission function
    msefilter = exp(-(wve-center)^2/(bandwidth^2/4/alog(2)))

    figure,1
    clf,/all
    ploti,result.wve,(result.bes_data.stokes_s0)/max(result.bes_data.stokes_s0),charsize=1.6,ptitle='Polarization',icolor=0
    oploti,result.wve,result.bes_data.poldegree,icolor=1
;    oploti,result.wve,result.bes_data.spsi,icolor=2
;    oploti,result.wve,result.bes_data.schi,icolor=3
    oploti,result.wve,msefilter,icolor=4
    legendi,label=['I/Imax','polarization','Filter'],charsize=1.6,icolor=[0,1,4]

    figure,2
    clf,/all
    ploti,result.wve,(result.bes_data.stokes_s0)/max(result.bes_data.stokes_s0)*msefilter,charsize=1.6,ptitle='Seen on spectrometer',icolor=0

    int_s0 = int_tabulated(result.wve,result.bes_data.stokes_s0*msefilter)
    int_s1 = int_tabulated(result.wve,result.bes_data.stokes_s1*msefilter)
    int_s2 = int_tabulated(result.wve,result.bes_data.stokes_s2*msefilter)
    int_s3 = int_tabulated(result.wve,result.bes_data.stokes_s3*msefilter)
    print,'integrated s0,s1,s2,s3:',int_s0,int_s1,int_s2,int_s3

    int_pol = sqrt(int_s1^2 + int_s2^2 + int_s3^2)/int_s0
    int_psi = atan(int_s2,int_s1)/2
    int_chi = atan(int_s3,sqrt(int_s1^2+int_s2^2))/2
    print,'polarization:',int_pol
    print,'psi (pol. angle):',int_psi
    print,'chi (ellipticity):',int_chi

    filtercenters = maken(center-8d,center+8d,200)
    all_pol = dblarr(200)
    all_psi = dblarr(200)
    all_chi = dblarr(200)
    for i=0,n_elements(filtercenters)-1 do begin
	filtercenter = filtercenters[i]
	msefilter = exp(-(wve-filtercenter)^2/(bandwidth^2/4/alog(2)))
	int_s0 = int_tabulated(result.wve,result.bes_data.stokes_s0*msefilter)
	int_s1 = int_tabulated(result.wve,result.bes_data.stokes_s1*msefilter)
	int_s2 = int_tabulated(result.wve,result.bes_data.stokes_s2*msefilter)
	int_s3 = int_tabulated(result.wve,result.bes_data.stokes_s3*msefilter)
	all_pol[i] = sqrt(int_s1^2 + int_s2^2 + int_s3^2)/int_s0
	all_psi[i] = atan(int_s2,int_s1)/2
	all_chi[i] = atan(int_s3,sqrt(int_s1^2+int_s2^2))/2
    endfor
    figure,3
    clf,/all
    ploti,filtercenters,all_pol,charsize=1.6,ptitle='Polarization vs filter center',xtitle='wavelength',ytitle='polarization',icolor=0
    figure,4
    clf,/all
    ploti,filtercenters,all_psi,charsize=1.6,ptitle='psi vs filter center',xtitle='wavelength',ytitle='psi (pol. angle)',icolor=0
    figure,5
    clf,/all
    ploti,filtercenters,all_chi,charsize=1.6,ptitle='chi vs filter center',xtitle='wavelength',ytitle='chi (ellipticity)',icolor=0
end
