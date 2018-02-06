@lib/get_alcbeam.pro
@lib/read_alcbeam.pro
@lib/vector.pro
@lib/extend.pro
@lib/outer.pro
@lib/mindex.pro
@lib/cmreplicate.pro
@lib/tic.pro
@lib/toc.pro
@lib/num2str.pro
@lib/savetomain.pro
@lib/defvalue.pro
@lib/figure.pro
@lib/sign.pro
@calc_geom_angles.pro
@calc_equil.pro
@get_beam_data.pro
@make_chord.pro
@get_profiles.pro
@synth_bes.pro
@synth_cxrs.pro
@synth_edge.pro
@synth_bremsstrahlung.pro
@get_parameters.pro

; This is the main level program
; Input:
;    run: name of a subdirectory in the runs directory. This subdirectory should contain a file
;         called parameters.pro which will be called automatically to load the run parameters
;    _EXTRA: extra parameters will be passed to the call of parameters.pro
pro nbass,run,_REF_EXTRA=extra
    @common_param.idl
    tic ;start calculation timer
    print,'Neutral Beam Active Spectroscopy Simulation (v 1.02)'
    print,'contact: Ken Liao <kenliao@physics.utexas.edu>'

    defvalue,run,'test'
    print,'Loading parameters for run '+run
    oldpath = !path
    !path = getenv('NBASS_RUNS_PATH')+'/'+run+':'+!path
    resolve_routine,'parameters',/compile_full_file,/either
    parameters,_EXTRA=extra
    !path = oldpath
    print,'done.'
    toc
    help,/memory
    print,''

    scalefactor = sens/10000 * detgain * int_time * abs(disp) ;conversion from ph/m2/s/sr/A to counts
    bes_spec = dblarr(n_elements(wve),num_beams>1)
    cxrs_spec = dblarr(n_elements(wve),num_beams>1)

    print,'Creating calculation grid'
    make_chord,chord
    print,'done.'
    toc
    help,/memory
    print,''

    print,'Calculating magnetic equilibrium'
    calc_equil,chord,field
    print,'done.'
    toc
    help,/memory
    print,''

    print,'Interpolating plasma profiles'
    get_profiles,chord,profiles
    print,'done.'
    toc
    help,/memory
    print,''

    for i=0,num_beams-1 do begin
	print,'Reading beam data'
	get_beam_data,i,chord,beam
	print,'done.'
	toc
	help,/memory
	print,''

	print,'Calculating geometric angles and Lorentz transform'
	calc_geom_angles,chord,field,beam,angle,lorentzE_norm,dangle_fg=dangle_fg,dangle_ap=dangle_ap,dangle_sp=dangle_sp
	print,'done.'
	toc
	help,/memory
	print,''

	;faraday,chord,field,angle,profiles ;Faraday rotation. Really quite negligible
	if keyword_set(enable_bes) then begin
	    print,'Calculating MSE spectrum'
	    bes_spec[*,i] = synth_bes(chord,field,beam,angle,lorentzE_norm,dangle_fg,dangle_ap,dangle_sp,profiles,bes_data) ;ph/m^2/s/sr/A
	    print,'done.'
	    toc
	    help,/memory
	    print,''
	endif else bes_data = 0

	if keyword_set(enable_cxrs) then begin
	    for j=0,n_elements(enable_cxrs)-1 do begin
		print,'Calculating CXRS spectrum',enable_cxrs[j]
		cxrs_spec[*,i] = synth_cxrs(chord,field,beam,angle,profiles,cxrs_data,enable_cxrs[j]) ;ph/m^2/s/sr/A
		print,'done.'
		toc
		help,/memory
		print,''
	    endfor
	endif else cxrs_data = 0
    endfor

    if n_elements(detstokes) eq 0 then detstokes = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]
    if keyword_set(enable_bes) then begin
	bes_spec_filt = detstokes[0,0]*bes_data.stokes_s0 + detstokes[1,0]*bes_data.stokes_s1 + detstokes[2,0]*bes_data.stokes_s2 + detstokes[3,0]*bes_data.stokes_s3
	bes_counts = bes_spec_filt * scalefactor
    endif else bes_counts = 0

    if keyword_set(enable_cxrs) then begin
	cxrs_spec_filt = detstokes[0,0]*cxrs_data.stokes_s0 + detstokes[1,0]*cxrs_data.stokes_s1 + detstokes[2,0]*cxrs_data.stokes_s2 + detstokes[3,0]*cxrs_data.stokes_s3
	cxrs_counts = cxrs_spec_filt * scalefactor
    endif else cxrs_counts = 0

    edge_spec = dblarr(n_elements(wve))
    if keyword_set(enable_edge) then begin
	for j=0,n_elements(enable_edge)-1 do begin
	    print,'Calculating edge spectrum for',enable_edge[j]
	    edge_spec += synth_edge(chord,field,beam,angle,profiles,cxrs_data,enable_edge[j]) ;ph/m^2/s/sr/A
	    print,'done.'
	    toc
	    help,/memory
	    print,''
	endfor
	edge_spec *= detstokes[0,0]
	edge_counts = edge_spec * scalefactor
    endif else begin
	edge_spec = 0
	edge_counts = 0
    endelse

    if keyword_set(enable_brem) then begin
	print,'Calculating bremsstrahlung spectrum'
	brem_spec = synth_bremsstrahlung(chord,profiles) ;ph/m^2/s/sr/A
	brem_spec *= detstokes[0,0]
	brem_counts = brem_spec * scalefactor
	print,'done.'
	toc
	help,/memory
	print,''
    endif else begin
	brem_spec = 0
	brem_counts = 0
    endelse

    print,'Calculating detector noise'
    if num_beams gt 0 then begin
	if beam.n_comp gt 1 then pure = total(bes_counts+cxrs_counts,2)+edge_counts+brem_counts $
	  else pure = bes_counts+cxrs_counts+edge_counts+brem_counts
    endif else begin
	pure = edge_counts+brem_counts
	bes_data = 0
	cxrs_data = 0
    endelse
    pure_bk = (edge_counts+brem_counts)
    noisy = dblarr(n_elements(wve),n_gen)
    noisy_bk = dblarr(n_elements(wve),n_gen)
    for i=0,n_elements(wve)-1 do begin
	noisy[i,*] = (randomu(seed,n_gen,poisson=pure[i]/detgain+darknoise,/double)-darknoise)*detgain
	noisy_bk[i,*] = (randomu(seed,n_gen,poisson=pure_bk[i]/detgain+darknoise,/double)-darknoise)*detgain
    endfor
    print,'done.'
    toc
    help,/memory
    print,''

    result = {wve:wve,noisy:noisy,noisy_bk:noisy_bk,pure:pure,pure_bk:pure_bk,sens:sens,scalefactor:scalefactor,$
                 bes_spec:bes_spec,cxrs_spec:cxrs_spec,edge_spec:edge_spec,brem_spec:brem_spec,$
                 bes_counts:bes_counts,cxrs_counts:cxrs_counts,edge_counts:edge_counts,brem_counts:brem_counts,bes_data:bes_data,cxrs_data:cxrs_data}
    parameters = get_parameters()
    data = {chord:chord,field:field,angle:angle,profiles:profiles,lorentzE_norm:lorentzE_norm,dangle_fg:dangle_fg,dangle_ap:dangle_ap,dangle_sp:dangle_sp}
    pathfilename = getenv('NBASS_RUNS_PATH')+'/'+run+'/'+filename
    save,parameters,data,result,filename=pathfilename,description='beam_spec generated synthetic spectra'
    print,'NBASS calculation complete.'
    print,'Results have been saved to '+pathfilename
    print,'type "plot_results,'''+run+''','''+filename+'''" to plot results'
end

pro plot_results,run,filename,_REF_EXTRA=extra
    if n_elements(filename) eq 0 then filename=run+'.sav'
    restore,getenv('NBASS_RUNS_PATH')+'/'+run+'/'+filename

    wve = result.wve
    bes_spec = result.bes_spec
    cxrs_spec = result.cxrs_spec
    edge_spec = result.edge_spec
    brem_spec = result.brem_spec
    bes_counts = result.bes_counts
    cxrs_counts = result.cxrs_counts
    edge_counts = result.edge_counts
    brem_counts = result.brem_counts
    noisy = result.noisy

    figure,1
    clf,/all
    yrange = max([bes_spec,cxrs_spec,edge_spec,brem_spec])
    ploti,charsize=1.5,xtitle='wavelength (A)',ytitle='emissivity ph/m^2/s/sr/A',ptitle='Synthetic Spectral Intensity',xrange=[min(wve),max(wve)],yrange=[0,yrange]
    if keyword_set(parameters.param.enable_bes) then oploti,wve,bes_spec,icolor=0
    if keyword_set(parameters.param.enable_cxrs) then oploti,wve,cxrs_spec,icolor=1
    if keyword_set(parameters.param.enable_edge) then oploti,wve,edge_spec,icolor=2
    if keyword_set(parameters.param.enable_brem) then oploti,wve,brem_spec,icolor=3
    enables = where([parameters.param.enable_bes,keyword_set(parameters.param.enable_cxrs),keyword_set(parameters.param.enable_edge),parameters.param.enable_brem])
    legendi,icolor=([0,1,2,3])[enables],label=(['BES','CXRS','edge','bremsstrahlung'])[enables],charsize=1.6
    ;plot_impurities

    figure,2
    clf,/all
    yrange = max([max(bes_counts),max(cxrs_counts),max(edge_counts),max(brem_counts)])
    ploti,charsize=1.5,xtitle='wavelength (A)',ytitle='detector counts',ptitle='Synthetic Detector Signal',xrange=[min(wve),max(wve)],yrange=[0,yrange]
    if keyword_set(parameters.param.enable_bes) then oploti,wve,bes_counts,icolor=0
    if keyword_set(parameters.param.enable_cxrs) then oploti,wve,cxrs_counts,icolor=1
    if keyword_set(parameters.param.enable_edge) then oploti,wve,edge_counts,icolor=2
    if keyword_set(parameters.param.enable_brem) then oploti,wve,brem_counts,icolor=3
    oploti,wve,reform(noisy[*,0]),icolor=4
    enables = [enables,4]
	;plot_impurities
    legendi,icolor=([0,1,2,3,4])[enables],label=(['BES','CXRS','edge','bremsstrahlung','total+noise'])[enables],charsize=1.6

    if size(result.bes_data,/type) eq 8 then begin
	;total emission stokes
	if n_elements(size(result.bes_data.stokes_s0,/dim)) gt 1 then begin
	    stokes_s0 = total(result.bes_data.stokes_s0,2) + total(result.cxrs_data.stokes_s0,2) + edge_spec + brem_spec
	    stokes_s1 = total(result.bes_data.stokes_s1,2) + total(result.cxrs_data.stokes_s1,2)
	    stokes_s2 = total(result.bes_data.stokes_s2,2) + total(result.cxrs_data.stokes_s2,2)
	    stokes_s3 = total(result.bes_data.stokes_s3,2) + total(result.cxrs_data.stokes_s3,2)
	endif else begin
	    stokes_s0 = result.bes_data.stokes_s0 + result.cxrs_data.stokes_s0 + edge_spec + brem_spec
	    stokes_s1 = result.bes_data.stokes_s1 + result.cxrs_data.stokes_s1
	    stokes_s2 = result.bes_data.stokes_s2 + result.cxrs_data.stokes_s2
	    stokes_s3 = result.bes_data.stokes_s3 + result.cxrs_data.stokes_s3
	endelse
	poldegree = sqrt(stokes_s1^2 + stokes_s2^2 + stokes_s3^2)/stokes_s0
	linpoldegree = sqrt(stokes_s1^2 + stokes_s2^2)/stokes_s0
	spsi = atan(stokes_s2,stokes_s1)/2
	schi = atan(stokes_s3,sqrt(stokes_s1^2+stokes_s2^2))/2

	figure,3
	clf,/all
	yrange = [min([stokes_s1,stokes_s2,stokes_s3]),max(stokes_s0)]
	ploti,wve,stokes_s0,charsize=1.5,ptitle='Stokes parameters',isection=0,xrange=[6550,6620],yrange=yrange,icolor=0,xtitle='wavelength (A)',ytitle='ph/m^2/s/sr/A'
	oploti,wve,stokes_s1,icolor=1,isection=0
	oploti,wve,stokes_s2,icolor=2,isection=0
	oploti,wve,stokes_s3,icolor=3,isection=0
	legendi,label=['Intensity','Linear +','Linear x','Circular'],charsize=1.6,icolor=[0,1,2,3]

	figure,4
	clf,/all
	ploti,wve,poldegree,charsize=1.5,ptitle='Degree of polarization',icolor=0,isection=1,xrange=[6550,6620],yrange=[0,1.2],xtitle='wavelength (A)'
	oploti,wve,linpoldegree,icolor=2,isection=1
	legendi,icolor=[0,2],label=['degree of polarization','deg of linear pol'],charsize=1.6

	figure,5
	clf,/all
	ploti,wve,spsi,charsize=1.5,ptitle='Polarization angle',ytitle='radians',isection=2,xrange=[6550,6620],yrange=[-!dpi/4,!dpi/4],icolor=0,xtitle='wavelength (A)'
	oploti,wve,spsi+!dpi/2,icolor=2,isection=2
	legendi,icolor=[0,2],label=['Polarization angle','Pol angle + \pi/2'],charsize=1.6
    endif
end
