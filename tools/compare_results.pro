@lib/figure.pro
pro compare_results,filenames
    figure,1
    clf,/all
    for i=0,n_elements(filenames)-1 do begin
	restore,getenv('NBASS_RUNS_PATH')+'/'+filenames[i]
	if i eq 0 then ploti,result.wve,result.pure,icolor=0,charsize=1.6,xtitle='wavelength',ytitle='counts' else oploti,result.wve,result.pure,icolor=i
    endfor
    legendi,label=filenames,charsize=1.6,icolor=indgen(n_elements(filenames))
end
