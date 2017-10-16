pro export_ncdf,run,filename,outfile
    if n_elements(filename) eq 0 then filename=run+'.sav'
    restore,getenv('NBASS_RUNS_PATH')+'/'+run+'/'+filename

;    id = ncdf_create(outfile)

;    wvedim = ncdf_dimdef(id,'wve',n_elements(result.wve))
;    gendim = ncdf_dimdef(id,'generation',parameters.param.n_gen)
;    sz = size(result.bes_counts,/dim)
;    energydim = ncdf_dimdef(id,'energycomp',sz[1])
;    np = data.chord.npoint
;    griddim = ncdf_dimdef(id,'gridpoint',np)
;
;    ncdf_varput,id,vid,
;    ncdf_close,id

    ;to hell with this
    data2 = {parameters:parameters,data:data,result:result}

    write_netCDF, data2, outfile, status, /clobber
end
