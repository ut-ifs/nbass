pro recurse_print,lun,struct,name
    printf,lun,name
    typ = size(struct,/type)
    case typ of
	8: begin ;structure
	    tags = tag_names(struct)
	    for i=0,n_elements(tags)-1 do begin
		recurse_print,lun,struct.(i),name+'.'+tags[i]
	    endfor
	end
	else: begin
	    printf,lun,struct
	end
    endcase
end

pro export_ascii,run,filename,outfile
    if n_elements(filename) eq 0 then filename=run+'.sav'
    restore,getenv('NBASS_RUNS_PATH')+'/'+run+'/'+filename

    openw,lun,outfile,/get_lun

    recurse_print,lun,parameters,'parameters'
    recurse_print,lun,data,'data'
    recurse_print,lun,result,'result'

    free_lun,lun
    info = file_info(outfile)
    print,strtrim(info.size,2)+' bytes written to '+outfile
end
