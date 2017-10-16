;Opens an EFIT G file and returns selected values in a structure
;Either the filename, or the shot and time (in ms) are needed as input
;The shot, time, and directory are only used to obtain the filename if it isn't given
;The return value is
;    gdata = {fpol:fpol,pres:pres,ffprim:ffprim,pprime:pprime,psirz:psirz,qpsi:qpsi,rgrid:rgrid,zgrid:zgrid,xlim:xlim,ylim:ylim,rmaxis:rmaxis,ssimag:ssimag,ssibdry:ssibdry,cpasma:cpasma,rbbbs:rbbbs,zbbbs:zbbbs}
;   fpol is the 
function read_geqdsk,shot,time,filename=filename,dir=dir,quiet=quiet
    if n_elements(filename) eq 0 then begin
	if n_elements(dir) eq 0 then dir = '.'
	if shot le 999999 then forshot='(i6.6,".",i5.5)' else forshot='(i10.10,".",i5.5)'
	filhd = dir+'/g'+string(shot,time,form=forshot)
    endif

    ;load geqdsk
    openr,lun,filename,/get_lun
    form2000 = '(6a8,3i4)'
    casee = strarr(6)
    idum=0l
    mw=0l ;# grid in R
    mh=0l ;# grid in z
    readf,lun,form=form2000,casee,idum,mw,mh
    if ~keyword_set(quiet) then begin
	print,casee
	print,idum,mw,mh
    endif
    form2020 = '(5d16.9)'
    readf,lun,form=form2020,xdim,zdim,rzero,rgrid1,zmid
    readf,lun,form=form2020,rmaxis,zmaxis,ssimag,ssibdry,bcentr
    readf,lun,form=form2020,cpasma,ssimag,xdum,rmaxis,xdum
    readf,lun,form=form2020,zmaxis,xdum,ssibry,xdum,xdum
    fpol = dblarr(mw)
    pres = dblarr(mw)
    readf,lun,form=form2020,fpol
    readf,lun,form=form2020,pres
    ffprim = dblarr(mw)
    readf,lun,form=form2020,ffprim
    pprime = dblarr(mw)
    readf,lun,form=form2020,pprime
    psirz = dblarr(mw,mh)
    readf,lun,form=form2020,psirz
    qpsi = dblarr(mw)
    readf,lun,form=form2020,qpsi

    nbbbs = 0l
    limitr = 0l
    readf,lun,nbbbs,limitr
    if (nbbbs ne 0) then begin
	rbbbs = dblarr(nbbbs)
	zbbbs = dblarr(nbbbs)
	rzbbbs = dblarr(2,nbbbs)
	readf,lun,form=form2020,rzbbbs
	rbbbs = rzbbbs[0,*]
	rbbbs = transpose(rbbbs)
	zbbbs = transpose(rzbbbs[1,*])
    endif
    xylim = dblarr(2,limitr)
    readf,lun,form=form2020,xylim
    xlim = transpose(xylim[0,*])
    ylim = transpose(xylim[1,*])
    rgrid = rgrid1 + xdim*dindgen(mw)/(mw-1)
    zgrid = zmid-0.5*zdim + zdim*dindgen(mh)/(mh-1)

    close,lun
    gdata = {fpol:fpol,pres:pres,ffprim:ffprim,pprime:pprime,psirz:psirz,qpsi:qpsi,rgrid:rgrid,zgrid:zgrid,xlim:xlim,ylim:ylim,rmaxis:rmaxis,ssimag:ssimag,ssibdry:ssibdry,cpasma:cpasma,rbbbs:rbbbs,zbbbs:zbbbs}
    return,gdata
end
