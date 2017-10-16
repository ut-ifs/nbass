;@lib/splinex.pro
@lib/spline2d.pro

;extract the B field from efit data
pro efit_psi2b,r,z,t1=t1,t2=t2,br=br,bz=bz,bt=bt,bmod=bmod,rgrid=rgrid,zgrid=zgrid,tgrid=tgrid,psirz=psirz,fpol=fpol,ssibry=ssibry,ssimag=ssimag
    nw = n_elements(rgrid)
    nh = n_elements(zgrid)

    if n_elements(t1) gt 0 then indx1 = mindex(tgrid,t1) else indx1 = indgen(n_elements(tgrid))
    if n_elements(t2) gt 0 then indx2 = mindex(tgrid,t2) else indx2 = indx1

    ntt = n_elements(indx1)
    np = n_elements(r)
    br = fltarr(np,ntt)
    bz = fltarr(np,ntt)
    bt = fltarr(np,ntt)
    xxx = findgen(nw)/(nw-1)
    ff = fltarr(nw,nh)
    psix = fltarr(nw,nh)
    psiz = fltarr(nw,nh)

    for itt=0,n_elements(indx1)-1 do begin
	nt = indx2[itt] - indx1[itt] + 1 ;number of efit time points to average
	for it=indx1[itt],indx2[itt] do begin
	    psin = (psirz[*,*,it]-ssimag[it])/(ssibry[it]-ssimag[it]) ;rho to r,z grid
	    inside = where(psin lt 1,nin,comp=outside,ncomp=nout)
	    ff[inside] = interpol(reform(fpol[*,it]),xxx,psin[inside],/spl) ;interpolate from rho grid to r,z grid
	    ff[outside]=fpol[nw-1,it]

	    ;sets2d,float(psirz[*,*,it]),c,float(rgrid),bkx,float(zgrid),bky,/IDL
	    ;seva2d,bkx,bky,c,r,z,pds,2
            ;
	    ;sets2d,float(ff),c,float(rgrid),bkx,float(zgrid),bky,/IDL
	    ;seva2d,bkx,bky,c,r,z,ff2,1
	    ;
	    ;if n_elements(indx1) eq 1 then begin
	    ;    br += -pds[2,*]/r/nt
	    ;    bz += pds[1,*]/r/nt
	    ;    bt += ff2[0,*]/r/nt
	    ;endif else begin
	    ;    br[0,itt] += reform(-pds[2,*]/r/nt)
	    ;    bz[0,itt] += reform(pds[1,*]/r/nt)
	    ;    bt[0,itt] += reform(ff2[0,*]/r/nt)
	    ;endelse
	    psiint = spline2d(psirz[*,*,it],rgrid,zgrid,r,z,dpsidr,dpsidz)
	    ffint = spline2d(ff,rgrid,zgrid,r,z)
	    if n_elements(indx1) eq 1 then begin
	        br += -dpsidz/r/nt
	        bz += dpsidr/r/nt
	        bt += ffint/r/nt
	    endif else begin
	        br[0,itt] += reform(-dpsidz/r/nt)
	        bz[0,itt] += reform(dpsidr/r/nt)
	        bt[0,itt] += reform(ffint/r/nt)
	    endelse
	endfor
    endfor
    bmod = sqrt(br^2+bz^2+bt^2)
end
