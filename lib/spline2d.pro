;This library provides 2d spline interpolation

;spl_interp_der is a replacement for spl_interp which also calculates the first derivative of the interpolating function
function spl_interp_der,x,y,dd,x2,der
    ind = value_locate(x,x2)
    ind = ((ind > 0) < (n_elements(x)-2))
    diff = x[ind+1] - x[ind]
    t = (x2 - x[ind])/diff
    b = (dd[ind]-dd[ind+1])/6d*diff^2
    a = -dd[ind]/2d*diff^2 + b
    result = y[ind]*(1-t) + y[ind+1]*t + t*(1-t)*(a+b*t)
    der = (y[ind+1]-y[ind]+a-2*(a-b)*t-3*b*t^2)/diff
    return,result
end

;spline2d performs a 2d cubic spline interpolation and optionally calculates derivatives
function spline2d,z,x,y,x2,y2,dx,dy
    nx = n_elements(x)
    ny = n_elements(y)
    nr = n_elements(x2)
    der = dblarr(nx,ny)
    for i=0,ny-1 do begin
	der[*,i] = spl_init(x,z[*,i],/double)
    end

    zi = dblarr(ny,nr)
    if arg_present(dx) then begin
	dzi = dblarr(ny,nr)
	for i=0,ny-1 do begin
	    zi[i,*] = spl_interp_der(x,z[*,i],der[*,i],x2,dtmp)
	    dzi[i,*] = dtmp
	endfor
    endif else begin
	for i=0,ny-1 do begin
	    zi[i,*] = spl_interp(x,z[*,i],der[*,i],x2,/double)
	endfor
    endelse

    result = dblarr(nr)
    if arg_present(dy) then dy = dblarr(nr)
    if arg_present(dx) then dx = dblarr(nr)
    for j=0,nr-1 do begin
	der2 = spl_init(y,zi[*,j])
	if arg_present(dy) then begin
	    result[j] = spl_interp_der(y,zi[*,j],der2,y2[j],ddy)
	    dy[j] = ddy
	endif else begin
	    result[j] = spl_interp(y,zi[*,j],der2,y2[j])
	endelse
	if arg_present(dx) then begin
	    derx = spl_init(y,dzi[*,j])
	    dx[j] = spl_interp(y,dzi[*,j],derx,y2[j])
	endif
    endfor

    return,result
end
