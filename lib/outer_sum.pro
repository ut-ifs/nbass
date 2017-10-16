;calculate the outer sum between two arrays.
;by Ken Liao
;
;Input:
;  A and B are two multidimensional numerical arrays.
;Output:
;  result has dimensions of n_dimension(A)+n_dimension(B), with type of A+B
;
;Mathematically, if A = A[i,j], B = B[k,l]
;result[i,j,k,l] = A[i,j]+B[k,l]
;in IDL use
;result = outer(A,B)

;dependencies
@cmproduct.pro

function outer_sum,A,B
    sa = size(A)
    a_flat = reform(A,cmproduct(sa[1:sa[0]])) ;flatten
    sb = size(B)
    b_flat = reform(B,cmproduct(sb[1:sb[0]])) ;flatten
    result = a_flat # make_array(n_elements(b_flat),value=1) + make_array(n_elements(a_flat),value=1) # b_flat
    result = reform(result,[sa[1:sa[0]],[sb[1:sb[0]]]],/overwrite) ;unflatten
    return,result
end

pro test
    ;a = indgen(4,4,4)
    ;b = indgen(4,4,4)+1
    a = randomu(seed,4,4,4)
    b = randomu(seed,4,4,4)

    starttime=systime(1)
    c = outer_sum(a,b)
    outertime=systime(1)-starttime
    starttime=systime(1)
    d = fltarr(4,4,4,4,4,4)
    for i1=0,3 do begin
    for i2=0,3 do begin
    for i3=0,3 do begin
    for i4=0,3 do begin
    for i5=0,3 do begin
    for i6=0,3 do begin
	d[i1,i2,i3,i4,i5,i6] += a[i1,i2,i3]+b[i4,i5,i6]
    endfor
    endfor
    endfor
    endfor
    endfor
    endfor
    testtime=systime(1)-starttime
    print, 'total error',total(abs(d - c))
    print, 'time function',outertime
    print, 'time test',testtime
end
