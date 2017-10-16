;calculate the outer product between two arrays.
;by Ken Liao
;
;Input:
;  A and B are two multidimensional numerical arrays.
;Output:
;  result has dimensions of n_dimension(A)+n_dimension(B).
;
;Mathematically, if A = A[i,j], B = B[k,l]
;result[i,j,k,l] = A[i,j]*B[k,l]
;in IDL use
;result = outer(A,B)
;example: if A and B are vectors, result = dyadic AB

;dependencies
@lib/cmproduct.pro

function outer,A,B
    sa = size(A)
    if sa[0] eq 0 then sa = [1,1] ;fix for scalar input
    a_flat = reform(A,cmproduct(sa[1:sa[0]])) ;flatten
    sb = size(B)
    if sb[0] eq 0 then sb = [1,1] ;fix for scalar input
    b_flat = reform(B,cmproduct(sb[1:sb[0]])) ;flatten
    result = a_flat # b_flat
    result = reform(result,[sa[1:sa[0]],[sb[1:sb[0]]]],/overwrite) ;unflatten

    return,result
end

;test code by comparing result to for loop implementation
pro test
    ;a = indgen(4,4,4)
    ;b = indgen(4,4,4)+1
    a = randomu(seed,4,4,4)
    b = randomu(seed,4,4,4)

    starttime=systime(1)
    c = outer(a,b)
    outertime=systime(1)-starttime
    starttime=systime(1)
    d = fltarr(4,4,4,4,4,4)
    for i1=0,3 do begin
    for i2=0,3 do begin
    for i3=0,3 do begin
    for i4=0,3 do begin
    for i5=0,3 do begin
    for i6=0,3 do begin
	d[i1,i2,i3,i4,i5,i6] += a[i1,i2,i3]*b[i4,i5,i6]
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
