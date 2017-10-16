function mindex, array, values, dist=dist

dist = fltarr(n_elements(values))

out = lonarr(n_elements(values))
for i = 0, n_elements(out)-1 do begin
    mindum = min(abs(array - values[i]), index, /nan)
    out[i] = index
dist[i] = values[i]-array[out[i]]
endfor

if n_elements(out) eq 1 then out = out[0]
if n_elements(out) eq 1 then dist = dist[0]

return, out
end

 
