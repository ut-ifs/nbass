;convert a number to a string, rounding at a chosen decimal place
;if decimalplaces is negative, round at a place before the radix point
;Ken Liao 10/25/2012
function num2str,value,decimalplaces
    if n_elements(decimalplaces) eq 0 then decimalplaces = 3
    if decimalplaces eq 0 then return,strtrim(round(value,/l64),2)
    if decimalplaces gt 0 then begin
	power10 = 10LL^decimalplaces
	intval = round(value*power10,/l64)
	intpart = intval/power10
	fracpart = intval - intpart*power10
	;if intpart eq 0 and fracpart lt 0 then return,'-0.'+string(abs(fracpart),format='(I0'+strtrim(decimalplaces,2)+')')
	;return,strtrim(intpart,2)+'.'+string(abs(fracpart),format='(I0'+strtrim(decimalplaces,2)+')')
	;Fix for vector input
	result = strtrim(intpart,2)+'.'+string(abs(fracpart),format='(I0'+strtrim(decimalplaces,2)+')')
	check = where(intpart eq 0 and fracpart lt 0)
	if check[0] ne -1 then result[check] = '-'+result[check]
	return,result
    endif else begin
	power10 = 10LL^(-decimalplaces)
	intval = round(value/power10,/l64)
	return,strtrim(intval,2)+string(0,format='(I0'+strtrim(-decimalplaces,2)+')')
    endelse
end
