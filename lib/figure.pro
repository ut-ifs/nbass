;This file contains routines which barely emulate the behavior of
;the far superior "Figure" library for IDL by Syun'ichi Shiraiwa
;and Naoto Tsujii, which is not included. Figure is a set of
;plotting routines which replace the clunky IDL plot routines.
;This file merely wraps the IDL plot routine to get something that runs.

pro ploti,x,y,ptitle=ptitle,icolor=icolor,charsize=charsize,_extra=ex
    if n_elements(x) eq 0 then x=[0]
    if n_elements(icolor) gt 0 then begin
	colarr = ['ff0000'xl, '7f00'xl, 'ff'xl, 'bfbf00'xl, 'bf00bf'xl, 'bfbf'xl, '3f3f3f'xl, 'ffffff'xl, 0l]
	color = colarr[icolor mod 9]
    endif
    !p.background='ffffff'xl
    if n_elements(charsize) gt 0 then !p.charsize=charsize
    sz = size(y)
    case sz[0] of
	0: begin
	    plot,x,_extra=ex,title=ptitle,color=0
	    oplot,x,_extra=ex,color=color
	    end
	1: begin
	    plot,x,y,_extra=ex,title=ptitle,color=0
	    oplot,x,y,_extra=ex,color=color
	    end
	2: begin
	    plot,x,y,_extra=ex,title=ptitle,color=0
	    for i=0,sz[2]-1 do begin
		oplot,x,y[*,i],_extra=ex,color=color
	    endfor
	    end
    endcase
end

pro oploti,x,y,icolor=icolor,charsize=charsize,_extra=ex
    if n_elements(icolor) gt 0 then begin
	colarr = ['ff0000'xl, '7f00'xl, 'ff'xl, 'bfbf00'xl, 'bf00bf'xl, 'bfbf'xl, '3f3f3f'xl, 'ffffff'xl, 0l]
	color = colarr[icolor mod 9]
    endif
    !p.background='ffffff'xl
    if n_elements(charsize) gt 0 then !p.charsize=charsize
    sz = size(y)
    case sz[0] of
	0: begin
	    oplot,x,_extra=ex,color=color
	    end
	1: begin
	    oplot,x,y,_extra=ex,color=color
	    end
	2: begin
	    for i=0,sz[2]-1 do begin
		oplot,x,y[*,i],_extra=ex,color=color
	    endfor
	    end
    endcase
end

pro clf,all=all
    erase
end

pro legendi,label=label,icolor=icolor,charsize=charsize
    colarr = ['ff0000'xl, '7f00'xl, 'ff'xl, 'bfbf00'xl, 'bf00bf'xl, 'bfbf'xl, '3f3f3f'xl, 'ffffff'xl, 0l]
    color = colarr[icolor mod 9]
    for i=0,n_elements(label)-1 do begin
	xyouts,0.18,0.86-i*0.05,'-'+label[i],charsize=charsize,color=color[i],/normal
    endfor
end

pro figure,id
    device, true_color = 24
    device, decomposed = 1
    device, retain = 2
    window,id,retain=2,xsize=650,ysize=500
end
