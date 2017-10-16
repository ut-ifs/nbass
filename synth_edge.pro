;This is essentially a placeholder for an edge emission calculation
;result in ph/(m^2*s*sr)
;input:
;  wve:      A      wavelength
;return:
;  spectrum  ph/(A m^2 s sr)
function synth_edge,chord,field,beam,angle,profiles,cxrs_data,species
    @param.idl
    @view_param.idl
    @detector_param.idl

    ;Data from von Hellermann PPCF 1995, with some tweaked wavelengths
    ; coef = [wavelength(A), Q_0(10^-15 m^3 s^-1), E_m(keV/amu), p', q']
    case species of
	'H': wvel = [4861.3298d,6562.7952d]
	'D': wvel = [4860.007d,6561.010d]
	'He3': wvel = [4685.9150d]
	'He4': wvel = [4685.7048d]
	'B': wvel = [4940.376d, 4944.6d]
	'C': wvel = [5290.5d,6578.05d,6582.88]
	else: message,species+' cxrs cross section not implemented'
    endcase

    case species of ;scaling in cm^3 s^-1, multiply by ne and ni
	;placeholder values... perhaps add C-R model results
	'H': scaling = [1d-5,1d-5]
	'D': scaling = [1d-5,1d-5]
	'He3': scaling = [1d-5]
	'He4': scaling = [1d-5]
	'B': scaling = [0,1d-5]
	'C': scaling = [0,1d-5,1d-5]
	else: message,species+' cxrs cross section not implemented'
    endcase

    case species of
	'H' : mc2kev = 938262.0813d ;mass of proton * c^2 in keV (2014 CODATA recommended values)
	'D' : mc2kev = 1865612.928d ;mass of deuteron * c^2 in keV (2014 CODATA recommended values)
	'He3': mc2kev = 2808920.906d ;mass of He3 nucleus * c^2 in keV
	'He4': mc2kev = 3727379.24082d ;mass of He4 nucleus * c^2 in keV
	'B'  : mc2kev = 10.81d * 931494.0954d
	'C'  : mc2kev = 12.01d * 931494.0954d
	else: message,species+' ion mass not implemented'
    endcase

    ind = (where(profiles.ions eq species))[0]
    rhind = mindex(chord.rho,1)
    ne_edge = profiles.n_e[rhind] ;fixme
    nz_edge = profiles.n_ion[rhind,ind] ;fixme
    ;Calculate zeeman spectrum
;    if keyword_set(calczeeman) then begin
;	print,'zeeman calculation'
;	normalize_zeeman,species,wvezeeman,ampzeeman,field.b_tot[fe],angle.psi[fe]
;	n_zeeman = (size(wvezeeman,/dimensions))[0]
;    endif

    ;A good model for the edge isn't implemented.
    ;Instead, this basic model is present as a filler.
    ;We can improve this model using collisional radiative to find edge excited population, but this requires many assumptions
    t_edge = 13.6d*2d/3d        ; eV       just a rough guess of the edge temperature
    w_edge = 1d             ; cm      just a rough guess of the edge "width", used to scale the density. Replace with better model

    pixboundl = [wve[0],(wve[0:n_elements(wve)-1]+wve[1:*])/2]
    pixboundr = [(wve[0:n_elements(wve)-1]+wve[1:*])/2,wve[n_elements(wve)-1]]

    d_edge = dblarr(n_elements(wve))
    for i=0,n_elements(wvel)-1 do begin
	width = sqrt((2*t_edge*wvel[i]^2/mc2kev/1000)^2 + l_instr^2)
	amplitude = ne_edge*nz_edge*scaling[i]*w_edge*100/4/!dpi ;ph/(m^2 s sr)
	d_edge += 0.5d*amplitude*(erf((pixboundl - wvel[i])/width) - erf((pixboundr - wvel[i])/width))
    endfor
    d_edge /= (pixboundl-pixboundr)
    return,d_edge
end
