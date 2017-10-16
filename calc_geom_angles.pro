;@/home/liao/idl_func/vector.pro

;a has dimensions [n,3]; b has dimensions [n,3]; axb has dimensions [n,3]
function crossp_2d2d,a,b
    return,[[a[*,1]*b[*,2]-a[*,2]*b[*,1]],[a[*,2]*b[*,0]-a[*,0]*b[*,2]],[a[*,0]*b[*,1]-a[*,1]*b[*,0]]]
end

function crossp_1d2d,a,b
    return,[[a[1]*b[*,2]-a[2]*b[*,1]],[a[2]*b[*,0]-a[0]*b[*,2]],[a[0]*b[*,1]-a[1]*b[*,0]]]
end

pro calc_geom_angles,chord,field,beam,angle,lorentzE_norm,dangle_fg=dangle_fg,dangle_ap=dangle_ap,dangle_sp=dangle_sp
    @view_param.idl
    @mesh_param.idl
    @equil_param.idl
;inputs
;  grid_rad:       grid radius (m)
;  grid_pos:       grid center in cartesian coordinates (m)
;  wall_pos:       beam line position in cartesian coordinates (m)
;outputs
;  angle (structure)
;    alpha: the angle between the beam and the viewing chord
;           the beam ray points into the plasma, and the viewing chord ray points into the dissector
;    psi:   the angle between the B field and the viewing chord
;    phi:   the angle between projection of viewing chord on to plane perpendicular to B,
;           and projection of beam on to plane perpendicular to B
;    phiL:  pi/2 - (the angle between E_L+E_r and the projection of the viewing
;                   chord on the plane perpendicular to B).
;    theta: the angle between the magnetic field and beam
;    theta1:the angle between the toroidal magnetic field and beam
;    gamma: the angle between the lorentz E field and viewing chord.
;           When E_lab=0, |cos(gamma)| = sin(psi)*sin(phi)
;    pitch: pitch angle, tan(pitch) = B_pol/B_tor
;  b_tot: magnetic field amplitude (T)
;  chord: structure containing chord step coordinates
;  dangle_fg, dangle_ap, dangle_sp: structures with spread in angles due to finitegrid, finiteaper, and finitespot
;         These are not combined into a single error because they have semicircular lineshapes

    ;define beam velocity direction using ray from center of grid to chord point
    ;todo:incorporate velocity distribution numerically
    ;This handles beam divergence but neglects finite beam grid size
    beam_vect = [[chord.x],[chord.y],[chord.z]] - rebin(transpose(beam.grid_pos),chord.npoint,3)
    beam_vect_norm = sqrt(total(beam_vect^2,2))
    beam_vect /= rebin(beam_vect_norm,chord.npoint,3) ; dim [npoint,3]: unit vectors
    ;FIXME TODO
;    beam_vect = extend([chord.npoint],norm_vect(beam.wall_pos - beam.grid_pos))

    ;calculate Lorentz E field, including motional and E_r contributions
    ;But this part does not include beam finite grid area effects
    c = 299792458d ;m/s
    b_vect = [[field.bx],[field.by],[field.bz]]   ;[npoint,3]  T
    vxB0 = crossp_2d2d(beam_vect,b_vect)*c ;[npoint,3]         (v/beta x B)
    vxB = outer(vxB0,beam.gamma*beam.beta) ;[npoint,3,n_comp]  gamma*(v x B)
    E_r_big = rebin([[field.Ex],[field.Ey],[field.Ez]],chord.npoint,3,beam.n_comp) ;[npoint,3,n_comp]
    beam_vect_big = rebin(beam_vect,chord.npoint,3,beam.n_comp)
    E_r_perp = E_r_big - beam_vect_big*rebin(reform(total(beam_vect_big*E_r_big,2),chord.npoint,1,beam.n_comp),chord.npoint,3,beam.n_comp)
    LorentzE = vxB + E_r_big + extend([chord.npoint,3],beam.gamma-1)*E_r_perp;[npoint,3,n_comp] V/m

    ;calculate angles between unit vectors A and B by taking acos(A dot B)
    alpha = acos(total(chord.vect*beam_vect,2)) ;[npoint]
    cos_psi = total(chord.vect*field.b_hat,2)        ;[npoint]
    psi = acos(cos_psi)                                   ;[npoint]
    cos_theta = total(beam_vect*field.b_hat,2)
    theta = acos(cos_theta)

    ;calculate two vectors perpendicular to the viewing direction
    vectH = crossp_1d2d([0,0,1],chord.vect)
    vectH /= rebin(sqrt(total(vectH^2,2)),chord.npoint,3)
    vectV = crossp_2d2d(chord.vect,vectH)
    vectV /= rebin(sqrt(total(vectV^2,2)),chord.npoint,3)

    ;calculate the Isler coordinate system used in starkzeeman.pro (z = b_vect, view in xz plane)
    xvect = chord.vect - rebin(cos_psi,chord.npoint,3)*field.b_hat
    xvect /= rebin(sqrt(total(xvect^2,2)),chord.npoint,3)
    yvect = crossp_2d2d(field.b_hat,xvect)
    yvect /= rebin(sqrt(total(yvect^2,2)),chord.npoint,3)

    ;transform vectH and vectV into the Isler coordinate system
    vectH = [[total(xvect*vectH,2)], [total(yvect*vectH,2)], [total(field.b_hat*vectH,2)]]
    vectV = [[total(xvect*vectV,2)], [total(yvect*vectV,2)], [total(field.b_hat*vectV,2)]]

    ;unit vector in purely toroidal direction with sign given by B field
    b_tor_vect = [[field.bx],[field.by],[dblarr(chord.npoint)]]
    b_tor_vect /= rebin(reform(sqrt(field.bx^2 + field.by^2),chord.npoint,1),chord.npoint,3)
    theta1 = acos(total(beam_vect*b_tor_vect,2))

    ;calculate phi, gamma
    beam_perp = beam_vect - field.b_hat*rebin(cos_theta,chord.npoint,3)
    beam_perp_norm = sqrt(total(beam_perp^2,2))
    beam_perp /= rebin(beam_perp_norm,chord.npoint,3) ; dim [npoint,3]: unit vectors
    view_perp = chord.vect - field.b_hat*rebin(cos_psi,chord.npoint,3)
    view_perp_norm = sqrt(total(view_perp^2,2))
    view_perp /= rebin(view_perp_norm,chord.npoint,3) ; dim [npoint,3]: unit vectors.  same thing as xvect.
    sign = sign(total(crossp_2d2d(chord.vect,beam_vect)*field.b_hat,2))
    phi = sign*acos(total(beam_perp*view_perp,2))
    gamma = acos(sin(psi)*sin(phi)) ;valid for E_lab=0

    ;calculate phiL
    LorentzE_norm = sqrt(total(LorentzE^2,2)) ;[npoint,n_comp]
    EL_hat = LorentzE/rebin(reform(LorentzE_norm,chord.npoint,1,beam.n_comp),chord.npoint,3,beam.n_comp) ;if Er != 0, different beam energy components will have different values here
    phiLcomp = acos(total(EL_hat*rebin(view_perp,chord.npoint,3,beam.n_comp),2)) ;dim [npoint,n_comp]
    phiL = dblarr(chord.npoint,beam.n_comp)
    for i=0,beam.n_comp-1 do begin
	sign = sign(total(crossp_2d2d(chord.vect,el_hat[*,*,i])*field.b_hat,2))
	phiL[*,i] = !dpi/2 + sign*phiLcomp[*,i]
    endfor

    pitch = atan(field.bz*sqrt(1+field.br^2/field.bz^2),abs(field.b_tor))   ; equation is rearranged so sign of bz determines the sign of pitch. sign of BT is ignored

    ;Calculate the effects of finite beam grid area using small angle approximations
    if (beam.grid_rad gt 0) then begin
	beam_axis = norm_vect(beam.wall_pos - beam.grid_pos)
	cos_axis_angle = beam_axis##beam_vect < 1d0  ;[nstep] dot product between beam_axis and beam_vect
	cos_axis_view = beam_axis##chord.vect < 1d0 ;[scalar] dot product between beam_axis and beam_vect
	cos_axis_b = beam_axis##b_vect < 1d0 ;[nstep] dot product between beam_axis and magnetic field
	;print,'Grid radius',beam.grid_rad
	delta = beam.grid_rad / beam_vect_norm
	;angle spread will have an approximately semicircular distribution (due to circular grid)
	dalpha = delta*(sqrt(1.0d0-cos_axis_view^2)-sqrt(1.0d0-cos_axis_angle^2>0))/sin(alpha)
	dtheta = delta*(sqrt(1.0d0-cos_axis_b^2)-sqrt(1.0d0-cos_axis_angle^2))/sin(theta)
	dphi = delta*(sqrt(1.0d0-total(chord.vect*beam_perp,2)^2)+sqrt(1.0d0-cos_axis_angle^2))/sin(phi)
	dgamma = dphi*cos(phi)*sin(psi)*sin(gamma)
	dangle_fg = {alpha:dalpha, psi:0, phi:dphi, theta:dtheta, gamma:dgamma}
    endif

    ;Calculate the effects of finite aperture using small angle approximations
    if (ap_rad gt 0) then begin
	;print,'Aperture radius',ap_rad
	view_dist = sqrt(total((optic_pos-spot_pos)^2))
	delap = ap_rad / view_dist
	dalpha = replicate(delap,chord.npoint)
	dpsi = replicate(delap,chord.npoint)
	dphi = delap*sqrt(1.0d0-(total(chord.vect*beam_perp,2) < 1d0)^2)/sin(phi)
	dgamma = (dphi*cos(phi)*sin(psi)+dpsi*cos(psi)*sin(phi))/sin(gamma) ;Not sure what happens to the lineshape here
	dangle_ap = {alpha:dalpha, psi:dpsi, phi:dphi, theta:0, gamma:dgamma}
    endif else dangle_ap = {alpha:0, psi:0, phi:0, theta:0, gamma:0}

    ;Calculate the effects of finite spot using small angle approximations
    if (sp_rad gt 0) then begin
	;print,'Spot focus radius',sp_rad
	view_dist = sqrt(total((optic_pos-spot_pos)^2))
	delsp = sp_rad / view_dist
	dalpha = replicate(delsp,chord.npoint)
	dpsi = replicate(delsp,chord.npoint)
	dphi = delsp*sqrt(1.0d0-(total(chord.vect*beam_perp,2) < 1d0)^2)/sin(phi)
	dgamma = (dphi*cos(phi)*sin(psi)+dpsi*cos(psi)*sin(phi))/sin(gamma) ;Not sure what happens to the lineshape here
	dangle_sp = {alpha:dalpha, psi:dpsi, phi:dphi, theta:0, gamma:dgamma}
    endif else dangle_sp = {alpha:0, psi:0, phi:0, theta:0, gamma:0}

    angle = {alpha:alpha, psi:psi, phi:phi, phiL:phiL, phiLcomp:phiLcomp, theta:theta, theta1:theta1, gamma:gamma, pitch:pitch, vectH:vectH, vectV:vectV, xvect:xvect, yvect:yvect}
end
