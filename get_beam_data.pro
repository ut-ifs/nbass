pro get_beam_data,beamid,chord,beam;,profiles
    @common_param.idl

    ;get ALCBEAM results
    print,'Reading ALCBEAM file '+alcbeam_file[beamid]
    read_alcbeam,alcbeam_file=alcbeam_file[beamid],alcbeam_input=alcbeam_input
    if alcbeam_input[0] eq ptr_new() then message,'invalid ALCBEAM file'
    alcbeam_input = *alcbeam_input[0]
    get_alcbeam, alcbeam_file=alcbeam_file[beamid],r_tor=chord.r,z_tor=chord.z,phi_tor=chord.phi,r0=chord.r[0],z0=chord.z[0],phi0=chord.phi[0],alcbeam_results=alcbeam_results

    E_beam = alcbeam_results.E ;keV
    if n_elements(energy_override) gt 0 then e_beam *= energy_override/E_beam[0]
    if min(E_beam) lt 0 then begin ;This indicates the presence of halo data
	n_halo = alcbeam_results[(where(E_beam lt 0))[0]].n_beam
	valid = where(E_beam ge 0)
	E_beam = E_beam[valid]
	alcbeam_results = alcbeam_results[valid]
    endif
    n_comp = n_elements(E_beam)
    beam_atom = alcbeam_input.beam_param.beam_atom
    case beam_atom of
	;'H': m_beam = 1.672621898e-24       ; g        proton mass
	;'D': m_beam = 3.343583719d-24       ; g        deuteron mass
	'H': mc2_beam = 938272.0813d          ; keV      proton rest energy
	'D': mc2_beam = 1875612.928d          ; keV      deuteron rest energy
    endcase
    E_amu = E_beam/mc2_beam*938272.0813d   ;approximate energy per nucleon
    gamma_beam = E_beam/mc2_beam + 1d0     ;[n_comp] relativistic gamma factor
    beta_beam = sqrt(1 - 1/gamma_beam^2)   ;[n_comp] v/c
    n_beam = alcbeam_results.n_beam ;1/cm^3 [npoint,n_comp]
    n_beam_n1 = n_beam*(1.0-alcbeam_results.exc_n2_frac-alcbeam_results.exc_n3_frac)
    n_beam_n2 = n_beam*alcbeam_results.exc_n2_frac
    n_beam_n3 = n_beam*alcbeam_results.exc_n3_frac

    grid_pos = cyl2cart([alcbeam_input.beam_geom.r_grid,alcbeam_input.beam_geom.phi_grid+alcbeam_input.beam_geom.beam_port_phi,alcbeam_input.beam_geom.z_grid])
    wall_pos = cyl2cart([alcbeam_input.beam_geom.r_wall,alcbeam_input.beam_geom.phi_wall+alcbeam_input.beam_geom.beam_port_phi,alcbeam_input.beam_geom.z_wall])
    ba = alcbeam_input.beam_geom.beam_apertures
    grid_rad = (max(ba[*,0]) - min(ba[*,0]))/2 + 0.005

    ;We call the holes in the acceleration grid "positive ion neutral injector" PINI although it refers equally to negative ion beams
    if num_pini gt 1 then begin
	if n_elements(size(alcbeam_results.angle_int,/dim)) lt 3 then message,'ALCBEAM output is missing velocity data! Required for num_pini>1'
	;decimate from every pini to num_pini to reduce computation
	all_pini = (size(ba,/dim))[0] ;This is the number of pini used in the ALCBEAM calculation
	if num_pini ge all_pini then begin
	    num_pini=all_pini
	    pini_x = ba[*,0]
	    pini_y = ba[*,1]
	    pini_w = alcbeam_results.angle_int[*,pini_index,*] ;will be an error if velocity distribution was not enabled in alcbeam
	    pini_w /= rebin(reform(total(pini_w,2),chord.npoint,1,n_comp),chord.npoint,num_pini,n_comp) ;normalize across dimension 2
	endif else begin
	    ;Todo: make a better beam source grid interpolation
	    pini_index = round(dindgen(num_pini)*(double(all_pini)-1)/(double(num_pini)-1))
	    pini_x = ba[pini_index,0]
	    pini_y = ba[pini_index,1]
	    pini_w = alcbeam_results.angle_int[*,pini_index,*] ;will be an error if velocity distribution was not enabled in alcbeam
	    pini_w /= rebin(reform(total(pini_w,2),chord.npoint,1,n_comp),chord.npoint,num_pini,n_comp)
	endelse
	badval = where(~finite(pini_w),badcount)
	if badcount gt 0 then pini_w[badval] = 0 ;hack for divide by zero problems

	;Note that alcbeam_results.angle corresponds to (pi-angle.alpha) because calc_geom_angles defines the viewing chord as pointing in the opposite direction

	;calculate the position of the pini in 3D
	beam_z = double(wall_pos) - double(grid_pos)
	beam_z /= norm(beam_z)
	beam_x = crossp([0,0,1],beam_z)
	beam_x /= norm(beam_x)
	beam_y = crossp(beam_z,beam_x)
	beam_y /= norm(beam_y)
	pini_pos = pini_x#beam_x + pini_y#beam_y + extend(num_pini,grid_pos)
    endif else begin
	pini_pos = grid_pos
	pini_w = 1
    endelse

    beam = {atom:beam_atom, E:E_beam, E_amu:E_amu, gamma:gamma_beam, beta:beta_beam, n_comp:n_comp, n:n_beam, n1:n_beam_n1, n2:n_beam_n2, n3:n_beam_n3, grid_pos:grid_pos, wall_pos:wall_pos, grid_rad:grid_rad, pini_pos:pini_pos, pini_w:pini_w}

;    if use_alcbeam_profiles
;	if alcbeam_input.n_e_prof.n_e_coord eq 0 then begin
;	    n_e = interpol(alcbeam_input.n_e_prof.n_e,alcbeam_input.n_e_prof.n_e_r,chord.r) > 0 ;1/cm^3
;	endif else begin
;	    n_e = interpol(alcbeam_input.n_e_prof.n_e,alcbeam_input.n_e_prof.n_e_r,chord.rho) > 0 ;1/cm^3
;	endelse
;	if alcbeam_input.t_e_prof.t_e_coord eq 0 then begin
;	    t_e = interpol(alcbeam_input.t_e_prof.t_e,alcbeam_input.t_e_prof.t_e_r,chord.r)*1000 > 0 ;eV
;	endif else begin
;	    t_e = interpol(alcbeam_input.t_e_prof.t_e,alcbeam_input.t_e_prof.t_e_r,chord.rho)*1000 > 0;eV
;	endelse
;	if alcbeam_input.z_eff_prof.z_eff_coord eq 0 then begin
;	    z_eff = interpol(alcbeam_input.z_eff_prof.z_eff,alcbeam_input.z_eff_prof.z_eff_r,chord.r) > 1
;	endif else begin
;	    z_eff = interpol(alcbeam_input.z_eff_prof.z_eff,alcbeam_input.z_eff_prof.z_eff_r,chord.rho) > 1
;	endelse
;	;Read ALCBEAM's table of impurities and convert this into density profiles of ions
;	;Unfortunately, ALCBEAM's method of specifying impurities requires some impurities with Z_imp != Z_main
;	;Cannot, for example, specify only impurity tritium in deuterium plasma
;	main_ion = alcbeam_input.plasma_param.main_ion
;	case main_ion of
;	    'H': Z_main = 1
;	    'D': Z_main = 1
;	    'He': Z_main = 2
;	endcase
;	impur_table = alcbeam_input.plasma_param.impur_table
;	imp_fr = double(impur_table[*,2])
;	imp_Z = double(impur_table[*,1]) ;allowed to be fractional to allow for partially ionized species
;	ions = [main_ion, impur_table[*,0]]
;	sumZimp = total(imp_fr*imp_Z)
;	sumZimp2 = total(imp_fr*imp_Z^2)
;	;in principle, Z_eff can be less than Z_main, e.g. He(H) plasma
;	if sumZimp2 - sumZimp*Z_main ne 0 then n_imp = n_e * (Z_eff - Z_main)/(sumZimp2 - sumZimp*Z_main) else n_imp = 0
;	n_main = (n_e - n_imp*sumZimp)/Z_main
;	n_ion = [[n_main], [n_imp # imp_fr]] ; [npoint, # impurities]
;	Z_ion = [Z_main, imp_Z]
;
;	;           cm^-3    eV       #            string
;	profiles = {n_e:n_e, T_e:T_e, Z_eff:Z_eff, main_ion:main_ion, ions:ions, n_ion:n_ion, Z_ion:Z_ion, T_i:T_e}
;    endif
end
