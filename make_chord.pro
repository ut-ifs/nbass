;Create structure 'chord' which contains the coordinates of the grid points which make up the chord
;Choice of grid structure depends on the chordmodel parameter
pro make_chord,chord
    @view_param.idl
    @mesh_param.idl

    ;Construct viewing chord mesh
    ;discretize the chord into steps starting from the optic_pos pointing toward the plasma
    view_vect = norm_vect(optic_pos - spot_pos) ;looking from plasma toward lens [3]
    case chordmodel of
	'ray': begin ;basic 1D ray, with analytical small-angle broadening applied
	    step = -view_vect * ds ;step from lens toward plasma
	    nstep = long(((cut_z1-optic_pos[2])/step[2] > (cut_z2-optic_pos[2])/step[2]) < ((-cut_r2-optic_pos[0])/step[0] > (cut_r2-optic_pos[0])/step[0]) < ((-cut_r2-optic_pos[1])/step[1] > (cut_r2-optic_pos[1])/step[1]))
	    chx = optic_pos[0]+dindgen(nstep)*step[0]
	    chy = optic_pos[1]+dindgen(nstep)*step[1]
	    chz = optic_pos[2]+dindgen(nstep)*step[2]
	    view_vect_all = view_vect
	    nray = 1
	end
	'grid': begin ;make a complete bipartite graph between aperture grid and spot grid
	    ;construct a ray which is perpendicular to view_vect
	    if abs(view_vect[2]) lt 0.99 then begin
		ray = norm_vect([0,0,1] - view_vect[2]*view_vect)  ;perpendicular to view, approximately "up"
	    endif else begin
		ray = norm_vect([1,0,0] - view_vect[0]*view_vect)  ;view is vertical, so look for a perpendicular ray oriented approximately "right"
	    endelse
	    ;construct rotation matrices
	    cos60 = 0.5d
	    sin60 = sqrt(0.75d)
	    rot60 = [[cos60+view_vect[0]^2*(1-cos60),view_vect[0]*view_vect[1]*(1-cos60)-view_vect[2]*sin60,view_vect[0]*view_vect[2]*(1-cos60)+view_vect[1]*sin60], $
		    [view_vect[1]*view_vect[0]*(1-cos60)+view_vect[2]*sin60,cos60+view_vect[1]^2*(1-cos60),view_vect[1]*view_vect[2]*(1-cos60)-view_vect[0]*sin60], $
		    [view_vect[2]*view_vect[0]*(1-cos60)-view_vect[1]*sin60,view_vect[2]*view_vect[1]*(1-cos60)+view_vect[0]*sin60,cos60+view_vect[2]^2*(1-cos60)]]
	    cos30 = sqrt(0.75d)
	    sin30 = 0.5d
	    rot30 = [[cos30+view_vect[0]^2*(1-cos30),view_vect[0]*view_vect[1]*(1-cos30)-view_vect[2]*sin30,view_vect[0]*view_vect[2]*(1-cos30)+view_vect[1]*sin30], $
		    [view_vect[1]*view_vect[0]*(1-cos30)+view_vect[2]*sin30,cos30+view_vect[1]^2*(1-cos30),view_vect[1]*view_vect[2]*(1-cos30)-view_vect[0]*sin30], $
		    [view_vect[2]*view_vect[0]*(1-cos30)-view_vect[1]*sin30,view_vect[2]*view_vect[1]*(1-cos30)+view_vect[0]*sin30,cos30+view_vect[2]^2*(1-cos30)]]
	    ;build hex grids in a plane perpendicular to view_vect
	    zerov = dblarr(1,3)
	    ray0 = reform(ray,1,3)
	    ray60 = rot60##ray0
	    ray90 = rot30##ray60
	    ray120 = rot60##ray60
	    array6 = [ray0,ray60,ray120,-ray0,-ray60,-ray120]
	    array7 = [zerov,array6*0.666666666666667d]  ;[7,3]
	    array19 = [zerov,array6*0.4d,array6*0.8d,rot30##array6*0.8d] ;[19,3]

	    case ap_shape of
		'point': optic_grid = transpose(optic_pos)
		'hex7':  optic_grid = ap_rad*array7 + extend(7,optic_pos)
		'hex19': optic_grid = ap_rad*array19 + extend(19,optic_pos)
		'rect': stop
		'grid': optic_grid = outer(ap_gridy,ray) + outer(ap_gridx,reform(ray90)) + extend(n_elements(ap_gridx),optic_pos)
	    endcase
	    n_optic = (size(optic_grid,/dim))[0]
	    case sp_shape of
		'point': spot_grid = transpose(spot_pos)
		'hex7':  spot_grid = sp_rad*array7 + extend(7,spot_pos)
		'hex19': spot_grid = sp_rad*array19 + extend(19,spot_pos)
		'rect': stop
		'grid': spot_grid = outer(sp_gridy,ray) + outer(sp_gridx,reform(ray90)) + extend(n_elements(sp_gridx),spot_pos)
	    endcase
	    n_spot = (size(spot_grid,/dim))[0]

	    step = -view_vect * ds ;step from lens toward plasma
	    nstep = long(((cut_z1-optic_pos[2])/step[2] > (cut_z2-optic_pos[2])/step[2]) < ((-cut_r2-optic_pos[0])/step[0] > (cut_r2-optic_pos[0])/step[0]) < ((-cut_r2-optic_pos[1])/step[1] > (cut_r2-optic_pos[1])/step[1]))
	    nray = n_optic*n_spot

	    chx = dblarr(nstep*nray)
	    chy = chx
	    chz = chx
	    view_vect_all = dblarr(3,nstep*nray)
	    count = 0
	    for i=0,n_optic-1 do begin
		for j=0,n_spot-1 do begin
		    view_vect2 = norm_vect(optic_grid[i,*] - spot_grid[j,*]) ;looking from plasma toward lens [3]
		    step = -view_vect2 * ds ;step from lens toward plasma
		    chx[count] = optic_grid[i,0]+dindgen(nstep)*step[0]
		    chy[count] = optic_grid[i,1]+dindgen(nstep)*step[1]
		    chz[count] = optic_grid[i,2]+dindgen(nstep)*step[2]
		    view_vect_all[0,count] = rebin(transpose(view_vect2),3,nstep)
		    count += nstep
		endfor
	    endfor
	end
	'montecarlo': begin ;make random rays between aperture and spot
	    stop ;not implemented
	end
    endcase

    chr = sqrt(chx^2+chy^2)
    chphi = atan(chy,chx)
    npoint = n_elements(chx)
    chrmid = dblarr(npoint)
    chrho = dblarr(npoint)

    chord = {npoint:npoint, x:chx, y:chy, z:chz, r: chr, phi: chphi, rmid:chrmid, rho:chrho, vect:transpose(view_vect_all), ds:ds, nray:nray} ;package chord data for output
end
