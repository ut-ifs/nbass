Function scalar_pr_xy,x1,y1,z1,x2,y2,z2
  ;calculates scalar product of two vectros in cartesian coordinates
  return,x1*x2+y1*y2+z1*z2
end

Function ang_2vect_xy,x1,y1,z1,x2,y2,z2
  ;calculates angle between two vectors in cartesian coordinates
  return,acos(scalar_pr_xy(x1,y1,z1,x2,y2,z2)/sqrt(scalar_pr_xy(x1,y1,z1,x1,y1,z1)*scalar_pr_xy(x2,y2,z2,x2,y2,z2)))
end


; The purpose of this programm to extracts ALCBEAM output results
; (beam density, excitation) for any point in plasma (1D,2D,3D arrays)
;Input R_tor (m),Z_tor(m),Phi_tor(rad)
;Either alcbeam_file or mdsplus_run value should be choosen
;template for the mdsplus_run is "BESPAM.RUN_1"
; 4 Aug 2012 (beam velocity angle is added)
;(R0,Phi0,Z0) - observation point (scalar or array like r_tor)
; if these are defined, calculate the angles
;Angle is calculated between chord_vect and beamlets
; 26 Aug 2012 KTL (output t1 and t2 via keywords)
; 31 Aug 2015 KTL fix to load x_bml and y_bml needed for velocity distribution n_bml
; 20 Jan 2016 KTL angles are now shifted by beam_port_phi. To get old behavior, set keyword "relative"
pro get_alcbeam, shot=shot, alcbeam_file=alcbeam_file, mdsplus_run=mdsplus_run, R_tor=R_tor, Z_tor=Z_tor, Phi_tor=Phi_tor, R0=R0, Phi0=phi0, Z0=Z0, alcbeam_results=alcbeam_results,alcbeam_geom=alcbeam_geom, t1=t1, t2=t2, quiet=quiet, relative=relative
  if keyword_set(alcbeam_file) then begin
    if (file_search(alcbeam_file))(0) eq "" then begin
      print, 'Error: Wrong filename'
      return
    endif
    openr,lun,alcbeam_file,/get_lun
    vel_vec_coef_arr=0
    vel_vec_x_arr=0
    vel_vec_y_arr=0

    val="template"
    while ~EOF(lun) do begin
      readf,lun,val

     ;load general info

      if val eq 'beam:' then begin
        readf,lun,val
        beam=strtrim(val,2)
        print,'Beam: '+beam
      endif
      if val eq 'shot:' then begin
        readf,lun,val
        shot=strtrim(val,2)
        print,'Shot: '+shot
      endif
      if val eq 'time_interv:' then begin
        readf,lun,val
        time_interv=float(strsplit(val,',',/extract))
        t1=time_interv(0)
        t2=time_interv(1)
        print,'Time interval: '+val+' sec'
      endif

      ;load beam geometry

      if val eq 'x_bml:' then begin
        readf,lun,val
        x_bml=float(strsplit(val,',',/extract))
      endif
      if val eq 'y_bml:' then begin
        readf,lun,val
        y_bml=float(strsplit(val,',',/extract))
      endif
      if val eq 'beam_port_phi:' then begin
        readf,lun,val
        beam_port_phi=float(val)
      endif
      if val eq 'tank_front_dist:' then begin
        readf,lun,val
        tank_front_dist=float(val)
      endif
      if val eq 'tank_size:' then begin
        readf,lun,val
        tank_size=float(val)
      endif
      if val eq 'neutr_size:' then begin
        readf,lun,val
        neutr_size=float(val)
      endif
      if val eq 'r_grid:' then begin
        readf,lun,val
        R_grid=float(val)
      endif
      if val eq 'z_grid:' then begin
        readf,lun,val
        z_grid=float(val)
      endif
      if val eq 'phi_grid:' then begin
        readf,lun,val
        phi_grid=float(val)
      endif
       if val eq 'r_wall:' then begin
        readf,lun,val
        R_wall=float(val)
      endif
      if val eq 'z_wall:' then begin
        readf,lun,val
        z_wall=float(val)
      endif
      if val eq 'phi_wall:' then begin
        readf,lun,val
        phi_wall=float(val)
      endif

     ;load beam parameters

      if val eq 'e_full:' then begin
        readf,lun,val
        e_full=float(val)
        print,'Beam full energy: '+val+' keV'
      endif
      if val eq 'e_frac:' then begin
        readf,lun,val
        e_frac=float(strsplit(val,',',/extract))
        print, 'Energy fractions: ','(E/'+strtrim(string(round(1.0/e_frac),format='(I2)'),2)+')'
      endif
     ;Load output data
      if val eq 'e_beam:' then begin
        readf,lun,val
        e_beam=float(strsplit(val,',',/extract))
      endif
      if val eq 'z_beam:' then begin
        readf,lun,val
        z_beam=float(strsplit(val,',',/extract))
      endif
      if val eq 'x_beam:' then begin
        readf,lun,val
        x_beam=float(strsplit(val,',',/extract))
      endif
      if val eq 'y_beam:' then begin
        readf,lun,val
        y_beam=float(strsplit(val,',',/extract))
      endif
      if val eq 'n_beam:' then begin
        n_ebeam=n_elements(e_beam)
        n_z=n_elements(z_beam)
        n_x=n_elements(x_beam)
        n_y=n_elements(y_beam)
        n_beam_arr=fltarr(n_ebeam,n_z,n_x,n_y)
        readu,lun,n_beam_arr
      endif
      if val eq 'exc_n2_frac:' then begin
        n_z=n_elements(z_beam)
        n_x=n_elements(x_beam)
        n_y=n_elements(y_beam)
        exc_n2_frac_arr=fltarr(n_ebeam,n_z,n_x,n_y)
        readu,lun,exc_n2_frac_arr
      endif
      if val eq 'exc_n3_frac:' then begin
        n_z=n_elements(z_beam)
        n_x=n_elements(x_beam)
        n_y=n_elements(y_beam)
        exc_n3_frac_arr=fltarr(n_ebeam,n_z,n_x,n_y)
        readu,lun,exc_n3_frac_arr
      endif
    if val eq 'vel_vec_x:' then begin
      n_bml=n_elements(x_bml)
      n_x=n_elements(x_beam)
      n_y=n_elements(y_beam)
      vel_vec_x_arr=fltarr(n_x,n_y,n_bml)
      readu,lun,vel_vec_x_arr
    endif
    if val eq 'vel_vec_y:' then begin
      n_bml=n_elements(x_bml)
      n_x=n_elements(x_beam)
      n_y=n_elements(y_beam)
      vel_vec_y_arr=fltarr(n_x,n_y,n_bml)
      readu,lun,vel_vec_y_arr
    endif
    if val eq 'vel_vec_coef:' then begin
      n_ebeam=n_elements(e_beam[where(e_beam gt 0)]) ;fix for data with halo
      n_z=n_elements(z_beam)
      n_x=n_elements(x_beam)
      n_y=n_elements(y_beam)
      vel_vec_coef_arr=fltarr(n_ebeam,n_z,n_x,n_y,9)
      readu,lun,vel_vec_coef_arr
    endif
    endwhile
    free_lun,lun
  endif


  if keyword_set(mdsplus_run) and keyword_set(shot) then begin
     MDSOPEN,'DNB',shot,/Quiet, status=st0
     if st0 then begin
       n_beam_arr=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run+'.RESULTS:n_beam',status=st1,/quiet)
       exc_n2_frac_arr=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run+'.RESULTS:exc_n2_frac',status=st2,/quiet)
       exc_n3_frac_arr=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run+'.RESULTS:exc_n3_frac',status=st3,/quiet)
       E_beam=MDSVALUE('Dim_Of(\DNB::TOP.ALCBEAM.'+mdsplus_run+'.RESULTS:n_beam,0)',status=st4,/quiet)
       z_beam=MDSVALUE('Dim_Of(\DNB::TOP.ALCBEAM.'+mdsplus_run+'.RESULTS:n_beam,1)',status=st5,/quiet)
       x_beam=MDSVALUE('Dim_Of(\DNB::TOP.ALCBEAM.'+mdsplus_run+'.RESULTS:n_beam,2)',status=st6,/quiet)
       y_beam=MDSVALUE('Dim_Of(\DNB::TOP.ALCBEAM.'+mdsplus_run+'.RESULTS:n_beam,3)',status=st7,/quiet)
       time_int=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run+':time_interv',status=st8,/quiet)
       vel_dis_type=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run+':vel_dis_type',status=st9,/quiet)
       vel_vec_x_arr=0 & vel_vec_y_arr=0 & vel_vec_coef_arr=0
       if vel_dis_type eq 'YES' then begin
         vel_vec_x_arr=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run+'.RESULTS:vel_vec_x',status=st9,/quiet)
         vel_vec_y_arr=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run+'.RESULTS:vel_vec_y',status=st9,/quiet)
         vel_vec_coef_arr=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run+'.RESULTS:vel_vec_coef',status=st9,/quiet)
       endif

       grid_focus=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run+'.INPUT.BEAM_GEOM:GRID_FOCUS', status=st9,/quiet)
       if not(st9) then grid_focus=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run+'.INPUT.BEAM_GEOM:X_GRID_FOCUS', status=st9,/quiet)
       beam_port_phi=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run+'.INPUT.BEAM_GEOM:BEAM_PORT_PHI')
       beam_port=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run+'.INPUT.BEAM_GEOM:BEAM_PORT', status=st10,/quiet)
       R_grid=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run+'.INPUT.BEAM_GEOM:R_grid', status=st11,/quiet)
       Z_grid=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run+'.INPUT.BEAM_GEOM:Z_grid', status=st12,/quiet)
       Phi_grid=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run+'.INPUT.BEAM_GEOM:Phi_grid', status=st13,/quiet)
       R_wall=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run+'.INPUT.BEAM_GEOM:R_wall', status=st14,/quiet)
       Z_wall=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run+'.INPUT.BEAM_GEOM:Z_wall', status=st15,/quiet)
       Phi_wall=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run+'.INPUT.BEAM_GEOM:Phi_wall', status=st16,/quiet)
       tank_front_dist=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run+'.INPUT.BEAM_GEOM:TANK_FRONT', status=st17,/quiet)
       tank_size=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run+'.INPUT.BEAM_GEOM:TANK_SIZE', status=st18,/quiet)
       neutr_size=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run+'.INPUT.BEAM_GEOM:NEUTR_SIZE', status=st19,/quiet)
       tank_magnet_dist=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run+'.INPUT.BEAM_GEOM:TANK_MAGNET', status=st20,/quiet)
       magnet_size=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run+'.INPUT.BEAM_GEOM:MAGNET_SIZE', status=st21,/quiet)
       tank_cal_dist=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run+'.INPUT.BEAM_GEOM:TANK_CAL', status=st22,/quiet)
       beam_apertures=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run+'.INPUT.BEAM_GEOM:BEAM_APERTUR', status=st23,/quiet)


       E_full=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run+'.INPUT.BEAM_PARAM:E_FULL', status=st24,/quiet)
       E_frac=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run+'.INPUT.BEAM_PARAM:E_FRAC', status=st25,/quiet)


       st_tot=st1*st2*st3*st4*st5*st6*st7*st8*st9*st10*st11*st12*st13*st14*st15*st16*st17*st18*st19*st20*st21*st22*st23*st24*st25
       if st_tot then begin
         t1=time_int(0)
         t2=time_int(1)
         if not keyword_set(quiet) then begin
           print,'Beam: DNBI_ALCATOR'
           print,'Shot: '+strtrim(shot,2)
           print,'Time interval: '+strtrim(string(t1,format='(F5.3)'),2)+', '+strtrim(string(t2,format='(F5.3)'),2)+' sec'
           print,'Beam full energy: '+strtrim(string(E_full,format='(F6.3)'),2)+' keV'
           print, 'Energy fractions: ','(E/'+strtrim(string(round(1.0/e_frac),format='(I2)'),2)+')'
         endif
       endif else begin
         print,'Error: Not sufficient input data in the MDSPLUS'
stop
         return
       endelse
    endif else begin
      print,'Error: Wrong shot number'
      return
    endelse
  endif

  if ~keyword_set(relative) then begin
    phi_grid += beam_port_phi
    phi_wall += beam_port_phi
  endif

  grid_cent_x=r_grid*cos(phi_grid)
  grid_cent_y=-r_grid*sin(phi_grid)
  wall_cent_x=r_wall*cos(phi_wall)
  wall_cent_y=-r_wall*sin(phi_wall)
  dist_all_XY=sqrt((wall_cent_x-grid_cent_x)^2.0+(wall_cent_y-grid_cent_y)^2.0)
  sin_pivot=-(wall_cent_y-grid_cent_y)/dist_all_XY;angle in tokamak XY plane (horizontal)
  cos_pivot=-(wall_cent_x-grid_cent_x)/dist_all_XY

  dist_all_XYZ=sqrt((wall_cent_x-grid_cent_x)^2.0+(wall_cent_y-grid_cent_y)^2.0+(z_wall-z_grid)^2.0)
  sin_alpha=(z_wall-z_grid)/dist_all_XYZ ; angle in tokamak XZ plane (vertical)
  cos_alpha=sqrt(1.0-sin_alpha^2.0)


  if n_elements(R_tor) eq 0 or n_elements(Z_tor) eq 0 or n_elements(phi_tor) eq 0 then begin
    print,'Error: Missing R_tor, Z_tor, or Phi_tor'
    return
  endif
  size_r=size(float(R_tor))
  size_Z=size(float(Z_tor))
  size_phi=size(float(Phi_tor))
  if total(size_r) ne total(size_z) or total(size_r) ne total(size_phi) or total(size_phi) ne total(size_z) then begin
    print,'Error: R_tor, Z_tor, and Phi_tor have inconsistent dimentions'
    return
  endif

  x_tor2=R_tor*cos(Phi_tor)
  y_tor2=-R_tor*sin(Phi_tor) ;why is there a negative sign here?
  z_tor2=z_tor
  z_beam1=-(x_tor2-r_grid*cos(phi_grid))*cos_pivot-(y_tor2+r_grid*sin(phi_grid))*sin_pivot
  x_beam1=(x_tor2-r_grid*cos(phi_grid))*sin_pivot-(y_tor2+r_grid*sin(phi_grid))*cos_pivot
  y_beam1=z_tor2
  z_b=z_beam1*cos_alpha+(y_beam1-z_grid)*sin_alpha
  y_b=-z_beam1*sin_alpha+(y_beam1-z_grid)*cos_alpha
  x_b=x_beam1

  n_e=n_elements(e_beam)
  n_x=n_elements(x_beam)
  n_y=n_elements(y_beam)
  n_z=n_elements(z_beam)

  x_ind=interpol(make_array(n_x,/index),x_beam,x_b)
  y_ind=interpol(make_array(n_y,/index),y_beam,y_b)
  z_ind=interpol(make_array(n_z,/index),z_beam,z_b)

  if n_elements(R0) gt 0 and n_elements(Phi0) gt 0 and n_elements(Z0) gt 0 then begin
    ;get chord vector
    x_tor2_vect=R_tor*cos(Phi_tor)-R0*cos(phi0)
    y_tor2_vect=-R_tor*sin(Phi_tor)+r0*sin(phi0) ;why is there a negative sign here?
    z_tor2_vect=z_tor-z0
    z_beam1_vect=-(x_tor2_vect)*cos_pivot-(y_tor2_vect)*sin_pivot
    x_beam1_vect=(x_tor2_vect)*sin_pivot-(y_tor2_vect)*cos_pivot
    y_beam1_vect=z_tor2_vect
    z_b_vect=z_beam1_vect*cos_alpha+(y_beam1_vect)*sin_alpha
    y_b_vect=-z_beam1_vect*sin_alpha+(y_beam1_vect)*cos_alpha
    x_b_vect=x_beam1_vect

  endif


  n_beam=interpolate(reform(n_beam_arr(0,*,*,*)),z_ind,x_ind,y_ind)
if total(n_beam lt 0) gt 0 then message,'negative beam density error',/continue
  n_beam>=0
  exc_n2_frac=interpolate(reform(exc_n2_frac_arr(0,*,*,*)),z_ind,x_ind,y_ind)
  exc_n3_frac=interpolate(reform(exc_n3_frac_arr(0,*,*,*)),z_ind,x_ind,y_ind)

  ; get velocity data
  angle=0
  angle_int=0
  if n_elements(vel_vec_coef_arr) gt 1 and n_elements(R0) gt 0 and n_elements(Phi0) gt 0 and n_elements(Z0) gt 0 then begin

     size_n=size(n_beam,/dimensions) ;number of points at which to evaluate beam
     size_coef=(size(vel_vec_coef_arr,/dimensions))(4) ;always 9
     size_x=(size(vel_vec_x_arr,/dimensions))(2) ;number of beam pores
    ; alpha=n_elements(size_x)
    ; alpha_int=n_elements(size_x)
    if size_n(0) eq 0 then vel_vec_coef=fltarr(size_coef) else vel_vec_coef=fltarr(size_n,size_coef)
    size_n_tot=n_elements(n_beam)
    unit_arr=make_array(size_n_tot,/index)
    for i=0, size_coef-1 do vel_vec_coef(unit_arr+size_n_tot*i)=interpolate(reform(vel_vec_coef_arr(0,*,*,*,i)),z_ind,x_ind,y_ind)
    ;if size_n(0) eq 0 then vel_vec_x=fltarr(size_x) else vel_vec_x=fltarr(size_n,size_x)
    ;vel_vec_y=vel_vec_x
    ;vel_vec_z=vel_vec_x
    if size_n(0) eq 0 then angle=fltarr(size_x) else angle=fltarr(size_n,size_x)
    angle_int=angle

    for i=0, size_x-1 do begin
      ;vel_vec_x(make_array(n_elements(n_beam),/index)+n_elements(n_beam)*i)=interpolate(reform(vel_vec_x_arr(*,*,i)),x_ind,y_ind)
      ;vel_vec_y(make_array(n_elements(n_beam),/index)+n_elements(n_beam)*i)=interpolate(reform(vel_vec_y_arr(*,*,i)),x_ind,y_ind)
      ;vel_vec_z(make_array(n_elements(n_beam),/index)+n_elements(n_beam)*i)=z_ind
      x_bml1=interpolate(reform(vel_vec_x_arr(*,*,i)),x_ind,y_ind,missing=!values.f_nan)
      y_bml1=interpolate(reform(vel_vec_y_arr(*,*,i)),x_ind,y_ind,missing=!values.f_nan)
      z_bml1=z_b
      angle(make_array(n_elements(n_beam),/index)+n_elements(n_beam)*i)=ang_2vect_xy(x_b_vect,y_b_vect,z_b_vect,x_bml1,y_bml1,z_bml1)
      angle_int(unit_arr+size_n_tot*i)=$
      vel_vec_coef(unit_arr+size_n_tot*0)+$
      vel_vec_coef(unit_arr+size_n_tot*1)*y_bml1+$
      vel_vec_coef(unit_arr+size_n_tot*2)*y_bml1^2.0+$
      vel_vec_coef(unit_arr+size_n_tot*3)*x_bml1+$
      vel_vec_coef(unit_arr+size_n_tot*4)*x_bml1*y_bml1+$
      vel_vec_coef(unit_arr+size_n_tot*5)*x_bml1*y_bml1^2.0+$
      vel_vec_coef(unit_arr+size_n_tot*6)*x_bml1^2+$
      vel_vec_coef(unit_arr+size_n_tot*7)*x_bml1^2*y_bml1+$
      vel_vec_coef(unit_arr+size_n_tot*8)*x_bml1^2*y_bml1^2.0
    ;This 2nd degree reconstruction has a problem! Extrapolation outside the ends of the grid can give
    ;high weights due to the parabola curving back up.
   endfor
 endif
  angle_int >= 0 ;any values < 0 are due to errors in the 2nd degree reconstruction and aren't real


  data={E:e_beam(0),n_beam:n_beam,exc_n2_frac:exc_n2_frac,exc_n3_frac:exc_n3_frac,angle:angle,angle_int:angle_int}
  alcbeam_results=replicate(data,n_elements(e_beam))
  alcbeam_geom={r_grid:r_grid,z_grid:z_grid,phi_grid:phi_grid,r_wall:r_wall,z_wall:z_wall,phi_wall:phi_wall}
  for j=1, n_elements(e_beam)-1 do begin ;changed from e_frac to e_beam to support halo data
    ;e_ind=make_array(n_elements(x_ind),value=i)
    ;n_beam=n_beam_arr(e_ind,z_ind,x_ind,y_ind)
    ;exc_n2_frac=exc_n2_frac_arr(e_ind,z_ind,x_ind,y_ind)
    ;exc_n3_frac=exc_n3_frac_arr(e_ind,z_ind,x_ind,y_ind)

    n_beam=interpolate(reform(n_beam_arr(j,*,*,*)),z_ind,x_ind,y_ind)
if total(n_beam lt 0) gt 0 then message,'negative beam density error',/continue
    n_beam >= 0
    exc_n2_frac=interpolate(reform(exc_n2_frac_arr(j,*,*,*)),z_ind,x_ind,y_ind)
    exc_n3_frac=interpolate(reform(exc_n3_frac_arr(j,*,*,*)),z_ind,x_ind,y_ind)

  ; get velocity data
  if n_elements(vel_vec_coef_arr) gt 1 and n_elements(R0) gt 0 and n_elements(Phi0) gt 0 and n_elements(Z0) gt 0 and e_beam[j] gt -1 then begin


    for i=0, size_coef-1 do vel_vec_coef(unit_arr+size_n_tot*i)=interpolate(reform(vel_vec_coef_arr(j,*,*,*,i)),z_ind,x_ind,y_ind)
    for i=0, size_x-1 do begin
      x_bml1=interpolate(reform(vel_vec_x_arr(*,*,i)),x_ind,y_ind)
      y_bml1=interpolate(reform(vel_vec_y_arr(*,*,i)),x_ind,y_ind)
      z_bml1=z_b
      angle_int(unit_arr+size_n_tot*i)=$
      vel_vec_coef(unit_arr+size_n_tot*0)+$
      vel_vec_coef(unit_arr+size_n_tot*1)*y_bml1+$
      vel_vec_coef(unit_arr+size_n_tot*2)*y_bml1^2.0+$
      vel_vec_coef(unit_arr+size_n_tot*3)*x_bml1+$
      vel_vec_coef(unit_arr+size_n_tot*4)*x_bml1*y_bml1+$
      vel_vec_coef(unit_arr+size_n_tot*5)*x_bml1*y_bml1^2.0+$
      vel_vec_coef(unit_arr+size_n_tot*6)*x_bml1^2+$
      vel_vec_coef(unit_arr+size_n_tot*7)*x_bml1^2*y_bml1+$
      vel_vec_coef(unit_arr+size_n_tot*8)*x_bml1^2*y_bml1^2.0
    endfor
   endif
    data={E:e_beam(j),n_beam:n_beam,exc_n2_frac:exc_n2_frac,exc_n3_frac:exc_n3_frac,angle:angle,angle_int:angle_int}
    alcbeam_results(j)=data
  endfor
  if not keyword_set(quiet) then help,alcbeam_results,/structure
end


pro get_alcbeam_test
  systime_1=systime(/seconds)
  R0=0.89
  Z0=-0.26
  phi0=-25.28*!Pi/180.0

  R_mid=0.8
  phi_mid=5*!Pi/180.0
  z_mid=0

  ;cartesian transformation
  x0=R0*sin(phi0)
  y0=R0*cos(phi0)

  x_mid=R_mid*sin(phi_mid)
  y_mid=R_mid*cos(phi_mid)

  n_points=500
  x_tor=interpol([x0,2.0*x_mid-x0],n_points)
  y_tor=interpol([y0,2.0*y_mid-y0],n_points)
  z_tor=interpol([z0,2.0*z_mid-z0],n_points)

  r_tor=sqrt(x_tor^2.0+y_tor^2.0)
  phi_tor=atan(x_tor/y_tor)

  ;get_alcbeam,file_name='DNBI_ALCATOR_2.abo',r_tor=r_tor,z_tor=z_tor,phi_tor=phi_tor,alcbeam_results=alcbeam_results
  ;get_alcbeam, shot='1070831028', mdsplus_run='bespam.run_4',r_tor=r_tor,z_tor=z_tor,phi_tor=phi_tor,alcbeam_results=alcbeam_results
  get_alcbeam, shot='1120621025', mdsplus_run='bespam.run_4',r_tor=r_tor,z_tor=z_tor,phi_tor=phi_tor,r0=r0,z0=z0,phi0=phi0,alcbeam_results=alcbeam_results
  ;get_alcbeam, shot='1070831028', alcbeam_file='DNBI_ALCATOR_12.abo',r_tor=r_tor,z_tor=z_tor,phi_tor=phi_tor,alcbeam_results=alcbeam_results

  if n_elements(alcbeam_results) ne 0 then begin
    window,0
    plot, x_tor,(alcbeam_results(0)).n_beam
  endif
  stop
end
