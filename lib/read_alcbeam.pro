; The purpose of this programm to extracts full set of ALCBEAM output results
;Either alcbeam_file or mdsplus_run value should be choosen
;template for the mdsplus_run is "BESPAM.RUN_1"
;
;modifications
; 8/19/2015 KTL, fix read limiters, remove status_wid lines
pro read_alcbeam, shot=shot, alcbeam_file=alcbeam_file, mdsplus_run=mdsplus_run, alcbeam_results=alcbeam_results,alcbeam_input=alcbeam_input

if keyword_set(alcbeam_file) then n_runs=n_elements(alcbeam_file)
if keyword_set(mdsplus_run) then n_runs=n_elements(mdsplus_run)

alcbeam_input = ptrarr(n_runs)
alcbeam_results = ptrarr(n_runs)
for i=0,n_runs-1 do begin
  
  if keyword_set(alcbeam_file) then begin
    if (file_search(alcbeam_file(i)))(0) eq "" then begin
      print, 'Error: Wrong filename'
      return
    endif
    
    beam_atom='H';default
    val='template'
    close,1
    openr,1,alcbeam_file(i)
    st=0
    ;put 0 for duct pressure
    duct_pressure=0.0
    duct_pressure_loc=0.0
    vel_dis_type='NO' ;default
    vel_vec_coef=0
    vel_vec_x=0
    vel_vec_y=0

    code_grid_arr={z:[0.0,0.0,0.0,0.0,0.0],x:[0.0,0.0,0.0],y:[0.0,0.0,0.0]}

    while ~EOF(1) do begin

    readf,1,val

    ;load general info

    if val eq 'beam:' then begin 
      readf,1,val
      beam=strtrim(val,2)
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'shot:' then begin 
      readf,1,val
      shot=strtrim(val,2)
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'time_interv:' then begin 
      readf,1,val
      time_interv=float(strsplit(val,', ',/extract))
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'div_type:' then begin 
      readf,1,val
      div_type_val=strtrim(val,2)
      div_type=div_type_val
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'atten_type:' then begin 
      readf,1,val
      atten_type_val=strtrim(val,2)
      atten_type=atten_type_val
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'vel_dis_type:' then begin 
      readf,1,val
      vel_dis_type_val=strtrim(val,2)
      vel_dis_type=vel_dis_type_val
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'stop_plsm_cs:' then begin 
      readf,1,val
      stop_plasma_type_val=strtrim(val,2)
      stop_plasma_type=stop_plasma_type_val
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'exc_plsm_cs:' then begin 
      readf,1,val
      exc_plasma_type_val=strtrim(val,2)
      exc_plasma_type=exc_plasma_type_val
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif

    ;load beam geometry
    if val eq 'x_bml:' then begin 
      readf,1,val
      x_bml=float(strsplit(val,', ',/extract))
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'y_bml:' then begin 
      readf,1,val
      y_bml=float(strsplit(val,', ',/extract))
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'grid_ap_diam:' then begin 
      readf,1,val
      grid_ap_diam=float(val)
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'beam_port_phi:' then begin 
      readf,1,val
      beam_port_phi=val
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'beam_port:' then begin 
      readf,1,val
      beam_port=val
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'grid_focus:' then begin 
      readf,1,val
      x_grid_focus=float(val)
      y_grid_focus=x_grid_focus
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'x_grid_focus:' then begin 
      readf,1,val
      x_grid_focus=float(val)
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'y_grid_focus:' then begin 
      readf,1,val
      y_grid_focus=float(val)
    endif
    if val eq 'tank_front_dist:' then begin 
      readf,1,val
      tank_front_dist=float(val)
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'tank_size:' then begin 
      readf,1,val
      tank_size=float(val)
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'neutr_size:' then begin 
      readf,1,val
      neutr_size=float(val)
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'tank_magnet_dist:' then begin 
      readf,1,val
      tank_magnet_dist=float(val)
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'magnet_size:' then begin 
      readf,1,val
      magnet_size=float(val)
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'tank_cal_dist:' then begin 
      readf,1,val
      tank_cal_dist=float(val)
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'r_grid:' then begin 
      readf,1,val
      R_grid=float(val)
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'z_grid:' then begin 
      readf,1,val
      z_grid=float(val)
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'phi_grid:' then begin 
      readf,1,val
      phi_grid=float(val)
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'r_wall:' then begin 
      readf,1,val
      R_wall=float(val)
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'z_wall:' then begin 
      readf,1,val
      z_wall=float(val)
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'phi_wall:' then begin 
      readf,1,val
      phi_wall=float(val)
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif

    ;load beam parameters
    if val eq 'beam_atom:' then begin 
      readf,1,val
      beam_atom=strtrim(val,2)
    endif
    if val eq 'e_full:' then begin 
      readf,1,val
      e_full=float(val)
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'e_frac:' then begin 
      readf,1,val
      e_frac=float(strsplit(val,', ',/extract))
      e_beam=e_full*e_frac
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'i_beam:' then begin 
      readf,1,val
      i_beam=float(val)
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'i_frac:' then begin 
      readf,1,val
      I_frac=float(strsplit(val,', ',/extract))
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'i_dens_par:' then begin 
      readf,1,val
      i_dens_par=float(val)
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'i_opt:' then begin 
      readf,1,val
      i_opt=float(val)
    endif
    if val eq 'div_bml_opt:' then begin 
      readf,1,val
      x_div_bml_opt=float(val)
      y_div_bml_opt=x_div_bml_opt
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'x_div_bml_opt:' then begin 
      readf,1,val
      x_div_bml_opt=float(val)
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'y_div_bml_opt:' then begin 
      readf,1,val
      y_div_bml_opt=float(val)
    endif
    if val eq 'div_dist_par:' then begin 
      readf,1,val
      div_dist_par=float(val)
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif

    ;Load beam limiters
    if val eq 'n_limiters:'then begin 
      readf,1,val
      n_limiters=fix(val)
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'limiters_table:'then begin 
      readf,1,val
      readf,1,val
      lim_l=n_elements(strsplit(val,':',/extract))
      readf,1,val
      if n_limiters eq 0 then limiters_table='' else begin
        limiters_table=strarr(7,n_limiters)
        limiters_table(*,*)='NAN'
        for j=0, n_limiters-1 do begin
          readf,1,val
          limiters_table(0:lim_l-1,j)=strtrim(strsplit(val,':',/extract,/regex),2)
        endfor
      endelse
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif

    ;Load neutral gas
    if val eq 'tank_pressure:' then begin 
      readf,1,val
      tank_pressure=float(val)
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'torus_pressure:' then begin 
      readf,1,val
      torus_pressure=float(val)
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'duct_pressure:' then begin 
      readf,1,val
      duct_pressure=float(val)
    endif
    if val eq 'duct_pressure_loc:' then begin 
      readf,1,val
      duct_pressure_loc=float(val)
    endif

    ;Load plasma geometry
    if val eq 'r_major:' then begin 
      readf,1,val
      r_major=float(val)
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'z_major:' then begin 
      readf,1,val
      z_major=float(val)
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'r_minor:' then begin 
      readf,1,val
      r_minor=float(val)
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'elong:' then begin 
      readf,1,val
      elong=float(val)
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'triang_upper:' then begin 
      readf,1,val
      triang_upper=float(val)
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'triang_lower:' then begin 
      readf,1,val
      triang_lower=float(val)
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif

    ;Load plasma parameters
    if val eq 'main_ion:' then begin 
      readf,1,val
      main_ion=strtrim(val,1)
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'n_impur:'then begin 
      readf,1,val
      n_impur=fix(val)
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'impur_table:'then begin 
      readf,1,val
      if n_impur gt 0 or n_impur le 7 then begin
        impur_table=strarr(n_impur,3)
        for j=0, 2 do begin
          readf,1,val
          impur_table(*,j)=strtrim(((strsplit(val,':',/extract,/regex)))(1:*),2)
        endfor
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
      endif
    endif

    ; load calculation grid info
    if val eq 'code_grid_arr.Z:' then begin 
      readf,1,val
      code_grid_arr.z=float(strsplit(val,', ',/extract))
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'code_grid_arr.X:' then begin 
      readf,1,val
      code_grid_arr.x=float(strsplit(val,', ',/extract))
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'code_grid_arr.Y:' then begin 
      readf,1,val
      code_grid_arr.y=float(strsplit(val,', ',/extract))
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif

    ;load used input profiles electron density
    if val eq 'n_e_coord:' then begin 
      readf,1,val
      n_e_coord=fix(val)
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'n_e_raw_r:' then begin 
      readf,1,val
      n_e_raw_r=float(strsplit(val,', ',/extract))
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'n_e_raw:' then begin 
      readf,1,val
      n_e_raw=float(strsplit(val,', ',/extract))
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'n_e_raw_err:' then begin 
      readf,1,val
      n_e_raw_err=float(strsplit(val,', ',/extract))
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'n_e_r:' then begin 
      readf,1,val
      n_e_r=float(strsplit(val,', ',/extract))
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'n_e:' then begin 
      readf,1,val
      n_e=float(strsplit(val,', ',/extract))
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'n_e_err:' then begin 
      readf,1,val
      n_e_err=float(strsplit(val,', ',/extract))
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif


    ;load used input profiles electron temperature
    if val eq 't_e_coord:' then begin 
      readf,1,val
      t_e_coord=fix(val)
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 't_e_raw_r:' then begin 
      readf,1,val
      t_e_raw_r=float(strsplit(val,', ',/extract))
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 't_e_raw:' then begin 
      readf,1,val
      t_e_raw=float(strsplit(val,', ',/extract))
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 't_e_raw_err:' then begin 
      readf,1,val
      t_e_raw_err=float(strsplit(val,', ',/extract))
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 't_e_r:' then begin 
      readf,1,val
      t_e_r=float(strsplit(val,', ',/extract))
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 't_e:' then begin 
      readf,1,val
      t_e=float(strsplit(val,', ',/extract))
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 't_e_err:' then begin 
      readf,1,val
      t_e_err=float(strsplit(val,', ',/extract))
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif

    ;load input profiles of Z_eff 
    if val eq 'z_eff_coord:' then begin 
      readf,1,val
      z_eff_coord=fix(val)
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'z_eff_raw_r:' then begin 
      readf,1,val
      z_eff_raw_r=float(strsplit(val,', ',/extract))
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'z_eff_raw:' then begin 
      readf,1,val
      z_eff_raw=float(strsplit(val,', ',/extract))
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'z_eff_raw_err:' then begin 
      readf,1,val
      z_eff_raw_err=float(strsplit(val,', ',/extract))
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'z_eff_r:' then begin 
      readf,1,val
      z_eff_r=float(strsplit(val,', ',/extract))
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'z_eff:' then begin 
      readf,1,val
      z_eff=float(strsplit(val,', ',/extract))
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'z_eff_err:' then begin 
      readf,1,val
      z_eff_err=float(strsplit(val,', ',/extract))
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif


    ;Load output data 
    if val eq 'e_beam:' then begin
      readf,1,val
      e_beam=float(strsplit(val,', ',/extract))
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'z_beam:' then begin 
      readf,1,val
      z_beam=float(strsplit(val,', ',/extract))
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'x_beam:' then begin 
      readf,1,val
      x_beam=float(strsplit(val,', ',/extract))
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'y_beam:' then begin 
      readf,1,val
      y_beam=float(strsplit(val,', ',/extract))
      if strlen(val) gt 0 and strmid(val,1,1) ne ';' then st=st+1
    endif
    if val eq 'n_beam:' then begin 
      n_ebeam=n_elements(e_beam) 
      n_z=n_elements(z_beam)
      n_x=n_elements(x_beam)
      n_y=n_elements(y_beam) 
      n_beam=fltarr(n_ebeam,n_z,n_x,n_y)
      readu,1,n_beam
      if n_elements(n_beam) gt 100 then st=st+1
    endif
    if val eq 'exc_n2_frac:' then begin  
      n_z=n_elements(z_beam)
      n_x=n_elements(x_beam)
      n_y=n_elements(y_beam) 
      exc_n2_frac=fltarr(n_ebeam,n_z,n_x,n_y)
      readu,1,exc_n2_frac
      if n_elements(exc_n2_frac) gt 100 then st=st+1
    endif
    if val eq 'exc_n3_frac:' then begin  
      n_z=n_elements(z_beam)
      n_x=n_elements(x_beam)
      n_y=n_elements(y_beam) 
      exc_n3_frac=fltarr(n_ebeam,n_z,n_x,n_y)
      readu,1,exc_n3_frac
      if n_elements(exc_n2_frac) gt 100 then st=st+1
    endif
   if val eq 'vel_vec_x:' then begin
      n_bml=n_elements(x_bml)  
      n_x=n_elements(x_beam)
      n_y=n_elements(y_beam) 
      vel_vec_x=fltarr(n_x,n_y,n_bml)
      readu,1,vel_vec_x
      if n_elements(vel_vec_x) gt 100 then st=st+1
    endif
    if val eq 'vel_vec_y:' then begin
      n_bml=n_elements(x_bml)  
      n_x=n_elements(x_beam)
      n_y=n_elements(y_beam) 
      vel_vec_y=fltarr(n_x,n_y,n_bml)
      readu,1,vel_vec_y
      if n_elements(vel_vec_y) gt 100 then st=st+1
    endif
    if val eq 'vel_vec_coef:' then begin
      n_ebeam=n_elements(e_beam[where(e_beam gt 0)]) ;fix for data with halo
      n_z=n_elements(z_beam)
      n_x=n_elements(x_beam)
      n_y=n_elements(y_beam) 
      vel_vec_coef=fltarr(n_ebeam,n_z,n_x,n_y,9)
      readu,1,vel_vec_coef
      if n_elements(vel_vec_coef) gt 100 then st=st+1
    endif
    
    endwhile
;    if vel_dis_type eq 1 and st lt 72 then begin
;      Widget_control, status_wid, Get_Value=status_tx
;      Widget_Control, status_wid,$
;      Set_Value=[status_tx,[strtrim(string(Fix(status_tx(n_elements(status_tx)-1))+1),1)+' : Full set of output data is missing in the file']], Set_text_top_line=n_elements(status_tx)-4
;      st_err=1
;    endif
;    if vel_dis_type eq 0 and st ne 76 then begin
;      Widget_control, status_wid, Get_Value=status_tx
;      Widget_Control, status_wid,$
;      Set_Value=[status_tx,[strtrim(string(Fix(status_tx(n_elements(status_tx)-1))+1),1)+' : Full set of output data is missing in the file']], Set_text_top_line=n_elements(status_tx)-4
;      st_err=1
;    endif
    
    close,1
    beam_apertures=[[x_bml],[y_bml]]
    if vel_dis_type eq 'NO' then if st eq 77 then st_tot=1 else st_tot=0
    if vel_dis_type eq 'YES' then if st eq 80 then st_tot=1 else st_tot=0
  endif
 
  
  if keyword_set(mdsplus_run) and keyword_set(shot) then begin
     MDSOPEN,'DNB',shot,/Quiet, status=st0
     if st0 then begin
       beam='DNBI_ALCATOR'
       n_beam=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.RESULTS:n_beam',status=st1,/quiet)
       exc_n2_frac=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.RESULTS:exc_n2_frac',status=st2,/quiet)
       exc_n3_frac=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.RESULTS:exc_n3_frac',status=st3,/quiet)
       E_beam=MDSVALUE('Dim_Of(\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.RESULTS:n_beam,0)',status=st4,/quiet)
       z_beam=MDSVALUE('Dim_Of(\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.RESULTS:n_beam,1)',status=st5,/quiet)
       x_beam=MDSVALUE('Dim_Of(\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.RESULTS:n_beam,2)',status=st6,/quiet)
       y_beam=MDSVALUE('Dim_Of(\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.RESULTS:n_beam,3)',status=st7,/quiet)
       vel_vec_x=0 & vel_vec_y=0 & vel_vec_coef=0
       time_interv=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+':time_interv',status=st8,/quiet)
       time_stamp=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+':time_stamp',status=st8,/quiet)
       div_type=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+':div_type',status=st8,/quiet)
       atten_type=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+':atten_type',status=st8,/quiet)
       vel_dis_type=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+':vel_dis_type',status=st8,/quiet)
       if vel_dis_type eq 'YES' then begin
         vel_vec_x=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.RESULTS:vel_vec_x',status=st8,/quiet)
         vel_vec_y=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.RESULTS:vel_vec_y',status=st8,/quiet)
         vel_vec_coef=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.RESULTS:vel_vec_coef',status=st8,/quiet)      
       endif
       stop_plasma_type=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+':stop_plsm_cs',status=st8,/quiet)
       exc_plasma_type=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+':exc_plsm_cs',status=st8,/quiet)
       
       z_min=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.CODE_GRID:Z_MIN', status=st9,/quiet)
       z_step1=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.CODE_GRID:Z_STEP1', status=st10,/quiet)
       z_mid=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.CODE_GRID:Z_MID', status=st11,/quiet)
       z_step2=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.CODE_GRID:Z_STEP2', status=st12,/quiet)
       z_max=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.CODE_GRID:Z_MAX', status=st13,/quiet)
       x_min=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.CODE_GRID:X_MIN', status=st14,/quiet)
       x_step=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.CODE_GRID:X_STEP', status=st15,/quiet)
       x_max=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.CODE_GRID:X_MAX', status=st16,/quiet)
       y_min=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.CODE_GRID:Y_MIN', status=st17,/quiet)
       y_step=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.CODE_GRID:Y_STEP', status=st18,/quiet)
       y_max=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.CODE_GRID:Y_MAX', status=st19,/quiet)
      
       if st9*st10*st11*st12*st13*st14*st15*st16*st17*st18*st19 then code_grid_arr={z:[z_min,z_step1,z_mid,z_step2,z_max],x:[x_min,x_step,x_max],y:[y_min,y_step,y_max]}
       x_grid_focus=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.BEAM_GEOM:X_GRID_FOCUS', status=st20,/quiet)
       y_grid_focus=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.BEAM_GEOM:Y_GRID_FOCUS', status=st20,/quiet)      
       if (not(st20)) then begin
         x_grid_focus=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.BEAM_GEOM:GRID_FOCUS', status=st20,/quiet)
         y_grid_focus=x_grid_focus
       endif
       beam_port_phi=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.BEAM_GEOM:BEAM_PORT_PHI')
       beam_port=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.BEAM_GEOM:BEAM_PORT', status=st21,/quiet)
       R_grid=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.BEAM_GEOM:R_GRID', status=st22,/quiet)
       Z_grid=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.BEAM_GEOM:Z_GRID', status=st23,/quiet)       
       Phi_grid=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.BEAM_GEOM:PHI_GRID', status=st24,/quiet)      
       R_wall=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.BEAM_GEOM:R_wall', status=st25,/quiet)
       Z_wall=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.BEAM_GEOM:Z_wall', status=st26,/quiet)       
       Phi_wall=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.BEAM_GEOM:PHI_wall', status=st27,/quiet) 
       tank_front_dist=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.BEAM_GEOM:TANK_FRONT', status=st28,/quiet)
       tank_size=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.BEAM_GEOM:TANK_SIZE', status=st29,/quiet)
       neutr_size=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.BEAM_GEOM:NEUTR_SIZE', status=st30,/quiet)
       tank_magnet_dist=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.BEAM_GEOM:TANK_MAGNET', status=st31,/quiet)
       magnet_size=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.BEAM_GEOM:MAGNET_SIZE', status=st32,/quiet)
       tank_cal_dist=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.BEAM_GEOM:TANK_CAL', status=st33,/quiet)
       beam_apertures=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.BEAM_GEOM:BEAM_APERTUR', status=st34,/quiet)
       grid_ap_diam=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.BEAM_GEOM:GRID_AP_DIAM', status=st35,/quiet)     


       tank_pressure=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.NEUTR_GAS:TANK_P', status=st36,/quiet)
       torus_pressure=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.NEUTR_GAS:TORUS_P', status=st37,/quiet)
       
       duct_pressure=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.NEUTR_GAS:DUCT_P', status=st361,/quiet)
       duct_pressure_loc=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.NEUTR_GAS:DUCT_Z', status=st361,/quiet)
       if (not(st361)) then duct_pressure=0.0
       if (not(st361)) then duct_pressure_loc=0.0

       beam_atom=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.BEAM_PARAM:BEAM_ATOM', status=st38,/quiet)
       if (not(st38)) then beam_atom='H'
       E_full=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.BEAM_PARAM:E_FULL', status=st38,/quiet)
       E_frac=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.BEAM_PARAM:E_FRAC', status=st39,/quiet)
       I_beam=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.BEAM_PARAM:I_BEAM', status=st40,/quiet)
       I_frac=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.BEAM_PARAM:I_FRAC', status=st41,/quiet)      
       I_opt=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.BEAM_PARAM:I_opt', status=st42,/quiet)
       I_dens_par=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.BEAM_PARAM:I_dens_par', status=st43,/quiet)     
       x_div_bml_opt=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.BEAM_PARAM:XDIV_BML_OPT', status=st44,/quiet)
       y_div_bml_opt=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.BEAM_PARAM:YDIV_BML_OPT', status=st44,/quiet)
       if (not(st44)) then begin
         x_div_bml_opt=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.BEAM_PARAM:DIV_BML_OPT', status=st44,/quiet)
         y_div_bml_opt=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.BEAM_PARAM:DIV_BML_OPT', status=st44,/quiet)  
       endif

       div_dist_par=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.BEAM_PARAM:div_dist_par', status=st45,/quiet)   
    
       r_major=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.PLASMA_GEOM:R_MAJOR', status=st46,/quiet)
       z_major=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.PLASMA_GEOM:Z_MAJOR', status=st47,/quiet)
       r_minor=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.PLASMA_GEOM:R_MINOR', status=st48,/quiet)
       elong=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.PLASMA_GEOM:ELONG', status=st49,/quiet)
       triang_upper=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.PLASMA_GEOM:TRIANG_U', status=st50,/quiet)
       triang_lower=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.PLASMA_GEOM:TRIANG_L', status=st51,/quiet)

       n_limiters=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.BEAM_LIMITER:N_LIMITERS', status=st52,/quiet)
       limiters_table_1=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.BEAM_LIMITER:LIM_TABLE', status=st53,/quiet)
       if st53 and n_elements(limiters_table_1) gt 2 then begin
         if (size(limiters_table_1))(1) eq 6 then begin
           limiters_table=reform(strarr(7,n_limiters))
           limiters_table(*)='NAN'
           limiters_table(0:5,*)=limiters_table_1
         endif else limiters_table=limiters_table_1
       endif 

       n_e_coord=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.N_E_PROF:N_E_COORD', /quiet)
       n_e_raw=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.N_E_PROF:N_E_RAW', status=st54,/quiet)
       n_e_raw_r=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.N_E_PROF:N_E_RAW_R', status=st55,/quiet)
       n_e_raw_err=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.N_E_PROF:N_E_RAW_ER', status=st56,/quiet)
       n_e=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.N_E_PROF:N_E', status=st57,/quiet)
       n_e_r=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.N_E_PROF:N_E_R', status=st58,/quiet)
       n_e_err=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.N_E_PROF:N_E_ER', status=st59,/quiet)
       t_e_coord=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.T_E_PROF:T_E_COORD', /quiet)
       t_e_raw=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.T_E_PROF:T_E_RAW', status=st60,/quiet)
       t_e_raw_r=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.T_E_PROF:T_E_RAW_R', status=st61,/quiet)
       t_e_raw_err=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.T_E_PROF:T_E_RAW_ER', status=st62,/quiet)
       t_e=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.T_E_PROF:T_E', status=st63,/quiet)
       t_e_r=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.T_E_PROF:T_E_R', status=st64,/quiet)
       t_e_err=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.T_E_PROF:T_E_ER', status=st65,/quiet)
       z_eff_coord=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.Z_Eff_PROF:Z_Eff_COORD', /quiet)
       Z_eff_raw=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.Z_Eff_PROF:Z_Eff_RAW', status=st66,/quiet)
       Z_eff_raw_r=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.Z_Eff_PROF:Z_Eff_RAW_R', status=st67,/quiet)
       Z_eff_raw_err=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.Z_Eff_PROF:Z_Eff_RAW_ER', status=st68,/quiet)
       Z_eff=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.Z_Eff_PROF:Z_Eff', status=st69,/quiet)
       Z_eff_r=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.Z_Eff_PROF:Z_Eff_R', status=st70,/quiet)
       Z_eff_err=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.Z_Eff_PROF:Z_Eff_ER', status=st71,/quiet)

       main_ion=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.PLASMA_PARAM:MAIN_ION', status=st72,/quiet)
       n_impur=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.PLASMA_PARAM:N_IMPUR', status=st73,/quiet)
       impur_table=MDSVALUE('\DNB::TOP.ALCBEAM.'+mdsplus_run(i)+'.INPUT.PLASMA_PARAM:IMPUR_TABLE', status=st74,/quiet)  
         
       st_tot=st1*st2*st3*st4*st5*st6*st7*st8*st9*st10*st11*st12*st13*st14*st15*st16*st17*st18*st19*st20*st21*st22*st23*st24*st25*st26*st27*st28*st29*st30*st31*st32*st33*st34*st35*st36*st37*st38*st39* $
       st40*st41*st42*st43*st44*st45*st46*st47*st48*st49*st50*st51*st52*st53*st54*st55*st56*st57*st58*st59*st60*st61*st62*st63*st64*st65*st66*st67*st68*st69*st70*st71*st72*st73*st74
        
    endif else begin
      print,'Error: Wrong shot number'
      return
    endelse
  endif
 
  if (not(st_tot)) then begin
         print,'Error: Not sufficient input data'
         return
  endif
 

  beam_geom={beam_port_phi:beam_port_phi,beam_port:beam_port,r_grid:r_grid,z_grid:z_grid,phi_grid:phi_grid,r_wall:r_wall,z_wall:z_wall,phi_wall:phi_wall,tank_front_dist:tank_front_dist,tank_size:tank_size,neutr_size:neutr_size,tank_magnet_dist:tank_magnet_dist,magnet_size:magnet_size,tank_cal_DIST:tank_cal_dist,beam_apertures:beam_apertures,grid_ap_diam:grid_ap_diam}
  neutr_gas={tank_pressure:tank_pressure,torus_pressure:torus_pressure,duct_pressure:duct_pressure,duct_pressure_loc:duct_pressure_loc}
  beam_param={beam_atom:beam_atom,e_full:e_full,e_frac:e_frac,i_beam:i_beam,i_frac:i_frac,i_opt:i_opt,i_dens_par:i_dens_par,x_div_bml_opt:x_div_bml_opt,y_div_bml_opt:y_div_bml_opt,div_dist_par:div_dist_par}
  plasma_geom={r_major:r_major,z_major:z_major,r_minor:r_minor,elong:elong,triang_u:triang_upper,triang_l:triang_lower}
  beam_limiter={n_limiters:n_limiters,limiters_table:limiters_table}
  n_e_prof={n_e_coord:n_e_coord,n_e_raw:n_e_raw,n_e_raw_r:n_e_raw_r,n_e_raw_err:n_e_raw_err,n_e:n_e,n_e_r:n_e_r,n_e_err:n_e_err}
  t_e_prof={t_e_coord:t_e_coord,t_e_raw:t_e_raw,t_e_raw_r:t_e_raw_r,t_e_raw_err:t_e_raw_err,t_e:t_e,t_e_r:t_e_r,t_e_err:t_e_err}
  z_eff_prof={z_eff_coord:z_eff_coord,z_eff_raw:z_eff_raw,z_eff_raw_r:z_eff_raw_r,z_eff_raw_err:z_eff_raw_err,z_eff:z_eff,z_eff_r:z_eff_r,z_eff_err:z_eff_err}
  plasma_param={main_ion:main_ion,n_impur:n_impur,impur_table:impur_table}

  input={beam:beam,shot:shot,time_interv:time_interv,div_type:div_type,atten_type:atten_type,stop_plasma_type:stop_plasma_type,exc_plasma_type:exc_plasma_type,code_grid:code_grid_arr,beam_geom:beam_geom,neutr_gas:neutr_gas,beam_param:beam_param,plasma_geom:plasma_geom,beam_limiter:beam_limiter,n_e_prof:n_e_prof,t_e_prof:t_e_prof,z_eff_prof:z_eff_prof,plasma_param:plasma_param}  
  results={e_beam:e_beam,z_beam:z_beam,x_beam:x_beam,y_beam:y_beam,n_beam:n_beam,exc_n2_frac:exc_n2_frac,exc_n3_frac:exc_n3_frac,vel_vec_x:vel_vec_x,vel_vec_y:vel_vec_y,vel_vec_coef:vel_vec_coef}

  alcbeam_input(i)=ptr_new(input)
  alcbeam_results(i)=ptr_new(results)
endfor
end 

