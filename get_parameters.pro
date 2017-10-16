function get_parameters
    @common_param.idl

    version='1.0'
    param = {filename:filename,enable_bes:enable_bes,enable_cxrs:enable_cxrs,enable_edge:enable_edge,enable_brem:enable_brem,linemodel:linemodel,cs_effects:cs_effects,calczeeman:calczeeman,nonstatistical:nonstatistical,nonstat_file:nonstat_file,autocenter:autocenter,n_gen:n_gen}
    detector_param = {l_instr:l_instr,disp:disp,npix:npix,wve:wve,sens:sens,int_time:int_time,detstokes:detstokes}
    beam_param = {num_beams:num_beams,alcbeam_file:alcbeam_file,beam_ripple:beam_ripple}
    view_param = {ap_shape:ap_shape,ap_rad:ap_rad,ap_gridx:ap_gridx,ap_gridy:ap_gridy,sp_shape:sp_shape,sp_rad:sp_rad,sp_gridx:sp_gridx,sp_gridy:sp_gridy,spot_pos:spot_pos,optic_pos:optic_pos}
    mesh_param = {ds:ds,cut_r1:cut_r1,cut_r2:cut_r2,cut_z1:cut_z1,cut_z2:cut_z2,chordmodel:chordmodel,num_pini:num_pini}
    equil_param = {Bmodel:Bmodel,Rmajor:Rmajor,aminor:aminor,B_tor0:B_tor0,upperelong:upperelong,lowerelong:lowerelong,uppertri:uppertri,lowertri:lowertri,z0:z0,shafranov0:shafranov0,qprof:qprof,qprof_rho:qprof_rho,efit_file:efit_file,efit_shot:efit_shot,efit_time:efit_time,Ermodel:Ermodel,Er_r:Er_r,Er_val:Er_val,Er_diamagnetic:Er_diamagnetic,Er_vtheta:Er_vtheta,Er_vphi:Er_vphi}
    prof_param = {main_ion:main_ion,impurities:impurities,imp_z:imp_z,imp_fr:imp_fr,ne_coord:ne_coord,ne_x:ne_x,ne_y:ne_y,te_coord:te_coord,te_x:te_x,te_y:te_y,zeff_coord:zeff_coord,zeff_x:zeff_x,zeff_y:zeff_y}
    parameters = {version:version,param:param,detector_param:detector_param,beam_param:beam_param,view_param:view_param,mesh_param:mesh_param,equil_param:equil_param,prof_param:prof_param}
    return,parameters
end
