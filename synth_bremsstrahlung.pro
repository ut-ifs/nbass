;input:
;  wve:      A      positions to calculate bremsstrahlung
;  ds:       m      chord step size
;  t_e:      eV     electron temperature at chord points (separated by ds)
;  n_e:      cm^-3  electron density at chord points (separated by ds)
;  z_eff:           z_eff at chord points (separated by ds)
;return:
;  spectrum  ph/(A m^2 s sr)
function synth_bremsstrahlung,chord,profiles
    common detector_param,l_instr,disp,npix,wve,sens,int_time
    h = 6.62608d-34 ;[J*sec]
    c_light = 299792458d ;[m/s]
    coef = 1d-10/h/c_light ;[1/(Å-J)] [N_photons/Energy]
    g_ff = sqrt(3d)/!dPi*alog(1.81d-4*wve#profiles.T_e) ;From Foord and Marmar Nucl Fusion 25 (2), 187-202 (1985)
    ;Bremsstrahlung equation, W. J. Karzas and R. Latter, The Astrophysical Journal Supplement Series 55, 167 (1961).
    ;dW/(dV dt dnu dOmega) = (2^5*pi*e^6)/(3mc^3)*sqrt(2pi/3mkT)*Z^2*n^2*exp(-hnu/kT)*g_ff            [CGS]
    ;dW/(dV dt dnu dOmega) = (e^6)/(6*pi^2*epsilon0^3*mc^3)*sqrt(2pi/3mkT)*Z^2*n^2*exp(-hnu/kT)*g_ff  [SI]
    ;dph/(dV dt dlambda dOmega) = (e^6)/(6*pi^2*epsilon0^3*mc^3)*sqrt(2pi/3mkT)*Z^2*n^2*exp(-hnu/kT)*g_ff/lamda/h
    ;evaluating (e^6)/(6*pi^2*epsilon0^3*mc^3)*sqrt(2pi/(3m*{electronvolt}))/h*{centimeter^-6} = 7.62642D-15/(Å cm^3 s sr)
    ;exp(-hc/lambda/T_e) ~= exp(-12398.4 Å*eV / (Te*lamda))
    brem_calc = (1D0/wve)#(7.62742D-9*profiles.n_e^2*profiles.Z_eff/sqrt(profiles.T_e))*exp(-12398.4/(wve#profiles.T_e))*g_ff ;(ph/(Å m^3 s sr))
    ;y/sqrt(t_e) is invalid where t_e = 0 (outside plasma)
    brem = total(brem_calc,2,/nan)*chord.ds/chord.nray ;ph/(A m^2 s sr)
    return,brem
end
