pro get_profiles,chord,profiles
    @prof_param.idl

    case ne_coord of
	'rhopsi': n_e = interpol(double(ne_y),ne_x,chord.rho)
	'rmid': n_e = interpol(double(ne_y),ne_x,chord.rmid)
    endcase
    case te_coord of
	'rhopsi': T_e = interpol(te_y,te_x,chord.rho)*1d3
	'rmid': T_e = interpol(te_y,te_x,chord.rmid)*1d3
    endcase
    case zeff_coord of
	'rhopsi': Z_eff = interpol(zeff_y,zeff_x,chord.rho)
	'rmid': Z_eff = interpol(zeff_y,zeff_x,chord.rmid)
    endcase
    main_ion = 'D'
    case main_ion of
	'H': Z_main = 1
	'D': Z_main = 1
	'T': Z_main = 1
	'He': Z_main = 2
    endcase
    ions = [main_ion,impurities]

    sumZimp = total(imp_fr*imp_Z)
    sumZimp2 = total(imp_fr*imp_Z^2)

    ;Calculate the density of each ion using Z_eff, imp_Z, and imp_fr
    if sumZimp2 - sumZimp*Z_main ne 0 then n_imp = n_e * (Z_eff - Z_main)/(sumZimp2 - sumZimp*Z_main) else n_imp = 0
    n_main = (n_e - n_imp*sumZimp)/Z_main
    n_ion = [[n_main], [n_imp # imp_fr]]
    Z_ion = [Z_main, imp_Z]

    profiles = {n_e:n_e, T_e:T_e, Z_eff:Z_eff, main_ion:main_ion, ions:ions, n_ion:n_ion, Z_ion:Z_ion, T_i:T_e}
end
