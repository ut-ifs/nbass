;Estimate Zeeman pattern based on Blom and Jupen. (2002) PPCF 44 1229
;for hydrogenic ions

;B in tesla
;kT in eV
;wve in angstrom
pro blom,species,wvezeeman,s0,s1,s2,s3,B,vectH,vectV,kT,widthfactor
    ;calculate the Stokes vectors for the zeeman pi and sigma components
    ;pi lines are linearly polarized in the z direction
    pi_s0 = vectH[*,2]^2 + vectV[*,2]^2 ;should be exactly equal to pi_s0
    pi_s1 = vectH[*,2]^2 - vectV[*,2]^2
    pi_s2 = 2*vectH[*,2]*vectV[*,2] ;not sure about sign
    ;pi_s0_1 = sqrt(pi_s1^2+pi_s2^2) ;should be exactly equal to pi_s0
    ;pi_s0_2 = sin(psi)^2 ;should be exactly equal to pi_s0

    ;sigma lines are circularly polarized in the xy (magnetic coordinates) plane
    ;The y component (magnetic coordinate) is phase delayed by +- pi/2 from the x component
    ;So, the y component should be understood as imaginary; the complex conjugate will flip the sign
    ;Stokes s2 = 2*Re(Ex Ey*), where Ex and Ey are in the (viewing coordinates)
    ;Stokes s3 = -2*Im(Ex Ey*), where Ex and Ey are in the (viewing coordinates)

    sigma_s0 = 0.5d*((vectH[*,0]^2+vectH[*,1]^2) + (vectV[*,0]^2+vectV[*,1]^2))
    sigma_s1 = 0.5d*((vectH[*,0]^2+vectH[*,1]^2) - (vectV[*,0]^2+vectV[*,1]^2))
    sigma_s2 = (vectH[*,0]*vectV[*,0] + vectH[*,1]*vectV[*,1])
    sigma_s3p = -(-vectH[*,0]*vectV[*,1] + vectH[*,1]*vectV[*,0])
    sigma_s3m = (-vectH[*,0]*vectV[*,1] + vectH[*,1]*vectV[*,0])
    ;sigma_s0_1 = sqrt(sigma_s1^2+sigma_s2^2+sigma_s3p^2) ;should be exactly equal to sigma_s0
    ;sigma_s0_2 = 0.5 + cos(psi)^2/2                      ;should be exactly equal to sigma_s0
    ;sigma_s3check = cos(psi) ;should be exactly equal to sigma_s3p

    case species of
	             ;  wve0       alpha    beta     gamma  nlines
	'H_3_2': P = [6562.76d, 0.402267d, 0.3415d, -0.5247d, 54]
	'H_4_2': P = [4861.30d, 0.220724d, 0.2837d, -0.5346d, 54]
	'D_3_2': P = [6560.99d, 0.402068d, 0.4384d, -0.5015d, 54]
	'D_4_2': P = [4859.99d, 0.220610d, 0.3702d, -0.5132d, 54]
	'He3_4_3': P = [4685.91d, 0.205200d, 1.4418d, -0.4892d, 146]
	'He3_5_3': P = [3203.24d, 0.095879d, 1.2576d, -0.5001d, 146]
	'He3_6_4': P = [6560.38d, 0.401980d, 0.8976d, -0.4971d, 286]
	'He3_7_4': P = [5411.75d, 0.273538d, 0.8529d, -0.5039d, 286]
	'He4_4_3': P = [4685.70d, 0.205206d, 1.6118d, -0.4838d, 146]
	'He4_5_3': P = [3203.10d, 0.095879d, 1.4294d, -0.4975d, 146]
	'He4_6_4': P = [6560.09d, 0.401955d, 1.0058d, -0.4918d, 286]
	'He4_7_4': P = [5411.51d, 0.273521d, 0.9563d, -0.4981d, 286]
	'Be_5_4': P = [2530.50d, 0.060354d, 2.1245d, -0.3190d, 286]
	'Be_6_5': P = [4658.63d, 0.202754d, 1.6538d, -0.3192d, 474]
	'Be_7_5': P = [2906.22d, 0.078966d, 1.7017d, -0.3348d, 474]
	'Be_8_6': P = [4685.32d, 0.205025d, 1.4581d, -0.3450d, 710]
	'B_6_5': P = [2981.38d, 0.083423d, 2.0519d, -0.2960d, 474]
	'B_7_6': P = [4944.67d, 0.228379d, 1.6546d, -0.2941d, 710]
	'B_8_6': P = [2998.49d, 0.084065d, 1.8041d, -0.3177d, 710]
	'B_8_7': P = [7618.59d, 0.541883d, 1.4128d, -0.2966d, 994]
	'B_9_7': P = [4519.80d, 0.190781d, 1.5440d, -0.3211d, 994]
	'B_10_8': P = [6478.44d, 0.391914d, 1.3569d, -0.3252d, 1326]
	'C_6_5': P = [2070.28d, 0.040900d, 2.4271d, -0.2818d, 474]
	'C_7_6': P = [3433.69d, 0.110398d, 1.9785d, -0.2816d, 710]
	'C_8_6': P = [2082.17d, 0.040747d, 2.1776d, -0.3035d, 710]
	'C_8_7': P = [5290.56d, 0.261405d, 1.6689d, -0.2815d, 994]
	'C_9_7': P = [3138.65d, 0.092096d, 1.8495d, -0.3049d, 994]
	'C_10_8': P = [4498.82d, 0.189020d, 1.6191d, -0.3078d, 1326]
	'C_11_8': P = [3438.02d, 0.110428d, 1.6600d, -0.3162d, 1326]
	'C_11_9': P = [6199.64d, 0.359009d, 1.4464d, -0.3104d, 1706]
	'N_7_6': P = [2522.61d, 0.060010d, 2.4789d, -0.2817d, 710]
	'N_8_7': P = [3886.83d, 0.141271d, 2.0249d, -0.2762d, 994]
	'N_9_8': P = [5669.39d, 0.300127d, 1.7415d, -0.2753d, 1326]
	'N_10_8': P = [3305.16d, 0.102089d, 1.9464d, -0.2975d, 1326]
	'N_11_9': P = [4555.48d, 0.193799d, 1.7133d, -0.2973d, 1706]
	'O_8_7': P = [2975.76d, 0.083081d, 2.4263d, -0.2747d, 994]
	'O_9_8': P = [4340.52d, 0.176049d, 2.0652d, -0.2721d, 1326]
	'O_10_8': P = [2530.44d, 0.059933d, 2.3445d, -0.2944d, 1326]
	'O_10_9': P = [6068.27d, 0.343805d, 1.8122d, -0.2718d, 1706]
	'O_11_9': P = [3487.70d, 0.113640d, 2.0268d, -0.2911d, 1706]
	'Ne_9_8': P = [2777.79d, 0.072488d, 2.8838d, -0.2758d, 1326]
	'Ne_10_9': P = [3883.54d, 0.141002d, 2.4755d, -0.2718d, 1706]
	'Ne_11_9': P = [2232.00d, 0.046673d, 2.8410d, -0.2917d, 1706]
	'Ne_11_10': P = [5248.93d, 0.257292d, 2.1890d, -0.2715d, 2134]
    endcase

    ;alpha*B is the splitting between lambda(sigma+) and lambda(sigma-)
    wvezeeman = replicate(P[0],3,n_elements(B)) + ([0d,-0.5d,0.5d] # (B*P[1]))
    widthfactor = P[2]*kT^P[3]

    ;divide by 2 so that total intensity = 1
    s0 = dblarr(3,n_elements(B))
    s0[0,*] = pi_s0/2
    s0[1,*] = sigma_s0/2
    s0[2,*] = sigma_s0/2
    s1 = dblarr(3,n_elements(B))
    s1[0,*] = pi_s1/2
    s1[1,*] = sigma_s1/2
    s1[2,*] = sigma_s1/2
    s2 = dblarr(3,n_elements(B))
    s2[0,*] = pi_s2/2
    s2[1,*] = sigma_s2/2
    s2[2,*] = sigma_s2/2
    s3 = dblarr(3,n_elements(B))
    s3[0,*] = 0
    s3[1,*] = sigma_s3m/2 ;not sure about sign
    s3[2,*] = sigma_s3p/2
end
