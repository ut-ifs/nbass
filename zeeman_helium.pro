;Calculate the Zeeman spectrum for helium
;Written by Ken Liao

pro calc_probabilities,Az_m,Ap_m,Am_m,verbose=verbose
;calculate the transition probabilites between l,j,m basis states
;of course, when field is turned on, need to take superposition over two j levels for given l and m

    alpha = 7.2973525376d-3 ;fine structure constant
    c = 29979245800d            ; cm/s     speed of light in vacuum
    lambda0 = 4685.9d-8         ; cm
    a0 = 0.52917720859d-8       ; cm       Bohr radius

    ;squared radial integrals {int(RnlRn'l'r^3 dr)^2} in bohr radius^2, from B & S, p264
    ;(not using because use adjusted NIST computed probabilities)
    RR=dblarr(4,3)
    ;RR[0,1]=6d       ;R4031 = R4s3p
    ;RR[1,0]=29.8d    ;R4130 = R4p3s
    ;RR[1,2]=1.7d     ;R4132 = R4p3d
    ;RR[2,1]=57d      ;R4231 = R4d3p
    ;RR[3,2]=104.7    ;R4332 = R4f3d
    ;calculated myself with help from mathematica. the factor of 4 comes from 1/Z^2 for helium
    RR[0,1]=(5750784d/5764801d)^2*6d/4d              ;R4031 = R4s3p
    RR[1,0]=(14100480d/5764801d)^2*5d/4d             ;R4130 = R4p3s
    RR[1,2]=(5308416d/5764801d)^2*2d/4d              ;R4132 = R4p3d
    RR[2,1]=(7962624d/5764801d)^2*30d/4d             ;R4231 = R4d3p
    RR[3,2]=(63700992d/5764801d)^2*6d/7d/4d          ;R4332 = R4f3d

    ;l resolved transition rates go like YY
    ; == Sum_m1[4*pi*(integrate{Y(l1,m1)* Y(1,0) Y(l2,m1)})^2] / (2*l1+1)
    ; == Sum_m1[4*pi*(integrate{Y(l1,m1)* Y(1,0) Y(l2,m1+1)})^2] / (2*l1+1)
    ; == Sum_m1[4*pi*(integrate{Y(l1,m1)* Y(1,0) Y(l2,m1-1)})^2] / (2*l1+1)
    ; Like magic, when everything is calculated and summed, YY = max(l1,l2)/(2*l1+1)

    YY=dblarr(4,3)
    YY[0,1]=1d/1d
    YY[1,0]=1d/3d
    YY[1,2]=2d/3d
    YY[2,1]=2d/5d
    YY[3,2]=3d/7d

    ;Here are the nl resolved rates
    rate_calc = 32d/3d*!dpi^3*alpha*c/lambda0^3*RR*YY*a0^2

    ;The nlj resolved rates are gotten by solving the angular momentum addition between L and S
    ; l  m_l  s   m_s   fup    fdown         avg fup   avg fdown
    ; 0   0  1/2  1/2    1       0              1          0

    ; 1  -1  1/2  1/2   1/3     2/3
    ; 1   0  1/2  1/2   2/3     1/3            2/3        1/3
    ; 1   1  1/2  1/2   3/3     0/3

    ; 2  -2  1/2  1/2   1/5     4/5
    ; 2  -1  1/2  1/2   2/5     3/5
    ; 2   0  1/2  1/2   3/5     2/5            3/5        2/5
    ; 2   1  1/2  1/2   4/5     1/5
    ; 2   2  1/2  1/2   5/5     0/5

    ; 3  -3  1/2  1/2   1/7     6/7
    ; 3  -2  1/2  1/2   2/7     5/7
    ; 3  -1  1/2  1/2   3/7     4/7
    ; 3   0  1/2  1/2   4/7     3/7            4/7        3/7
    ; 3   1  1/2  1/2   5/7     2/7
    ; 3   2  1/2  1/2   6/7     1/7
    ; 3   3  1/2  1/2   7/7     0/7
    ;fup is fraction of additions that result in j=l+1/2      (clebsh-gordon coefficient squared)
    ;fdown is fraction of additions that result in j=l-1/2

    ;the pattern is clear: m resolved rates are calculated from 

    rate_calc2 = 32d/3d*!dpi^3*alpha*c/lambda0^3*RR*a0^2

    ;transitions are enumerated by 4 dimensions
    ;first is li
    ;second is 1 if ji=li+1/2, 0 if ji=li-1/2
    ;third is 1 if lf = li+1, 0 if lf = li-1
    ;fourth is 1 if jf=lf+1/2, 0 if jf=lf-1/2
    ;if transition is not allowed, transition probability is 0
    ;data from NIST tables
    gA=make_array(4,2,2,2,/double,value=0d) ;gA is degeneracy*transition amplitude for a LS transition
    gA[2,0,0,0]=3.755800D+08
    gA[1,1,0,1]=1.962800D+08
    gA[0,1,1,0]=1.958900D+07
    gA[1,0,0,1]=9.814000D+07
    gA[3,0,0,0]=1.236000D+09
    gA[2,1,0,1]=6.759600D+08
    gA[1,1,1,0]=2.225200D+06
    gA[2,0,0,1]=7.510400D+07
    gA[3,1,0,1]=1.765600D+09
    gA[3,0,0,1]=8.827800D+07
    gA[1,1,1,1]=2.002700D+07
    gA[0,1,1,1]=3.917600D+07
    gA[1,0,1,0]=1.112700D+07

    ;transitions are enumerated by 4 dimensions
    ;first is li
    ;second is 1 if ji=li+1/2, 0 if ji=li-1/2
    ;third is 1 if lf = li+1, 0 if lf = li-1
    ;fourth is 1 if jf=lf+1/2, 0 if jf=lf-1/2
    ;fifth is imi
	;imi :  0    1    2    3    4    5    6    7
	;mi  :-7/2 -5/2 -3/2 -1/2  1/2  3/2  5/2  7/2
    ;if transition is not allowed, transition probability is 0
    ;CG coefficients are for angular momentum addition between photon state (j=1,m_photon) and final state (jf,mf) to give initial state (ji,mi)
    CG0 =make_array(4,2,2,2,8,/double,value=0d) ;Clebsch-Gordan coefficients for m_photon=0
    CGp =make_array(4,2,2,2,8,/double,value=0d) ;Clebsch-Gordan coefficients for m_photon=-1
    CGm =make_array(4,2,2,2,8,/double,value=0d) ;Clebsch-Gordan coefficients for m_photon=+1
    Az_m=make_array(4,2,2,2,8,/double,value=0d) ;transition probability for z polarization
    Ap_m=make_array(4,2,2,2,8,/double,value=0d) ;transition probability for (x+iy)/sqrt(2) polarization
    Am_m=make_array(4,2,2,2,8,/double,value=0d) ;transition probability for (x-iy)/sqrt(2) polarization

    ;loop over all states with transitions and calculate the m weighting of transition probabilities
    for li=0,3 do begin
	for ilf=0,1 do begin
	    lf=li+2*ilf-1
	    if lf lt 0 or lf gt 2 then continue
	    for iji=(1-li)>0,1 do begin
		ji=li+iji-0.5d
		for ijf=(1-lf)>0,1 do begin
		    jf=lf+ijf-0.5d
		    if abs(jf - ji) gt 1 then continue
		    for imi=3-li,4+li do begin
			mi=imi-3.5d
			;z polarization case, see B&S 60.7, 60.11, 64.14, 64.15
			;initial means the upper state
			mf=mi ;delta m = 0 case
			if mf ge -jf and mf le jf then begin
			    case ji-jf of
				 1: CG0[li,iji,ilf,ijf,imi] = (ji^2-mi^2)/(ji*(1d0+2d0*jf))
				 0: CG0[li,iji,ilf,ijf,imi] = (mi^2)/(jf*(1d0+jf))
				-1: CG0[li,iji,ilf,ijf,imi] = (jf^2-mi^2)/(jf*(1d0+2d0*jf)) ;minus sign in front of root
;				 1: CG0[li,iji,ilf,ijf,imi] = (ji^2-mi^2)/((2d*ji+1d)*(2d*ji-1d))
;				 0: CG0[li,iji,ilf,ijf,imi] = (jf^2-mf^2)/((2d*jf+1d)*(2d*jf-1d))
;				-1: CG0[li,iji,ilf,ijf,imi] = (mi^2)/(2d*ji*(ji+1d))
			    endcase
			endif
			;x+iy polarization case
			mf=mi+1d
			if mf ge -jf and mf le jf then begin
			    case ji-jf of
				 1: CGp[li,iji,ilf,ijf,imi] = (jf-mi)*(ji-mi)/(2d0*ji*(1d0+2d0*jf))
				 0: CGp[li,iji,ilf,ijf,imi] = (1d0+jf+mi)*(jf-mi)/(2d0*jf*(1d0+jf))
				-1: CGp[li,iji,ilf,ijf,imi] = (jf+mi)*(1d0+jf+mi)/(2d0*jf*(1d0+2d0*jf))
;				 1: CGp[li,iji,ilf,ijf,imi] = (ji-mi)*(ji-mi-1d)/((2*ji+1)*(2*ji-1))
;				 0: CGp[li,iji,ilf,ijf,imi] = (ji+mi+2d)*(ji+mi+1d)/((2*ji+3)*(2*ji+1))
;				-1: CGp[li,iji,ilf,ijf,imi] = (ji+mi+1d)*(ji-mi)/(2d*ji*(ji+1d))
			    endcase
			endif
			;x-iy polarization case
			mf=mi-1d
			if mf ge -jf and mf le jf then begin
			    case ji-jf of
				 1: CGm[li,iji,ilf,ijf,imi] = (jf+mi)*(ji+mi)/(2d0*ji*(1d0+2d0*jf))
				 0: CGm[li,iji,ilf,ijf,imi] = (1d0+jf-mi)*(ji+mi)/(2d0*jf*(1d0+jf)) ;minus sign in front of root
				-1: CGm[li,iji,ilf,ijf,imi] = (jf-mi)*(1d0+jf-mi)/(2d0*jf*(1d0+2d0*jf))
;				 1: CGm[li,iji,ilf,ijf,imi] = (ji+mi)*(ji+mi-1d)/((2d*ji+1d)*(2d*ji-1d))
;				 0: CGm[li,iji,ilf,ijf,imi] = (ji-mi+2d)*(ji-mi+1d)/((2d*ji+3d)*(2d*ji+1d))
;				-1: CGm[li,iji,ilf,ijf,imi] = (ji+mi)*(ji-mi+1d)/(2d*ji*(ji+1d))
			    endcase
			endif
		    endfor ;imi
		    coef = rate_calc2[li,lf] * (2d0*jf + 1d0) / (2d0)
		    ;coef = rate_calc[li,lf] * (2*ji + 1d)
		    if keyword_set(verbose) then begin
;			print,'(li,ilf,iji,ijf)=',li,ilf,iji,ijf
			print,'For (li,ji,lf,jf)=',li,ji,lf,jf,' compare calculated sum over m to NIST'
			print,'calc,NIST,ratio:',total(CG0[li,iji,ilf,ijf,*])*coef,gA[li,iji,ilf,ijf],total(CG0[li,iji,ilf,ijf,*])*coef/gA[li,iji,ilf,ijf]
			print,'calc,NIST,ratio:',total(CGp[li,iji,ilf,ijf,*])*coef,gA[li,iji,ilf,ijf],total(CGp[li,iji,ilf,ijf,*])*coef/gA[li,iji,ilf,ijf]
			print,'calc,NIST,ratio:',total(CGm[li,iji,ilf,ijf,*])*coef,gA[li,iji,ilf,ijf],total(CGm[li,iji,ilf,ijf,*])*coef/gA[li,iji,ilf,ijf]
			print,total(CG0[li,iji,ilf,ijf,*]),total(CGp[li,iji,ilf,ijf,*]),total(CGm[li,iji,ilf,ijf,*])
		    endif
		    ;scale the transition probabilities such that the sum over m equals to NIST values
;		    coef = gA[li,iji,ilf,ijf]/total(Az_m[li,iji,ilf,ijf,*])
		    if gA[li,iji,ilf,ijf] gt 0d then begin
			;Az_m[li,iji,ilf,ijf,*] = CG0[li,iji,ilf,ijf,*]*gA[li,iji,ilf,ijf]/total(CG0[li,iji,ilf,ijf,*])
			;Ap_m[li,iji,ilf,ijf,*] = CGp[li,iji,ilf,ijf,*]*gA[li,iji,ilf,ijf]/total(CGp[li,iji,ilf,ijf,*])
			;Am_m[li,iji,ilf,ijf,*] = CGm[li,iji,ilf,ijf,*]*gA[li,iji,ilf,ijf]/total(CGm[li,iji,ilf,ijf,*])
			Az_m[li,iji,ilf,ijf,*] = CG0[li,iji,ilf,ijf,*]*coef
			Ap_m[li,iji,ilf,ijf,*] = CGp[li,iji,ilf,ijf,*]*coef
			Am_m[li,iji,ilf,ijf,*] = CGm[li,iji,ilf,ijf,*]*coef
		    endif
		    ;note: total(Ap) = total(Am) = 2*total(Az) due to symmetry
		endfor ;ijf
	    endfor ;iji
	endfor ;ilf
    endfor ;li
end

pro zeeman_helium,B,invwve,pol,ampz,ampp,ampm,dbg=dbg,species=species,nist=nist
defvalue,nist,1
;calculate Zeeman splitting of hydrogenic helium n=4->3 line
if n_elements(B) eq 0 then B = 5.3
if n_elements(dbg) eq 0 then dbg = 0

;value of Bohr magneton/hc in inverse meters per Tesla, from NIST database of constants (1/19/2011)
bohr = 46.6864515D ;m^-1 T^-1
rydbergC = 10973731.568539 ;m^-1 (updated 6/10/2013)
alpha = 7.2973525698d-3 ;fine structure constant (updated 6/10/2013)

Z=2

case species of
    'He3': mhemeratio = 5495.885275450d ;ratio of helion mass to electron mass, NIST Jun 18 2012
    'He4': mhemeratio = 7294.299536129d ;ratio of alpha mass to electron mass, NIST Jun 18 2012
endcase
nistadjust = (1+1/7294.299536129d)/(1+1/mhemeratio) ;isotope adjustment to NIST energy levels

rydbergreduced = rydbergC/(1d + 1d/mhemeratio)

;value of Landé g-factor is
;g_L = 1-1/mhemeratio ;reduced mass effect on g_L, Bethe and Salpeter p214
g_L = 1d/(1d + 1d/mhemeratio)
g_S = 2.0023193043622d ;electron g-factor, NIST
;g_J = g_L*(J*(J+1)-S*(S+1)+L*(L+1))/(2*J*(J+1)) + g_S*(J*(J+1)+S*(S+1)-L*(L+1))/(2*J*(J+1))
;g_J ~= (j+1/2)/(l+1/2)

;the weak field case
;ds_z = Bmag*bohr*m_j*lande
;the strong field case
;ds_pb = Bmag*bohr*(m_l+g_S*m_s)
;the separation between weak and strong field occurs at H_B ~= H_LS
;Bmag*bohr*(m_l+g_S*m_s) ~=
;H_LS ~ (Z/n^3)*3.6e-4 eV
;B >>? Z^4/(n^3*l*(l+1))*12.5T
;for B ~ 5.3T, Z=2, n=4
;1.696 >>? 1/(l*(l+1))
;unfortunately, for l=1, strong field approximation is not so good, and a full model would be more accurate.

;what about stark effect? Let's compare stark effect splitting magnitude to zeeman effect.
;    c = 29979245800d            ; cm/s     speed of light in vacuum
;    ec = 4.8032068d-10          ; esu
;    a0 = 0.52917720859d-8       ; cm       Bohr radius
;    eps = 3*ec*a0/Z*vperp/c*Bmag       ;epsilon = Stark energy splitting coefficient
;  for Bmag = 5.3T, vperp = vthermal = sqrt(2T/m_he3), T = 1keV,
; Stark splitting ~= 9E-5 eV
; LS splitting ~= m_e*c^2*(Z*alpha)^4/(2*n^3)*(2/3) = 3E-4 eV
; B splitting ~= bohr*B = 3E-4 eV
; So, Stark splitting (for thermal ions) is probably small enough to treat as a perturbation (just some additional broadening)
; (for fast ions, the velocity can be much higher, but, it gets complicated)

;n=4 states, l=0,1,2,3. s=+-1/2

;Designate states with the following quantum numbers: n, l, m, o
; where m is m_j and o is a new quantum number which distinguishes the two spin states, label them l and u
; In the limit of no B field, o=0 corresponds to j=l-1/2, o=1 to j=l+1/2
;n=3
;index, n, l, m, o
; 0, 3, 0, -1/2, xl
; 1, 3, 0, +1/2, xu
; 2, 3, 1, -3/2, xl
; 3, 3, 1, +3/2, xu
; 4, 3, 2, -5/2, xl
; 5, 3, 2, +5/2, xu

; 6, 3, 1, -1/2, l
; 7, 3, 1, -1/2, u
; 8, 3, 1, +1/2, l
; 9, 3, 1, +1/2, u
;10, 3, 2, -3/2, l
;11, 3, 2, -3/2, u
;12, 3, 2, -1/2, l
;13, 3, 2, -1/2, u
;14, 3, 2, +1/2, l
;15, 3, 2, +1/2, u
;16, 3, 2, +3/2, l
;17, 3, 2, +3/2, u
l3=[0,0,1,1,2,2,1,1,1,1,2,2,2,2,2,2,2,2]
m3=[-.5d,.5d,-1.5d,1.5d,-2.5d,2.5d,-.5d,-.5d,.5d,.5d,-1.5d,-1.5d,-.5d,-.5d,.5d,.5d,1.5d,1.5d]
o3=[0,0,0,0,0,0,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1]
W3=make_array(3,6,2,/double,value=!values.d_nan) ;field-free energies (m^-1), 36 space for 18 states
E3=make_array(3,6,2,/double,value=!values.d_nan) ;state energies (m^-1), 36 space for 18 states
a3=make_array(3,6,2,/double,value=0d)
b3=make_array(3,6,2,/double,value=0d)

;n=4
;index, n, l, m, o, allowed transitions
; 0, 4, 0, -1/2, xl,
; 1, 4, 0, +1/2, xu,
; 2, 4, 1, -3/2, xl,
; 3, 4, 1, +3/2, xu,
; 4, 4, 2, -5/2, xl,
; 5, 4, 2, +5/2, xu,
; 6, 4, 3, -7/2, xl,
; 7, 4, 3, +7/2, xu,

; 8, 4, 1, -1/2, l ,
; 9, 4, 1, -1/2, u ,
;10, 4, 1, +1/2, l ,
;11, 4, 1, +1/2, u ,
;12, 4, 2, -3/2, l ,
;13, 4, 2, -3/2, u ,
;14, 4, 2, -1/2, l ,
;15, 4, 2, -1/2, u ,
;16, 4, 2, +1/2, l ,
;17, 4, 2, +1/2, u ,
;18, 4, 2, +3/2, l ,
;19, 4, 2, +3/2, u ,
;20, 4, 3, -5/2, l ,
;21, 4, 3, -5/2, u ,
;22, 4, 3, -3/2, l ,
;23, 4, 3, -3/2, u ,
;24, 4, 3, -1/2, l ,
;25, 4, 3, -1/2, u ,
;26, 4, 3, +1/2, l ,
;27, 4, 3, +1/2, u ,
;28, 4, 3, +3/2, l ,
;29, 4, 3, +3/2, u ,
;30, 4, 3, +5/2, l ,
;31, 4, 3, +5/2, u ,
;states 0-7 are diagonal in nlmj, and do not need o to distinguish anything
;x is prepended to remind me to calculate this differently
;as expected, there are 32 states
l4=[0,0,1,1,2,2,3,3,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3]
m4=[-.5d,.5d,-1.5d,1.5d,-2.5d,2.5d,-3.5d,3.5d,-.5d,-.5d,.5d,.5d,-1.5d,-1.5d,-.5d,-.5d,.5d,.5d,1.5d,1.5d,-2.5d,-2.5d,-1.5d,-1.5d,-.5d,-.5d,.5d,.5d,1.5d,1.5d,2.5d,2.5d]
o4=[0,0,0,0,0,0,0,0,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1]
W4=make_array(4,8,2,/double,value=!values.d_nan) ;field-free energies (m^-1), 64 spaces for 32 states
E4=make_array(4,8,2,/double,value=!values.d_nan) ;state energies (m^-1), 64 spaces for 32 states
a4=make_array(4,8,2,/double,value=0d)
b4=make_array(4,8,2,/double,value=0d)

n=3
if keyword_set(nist) then begin
    w3 = make_array(3,6,2,value=!values.d_nan)
    w3[0,2:3,1] = nistadjust*390140.964175d2 ;3S1/2
    w3[1,2:3,0] = nistadjust*390140.824628d2 ;3P1/2
    w3[1,1:4,1] = nistadjust*390142.560089d2 ;3P3/2
    w3[2,1:4,0] = nistadjust*390142.55723373d2 ;3D3/2
    w3[2,0:5,1] = nistadjust*390143.13566549d2 ;3D5/2
endif else begin
    for l=0,2 do begin
	for mi=2-l,3+l do begin
	    ;mi :  0    1    2    3    4    5
	    ;m  :-5/2 -3/2 -1/2  1/2  3/2  5/2
	    if mi gt 2-l and mi lt 3+l then W3[l,mi,0] = -rydbergreduced*alpha^2*Z^4/n^3*(1d/l-3d/(4d*n)) - rydbergreduced*Z^2/3d^2
					    W3[l,mi,1] = -rydbergreduced*alpha^2*Z^4/n^3*(1d/(l+1d)-3d/(4d*n)) - rydbergreduced*Z^2/3d^2
	endfor
    endfor
endelse
for l=0,2 do begin
    for mi=2-l,3+l do begin
	;mi :  0    1    2    3    4    5
	;m  :-5/2 -3/2 -1/2  1/2  3/2  5/2
	m = mi-2.5d
	if mi gt 2-l and mi lt 3+l then begin ;in this case, there is mixing between j+1/2 and j-1/2 states
	    deltaW = W3[l,mi,1]-W3[l,mi,0] ;m^-1
	    xi = B*bohr/deltaW ; dimensionless
	    avgW = (W3[l,mi,1]+W3[l,mi,0])/2
	    commode = avgW + deltaW*(xi*m)
	    diffmode = 0.5d*deltaW*sqrt(1d +xi*4d*m/(2*l+1d)+xi^2)
	    E3[l,mi,0] = commode-diffmode
	    E3[l,mi,1] = commode+diffmode
	    gam=(1d +xi*2d*m/(2d*l+1))/sqrt(1+ xi*4d*m/(2d*l+1)+xi^2)
	    a3[l,mi,0] = -sqrt(0.5d*(1d -gam))
	    b3[l,mi,0] = sqrt(0.5d*(1d +gam))
	    a3[l,mi,1] = sqrt(0.5d*(1d +gam))
	    b3[l,mi,1] = sqrt(0.5d*(1d -gam))
	    ;squared values
;	    a3[l,mi,0] = 0.5d*(1d -gam)
;	    b3[l,mi,0] = 0.5d*(1d +gam)
;	    a3[l,mi,1] = 0.5d*(1d +gam)
;	    b3[l,mi,1] = 0.5d*(1d -gam)
	endif else begin ;else, the element is already diagonal, so add anomolous zeeman effect
	    j = abs(m)
	    g = g_L*(J*(J+1d)-3d/4d +L*(L+1d))/(2d*J*(J+1d)) + g_S*(J*(J+1d)+3d/4d -L*(L+1d))/(2d*J*(J+1d))
	    E3[l,mi,1] = W3[l,mi,1]+m*g*bohr*B
	    a3[l,mi,1] = 1d
	    b3[l,mi,1] = 0d
	endelse
    endfor
endfor

n=4
if keyword_set(nist) then begin
    w4 = make_array(4,8,2,value=!values.d_nan)
    w4[0,3:4,1] = nistadjust*411477.181412d2   ;4S1/2
    w4[1,3:4,0] = nistadjust*411477.1224141d2  ;4P1/2
    w4[1,2:5,1] = nistadjust*411477.8545580d2  ;4P3/2
    w4[2,2:5,0] = nistadjust*411477.853339d2   ;4D3/2
    w4[2,1:6,1] = nistadjust*411478.097366d2   ;4D5/2
    w4[3,1:6,0] = nistadjust*411478.096926d2   ;4F5/2
    w4[3,0:7,1] = nistadjust*411478.218936d2   ;4F7/2
endif else begin
    for l=0,3 do begin
	for mi=3-l,4+l do begin
	    ;mi :  0    1    2    3    4    5    6    7
	    ;m  :-7/2 -5/2 -3/2 -1/2  1/2  3/2  5/2  7/2
	    if mi gt 3-l and mi lt 4+l then W4[l,mi,0] = -rydbergreduced*alpha^2*Z^4/n^3*(1d/l-3d/(4d*n)) - rydbergreduced*Z^2/4d^2
	    W4[l,mi,1] = -rydbergreduced*alpha^2*Z^4/n^3*(1/(l+1d)-3d/(4d*n)) - rydbergreduced*Z^2/4d^2
	endfor
    endfor
endelse
for l=0,3 do begin
    for mi=3-l,4+l do begin
	;mi :  0    1    2    3    4    5    6    7
	;m  :-7/2 -5/2 -3/2 -1/2  1/2  3/2  5/2  7/2
	m = mi-3.5d
	if mi gt 3-l and mi lt 4+l then begin ;in this case, there is mixing between j+1/2 and j-1/2 states
	    deltaW = W4[l,mi,1]-W4[l,mi,0] ;m^-1
	    xi = B*bohr/deltaW ; dimensionless
	    avgW = (W4[l,mi,1]+W4[l,mi,0])/2
	    commode = avgW + deltaW*(xi*m)
	    diffmode = 0.5d*deltaW*sqrt(1d +xi*4d*m/(2*l+1d)+xi^2)
	    E4[l,mi,0] = commode-diffmode
	    E4[l,mi,1] = commode+diffmode
	    gam=(1d +xi*2d*m/(2d*l+1))/sqrt(1+ xi*4d*m/(2d*l+1)+xi^2)
	    a4[l,mi,0] = -sqrt(0.5d*(1d -gam))
	    b4[l,mi,0] = sqrt(0.5d*(1d +gam))
	    a4[l,mi,1] = sqrt(0.5d*(1d +gam))
	    b4[l,mi,1] = sqrt(0.5d*(1d -gam))
	    ;squared values
;	    a4[l,mi,0] = 0.5d*(1d -gam)
;	    b4[l,mi,0] = 0.5d*(1d +gam)
;	    a4[l,mi,1] = 0.5d*(1d +gam)
;	    b4[l,mi,1] = 0.5d*(1d -gam)
	endif else begin ;else, the element is already diagonal, so add anomolous zeeman effect
	    j = abs(m)
	    g = g_L*(J*(J+1d)-3d/4d +L*(L+1d))/(2d*J*(J+1d)) + g_S*(J*(J+1d)+3d/4d -L*(L+1d))/(2d*J*(J+1d))
	    E4[l,mi,1] = W4[l,mi,1]+m*g*bohr*B
	    a4[l,mi,1] = 1d
	    b4[l,mi,1] = 0d
	endelse
    endfor
endfor

;now need to handle transitions between these states
calc_probabilities,Az_m,Ap_m,Am_m

count=0
invwve=dblarr(146)
pol=intarr(146)
ampz=dblarr(146) ;z polarized
ampp=dblarr(146) ;(x+iy)/sqrt(2) polarized
ampm=dblarr(146) ;(x-iy)/sqrt(2) polarized
if dbg then print,'transition: li, mi, ji, lf, mf, jf, deltam, deltal, deltaj, energy'
for l=0,3 do begin
    for mi=3-l,4+l do begin
	;mi :  0    1    2    3    4    5    6    7
	;mf :       0    1    2    3    4    5
	;m  :-7/2 -5/2 -3/2 -1/2  1/2  3/2  5/2  7/2
	m = mi-3.5d
	for lf=abs(l-1),l+1<2,2 do begin ;allowed final l state
	    ilf=(lf-l)>0
	    for mf=(mi-2)>(2-lf),mi<(3+lf) do begin
		if mi gt 3-l and mi lt 4+l then iup=1 else iup=0 ; is there more than one initial spin state?
		if mf gt 2-lf and mf lt 3+lf then fup=1 else fup=0 ; is there more than one final spin state?
		for oi=1-iup,1 do begin
		    for o_f=1-fup,1 do begin
			;if abs(l+oi-lf-o_f) gt 1 then continue ;selection, delta j = -1, 0, 1
			invwve[count]=E4[l,mi,oi]-E3[lf,mf,o_f]
			pol[count]=mf-mi+1 ; pol={-1,0,1} for sigma-, pi, sigma+
			if dbg then print,'transition:',l,mi-3.5,l+oi-0.5,lf,mf-2.5,lf+o_f-0.5,pol[count],l-lf,l+oi-lf-o_f,invwve[count]
			case pol[count] of
;			  -1: ampm[count]=b4[l,mi,oi]*b3[lf,mf,o_f]*Am_m[l,0,ilf,0,mi] + b4[l,mi,oi]*a3[lf,mf,o_f]*Am_m[l,0,ilf,1,mi] + a4[l,mi,oi]*b3[lf,mf,o_f]*Am_m[l,1,ilf,0,mi] + a4[l,mi,oi]*a3[lf,mf,o_f]*Am_m[l,1,ilf,1,mi]
;			   0: ampz[count]=b4[l,mi,oi]*b3[lf,mf,o_f]*Az_m[l,0,ilf,0,mi] + b4[l,mi,oi]*a3[lf,mf,o_f]*Az_m[l,0,ilf,1,mi] + a4[l,mi,oi]*b3[lf,mf,o_f]*Az_m[l,1,ilf,0,mi] + a4[l,mi,oi]*a3[lf,mf,o_f]*Az_m[l,1,ilf,1,mi]
;			   1: ampp[count]=b4[l,mi,oi]*b3[lf,mf,o_f]*Ap_m[l,0,ilf,0,mi] + b4[l,mi,oi]*a3[lf,mf,o_f]*Ap_m[l,0,ilf,1,mi] + a4[l,mi,oi]*b3[lf,mf,o_f]*Ap_m[l,1,ilf,0,mi] + a4[l,mi,oi]*a3[lf,mf,o_f]*Ap_m[l,1,ilf,1,mi]
			  -1: ampm[count]=(b4[l,mi,oi]*b3[lf,mf,o_f]*sqrt(Am_m[l,0,ilf,0,mi]) + b4[l,mi,oi]*a3[lf,mf,o_f]*sqrt(Am_m[l,0,ilf,1,mi]) + a4[l,mi,oi]*b3[lf,mf,o_f]*sqrt(Am_m[l,1,ilf,0,mi]) + a4[l,mi,oi]*a3[lf,mf,o_f]*sqrt(Am_m[l,1,ilf,1,mi]))^2
			   0: ampz[count]=(b4[l,mi,oi]*b3[lf,mf,o_f]*sqrt(Az_m[l,0,ilf,0,mi]) + b4[l,mi,oi]*a3[lf,mf,o_f]*sqrt(Az_m[l,0,ilf,1,mi]) + a4[l,mi,oi]*b3[lf,mf,o_f]*sqrt(Az_m[l,1,ilf,0,mi]) + a4[l,mi,oi]*a3[lf,mf,o_f]*sqrt(Az_m[l,1,ilf,1,mi]))^2
			   1: ampp[count]=(b4[l,mi,oi]*b3[lf,mf,o_f]*sqrt(Ap_m[l,0,ilf,0,mi]) + b4[l,mi,oi]*a3[lf,mf,o_f]*sqrt(Ap_m[l,0,ilf,1,mi]) + a4[l,mi,oi]*b3[lf,mf,o_f]*sqrt(Ap_m[l,1,ilf,0,mi]) + a4[l,mi,oi]*a3[lf,mf,o_f]*sqrt(Ap_m[l,1,ilf,1,mi]))^2
			   else: stop
			endcase
			if dbg then print,'amplitude (z,p,m)',ampz[count],ampp[count],ampm[count]
			count++
		    endfor
		endfor
	    endfor
	endfor
    endfor
endfor
;print,'number of transitions:',count

if dbg then stop
end

pro testBdependence,Bmin,Bmax,n_B,species=species,angle=angle,nist=nist
;plot the transition line energy shifts versus B field
    if n_elements(Bmin) eq 0 then Bmin = 0d
    if n_elements(Bmax) eq 0 then Bmax = 10d
    if n_elements(n_B) eq 0 then n_B = 100
    if n_elements(species) eq 0 then species='He3'
    B=interpol([Bmin,Bmax],n_B)
    lines=dblarr(n_B,146)
    centroidz=dblarr(n_B)
    centroidp=dblarr(n_B)
    centroidm=dblarr(n_B)
    sze = dblarr(n_B,146)
    for i=0,n_B-1 do begin
	zeeman_helium,B[i],invwve,pol,ampz,ampp,ampm,species=species,nist=nist
	lines[i,*]=invwve
	centroidz[i] = total(ampz*invwve)/total(ampz)
	centroidp[i] = total(ampp*invwve)/total(ampp)
	centroidm[i] = total(ampm*invwve)/total(ampm)
	sze[i,*] = (ampz+ampp+ampm)
    endfor
    figure,1
    clf,/all
    ploti,xtitle='Magnetic field (T)',ytitle='Energy shift (meV)'
    colchoose=['ff0000'x,'00ff00'x,'0000ff'x]
    invmtomeV = 1.23984d-3
    for i=0,145 do begin
	oploti,B,lines[*,i]*invmtomev,color=colchoose[pol[i]+1]
    endfor
    ;invm0 = 2134149.8d
    invm0 = 0d

    figure,4
    clf,/all
    ploti,xtitle='Magnetic field (T)',ytitle='Wavelength (A)'
    colchoose=['ff0000'x,'00ff00'x,'0000ff'x]
    lambdas = 1d10/(lines+invm0)
    sze /= max(sze)
    threshold = 0.01
    for i=0,145 do begin
	;oploti,B,lambdas[*,i],color=colchoose[pol[i]+1]
	visible = where(sze[*,i] gt threshold, complement=invis)
	if visible[0] ne -1 then begin
	    lamvis = lambdas[*,i]
	    if invis[0] ne -1 then lamvis[invis] = !values.d_nan
	    ;oploti,B[visible],lambdas[visible,i],color=colchoose[pol[i]+1]
	    oploti,B,lamvis,color=colchoose[pol[i]+1],thick=0.4
	endif
	if invis[0] ne -1 then begin
	    laminv = lambdas[*,i]
	    if visible[0] ne -1 then lamvis[visible] = !values.d_nan
	    ;oploti,B[invis],lambdas[invis,i],color=colchoose[pol[i]+1],linestyle=2
	    oploti,B,laminv,color=colchoose[pol[i]+1],linestyle=2,thick=0.2
	endif
    endfor

    figure,2
    clf,/all
    ploti,xrange=[Bmin,Bmax],yrange=[min(invwve),max(invwve)],xtitle='Magnetic field (T)',ytitle='Energy shift (meV)' 
    oploti,B,centroidz*invmtomev,color=colchoose[0]
    oploti,B,centroidp*invmtomev,color=colchoose[1]
    oploti,B,centroidm*invmtomev,color=colchoose[2]
    if n_elements(angle) eq 0 then angle=1.2
    figure,3
    clf,/all
    centroidz *= sin(angle)^2
    centroidp *= (1d + cos(angle)^2)/2
    centroidm *= (1d + cos(angle)^2)/2
    centroid = (centroidz + centroidp + centroidm) / (sin(angle)^2 + (1d + cos(angle)^2)/2 + (1d + cos(angle)^2)/2)
    lambda = 1d10/(centroid+invm0)
    ploti,B,lambda,yrange=[min(lambda),max(lambda)],xtitle='Magnetic field (T)',ytitle='Wavelength (A)'
end

pro testpattern,B,nist=nist
    zeeman_helium,B,invwve,pol,ampz,ampp,ampm,species='He4',nist=nist

    lambda=1d10/invwve
    amp=ampz+ampp+ampm
    figure,1
    clf,/all
    ploti,xtitle='Wavelength (A)', ytitle='Intensity'
    colchoose=['ff0000'x,'00ff00'x,'0000ff'x]
    for i=0,145 do begin
	oploti,[lambda[i],lambda[i]],[0d,amp[i]],color=colchoose[pol[i]+1]
    endfor
end

pro testpattern2,B,width,wve,angle
    zeeman_helium,B,invwve,pol,ampz,ampp,ampm,species='He4'
    if n_elements(width) eq 0 then width=1.39d  ;approx instrument function
    if n_elements(angle) ne 0 then begin
	ampz *= sin(angle)^2
	ampp *= (1d + cos(angle)^2)/2
	ampm *= (1d + cos(angle)^2)/2
    end

    invm0 = 2134149.8d
    lambda=1d10/(invwve+invm0)
    amp=ampz+ampp+ampm
    low=min(lambda)
    high=max(lambda)
    delta=high-low
    np=2000
    if n_elements(wve) eq 0 then wve=interpol([low-delta*0.1d -1,high+delta*0.1d +1],np)
    intz=dblarr(np)
    intp=dblarr(np)
    intm=dblarr(np)
    ;window,xsize=1000,ysize=800
    ;cgplot,[0],[0],xrange=[min(wve),max(wve)],yrange=[0d,max(amp)],/nodata,title='!U3!NHe!U1+!N (n=4->3) Zeeman pattern',xtitle='wavelength (angstroms)',ytitle='transition rates (s!U-1!N)'
    figure,1
    clf,/all
    ploti,xtitle='wavelength (angstroms)',ytitle='transition rates (s!U-1!N)'
    colchoose=['ff0000'x,'00ff00'x,'0000ff'x]
    for i=0,145 do begin
	case pol[i] of
	  -1: intm += ampm[i]*exp(-(wve-lambda[i])^2/width^2)/sqrt(!dpi*width^2)
	   0: intz += ampz[i]*exp(-(wve-lambda[i])^2/width^2)/sqrt(!dpi*width^2)
	   1: intp += ampp[i]*exp(-(wve-lambda[i])^2/width^2)/sqrt(!dpi*width^2)
	   else: stop
	endcase
	;cgplot,[lambda[i],lambda[i]],[0d,amp[i]],color=colchoose[pol[i]+1],/overplot
	oploti,[lambda[i],lambda[i]],[0d,amp[i]],color=colchoose[pol[i]+1]
    endfor
    scale = max(amp)/max(intm+intz+intp)*0.8
    ;cgplot,wve,intm*scale,color=colchoose[0],/overplot
    ;cgplot,wve,intz*scale,color=colchoose[1],/overplot
    ;cgplot,wve,intp*scale,color=colchoose[2],/overplot
    ;cgplot,wve,(intm+intz+intp)*scale,/overplot
    oploti,wve,intm*scale,color=colchoose[0]
    oploti,wve,intz*scale,color=colchoose[1]
    oploti,wve,intp*scale,color=colchoose[2]
    oploti,wve,(intm+intz+intp)*scale,color='000000'x
end

pro blom,B,T,wve,width,he3=he3,sigmadopp=sigmadopp
;A Blom and C Jupen, Parameterization of the Zeeman effect for hydrogen-like spectra in high-temperature plasma, Plasma Phys. Control. Fusion 44 (2002) 1229-1241
;  Ion     n-n'  lambda   alpha       beta     gamma     lines
;  H I     2-3   6562.76  0.402267    0.3415   -0.5247   54
;  H I     2-4   4861.30  0.220724    0.2837   -0.5346   54
;  D I     2-3   6560.99  0.402068    0.4384   -0.5015   54
;  D I     2-4   4859.99  0.220610    0.3702   -0.5132   54
;  3He II  3-4   4685.91  0.205200    1.4418   -0.4892   146
;  3He II  3-5   3203.24  0.095879    1.2576   -0.5001   146
;  3He II  4-6   6560.38  0.401980    0.8976   -0.4971   286
;  3He II  4-7   5411.75  0.273538    0.8529   -0.5039   286
;  4He II  3-4   4685.70  0.205206    1.6118   -0.4838   146
;  4He II  3-5   3203.10  0.095879    1.4294   -0.4975   146
;  4He II  4-6   6560.09  0.401955    1.0058   -0.4918   286
;  4He II  4-7   5411.51  0.273521    0.9563   -0.4981   286
    if keyword_set(he3) then begin
	lambda = 4685.91d
	alpha = 0.205200d
	beta = 1.4418d
	gamma = -0.4892d
    endif else begin
	lambda = 4685.70d
	alpha = 0.205206d
	beta = 1.6118d
	gamma = -0.4838d
    endelse
    dlambda = alpha*B
    widthratio = beta*T^gamma
    ;width = sqrt(widthzeeman^2+widthdoppler^2) = widthdoppler*sqrt(widthratio^2+1)
    width=sigmadopp*sqrt(widthratio^2+1)
    ;pi2sigmaratio = 2d*sin(theta)^2/(1d +cos(theta)^2)
    wve = [lambda-dlambda/2,lambda,lambda+dlambda/2]
end

pro compareblomfull,B,T,theta,nudge=nudge
    mc2_He = 3727.379109d6 ; eV
    lambda0 = 4685.70d
    if n_elements(B) eq 0 then B=5.3d
    if n_elements(T) eq 0 then T=1500d
    if n_elements(nudge) eq 0 then nudge=2134079.0d

    ;doppler broadening width (standard deviation)
    sigmadopp = sqrt(T/mc2_He)*lambda0

    zeeman_helium,B,invwve,pol,ampz,ampp,ampm,species='He4'

    ;lambda=1d10/(invwve+2134083.5)
    lambda=1d10/(invwve+nudge)
    amp=ampz+ampp+ampm
    low=min(lambda)
    high=max(lambda)
    delta=high-low
    np=4000
    ;wve=interpol([low-delta*0.1d,high+delta*0.1d],np)
    wve=interpol([4666,4707],np)
    intz=dblarr(np)
    intp=dblarr(np)
    intm=dblarr(np)
    for i=0,145 do begin
	case pol[i] of
	  -1: intm += ampm[i]*exp(-(wve-lambda[i])^2/(2*sigmadopp^2))/sqrt(!dpi*2*sigmadopp^2)
	   0: intz += ampz[i]*exp(-(wve-lambda[i])^2/(2*sigmadopp^2))/sqrt(!dpi*2*sigmadopp^2)
	   1: intp += ampp[i]*exp(-(wve-lambda[i])^2/(2*sigmadopp^2))/sqrt(!dpi*2*sigmadopp^2)
	   else: stop
	endcase
    endfor
    window,0
    plot,[0],[0],xrange=[min(wve),max(wve)],yrange=[0d,max(intm+intz+intp)],charsize=2
    colchoose=['ff0000'x,'00ff00'x,'0000ff'x]
    oplot,wve,intm,color=colchoose[0]
    oplot,wve,intz,color=colchoose[1]
    oplot,wve,intp,color=colchoose[2]
    oplot,wve,intm+intz+intp

    scl = total(ampz)
    blom,B,T,wvel,width,sigmadopp=sigmadopp
    ampblom0 = exp(-(wve-wvel[0])^2/(2*width^2))/sqrt(!dpi*2*width^2)
    ampblom1 = exp(-(wve-wvel[1])^2/(2*width^2))/sqrt(!dpi*2*width^2)
    ampblom2 = exp(-(wve-wvel[2])^2/(2*width^2))/sqrt(!dpi*2*width^2)
    ampblomt = ampblom0+ampblom1+ampblom2
    colchoose=['00ffff'x,'ff00ff'x,'ffff00'x]
    oplot,wve,ampblom0*scl,color=colchoose[0]
    oplot,wve,ampblom1*scl,color=colchoose[1]
    oplot,wve,ampblom2*scl,color=colchoose[2]
    oplot,wve,ampblomt*scl

    window,1
    discrepancy = intm+intz+intp-ampblomt*scl
    plot,wve,discrepancy,title='difference between full zeeman and blom',charsize=2

end

pro finestruct,invwve
bohr = 46.6864515D ;m^-1 T^-1
rydbergC = 10973731.568527 ;m^-1
alpha = 7.2973525376d-3 ;fine structure constant
Z=2
mhe3meratio = 5495.8852765d ;ratio of masses, NIST
g_L = 1-1/mhe3meratio ;reduced mass effect on g_L, Bethe and Salpeter p214
g_S = 2.0023193043622d ;electron g-factor, NIST

W3=make_array(3,2,/double,value=!values.d_nan)
W4=make_array(4,2,/double,value=!values.d_nan)
n=3
for l=0,2 do begin
    for o=0,1 do begin
	W3[l,o] = -rydbergC*alpha^2*Z^4/n^3*(1d/(l+o)-3d/(4d*n))
    endfor
endfor
n=4
for l=0,3 do begin
    for o=0,1 do begin
	W4[l,o] = -rydbergC*alpha^2*Z^4/n^3*(1d/(l+o)-3d/(4d*n))
    endfor
endfor
;now need to handle transitions between these states
count=0
invwve=dblarr(16)
pol=intarr(16)
print,'transition: li, mi, ji, lf, mf, jf, deltam, deltal, deltaj, energy'
for l=0,3 do begin
    for oi=(1-l)>0,1 do begin
	for lf=abs(l-1),l+1<2,2 do begin ;allowed final l state
	    for o_f=(1-lf)>0,1 do begin
		if abs(l+oi-lf-o_f) gt 1 then continue ;selection, delta j = -1, 0, 1
		print,'transition:',l,l+oi-.5,lf,lf+o_f-.5,W4[l,oi]-W3[lf,o_f]
		invwve[count]=W4[l,oi]-W3[lf,o_f]
		count++
	    endfor
	endfor
    endfor
endfor
print,'number of transitions:',count
lambdas = 1d10/(invwve+2134149.8d)
figure,4
    for i=0,count-1 do begin
	plothl,lambdas[i],linestyle=2
    endfor
end
