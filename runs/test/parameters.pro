forward_function get_spot_position

;This procedure loads all parameters for an NBASS run into "common" variables
;Parameters are organized into groups:
;   param (general parameters)
;   detector_param (parameters relating to the spectrometer)
;   beam_param (parameters relating to the beam)
;   view_param (parameters relating to the view geometry)
;   mesh_param (parameters relating to the mesh)
;   equil_param (parameters relating to the magnetic equilibrium)
;   prof_param (parameters relating to plasma profiles)
;   hack_param (undocumented parameters used for specialized calculations)

;It is permissable to add keywords here. They will carry to the nbass call. e.g.
;nbass,'test',channel=5
;This can help in programmatically looping over a parameter value
;Since the parameter file is an idl procedure, calculations can be performed
pro parameters,channel=channel,filename=filename_
    if n_elements(channel) eq 0 then channel=4
    @../../common_param.idl ;This line sets up the common blocks

; variables in group param
;   filename:    string. filename for results save file
;   enable_bes:  boolean. output Stark BES spectrum
;   enable_cxrs: string. output CXRS spectrum
;   enable_edge: boolean. output emission from cold edge (placeholder; poorly modeled)
;   enable_brem: boolean. output bremsstrahlung
;   linemodel:   string.
;       'gaussian' use gaussian
;       'erf'      use error function in calculating the spectral lineshape
;       'delta'    use a dirac delta function
;   cs_effects:  boolean. make adjustments in CX cross section and width for finite plasma temperature
;   calczeeman: string.
;       '1'    use full Zeeman calculation for CXRS lines
;       'Blom' use Blom Zeeman calculation for CXRS lines
;       else   ignore Zeeman effect for CXRS lines
;   nonstatistical: boolean. use non-statistically populated beam excited states in Stark model
;   autocenter:  boolean. automatically center the wavelength range of the detector on the BES spectrum, overriding the wve parameter in detector_param
;   n_gen:  number of noisy spectra to generate from a Poisson distribution
    filename = 'test_'+strtrim(channel,2)+'.sav'
    enable_bes = 1
    enable_cxrs = 'D_3_2'
    enable_edge = 0 ;['D','C']
    enable_brem = 1
    linemodel = 'erf'
    cs_effects = 1  ;for CXRS
    calczeeman = 'Blom'  ;for CXRS
    nonstatistical = 1
    nonstat_file = 'runs/test/timemagneticlines.txt'
    autocenter = 0 ;should not be set when more than one beam in calculation
    n_gen = 10

; variables in group detector_param
;   l_instr:     double. instrument width in angstroms of spectrometer/detector, assumed to be Gaussian
;   disp:        double or double[npix]. dispersion in angstroms/pixel
;   npix:        int. number of pixels
;   wve:         double[npix]. wavelength of each pixel (arranged monotonically)
;   sens:        double or double[npix]. transmission * QE * etendue in cm^2/steradian
;   int_time:    double. integration time in seconds
;   detstokes:   4x4 matrix. multiply results by an arbitrary stokes matrix [optional]
;   detgain:     double. Number of detector counts per photon detected
;   darknoise:   double. Detector dark noise counts
    defvalue,l_instr, 0.6d  ;A     instrument function HW@1/e
    disp = 0.313d           ;A     much better than we get
    npix = 512
      wve_center = 6590       ;A
    wve = (dindgen(npix)-(npix-1d)/2d)*disp + wve_center
      ;Make some calculations for sens
      fiberd = 0.001d   ;m
      num_fibers = 18
      ;etendue of a fiber is pi^2/4*diameter^2*NA^2
      NA = 0.22d   ;numerical aperture of the fiber
      etendue = num_fibers*!dpi^2*(fiberd*100)^2*NA^2/4 ;cm^2*sr
      slitw = 50d-6     ;m
      slitfactor = (2*asin(slitw/fiberd)+2*slitw/fiberd*sqrt(1-(slitw/fiberd)^2))/!dpi
      qe = 0.4  ; guess quantum efficiency
      transmission = 0.2*slitfactor ;rough guess
    sens = qe*transmission*etendue
    int_time = 0.003d
    detstokes = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]
    detgain = 1.2d
    darknoise = 20d

; variables in group beam_param
;   num_beams:    number of neutral beams
;   alcbeam_file: string or string[num_beams]. filename or array of filenames of ALCBEAM output files for each neutral beam
;   beam_ripple:  double. beam_ripple = stddev(E_beam)/E_beam
    num_beams = 1
    ;alcbeam should be run beforehand to generate an alcbeam file
    alcbeam_file = 'runs/test/alcbeam_test.abo'
    beam_ripple = 0.005d ;raised cosine peak-peak amplitude with same stddev beam ripple: = 0.005*60kV/sqrt(1/3 - 2/pi^2)/2 = 0.4149kV

; variables in group view_param
;   ap_shape:    string. shape of aperture used to generate mesh points
;       'point'  single point
;       'hex7'   7 point hexagonal grid
;       'hex19'  19 point hexagonal grid
;       'grid'   arbitrarily specified grid
;   ap_rad:      double. for 'point','hex7','hex19' this is aperture radius (meters). for 'grid' this is size of grid point
;   ap_gridx:    double array. positions of each point in aperture grid in H,V coordinates (meter)
;   ap_gridy:    double array. positions of each point in aperture grid in H,V coordinates (meter)
;   sp_shape:    string. shape of spot used to generate mesh points
;       'point'  single point
;       'hex7'   7 point hexagonal grid
;       'hex19'  19 point hexagonal grid
;       'grid'   arbitrarily specified grid
;   sp_rad:      double. for 'point','hex7','hex19' this is spot radius (meters). for 'grid' this is size of grid point
;   sp_gridx:    double array. positions of each point in spot grid in H,V coordinates (meter)
;   sp_gridy:    double array. positions of each point in spot grid in H,V coordinates (meter)
;   spot_pos:    double[3]. XYZ position of spot (meters) in the plasma where view is focused
;   optic_pos:   double[3]. effective XYZ position of the optics (meters). Approximated by lens center position
    ap_shape = 'hex7'
    ap_gridx = [0] ;not used
    ap_gridy = [0] ;not used
    ap_rad = 0.070/2 ;meter
    sp_shape = 'grid'
    sp_gridx = [0.0,0.0,0.0,0.0,0.0,0.0, 9.51,9.51,9.51,9.51,9.51,9.51, -9.51,-9.51,-9.51,-9.51,-9.51,-9.51]/1000d
    sp_gridy = [24.705,13.725,2.745,-8.235,-19.215,-30.195, 30.195,19.215,8.235,-2.745,-13.725,-24.705, 30.195,19.215,8.235,-2.745,-13.725,-24.705]/1000d
    sp_rad = 0.01098/2 ;meter
    spot_pos = get_spot_position(channel)
    optic_pos = [-2.1452962d, 1.7867581d, -0.033d]

; variables in group mesh_param
;   ds:          double. step length (meters) along viewing chord for calculation mesh
;   cut_r1:      double. calculation boundary inner wall (meters)
;   cut_r2:      double. calculation boundary outer wall (meters)
;   cut_z1:      double. calculation boundary bottom (meters)
;   cut_z2:      double. calculation boundary top (meters)
;   chordmodel:  string. options for handling finite chord volume
;       'ray'        basic 1D ray, with analytical small-angle broadening applied. If set, ap_shape and sp_shape are ignored
;       'grid'       use ap_shape and sp_shape to determine chord mesh
;   num_pini:    integer. number of points to use in the beam source to model beam divergence
    ds = 0.01d ;default mesh step size, 1cm
    cut_r1 = 1.35  ;m
    cut_r2 = 2.35  ;m
    cut_z1 = -1.20 ;m
    cut_z2 = 1.16  ;m
    ;chordmodel = 'ray'  ;quick testing
    chordmodel = 'grid' ;full
    num_pini = 20

; variables in group equil_param
;   Bmodel:      string.
;       'miller'    use a Miller equilibrium with poloidal field generated from a known q profile
;       'efit_file' read efit data from a g file
;       'efit_mds'  read efit data from mdsplus. Requires some tokamak specific implementation
;  the following parameters are used with Bmodel='miller'
;   Rmajor:      double. Miller major radius (meters)
;   aminor:      double. Miller minor radius (meters)
;   B_tor0:      double. magnetic field at geometric center Rmajor (Tesla)
;   upperelong:  double. Miller elongation for upper half of plasma, assumed to be constant across flux surfaces
;   lowerelong:  double. Miller elongation for lower half of plasma, assumed to be constant across flux surfaces
;   uppertri:    double. Miller triangulation for upper half of plasma at r/a=0.95. model: delta(r) = delta95/0.95^2 * (r/a)^2
;   lowertri:    double. Miller triangulation for lower half of plasma at r/a=0.95. model: delta(r) = delta95/0.95^2 * (r/a)^2
;   z0:          double. z coordinate position of midplane (meters).
;   shafranov0:  double. Shafranov shift in center. model: Delta(r) = shafranov0*(1-(r/a)^2)
;   qprof:       double array. Table of q-values to use to generate poloidal field
;   qprof_rho:   double array. normalized poloidal flux associated with q values
;  the following parameters are used with Bmodel='efit_file'
;   efit_file:   string. name of a g-file to open
;  the following parameters are used with Bmodel='efit_mds'
;   efit_shot:   double. milliseconds
;   efit_time:   double. seconds
;
;   Ermodel:     string. 'none': Er=0. 'tabulated': given by Er_r,Er_val
;todo	'forcebalance': given by 1/(n*Z*e) * dp_i/dr - v_{theta,i}*B_phi + v_{phi,i}*B_theta  (c.f. Wesson 2004 4.19.4)
;todo	'romannikov': E_r(r) = -B_p(r)/c * int{xi:r->a}{[v_p(xi)*Bt/Bp(xi) - c/(e*n_e(xi)*B_p(xi)*dp_i/dxi) - j/(2*e*n_e(xi)) + e*(n_i-n_e)*c^2/j(xi)]*(1/xi+1/B_p(xi)*dB_p/dxi) dxi} + [v_p(r)*B_t/c - 1/(e*n_e(r))*dP_i(r)/dr]{r=a}
;   Er_r:        double[n]. array of Rmid locations where Er is given
;   Er_val:      double[n]. array of Er values (V/m)
;   Er_diamagnetic: double[n]. V/m = kg*m*s^-3*A-1 = m^3*C^-1*Pa/m
;   Er_vtheta:   double[n]. poloidal ion velocity used for calculating Er. m/s
;   Er_vphi:     double[n]. toroidal ion velocity used for calculating Er. m/s
    Bmodel = 'efit_file' ;select from 'miller', 'efit_mds', 'efit_file'
    ;parameters for option 'miller'
    Rmajor = 1.868 ;m
    aminor = 0.454 ;m
    B_tor0 = -2.256 ;m
    upperelong = 1.655 ;upper elongation, assumed to be a constant function of rho
    lowerelong = 1.655 ;lower elongation, assumed to be a constant function of rho
    uppertri = 0.523 ;upper triangularity, assumed to be a quadratic function of rho
    lowertri = 0.280 ;lower triangularity, assumed to be a quadratic function of rho
    z0 = 0.016 ;plasma center
    shafranov0 = 0.024 ;m (difference Rm-Rc), shift assumed to be a quadratic function of rho
    q95 = 6.534
    qprof_rho = interpol([0,1],100)
    qprof = 1.909 + qprof_rho^8.886*(q95-1.909)/0.95^8.886 ;for simplicity, use this basic shape
    ;parameters for option 'efit_mds'
    ;procedure called 'connect_efit_mdsplus' should also be defined
    efit_shot = 70079 ;only used for MDSplus EFIT retrieval
    efit_time = 3.3  ;(s), only used for MDSplus EFIT retrieval
    ;parameters for option 'efit_file'
    efit_file = 'runs/test/g070079.03300'

    Ermodel = 'none'
    Er_r = [1.4,3]
    Er_val = [0,0]
    Er_diamagnetic = [0,0]
    Er_vtheta = [0,0]
    Er_vphi = [0,0]

; variables in group prof_param
;   main_ion:    string.       should be 'H', 'D', or 'He'
;   impurities:  string[nimp]. array of atomic symbols for impurities to include
;   imp_z:       double[nimp]. charge number of each impurity. allowed to be non-integer for partially ionized species
;   imp_fr:      double[nimp]. fraction of total impurity density for each impurity
;   ne_coord:    string. 'rhopsi' or 'rmid'. Coordinate to use for ne.
;   ne_x:        double[].     array of ne measurement positions
;   ne_y:        double[].     array of ne measurements in cm^-3
;   te_coord:    string. 'rhopsi' or 'rmid'. Coordinate to use for Te.
;   te_x:        double[].     array of Te measurement positions
;   te_y:        double[].     array of Te measurements in keV
;   zeff_coord:  string. 'rhopsi' or 'rmid'. Coordinate to use for Zeff.
;   zeff_x:      double[].     array of Zeff measurement positions
;   zeff_y:      double[].     array of Zeff measurements
    main_ion = 'D'
    impurities = ['C'] ;array of all impurities present
    imp_z = [6.0] ;charge number of each impurity. allowed to be non-integer for partially ionized species
    imp_fr = [1.0] ;fraction of total impurity density for each impurity
    ;ne in cm^-3
    ne_coord = 'rhopsi'
    ne_x = [0.0000, 0.0100, 0.0200, 0.0300, 0.0400, 0.0500, 0.0600, 0.0700, 0.0800, 0.0900, 0.1000, 0.1100, 0.1200, 0.1300, 0.1400, 0.1500, 0.1600, 0.1700, 0.1800, 0.1900, 0.2000, 0.2100, 0.2200, 0.2300, 0.2400, 0.2500, 0.2600, 0.2700, 0.2800, 0.2900, 0.3000, 0.3100, 0.3200, 0.3300, 0.3400, 0.3500, 0.3600, 0.3700, 0.3800, 0.3900, 0.4000, 0.4100, 0.4200, 0.4300, 0.4400, 0.4500, 0.4600, 0.4700, 0.4800, 0.4900, 0.5000, 0.5100, 0.5200, 0.5300, 0.5400, 0.5500, 0.5600, 0.5700, 0.5800, 0.5900, 0.6000, 0.6100, 0.6200, 0.6300, 0.6400, 0.6500, 0.6600, 0.6700, 0.6800, 0.6900, 0.7000, 0.7100, 0.7200, 0.7300, 0.7400, 0.7500, 0.7600, 0.7700, 0.7800, 0.7900, 0.8000, 0.8100, 0.8200, 0.8300, 0.8400, 0.8500, 0.8600, 0.8700, 0.8800, 0.8900, 0.9000, 0.9100, 0.9200, 0.9300, 0.9400, 0.9500, 0.9600, 0.9700, 0.9800, 0.9900, 1.0000]
    ne_y = [1.600E+13,   1.600E+13,   1.599E+13,   1.599E+13,   1.597E+13,   1.596E+13,   1.594E+13,   1.592E+13,   1.590E+13,   1.587E+13,   1.584E+13,   1.581E+13,   1.577E+13,   1.573E+13,   1.569E+13,   1.564E+13,   1.559E+13,   1.554E+13,   1.548E+13,   1.542E+13,   1.536E+13,   1.529E+13,   1.523E+13,   1.515E+13,   1.508E+13,   1.500E+13,   1.492E+13,   1.483E+13,   1.475E+13,   1.465E+13,   1.456E+13,   1.446E+13,   1.436E+13,   1.426E+13,   1.415E+13,   1.404E+13,   1.393E+13,   1.381E+13,   1.369E+13,   1.357E+13,   1.344E+13,   1.331E+13,   1.318E+13,   1.304E+13,   1.290E+13,   1.276E+13,   1.261E+13,   1.247E+13,   1.231E+13,   1.216E+13,   1.200E+13,   1.184E+13,   1.167E+13,   1.151E+13,   1.133E+13,   1.116E+13,   1.098E+13,   1.080E+13,   1.062E+13,   1.043E+13,   1.024E+13,   1.005E+13,   9.850E+12,   9.650E+12,   9.446E+12,   9.240E+12,   9.030E+12,   8.818E+12,   8.602E+12,   8.382E+12,   8.160E+12,   7.934E+12,   7.706E+12,   7.474E+12,   7.238E+12,   7.000E+12,   6.758E+12,   6.514E+12,   6.266E+12,   6.014E+12,   5.760E+12,   5.502E+12,   5.242E+12,   4.978E+12,   4.710E+12,   4.440E+12,   4.166E+12,   3.890E+12,   3.610E+12,   3.326E+12,   3.040E+12,   2.750E+12,   2.458E+12,   2.162E+12,   1.862E+12,   1.560E+12,   1.254E+12,   9.456E+11,   6.336E+11,   3.184E+11,   0.000E+00]
    ;te in keV
    te_coord = 'rhopsi'
    te_x = [0.0000, 0.0100, 0.0200, 0.0300, 0.0400, 0.0500, 0.0600, 0.0700, 0.0800, 0.0900, 0.1000, 0.1100, 0.1200, 0.1300, 0.1400, 0.1500, 0.1600, 0.1700, 0.1800, 0.1900, 0.2000, 0.2100, 0.2200, 0.2300, 0.2400, 0.2500, 0.2600, 0.2700, 0.2800, 0.2900, 0.3000, 0.3100, 0.3200, 0.3300, 0.3400, 0.3500, 0.3600, 0.3700, 0.3800, 0.3900, 0.4000, 0.4100, 0.4200, 0.4300, 0.4400, 0.4500, 0.4600, 0.4700, 0.4800, 0.4900, 0.5000, 0.5100, 0.5200, 0.5300, 0.5400, 0.5500, 0.5600, 0.5700, 0.5800, 0.5900, 0.6000, 0.6100, 0.6200, 0.6300, 0.6400, 0.6500, 0.6600, 0.6700, 0.6800, 0.6900, 0.7000, 0.7100, 0.7200, 0.7300, 0.7400, 0.7500, 0.7600, 0.7700, 0.7800, 0.7900, 0.8000, 0.8100, 0.8200, 0.8300, 0.8400, 0.8500, 0.8600, 0.8700, 0.8800, 0.8900, 0.9000, 0.9100, 0.9200, 0.9300, 0.9400, 0.9500, 0.9600, 0.9700, 0.9800, 0.9900, 1.0000]
    te_y = [2.2000, 2.1998, 2.1991, 2.1980, 2.1965, 2.1945, 2.1921, 2.1892, 2.1859, 2.1822, 2.1780, 2.1734, 2.1683, 2.1628, 2.1569, 2.1505, 2.1437, 2.1364, 2.1287, 2.1206, 2.1120, 2.1030, 2.0935, 2.0836, 2.0733, 2.0625, 2.0513, 2.0396, 2.0275, 2.0150, 2.0020, 1.9886, 1.9747, 1.9604, 1.9457, 1.9305, 1.9149, 1.8988, 1.8823, 1.8654, 1.8480, 1.8302, 1.8119, 1.7932, 1.7741, 1.7545, 1.7345, 1.7140, 1.6931, 1.6718, 1.6500, 1.6278, 1.6051, 1.5820, 1.5585, 1.5345, 1.5101, 1.4852, 1.4599, 1.4342, 1.4080, 1.3814, 1.3543, 1.3268, 1.2989, 1.2705, 1.2417, 1.2124, 1.1827, 1.1526, 1.1220, 1.0910, 1.0595, 1.0276, 0.9953, 0.9625, 0.9293, 0.8956, 0.8615, 0.8270, 0.7920, 0.7566, 0.7207, 0.6844, 0.6477, 0.6105, 0.5729, 0.5348, 0.4963, 0.4574, 0.4180, 0.3782, 0.3379, 0.2972, 0.2561, 0.2145, 0.1725, 0.1300, 0.0871, 0.0438, 0.0000]
    ;zeff
    zeff_coord = 'rhopsi'
    zeff_x = [0.0,1.0] ;value will be interpolated over range
    zeff_y = [1.55,1.55]
    ;The impurity density is determined by using ne, zeff, imp_z, and imp_fr
end

pro connect_efit_mdsplus
    mdsconnect,'mds.ipp.ac.cn' ;= 202.127.204.12
    mdsopen,'efit_east',61973
end

; calculate the spot position for a channel
; This is an example where custom calculations are used within the parameter file
function get_spot_position,channel
    channels = ['1'    , '2'    , '3'    , '4'    , '5'    , '6'    , '7'    , '8'    , '9'    , '10'   , '11'   , '12'   , '13'   , '14'   , '15'   , '16'   , '17'   , '18'   ]
    radii =    [1.81500, 1.83765, 1.86029, 1.88294, 1.90559, 1.92824, 1.95088, 1.97353, 1.99618, 2.01882, 2.04147, 2.06412, 2.08676, 2.10941, 2.13206, 2.15471, 2.17735, 2.20000]
    chind = (where(channels eq channel))[0]
    spot_radius = radii[chind]

    optic_pos = [-2.1452962d, 1.7867581d, -0.033d]
    grid_pos = [-5.4613283d, 7.3322029d, 0]
    wall_pos = [-3.5115810d, 5.3011930d, 0]
    ;find the intersection of line and cylindar to get the spot position from a desired radius
    dx = wall_pos[0]-grid_pos[0]
    dy = wall_pos[1]-grid_pos[1]
    dr2 = dx^2 + dy^2
    r2 = spot_radius^2
    dd = grid_pos[0]*wall_pos[1]-wall_pos[0]*grid_pos[1]
    plus = dblarr(3)
    plus[0] = (dd*dy + sign(dy)*dx*sqrt(r2*dr2 - dd^2))/dr2
    plus[1] = (-dd*dx + abs(dy)*sqrt(r2*dr2 - dd^2))/dr2
    dplus = (plus[0] - optic_pos[0])^2 + (plus[1] - optic_pos[1])^2
    minus = dblarr(3)
    minus[0] = (dd*dy - sign(dy)*dx*sqrt(r2*dr2 - dd^2))/dr2
    minus[1] = (-dd*dx - abs(dy)*sqrt(r2*dr2 - dd^2))/dr2
    dminus = (minus[0] - optic_pos[0])^2 + (minus[1] - optic_pos[1])^2
    if dplus lt dminus then spot_pos=plus else spot_pos=minus
    print,'spot position:',spot_pos
    return,spot_pos
end
