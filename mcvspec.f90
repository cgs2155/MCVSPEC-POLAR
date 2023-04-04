!
!  MCVSPEC
!  ----
!
!  May 2022: Accretion column temperature/density profile calculator for magnetic CVs.
!
!  INPUT PARAMETERS:
!  PARAM(1) - B-field
!  PARAM(2) - Specific Accretion Rate
!  PARAM(3) - WD Mass in Solar Masses
!  PARAM(4) - Abundance
!  PARAM(5) - SIGMA_S, ratio of partial electron-ion pressures

! N.B: Large eqautions are often broken up into "parts" for readability

SUBROUTINE MCVSPEC(EAR,NE,PARAM,IFL,PHOTAR,PHOTER)

  IMPLICIT NONE

  INTEGER    NE,IFL

  INTEGER    VGRIDMAX,RGRIDMAX
! INTEGER    STABLE !global boolean to check for plasma instability, this constrains the lower bound of SIGMA

  PARAMETER  (VGRIDMAX=15000)
  PARAMETER  (RGRIDMAX=5)

  REAL     EAR(0:NE),PARAM(*),PHOTAR(NE),PHOTER(NE)
  REAL     METABUN

  REAL     M, B, SIGMA
  REAL     DISTNORM
  REAL     RHO(VGRIDMAX),P(VGRIDMAX),TK(VGRIDMAX),                &
       &               X(VGRIDMAX),NELEC(VGRIDMAX), SOLN(VGRIDMAX), &
       &   TAU(VGRIDMAX), PI_E(VGRIDMAX)

  INTEGER    J,VGRID,RGRID
  REAL, EXTERNAL :: func, func_prime


  REAL Msun, Rsun, G, mue, mmw, A, mu, k, mH, alpha, beta, gamma
  REAL M_3 , M_wd, R_wd, shock_height, vff, n_elec_shock, Tshock, ES0, MDOT0, R_m, Bm, L

  !Intermediate Variables
  REAL Lpart1, Lpart2

  REAL m_e, pi, h, e, E0, X__, chi, p_a
  REAL mbar, zbar, zbarsqr, g_b, R_

  common  /wdparams/ M_3, M_wd, R_wd, shock_height, vff, n_elec_shock, Tshock, ES0, MDOT0, R_m, Bm, L, p_a, SIGMA
  common /constants/ Msun, Rsun, G, mue, mmw, A, mu, k, mH, alpha, beta, gamma, X__, chi, zbar, R_, m_e, mbar
  !common /bool/ STABLE


  Msun     = 1.989100e+33 ! solar mass [gram]
  Rsun     = 6.9599e10 ! solar radius [cm]
  G        = 6.672590e-8 ! Gravitational constant
  mue      = 1.173 ! - mean number of electrons per baryons
  mmw      = 2.0    !- Mean molecular mass in white dwarf-
  mu       = 0.615  ! - mean molecular weight in accretion column (some papers use mu = 0.5 and Cropper 1999 did comparison)
  A        = 6.99e16 !used by Saxton
  k        = 1.380658e-16 ! -Boltzmann constant-
  mH       = 1.672623e-24 ! -Mass of hydrogen-
  alpha    = 2.0 ! -Thermodynamic Constant-
  beta     = 3.85 ! -Thermodynamic Constant-
  gamma    = 5.0/3.0  ! -Adidabatic Constant-
  e        = 1.60218e-19 ! -Charge of an Electron-
  E0       = 8.75419e-12 ! -Permissivity of Free Space-
  h        = 6.62607e-34 ! Pplanck's Constant-
  m_e      = 9.11e-28 ! -Electron Mass-
  pi       = 3.141592654
  X__      = 2.7e34 ! -Constant Based on Plasma Composition-
  chi      = 1.9091 ! -Constant that Depends on the Abundance-Weighted-Mean Charge of Ions-
  R_       = 8.314e7 ! -Gas Constant-
  mbar     = 1.2886*mH
  zbar     = 1.09987 ! -Mean Molecular Mass in Column-
  zbarsqr  = 1.391
  g_b      = 1.25 ! -Gaunt Factor-

!
! RGRID is the number of radial zones in the accretion column.
!
  RGRID = 1
!
! These are the input parameters:
!  Mdot(0) is the specific accretion rate on the central axis of the col
!  Rc is the radius of the column in 10^9 cm
!  M is the mass of the primary in solar masses
!  METABUN is the metal abundance: solar = 1.0
!  VGRID is the number of vertical layers in the accretion column
!  ANG is the viewing angle used to modify the albedo in degrees
!
  B       = PARAM(1) ! Surface B-field [10^6 G]
  MDOT0   = PARAM(2) ! mass accretion rate [g/cm2/s]
  M       = PARAM(3) ! WD mass in solar mass
  METABUN = PARAM(4) ! Abundance
  SIGMA   = PARAM(5) ! Ratio of Partial Electron and Ion Pressures
  VGRID   = 250 ! Number of vertical grids fixed to 50
!
! Flux normalization factor for APEC model as we define our flux norm as (R [km] / d [kpc])^2
! Note: PI*R_C^2/(4*PI*Distance^2) where R_C = accretion column radius [cm] and distance [cm]
!
  DISTNORM = 2.62511E-34

  DO J=1,NE
     PHOTER(J)=0.0
  ENDDO

! **This checks to make sure that there aren't too many vertical layers
! Checks that VGRID is not too large

  IF(VGRID.GT.VGRIDMAX)THEN
     PRINT *,'Too many vertical elements, max = ',VGRIDMAX
     VGRID = VGRIDMAX
  ENDIF

!           Calculating Neccesary WD Characteristics
  M_3  = (5.816*Msun)/(mmw**2.0) ! -Chandresekar Mass [gram]
  M_wd = M*Msun ! WD mass [grams]

  R_wd = Rsun*(0.0225/mmw)*SQRT(1.0-(M_wd/M_3)**(4.0/3.0))/((M_wd/M_3)**(1.0/3.0))  ! WD radius [cm]
  Bm   = (B*1e6)*((R_wd)**3) !Magnetic Moment

  Lpart1 = G*M_wd/R_wd !Intermediate Luminosity Component
  Lpart2 = (8.6e17*(MDOT0**(7./5.))*(M**(1./7.))*(((R_wd/(1.e9))**(9./5.)))*(B**(-4./5.))) !Intermediate Luminosity Component

  L    = Lpart1*Lpart2 !Accretion Luminosity

  R_m  = (2.75e8)*((M)**(1.0/7.0))*(((R_wd/(1.e9))**(-2.0/7.0)))*((L/1.0e33)**(-2.0/7.0))*((Bm/1.0e30)**(4.0/7.0)) !Magnetospheric Radius


  IF (R_wd.LT.(R_m)) THEN
    vff = ((2.*G*M_wd)*((1./(R_wd))-(1./R_m)))**(1./2.) ! Free Fall Velocity with R_m correction, Suleimanov (2016)
  ELSE
    vff = ((2.*G*M_wd)/(R_wd))**(1./2.) ! -Free Fall Velocity- at the shock height, Mukai 2017
  ENDIF


  p_a = MDOT0/vff ! -Density of Pre-Shock Flow- used in the normalization of density, Saxton 2005

  n_elec_shock = 4.0*7.013e23*MDOT0/vff !electron number density at shock

  TShock = (3./16.)* (m_e*(vff**2)/k)*(1/(1/SIGMA+1))*((zbar+(mbar/m_e))/zbar) ! -Temperature at Shock- Saxton 2005
  
  !the efficiency of cyclotron cooling relative to bremsstrahlung cooling, Saxton 2005
  ES0 =  ((vff/1e8)**4)*((p_a/1e-8)**(-1.85))* &
      & ((B/10)**2.85)*((1)**(-.425))*( 2.13e-16)* &
      & ((zbar+ (mbar/m_e))**3.85)/(g_b*(zbar**2.85)*zbarsqr*((1+1/SIGMA)**2))

  !finds the shock height of the WD
  CALL SHOCK_APPROX( (3./(4. *( 1.+1./SIGMA))), VGRIDMAX, TAU, PI_E)

! Checking WD characteristics calculation
!  write(*,*) M_3, M_wd, R_wd, shock_height, vff, MDOT0
!  write(*,*) "B [MG] = ", B, "  ES0 = ", ES0


  CALL MCVSPEC_SHOOTING(VGRID,VGRIDMAX,RHO,P,TK,X,NELEC, SOLN, TAU, PI_E)


! Calculate the magnetic field. This is done using eqn(10) of Wu, Chanmu
! and Shaviv (1994) ApJ 426, 664, inverted to determine B from the other
! parameters. KM note: I think B is in unit of MegaGauss.

!  PRINT '(A, F10.2)', 'White dwarf radius [10^7 cm] = ', R_wd*1e-7
!  PRINT '(A, F10.2)', 'Shock height [10^7 cm] (numerical) = ', X(VGRID)*1e-7
!  PRINT '(A, F10.2)', 'Shock height [10^7 cm] (approximate) = ', shock_height*1e-7
!  PRINT '(A, F10.2)', 'Tshock [keV] = ', TShock*8.618e-8
!  PRINT *,'T [keV] (numerical) = ', TK(VGRID)*8.618e-8
!  PRINT *,'Base T [kEV] = ', TK(1)*8.618e-8
 !  PRINT '(A, F10.2)', 'v_ff [10^8 cm/s] = ', vff*1e-8

! write(*,*) ES0, TK(VGRID), TSHOCK, NELEC(VGRID), X(VGRID), shock_height
  
  IF(X(VGRID).GT.(R_wd*.8)) THEN
    PRINT *,'WARNING WARNING WARNING'
    PRINT *,'SHOCK HEIGHT IS COMPARABLE TO WHITE DWARF RADIUS'
    PRINT *,'1D ACCRETION COLUMN IS NO LONGER VALID'
  ENDIF

!PRINT '(A, F10.2)',"R_m/R_wd = ", R_m/R_wd


 IF(R_wd.GT.(R_m)) THEN
  PRINT *,'WARNING WARNING WARNING'
  PRINT *,'Magnetospheric Radius is Less than White Dwarf Radius'
  PRINT *,'Must Adjust for Higher B or Lower MDOT'
  PRINT *,'CANNOT FIT WITH CURRENT PARAMETERS, IGNORING R_m'
ENDIF



  CALL MCVSPEC_APEC(VGRID,X,TK,NELEC,METABUN,    &
       &              DISTNORM,EAR,NE,IFL,PHOTAR)

  RETURN
END SUBROUTINE MCVSPEC

SUBROUTINE MCVSPEC_APEC(VGRID,X,TK,NELEC,METABUN,    &
     &                        DISTNORM,EAR,NE,IFL,PHOTAR)

  IMPLICIT NONE

  INTEGER    IFL,NE
  INTEGER    VGRID,J,L
  REAL     EAR(0:NE),PHOTAR(NE)!,ANG!,FACT
  REAL     PARAM1(3)
  REAL     METABUN,FLX(NE),FLXERR(NE)
  REAL     TK(VGRID),X(VGRID),NELEC(VGRID)
  REAL     KK,PI,NENH
  REAL     DISTNORM!,DIST,NORM_FACTOR

  PI = 3.141592654
  !
  ! Constant to convert T(Kelvin) to kT
  !
  KK  = 11.6048E6
  !
  ! Constant to convert Ne to Nh
  !
  NENH = 1.21
  !
  ! First zero the flux array
  !
  DO L=1,NE
     PHOTAR(L) = 0.0
  ENDDO
  !
  ! Main loop to label 100.
  ! Loop over each vertical element to calculate and add the flux into the
  ! array
  !
  DO 100 J=1,VGRID
     !
     ! Calculates the APEC spectrum for each vertical element on the energy
     ! grid passed into it
     ! APEC parameter set
     PARAM1(1) = TK(J)/KK
     PARAM1(2) = METABUN
     PARAM1(3) = 0.0

     IF(TK(J)/KK.LT.86.)THEN
        CALL APEC(EAR,NE,PARAM1,IFL,FLX,FLXERR)
     ELSE
        CALL BREMSS(EAR,NE,PARAM1,IFL,FLX,FLXERR) ! added by KM if kT > 86 keV as it exceeds APEC's temperature upper limit.
     ENDIF

     DO L=1,NE
        IF(J.EQ.1)THEN
            FLX(L) = DISTNORM*X(J)*((NELEC(J)**2)/NENH)*1.0E-14*FLX(L)
        ELSE
           FLX(L) = DISTNORM*(X(VGRID)-X(J-1))*((NELEC(J)**2)/NENH)*1.0E-14*FLX(L)
        ENDIF

     ENDDO
     !
     ! adds the flux for this vertical element to the existing flux in each
     ! energy bin
     !
     DO L=1,NE
        PHOTAR(L) = PHOTAR(L) + FLX(L)
     ENDDO

100 END DO

  RETURN
END SUBROUTINE MCVSPEC_APEC

!! Two Temperature Equations
SUBROUTINE MCVSPEC_SHOOTING(VGRID,VGRIDMAX, RHO,P,TK,X,NELEC,SOLN, TAU_, PI_E_)
  !shock temp is temperature at shock which is the same from ions and electrons
  common  /wdparams/ M_3, M_wd, R_wd, shock_height, vff, n_elec_shock, Tshock, ES0, MDOT0, R_m, Bm, L, p_a, SIGMA
  common /constants/ Msun, Rsun, G, mue, mmw, A, mu, k, mH, alpha, beta, gamma, X__, chi, zbar, R_, m_e, mbar

  INTEGER steps, r, VGRID, VGRIDMAX, JUMP

  REAL    soln(VGRIDMAX), taufinal(VGRIDMAX), piefinal(VGRIDMAX), MDOT0
  REAL    RHO(VGRID), X(VGRID), NELEC(VGRID)
  REAL    P(VGRID),TK(VGRID), TkeV(VGRID)
  REAL    TAU_(VGRIDMAX), PI_E_(VGRIDMAX)

  REAL k, m_e, mbar !constants need to be redeclared at times

  REAL xi, xi_init, n
  REAL tau_init, tau_f

  xi_init = 0
  xi = xi_init

  tau_init = TAU_(VGRIDMAX)
  tau_f    = .25

  n = (tau_f - tau_init)/float(VGRIDMAX)

  taufinal(1)  = TAU_(VGRIDMAX)
  piefinal(1)  = PI_E_(VGRIDMAX)
  soln(1)      = 0

  k        = 1.380658e-16
  m_e      = 9.11e-28
  mbar     = 2.1553419978e-24

  DO steps = 2,VGRIDMAX
     taufinal(steps) = TAU_(VGRIDMAX - steps+1)
     piefinal(steps) = PI_E_(VGRIDMAX - steps+1)
     xi   = xi+ n*(dxidtau(taufinal(steps), piefinal(steps)))
     soln(steps)     = xi ! -This is the 'non-normalized position grid-
  ENDDO

  taufinal(VGRIDMAX) = tau_f
  piefinal(VGRIDMAX) = PI_E_(1)

  JUMP = VGRIDMAX/VGRID

  DO r = 1, VGRID
    X(r)         = soln(r*JUMP)*shock_height ! -This is the position grid-
    RHO(r)       = p_a/taufinal(r*JUMP) ! -This is the density grid- in terms of mass per cubic centimeter
    P(r)         = piefinal(r*JUMP) * vff !electron pressure
    NELEC(r)     = RHO(r) * 7.01e23!electron number density
    TK(r)        = piefinal(r*JUMP) * taufinal(r*JUMP)  * (vff**2.) * m_e * (zbar + mbar/m_e)/(zbar*k)!-Electron temperature grid-
    TkeV(r)      = (8.6173e-8)*TK(r) !temprature in KeV
  ENDDO
END SUBROUTINE MCVSPEC_SHOOTING

!Find the shock_height of a given WD and calculates if the plasma is stable based on SIGMA_S
!and the normalized velocity and electron pressures
SUBROUTINE SHOCK_APPROX(p, tnumber, t, pi)
  common  /wdparams/ M_3, M_wd, R_wd, shock_height, vff, n_elec_shock, Tshock, ES0, MDOT0, R_m, Bm, L, p_a, SIGMA
  common /constants/ Msun, Rsun, G, mue, mmw, A, mu, k, mH, alpha, beta, gamma, X__, chi, zbar, R_, m_e, mbar
  common /bool/ STABLE

  REAL pi_e_integral, step_size, pi_e, int_coeff1, int_coeff2, p
  REAL rt_coeff1, rt_coeff2
  REAL limit, m_s
  INTEGER tnumber, steps, local_stable
  REAL pi(tnumber), t(tnumber)
  REAL M_wd, Msun
  REAL arg1, arg2

  m_s = M_wd/Msun ! mass of the white dwarf in solar masses

  !this is a quintic regression fit for what the lower limit of tau
  !should be such that the ODE does not encounter stiffness
  limit = 0.005696503 + (0.06507075*m_s) - 0.1599213*(m_s**2) + 0.1441434*(m_s**3) - 0.05448718*(m_s**4) + 0.006410256*(m_s**5)
  !"safety padding"
  limit = limit + .002

  step_size = (.25-limit)/tnumber
  t_init = .25
  pi_e = p
  pi(1) = pi_e
  t(1) = t_init
  pi_e_integral = 0.
  local_stable = 1


  DO steps = 1, tnumber
      t(steps+1) = t(steps)-step_size
      pi(steps+1) = pi(steps) - dpiedtau(t(steps),pi(steps))*(step_size)
      int_coeff1 = (gamma*(1.-t(steps)) - t(steps)) /(1. + ES0*f_(t(steps),pi(steps)))
      rt_coeff1 = SQRT((t(steps)**3.)/pi(steps))
      int_coeff2 = (gamma*(1.-t(steps+1)) - t(steps+1))/(1. + ES0*f_(t(steps+1),pi(steps+1)))
      rt_coeff2 = SQRT((t(steps+1)**3.)/pi(steps+1))

      !checking stability condition, if not iterate to a higher sigma_s
      arg1 = dpiedtau(t(steps),pi(steps))
      arg2 = -1.*gamma*pi(steps)/t(steps)

      IF (arg1.LT.arg2) THEN
          local_stable = 0
      ENDIF
      !Trapezoidal Riemmann Sum for Calculating the Integral to Find Shock Height
      pi_e_integral = pi_e_integral + (step_size*(((rt_coeff1*int_coeff1) + (rt_coeff2*int_coeff2))/2.))
  ENDDO

  IF(local_stable.EQ.0) THEN
    STABLE = 0
  ENDIF

  rhs_coeff = (gamma-1.)*p_a*A/(vff**2.)

  shock_height = (pi_e_integral/rhs_coeff)
END SUBROUTINE SHOCK_APPROX

!Differential Equation for Normalized Height with respect to Normalized Velocity
REAL FUNCTION dxidtau(t,p)
  common  /wdparams/ M_3, M_wd, R_wd, shock_height, vff, n_elec_shock, Tshock, ES0, MDOT0, R_m, Bm, L, p_a, SIGMA
  common /constants/ Msun, Rsun, G, mue, mmw, A, mu, k, mH, alpha, beta, gamma, X__, chi, zbar, R_, m_e, mbar

  REAL lambda, t, p

  lambda    = (gamma - 1.)*shock_height*(p_a)*(vff**(-2.))*A*SQRT(p/(t**3.)) * (1. + (ES0*f_(t,p)))

  dxidtau = (gamma*(1-t) - t)/lambda

  RETURN
END FUNCTION dxidtau

!Differential Equation for Normalized Electron Pressure with respect to Normalized Velocity
!Used for Shock Height Approximation
REAL FUNCTION dpiedtau(t,p)
  common  /wdparams/ M_3, M_wd, R_wd, shock_height, vff, n_elec_shock, Tshock, ES0, MDOT0, R_m, Bm, L, p_a, SIGMA
  common /constants/ Msun, Rsun, G, mue, mmw, A, mu, k, mH, alpha, beta, gamma, X__, chi, zbar, R_, m_e, mbar

  REAL outterm, midterm, gamlam1, gamlam2

  outterm = (gamma*(1.-t)-t)/t
  midterm = (gamma*p)/(gamma*(1.-t)-t)
  gamlam1 =  (X__/(A*(vff**(2.))))
  gamlam2 = (1.-t-(chi*p))/(t*(p**2)*(1.+ (ES0*f_(t,p))))

  dpiedtau = outterm*(1. - midterm - (gamlam1*gamlam2))

  RETURN
END FUNCTION dpiedtau

!Uses ES0 as a coefficient for describing efficiency of a secondary cooling process
!relative to thermal bremsstrahlung cooling as a function of pressure and velocity
REAL FUNCTION f_(t,p)
  common  /wdparams/ M_3, M_wd, R_wd, shock_height, vff, n_elec_shock, Tshock, ES0, MDOT0, R_m, Bm, L, p_a, SIGMA
  common /constants/ Msun, Rsun, G, mue, mmw, A, mu, k, mH, alpha, beta, gamma, X__, chi, zbar, R_, m_e, mbar
  f_ = (4.**(alpha+beta))/(3.**alpha)*(((1.+sigma)/(sigma))**alpha)*(p**alpha)*(t**beta)
  RETURN
END FUNCTION f_
