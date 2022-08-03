!
Module Wave_Propagator
!
! Wave optics propagator.
!----------------------------------------------------------
! (C) Copyright 1999, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 27 Feb 1999 | Original version.
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Parameters:
    Double, Pi,     &
    C_Light
!
Use Externals, only: &
! Imported Routines:
    CPrintf
!
Use IO, only: &
! Imported Routines:
    PutXY
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
!
Contains


!==========================================================
Subroutine Wave_Propagate &
  (Year,      & ! <-- Occultation year
   Month,     & ! <-- Occultation month
   Day,       & ! <-- Occultation day
   Hour,      & ! <-- Occultation hour
   Minute,    & ! <-- Occultation begin minute
   Second,    & ! <-- Occultation begin second
   TR,        & ! <-- Relative time of samples [sec]
   Freq,      & ! <-- Frequencies [Hz]
   AFreq,     & ! <-- Frequencies for absorption [Hz]
   RLEO,      & ! <-- LEO coordinates (J2000/ECEF)
   RGPS,      & ! <-- GPS coordinates (J2000/ECEF)
   Hmax,      & ! <-- Maximum height
   Hmin,      & ! <-- Minimum height
   DYN,       & ! <-- Minimum vertical scale of N
   DX,        & ! <-- Step between phase screens
   XLS,       & ! <-- Additional screen position
   opt_LS,    & ! <-- Propagate to additional screen
   FZ,        & ! <-- Initial last Fresnel zone size [rad]
   SR,        & ! <-- Sampling rate [Hz]
   opt_TS,    & ! <-- Time scaling
   opt_AS,    & ! <-- Asymptotic solution
   opt_AB,    & ! <-- Modeling of absorption
   opt_ECEF,  & ! <-- Coordinates in ECEF
   opt_FIO,   & ! <-- FIO forward propagation
   WT,        & ! --> WP data time
   WA,        & ! --> WP amplitudes (time,channel)
   WS,        & ! --> WP phase excess (time,channel) [m]
   WRLEO,     & ! --> WP LEO coordinates (J2000/ECEF)
   WVLEO,     & ! --> WP LEO velocity (J2000/ECEF)
   WRGPS,     & ! --> WP GPS coordinates (J2000/ECEF)
   WVGPS,     & ! --> WP GPS velocity (J2000/ECEF)
   PCT,       & ! --> CT impact parameter 
   ECT,       & ! --> CT refraction angle
   ACT,       & ! --> CT amplitude
   Stat,      & ! ~~> Error status
   Vrb)         ! <~~ Verbosity level
!
! Calculation of wave propagation in atmosphere, simulation
! refractometric soundings.
!----------------------------------------------------------
! Method:
!   Multiple phase screens or asymptotic solution based
!   on the inverse canonical transform.
!----------------------------------------------------------
! (C) Copyright 1999-2009, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 27 Feb 1999 | Original version.
!   2.0   | 07 Mar 1999 | Hmax, better apodization.
!   3.0   | 14 Mar 1999 | Dynamical determination of
!         |             | mumber of output data.
!   4.0   | 17 Mar 1999 | Error status.
!   4.1   | 14 May 1999 | Use of Modulo.
!   5.0   | 08 Jul 1999 | Verbosity parameter.
!   5.1   | 17 Dec 2000 | FZ as argument.
!   5.2   | 21 Mar 2001 | SR.
!   6.0   | 15 May 2001 | Absorption.
!   6.1   | 16 May 2001 | Asymptotic solution.
!   6.2   | 13 Oct 2001 | Compute absorption only if opt_AB set;
!         |             | check for 0-amplitude.
!   6.3   | 12 Dec 2001 | Correction for Earth rotation.
!   6.4   | 22 Dec 2001 | Limitation of ray direction angle.
!   6.5   | 09 Aug 2002 | Merged with WaveL_Propagator.f90.
!   6.6   | 10 Aug 2002 | DYN as argument.
!   7.0   | 13 Aug 2002 | Accurate determination of
!         |             | propagation area.
!   7.1   | 21 Aug 2002 | CT in last phase screen.
!   7.2   | 01 Feb 2003 | Variable number of channels.
!   8.0   | 20 Aug 2003 | Verbosity level.
!   8.1   | 29 Oct 2003 | Corrected mistake: IC=1,2
!   8.2   | 24 Mar 2004 | Accurate check for Verbosity level.
!   9.0   | 08 May 2004 | Choice of FP_FIO and FP_Fresnel.
!   9.1   | 27 Aug 2004 | Save/read field before propagation
!         |             | to additional screen.
!   9.2   | 03 Sen 2004 | Variable integration step
!         |             | and grid reduction.
!   9.3   | 04 Sen 2004 | Save and restore checkpoint.
!   9.4   | 16 Sen 2004 | Corrected invoke of Ray_Trace,
!         |             | limitation of DYN for some estimates.
!  10.0   | 07 Nov 2005 | opt_J2000, opt_FIO, Hmin.
!  10.1   | 13 Dec 2005 | opt_ECEF.
!  10.2   | 25 Dec 2005 | Accurate dimension handling.
!  10.3   | 07 Feb 2007 | Correction of phase for 0 amplitude.
!  10.4   | 23 Feb 2007 | opt_TS.
!  10.5   | 18 Apr 2008 | DXR >= DX.
!  10.6   | 04 Jul 2008 | Corrected initialization of Xmin.
!  10.7   | 05 Sep 2008 | Empirical definition of effective
!         |             | phase screen position, HmaxGO.
!  10.8   | 17 Jun 2009 | Iocc defined as index of stationary
!         |             | GPS position.
!  10.9   | 09 Sep 2009 | More screen output.
!  11.0   | 23 Sep 2009 | Index order changed to (time,channel).
!  11.1   | 25 Sep 2009 | Optimized code.
!----------------------------------------------------------
! Modules used:
!
Use Time, only: &
! Imported Routines:
    Elapsed_Time
!
Use Coordinates, only: &
! Imported Type Definitions:
    Cartesian,      &
! Imported Routines:
    Vector_Normed,  &
    Vector_Angle,   &
    Vector_Norm,    &
    Rotate,         &
! Imported Operators:
    Operator(+),    &
    Operator(-),    &
    Operator(*),    &
    Operator(.xx.)
!
Use Earth, only: &
! Imported Parameters:
    R_Earth, H_atm,   &
! Imported Type Definitions:
    Geodetic,         &
! Imported Routines:
    GAST
!
Use Occ_Coordinates, only: &
! Imported Routines:
    Occ_Geometry,            &
    Occ_Geometry_ECEF,       &
    Occ_Point,               &
    Satellite_Velocities,    &
    Plane_Coordinates,       &
    Plane_Basis
!
Use Occ_Diffraction, only: &
! Imported Routines:
    Accumulate_Phase,        &
    Propagate
!
Use Vacuum_Propagator, only: &
! Imported Routines:
    Forward_Propagation_FIO,    &
    Forward_Propagation_Fresnel
!
Use Atmosphere, only: &
! Imported Routines:
    Atmosphere_NGradN
!
Use Atmosphere_Rays, only: &
! Imported Routines:
    Ray_Trace
!
Use Signal, only: &
! Imported Routines:
    Monotonize,   &
    Apodize
!
Use FIO, only: &
! Imported Routines:
    Canonical_Transform
!
Use FFTW, only: &
! Imported Routines:
    Fourier_Filter
!
Use Interpolation, only: &
! Imported Routines:
    Linear,          &
    Nearest_Power2
!
Use IO, only: &
! Imported Routines:
    GetFreeUnit
!
Use Debug, only: &
! Imported Scalars:
    IDbg   ! Mode: 0 - normal run
           !       1 - save field
           !       2 - read field
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Integer, Intent(In)          :: &
   Year       ! Occultation year
!
Integer, Intent(In)          :: &
   Month      ! Occultation month
!
Integer, Intent(In)          :: &
   Day        ! Occultation day
!
Integer, Intent(In)          :: &
   Hour       ! Occultation hour
!
Integer, Intent(In)          :: &
   Minute     ! Occultation begin minute
!
Real(Double), Intent(In)     :: &
   Second     ! Occultation begin second
!
Real(Double), Intent(In)     :: &
   TR(:)      ! Relative time of samples [sec]
!
Real(Double), Intent(In) :: &
   Freq(:)    ! Frequencies [Hz]
!
Real(Double), Intent(In) :: &
   AFreq(:)   ! Frequencies for absorption [Hz]
!
Type(Cartesian), Intent(In) :: &
   RLEO(:)    ! LEO coordinates (J2000/ECEF)
!
Type(Cartesian), Intent(In) :: &
   RGPS(:)    ! GPS coordinates (J2000/ECEF)
!
Real(Double), Intent(In) :: &
   Hmax       ! Maximum height
!
Real(Double), Intent(In) :: &
   Hmin       ! Minimum height
!
Real(Double)       :: &
   DYN        ! Minimum vertical scale of N
!
Real(Double), Intent(In) :: &
   DX         ! Step between phase screens
!
Real(Double), Intent(In) :: &
   XLS        ! Additional screen position
!
Logical, Intent(In)      :: &
   opt_LS     ! Propagate to additional screen
!
Real(Double), Intent(In) :: &
   FZ         ! Initial last Fresnel zone size [rad]
!
Real(Double), Intent(In) :: &
   SR         ! Sampling rate [Hz]
!
Logical, Intent(In)      :: &
   opt_TS     ! Time scaling
              !   True  - scale time to conserve trajectory start/finish time
              !   False - no time scaling
!
Logical, Intent(In)      :: &
   opt_AS     ! Asymptotic solution
!
Logical, Intent(In)      :: &
   opt_AB     ! Modeling of absorption
!
Logical, Intent(In)      :: &
   opt_ECEF   ! Coordinates in ECEF
!
Logical, Intent(In)      :: &
   opt_FIO    ! FIO forward propagation
!
! Output arguments:
!
Real(Double), Pointer :: &
   WT(:)      ! WP data time
!
Real(Double), Pointer :: &
   WA(:,:)    ! WP Amplitudes (time,channel)
!
Real(Double), Pointer :: &
   WS(:,:)    ! WP Phase excess (time,channel) [m]
!
Type(Cartesian), Pointer :: &
   WRLEO(:)   ! WP LEO coordinates (J2000/ECEF)
!
Type(Cartesian), Pointer :: &
   WVLEO(:)   ! WP LEO velocities (J2000/ECEF)
!
Type(Cartesian), Pointer :: &
   WRGPS(:)   ! WP GPS coordinates (J2000/ECEF)
!
Type(Cartesian), Pointer :: &
   WVGPS(:)   ! WP GPS velocities (J2000/ECEF)
!
Real(Double), Pointer :: &
   PCT(:)     ! CT impact parameter 
!
Real(Double), Pointer :: &
   ECT(:,:)   ! CT refraction angle (p,channel)
!
Real(Double), Pointer :: &
   ACT(:,:)   ! CT amplitude (p,channel)
!
Integer, Intent(Out)  :: &
   Stat       ! Error status
!
! Input optional arguments:
!
Integer, Optional, Intent(In) :: &
   Vrb        ! Verbosity mode.
              ! True by default.
!----------------------------------------------------------
! Local Parameters:
!
Character(Len=*), Parameter :: &
   CkptName = 'saved-wop.dat' ! Chekpoint file name
!
Real(Double), Parameter :: &
   DYNmin = 0.03_Double       ! Minimum DYN for some estimates [km]
!
Real(Double), Parameter :: &
   HmaxGO = 5.0_Double        ! Maximum height for GO simulation [km]
!
Real(Double), Parameter :: &
   A0 = 750.0_Double          ! Vacuum amplitude
!
Real(Double), Parameter :: &
   ELim = 0.1_Double          ! Limit ray direction angle for
                              ! spatial spectra
!
Complex(Double), Parameter ::  &
   Ci = (0.0_Double, 1.0_Double)  ! I = Sqrt(-1)
!
Real(Double), Parameter :: &
   DYA = 5.0_Double           ! Apodization border width
!
Real(Double), Parameter :: &
   Ain = 0.95                 ! Inhomogeneity parameter of grid of p
                              ! for GO propagation
!
Real(Double), Parameter :: &
   XB  = 400.0_Double         ! BP plane position for inverse
                              ! canonical transforms
!
Real(Double), Parameter :: &
   YB    = 10.0_Double,     & ! Lenght of safety border area
   Pmax  = 10.0_Double        ! Maximum ray height
!
Type(Cartesian), Parameter :: &
   PA =  &  ! Polar axis
      Cartesian((/0,0,1/))
!
Real(Double), Parameter :: &
   DXA(6) = (/  0.5,   1.0,   2.0,    5.0,  10.0,   25.0   /),  &
   KA(6)  = (/ -0.25, -0.25, -0.235, -0.12, -0.062, -0.006 /)
!
! Local Scalars:
!
! --- Geometry of occultation
!
Integer           :: N        ! Number of data read
Type(Geodetic)    :: GP       ! Occultation point
Type(Cartesian)   :: ERLC     ! Local curvature center (ECEF)
Type(Cartesian)   :: RLC      ! Local curvature center (J2000)
Real(Double)      :: RE       ! Local curvature radius
Type(Cartesian)   :: AX, AY   ! Occultation plane XY-basis
Type(Cartesian)   :: EX, EY   ! Occultation plane basis (ECEF)
Real(Double)      :: PhiN     ! Last Earth rotation angle
Real(Double)      :: PhiT     ! Current Earth rotation angle
Integer           :: NW       ! Number of simulated data
Integer           :: Iocc     ! Occultation point index
!
! --- Geometry of propagation area
!
Real(Double)      :: Xmin     ! Minimum X for phase screens
Real(Double)      :: Xmax     ! Maximum X for phase screens
Real(Double)      :: Ymin     ! Minimum Y for phase screens
Real(Double)      :: Ymax     ! Maximum Y for phase screens
Real(Double)      :: XImin    ! Minimum momentum
Real(Double)      :: XImax    ! Maximum momentum
Real(Double)      :: DXR      ! Recalculated step between screens
Integer           :: NUmax    ! Maximum number of complex field samples
Integer           :: NU       ! Current number of samples
Integer           :: NN       ! Number of samples of N
Real(Double)      :: E0       ! Divergence angle of spheric wave
Real(Double)      :: DY       ! Sample step in phase screen
Real(Double)      :: DYL      ! Sample step in last phase screen
Real(Double)      :: AS       ! Synthesized aperture size
Integer           :: NL       ! Number of samples in last screen
Real(Double)      :: R        ! Distance from GPS to point in phase screen
Real(Double)      :: YE       ! Position of Earth surface
Integer           :: IR       ! Reduced grid index
Real(Double)      :: YBC      ! Boundary condition location
Integer           :: IL       ! Index of lowest point out of apodization
Integer           :: IH       ! Index of heighest point out of apodization
Integer           :: DI       ! Index step for field for forward propagation
!
! --- Variables for ray-tracing
!
Integer           :: NGO      ! Number of rays for GO propagation
Type(Cartesian)   :: XT       ! Transmitter position
Type(Cartesian)   :: XN       ! Final ray point
Type(Cartesian)   :: UT       ! Ray direction at transmitter
Type(Cartesian)   :: UN       ! Ray direction at XN
Real(Double)      :: t        ! Parameter for scaling grid of p
Real(Double)      :: YP       ! Ray leveling height
Real(Double)      :: YPmin    ! Minimum leveling height
Real(Double)      :: Hper     ! Perigee height
Real(Double)      :: Y1, Y2   ! Perigee limits for dichotomy
Real(Double)      :: Eps      ! Refraction angle
Integer           :: RStat    ! Ray tracer status
!
! --- Variables for canonical transform
!
Integer           :: OCD      ! GO impact parameter direction
Integer           :: NY       ! Number of high-res field samples
Real(Double)      :: YHmin    ! Minimum YH for CT
Real(Double)      :: YHmax    ! Maximum YH for CT
Real(Double)      :: DYH      ! Sample step for CT
Real(Double)      :: DYS      ! Standard step of high-res Y-grid
Real(Double)      :: AHI      ! Interpolated amplitude
Real(Double)      :: PHI      ! Interpolated phase
Real(Double)      :: PHHI     ! Interpolated hi-res phase
Real(Double)      :: R0       ! Shift of Y-coordinate
Integer           :: CTSgn    ! Canonical transform direction:
Integer           :: IC       ! Channel number
Real(Double)      :: YBmax    ! Maximum usable Y in BP plane
Integer           :: IBmax    ! Maximum usable index in BP plane
Real(Double)      :: Pmin     ! Minimum ray height estimate
Real(Double)      :: XG       ! X-coordinate of GPS satellite
Real(Double)      :: WF       ! Phase filtering width
Integer           :: NCT      ! Number of output CT data
!
! --- Variables for wave propagation
!
Integer           :: NC       ! Number of channels
Real(Double)      :: HE       ! Estimated altitude of lowest visible ray
Real(Double)      :: DPh      ! Phase rate
Real(Double)      :: DYF      ! Y step of field for forward propagation
Real(Double)      :: KS       ! Coefficient for effective screen position
!
! --- Work variables
!
Integer           :: LVrb     ! Verbosity level
Integer           :: i        ! Array indices
Character(Len=80) :: Line     ! Terminal line
!
!
! Local Arrays:
!
! --- Geometry of occultation
!
Type(Cartesian), Allocatable :: &
   ERLEO(:),  & ! LEO coordinates (ECEF)
   ERGPS(:)     ! GPS coordinates (ECEF)
Real(Double), Allocatable    :: &
   XLEO(:),   & ! X coordinates of LEO
   YLEO(:),   & ! Y coordinates of LEO
   XGPS(:),   & ! X coordinates of GPS
   YGPS(:)      ! Y coordinates of GPS
Real(Double), Pointer        :: &
   WX(:),     & ! X coordinates of observation curve
   WY(:)        ! Y coordinates of observation curve
Type(Cartesian), Allocatable :: &
   PRLEO(:),  & ! Plane LEO coordinates
   PRGPS(:)     ! Plane GPS coordinates
!
! --- Wave field
!
Real(Double), Allocatable :: &
   k(:),      & ! Wave vectors (channel)
   Ak(:),     & ! Wave vectors (channel)
   Lam(:),    & ! Wavelengths (channel)
   At(:)        ! Attenuation [exp factor]
Real(Double), Allocatable, Target    :: &
   Y(:),      & ! Y-coordinate in phase screen
   X(:)         ! X-coordinate in phase screen
Complex(Double), Allocatable, Target :: &
   U0(:,:),   & ! Complex field in previous screen (Y,channel)
   U(:,:)       ! Propagated complex field in screen (Y,channel)
Real(Double), Allocatable, Target    :: &
   A(:,:),    & ! Amplitude of complex field in screen (Y,channel)
   Ph(:,:)      ! Phase of complex field in screen (Y,channel)
Real(Double), Pointer        :: &
   YL(:),     & ! Y-coordinate in last phase screen
   XL(:)        ! X-coordinate in last phase screen
Complex(Double), Pointer     :: &
   UL0(:,:),  & ! Complex field in last screen (Y,channel)
   UL(:,:)      ! Propagated complex field in last screen (Y,channel)
Real(Double), Pointer        :: &
   AL(:,:),   & ! Amplitude of complex field in last screen (Y,channel)
   PL(:,:)      ! Phase of complex field in last screen (Y,channel)
Real(Double), Pointer        :: &
   WP(:,:)      ! Accumulated phase of observed field (time,channel)
Real(Double), Allocatable    :: &
   RN(:),     & ! Refractive index in phase screen (Y)
   RNI(:,:)     ! Imaginary part of refractive index
                ! in phase screen (Y,channel)
!
! --- Variables for Geometric_Optics ray tracer
!
Real(Double), Allocatable    :: &
   PGO(:),    & ! GO impact parameter
   EGO(:),    & ! P-momentum
   AGO(:,:),  & ! GO amplitude (time,channel)
   AP(:,:)      ! Amplitude in (p,xi)-representation (p,channel)
Real(Double)                 :: &
   Alf(3),    & ! Sequential ingoing ray andles
   YN(3),     & ! Sequential outgoing ray positions
   PN(3)        ! Sequential outgoing ray impact parameters
!
! --- Variables for canonical transform
!
Real(Double), Allocatable    :: &
   YGO(:),    & ! GO Y-coordinate
   HGO(:),    & ! Y-momentum
   HY(:)        ! Interpolated momentum
Real(Double), Allocatable :: &
   YH(:),         & ! High-res Y-grid
   PHH(:,:),      & ! Accumulated phase of wave function (p,channel)
   AH(:,:),       & ! Amplitude of wave function (p,channel)
   AF(:,:),       & ! Filtered amplitude (p,channel)
   EH(:,:)          ! Refraction angles (p,channel)
Real(Double), Allocatable :: &
   PminC(:)         !  Estimate of shadow border for each channel.
!
Complex(Double), Allocatable :: &
   UH(:,:)          ! High-res interpolated complex field (p,channel)
!----------------------------------------------------------


!----------------------------------------------------------
! 0. INITIALIZATION
!----------------------------------------------------------


!--- 0.1. Determination of number of input data and channels

N  = Size(TR)
NC = Size(Freq)


!--- 0.2. Array size check

If (.not. All((/         &
      Size(RLEO) == N,   &
      Size(RGPS) == N    &
      /))) then
   Stat = -1
   Return
Else
   Stat = 0
End If


!--- 0.3. Array allocation

Allocate(ERLEO(1:N))
Allocate(ERGPS(1:N))
Allocate(XLEO(1:N))
Allocate(YLEO(1:N))
Allocate(XGPS(1:N))
Allocate(YGPS(1:N))
Allocate(k(NC))
Allocate(Ak(NC))
Allocate(Lam(NC))
Allocate(At(NC))
Allocate(PminC(NC))


!--- 0.4. Setting verbosity mode

If (Present(Vrb)) then
   LVrb = Vrb
Else
   LVrb = 0
End If


!--- 0.6. Calculation of wave vectors

k(:)   = 2*Pi*Freq(:)/C_Light
Lam(:) = 2*Pi/k(:)
Ak(:)  = 2*Pi*AFreq(:)/C_Light


!----------------------------------------------------------
! 1. COORDINATE TRANSFORMS
!----------------------------------------------------------


!--- 1.1. Calculation of curvature center and
!---      coordinates in ECEF frame

If (opt_ECEF) then
   ERLEO(:) = RLEO(:)
   ERGPS(:) = RGPS(:)
   Call Occ_Geometry_ECEF &
     (ERLEO,     & ! <-- LEO coordinates (ECEF)
      ERGPS,     & ! <-- GPS coordinates (ECEF)
      GP,        & ! --> Geodetic coordinates of occultation point
      ERLC,      & ! --> Curvature center (ECEF)
      RE)          ! --> Local curvature radius
   RLC = ERLC
Else
   Call Occ_Geometry &
     (Year,      & ! <-- Occultation year
      Month,     & ! <-- Occultation month
      Day,       & ! <-- Occultation day
      Hour,      & ! <-- Occultation hour
      Minute,    & ! <-- Occultation begin minute
      Second,    & ! <-- Occultation begin second
      TR,        & ! <-- Relative time of samples [sec]
      RLEO,      & ! <-- LEO coordinates (J2000)
      RGPS,      & ! <-- GPS coordinates (J2000)
      ERLEO,     & ! --> LEO coordinates (ECEF)
      ERGPS,     & ! --> GPS coordinates (ECEF)
      GP,        & ! --> Geodetic coordinates of occultation point
      ERLC,      & ! --> Curvature center (ECEF)
      RLC,       & ! --> Curvature center (J2000)
      RE,        & ! --> Local curvature radius
      Stat=Stat)   ! ~~> Error status
End If

If (Stat /= 0) then
   Return
End If


!--- 1.2. Calculation of coordinates in occultation plane

Call Plane_Coordinates &
  (RLEO,      & ! <-- LEO coordinates
   RGPS,      & ! <-- GPS coordinates
   RLC,       & ! <-- Curvature center
   RE,        & ! <-- Local curvature radius [km]
   XLEO,      & ! --> X coordinates of LEO
   YLEO,      & ! --> Y coordinates of LEO
   XGPS,      & ! --> X coordinates of GPS
   YGPS,      & ! --> Y coordinates of GPS
   AX,        & ! --> Occultation plane X basis vector
   AY)          ! --> Occultation plane Y basis vector

YLEO(:)  = YLEO(:) - RE
YGPS(:)  = YGPS(:) - RE


!--- 1.3. Calculation of plane basis (ECEF)

Call Plane_Basis &
  (ERLEO,     & ! <-- LEO coordinates
   ERGPS,     & ! <-- GPS coordinates
   ERLC,      & ! <-- Curvature center
   RE,        & ! <-- Local curvature radius [km]
   EX,        & ! --> Occultation plane X basis vector
   EY)          ! --> Occultation plane Y basis vector


!--- 1.4. Occultation point determination

Call Occ_Point &
  (ERLEO,     & ! <-- LEO coordinates (ECEF)
   ERGPS,     & ! <-- GPS coordinates (ECEF)
   GP,        & ! --> Occultation point (Geodetic)
   Iocc)        ! --> Occultation point index


!----------------------------------------------------------
! 2. DETERMINATION OF PROPAGATION AREA
!----------------------------------------------------------


!--- 2.1. Determination of shadow border

If (LVrb >= 1) then
   Write(*,'(2X,A)') 'Determination of shadow border.'
End If

Y1  = -20.0_Double
Y2  =  20.0_Double
Eps = 0

XT = ERLC + EX*XGPS(Iocc) + EY*(RE+YGPS(Iocc))

If (LVrb >= 2) then
   Write(*,'(2X,A,F10.3)') 'XGPS = ', XGPS(Iocc)
End If

Dichotomy: Do
   YP = (Y1 + Y2)/2
   UT = Vector_Normed(ERLC + EY*(RE+YP) - XT)
!   print*,'Dichotomy: calling Ray_Trace'
   Call Ray_Trace &
     (XT,          & ! <-- Transmitter position
      UT,          & ! <-- Transmitter ray direction
      DS   = DYN,  & ! <-- Integration step parameter
      XN   = XN,   & ! --> Ray point nearest to receiver
      Hper = Hper, & ! ~~> Perigee altitude [km]
      UN   = UN,   & ! --> Ray direction at XN
      Stat = RStat)  ! --> Error status
!   print*,'Dichotomy: after Ray_Trace'
   If (RStat == 0 .and. Hper >= Hmin) then
      Y2   = YP
   Else
      Y1 = YP
   End If
   If (LVrb >= 3) then
      Write (Line,'(2X,A,F9.4,A,F9.4,A,I1,A1)')  &
         'Y1 = ', Y1,  ', Y2 = ', Y2,  ', RStat = ', RStat, &
         Char(0)
      Call CPrintf(Line)
   End If
   If (Abs(Y2-Y1) < 0.001) then
      YPmin = Y2
      Exit Dichotomy
   End If
End Do Dichotomy


If (LVrb >=3) then
   Write(*,'()')
End If


!--- 2.2. GO propagation

If (LVrb >= 1) then
   Write(*,'(2X,A)') 'GO propagation.'
End If

NGO = 2*(1-Ain)*Hmax/Max(DYN,DYN)
If (Vrb >= 2) then
   Write(*,'(2X,A,I5)') 'NGO = ', NGO
End If

Allocate(PGO(NGO))
Allocate(EGO(NGO))
Allocate(AGO(NGO,NC))
Allocate(AP(NGO,NC))

Alf(:) = 0
YN(:)  = 0

Xmax =   Sqrt((R_Earth + H_atm)**2 - R_Earth**2)
Xmin = - Sqrt((R_Earth + H_atm)**2 - R_Earth**2)

Do i=1,NGO
   t  = Real(NGO - i)/Real(NGO - 1)
   YP = YPmin + ((1-Ain)*t + Ain*t*t)*Hmax
   UT = Vector_Normed(ERLC + EY*(RE+YP) - XT)
   If (opt_AB) then
      Call Ray_Trace &
         (XT,               & ! <-- Transmitter position
          UT,               & ! <-- Transmitter ray direction
          DS   = DYN,       & ! <-- Integration step parameter
          XN   = XN,        & ! --> Ray point nearest to receiver
          UN   = UN,        & ! --> Ray direction at XN
          At   = At,        & ! ~~> Attenuation [exponential factor]
          Stat = RStat)       ! --> Error status
   Else
      Call Ray_Trace &
         (XT,               & ! <-- Transmitter position
          UT,               & ! <-- Transmitter ray direction
          DS   = DYN,       & ! <-- Integration step parameter
          XN   = XN,        & ! --> Ray point nearest to receiver
          UN   = UN,        & ! --> Ray direction at XN
          Stat = RStat)       ! --> Error status
      At(:) = 0.0
   End If
   If (Vector_Angle(EX, UN) < ELim) then
      Xmax = Max(Xmax, EX*(XN - ERLC))
   End If
   PGO(i) = Vector_Norm((XN-ERLC) .xx. UN)
   EGO(i) = Vector_Angle(EX, UN, (EX).xx.(EY))
   Alf(1) = Alf(2)
   Alf(2) = Alf(3)
   Alf(3) = Vector_Angle(EX, UT, (EX).xx.(EY))
   YN(1)  = YN(2)
   YN(2)  = YN(3)
   YN(3)  = (PGO(i) + XB*Sin(EGO(i)))/Cos(EGO(i))
   PN(1)  = PN(2)
   PN(2)  = PN(3)
   PN(3)  = PGO(i)
   If (i > 2) then
      AGO(i-1,1:NC) = 100*Exp(-Ak(1:NC)*At(1:NC))*                     &
                   Sqrt(Abs((Xmin-XGPS(Iocc))*(Alf(3) - Alf(1))/       &
                               ((YN(3) - YN(1))*Cos(EGO(i-1)))))
      AP(i-1,1:NC)  = 100*Exp(-Ak(1:NC)*At(1:NC))*                     &
                   Sqrt(Abs((Xmin-XGPS(Iocc))*(Alf(3) - Alf(1))/       &
                              (PN(3) - PN(1))))
   End If
   If (LVrb >= 3) then
      Write (Line,'(2X,A,I5,A,F9.3,A,ES10.3,A,F6.0,A1)')  &
         'i = ', i,                                    &
         ', P = ', PGO(i)-RE,  ', Xi = ', EGO(i),      &
         ', Xmax = ', Xmax,                            &
         Char(0)
      Call CPrintf(Line)
   End If
End Do

AGO(1,1:NC)   = AGO(2,1:NC)
AGO(NGO,1:NC) = AGO(NGO-1,1:NC)
AP(1,1:NC)    = AP(2,1:NC)
AP(NGO,1:NC)  = AP(NGO-1,1:NC)

If (LVrb >= 3) then
   Write(*,'()')
End If

!!!DEBUG
! Call PutXY('wop-ap1.dat', PGO(:), AP(:,1), Stat=RStat)
! Call PutXY('wop-ap2.dat', PGO(:), AP(:,2), Stat=RStat)


!--- 2.3. Determination of Ymin, Ymax, Xmin, and Eps

Ymin = MinVal((PGO(:) + Xmax*Sin(EGO(:)))/Cos(EGO(:)) - RE,   &
              Mask = (-EGO(:) < ELim))
Ymax = MaxVal((PGO(:) + Xmax*Sin(EGO(:)))/Cos(EGO(:)) - RE)

AS   = 0.25*MaxVal(Xmax*2*Pi/(k(:)*Max(DYNmin,DYN)))

Ymin = Ymin - DYA - AS
Ymax = Ymax + DYA

Eps  = MaxVal(Abs(EGO(:)), Mask = (-EGO(:) < ELim))

If (LVrb >= 2) then
   Write(*,'(2X,A,F6.0,A,F6.0,A,F6.0)') &
      'Ymin = ', Ymin, ', Ymax = ', Ymax, ', Xmin = ', Xmin
End If


!----------------------------------------------------------
! 3. DETERMINATION OF RESOLUTION
!----------------------------------------------------------


!--- 3.1. Determination of initial resolution

E0 = Max(Abs(Ymax/XGPS(Iocc)), Abs(Ymin/XGPS(Iocc)))
NU = Ceiling(MaxVal((Ymax-Ymin)*k(:)*E0*2/Pi))
NN = Ceiling((Ymax-Ymin)/Max(DYNmin,DYN))
NU = Max(NU, NN)
NU = Nearest_Power2(NU)


!--- 3.2. Determination of maximum resolution

NUmax = Ceiling(MaxVal((Ymax-Ymin)*k(:)*Eps*2/Pi))
NUmax = Max(NU, NUmax)
NUmax = Nearest_Power2(NUmax)

If (LVrb >= 2) then
   Write(*,'(2X,A,I10)') 'NUmax = ', NUmax
End If

Allocate(U0(NUmax,NC))
Allocate(U(NUmax,NC))
Allocate(Ph(NUmax,NC))
Allocate(A(NUmax,NC))
Allocate(Y(NUmax))
Allocate(X(NUmax))
Allocate(RN(NUmax))
Allocate(RNI(NUmax,NC))


!--- 3.3. Determination of max/min momentum

XImax = Abs(Ymax/XGPS(Iocc))
XImin = -Eps


!----------------------------------------------------------
! 4. MULTIPLE PHASE SCREEN PROPAGATION
!----------------------------------------------------------

If (.not. opt_AS) then

   If (LVrb >= 1) then
      Write(*,'(2X,A)') 'Multiple phase screens.'
   End If


   !--- 4.1. Setting initial condition in first screen

   X(:) = Xmin
   Do i=1,NU
      Y(i)    = (Ymin*(NU-i) + Ymax*(i-1))/(NU-1)
   End Do
   DY   = Y(2) - Y(1)

   YE   = Ymin

   Do IC=1,NC

      Do i=1,NU
         R       = Sqrt((Xmin-XGPS(Iocc))**2 + Y(i)**2)
         U(i,IC) = A0*Sqrt(Xmin-XGPS(Iocc))*Exp(Ci*k(IC)*R)/Sqrt(R)
      End Do

      Ph(1:NU,IC) = Atan2(AImag(U(1:NU,IC)),Real(U(1:NU,IC)))
      Call Accumulate_Phase(Ph(1:NU,IC))

      Call Apodize(Y(1:NU), U(1:NU,IC), DYA)
      U0(1:NU,IC) = U(1:NU,IC)

   End Do


   !--- 4.2. Multiple phase screens

   If (IDbg == 2) then
      Call Restore_Ckpt
   End If

   Phase_Screens: Do While (X(1) < Xmax)


      !--- 4.2.1. Boundary condition

      Call Boundary_Condition &
        (Y(1:NU),     & ! <-- Grid
         YE,          & ! <-- Earth's boundary
         X(1),        & ! <-- X-coordinate of phase screen
         DYA,         & ! <-- Apodization border width
         U0(1:NU,:),  & ! <-> Complex field (Y,channel)
         YBC,         & ! --> Boundary condition location
         IR)            ! --> Reduced grid index


      !--- 4.2.2. Vacuum propagation

      HE = Sqrt(RE**2 + X(1)**2) - RE

      If (HE < 25.0) then
         DXR = DX
      Else
         DXR = Max(DX, Min(20.0_Double, 10.*DX))
      End If

      DY = Y(2) - Y(1)

      Do IC=1,NC
         U(:,IC) = 0
         Call Propagate &
           (U0(IR:NU,IC), & ! <-- Complex field in source plane
            k(IC),        & ! <-- Wave vector
            DY,           & ! <-- Discretization step
            DXR,          & ! <-- Distance to observation plane
            U(IR:NU,IC))    ! --> Propagated field in observation plane
      End Do

      X(:) = X(:) + DXR


      !--- 4.2.3. Refractivity calculation

      Call Linear(DXA, KA, DXR, KS)

!      print*,'calling make_refractivity'

      Call Make_Refractivity(X(1)+KS*DXR, Y(1:NU), RN(1:NU), RNI(1:NU,:), YE)


      !--- 4.2.4. Phase correction

      Do IC=1,NC
         Do i=1,NU
            If (Y(i) > YE) then
               U(i,IC) = U(i,IC)*                                 &
                        Exp(Ci*Modulo(k(IC)*DXR*RN(i), 2*Pi))*    &
                        Exp(-Abs(Ak(IC)*DXR*RNI(i,IC)))
            End If
         End Do
      End Do


      !--- 4.2.5. Apodization

      Call Apodize(Y(1:NU), U(1:NU,:), DYA)


      !--- 4.2.6. Calculation of amplitude and phase

      Do IC=1,NC

         A(1:NU,IC)  = Abs(U(1:NU,IC))

         Where (A(1:NU,IC) /= 0.0_Double)
            Ph(1:NU,IC) = Atan2(AImag(U(1:NU,IC)),Real(U(1:NU,IC)))
         Elsewhere
            Ph(1:NU,IC) = 0.0
         End Where

         Call Accumulate_Phase(Ph(1:NU,IC))

      End Do

      DPh = MaxVal(Abs(Ph(2:NU,1)-Ph(1:NU-1,1)), &
                   Mask = (Y(1:NU-1) > Max(YBC, Ymin+DYA)) .and. &
                          (Y(2:NU)   < Ymax-DYA))


      !--- 4.2.7. Doubling resolution if necessary

      If ((DPh > Pi/2) .and. (NU < NUmax)) then
         Call Double_Resolution(Y, Ph, U, RN, DY, NU)
      End If


      !--- 4.2.8. Setting conditions for next step

      U0(1:NU,1:NC) = U(1:NU,1:NC)


      !--- 4.2.9. Printing progress message

      If (LVrb >= 3) then
         Write (Line,'(2X,2(A,F7.1),A,F4.1,A,I8,A,I8,A1)')   &
            'X = ', X(1),  ', Xmax = ', Xmax,                &
            ', DPh = ', DPh, ', IR:NU = ', IR, ':', NU,      &
            Char(0)
         Call CPrintf(Line)
      End If


      !--- 4.2.10. Saving checkpoint

      If (IDbg == 1 .or. IDbg == 2) then
         Call Save_Ckpt(20)
      End If

   End Do Phase_Screens

   If (LVrb >= 3) then
      Write(*,'()')
   End If


   !--- 4.3. Recomputing Xmax

   Xmax = X(1)


!----------------------------------------------------------
! 5. ASYMPTOTIC SOLUTION
!----------------------------------------------------------

Else

   If (LVrb >= 1) then
      Write(*,'(2X,A)') 'Asymptotic solution.'
   End If

   !--- 5.1. GO canonical transform

   Allocate(YGO(NGO))
   Allocate(HGO(NGO))

   Do i=1,NGO
      YGO(i) = (PGO(i) + XB*Sin(EGO(i)))/Cos(EGO(i)) - RE
      HGO(i) = Sin(EGO(i))
   End Do

   OCD = Nint(Sign(1.0_Double, YGO(NGO)-YGO(1)))
   Call Monotonize(OCD, 0.05_Double, YGO(:))


   !--- 5.2. Computation of single-ray wave field

   Allocate(HY(NUmax))

   NU      = NUmax
   X(:)    = Xmin

   DY      = (Ymax-Ymin)/(NU-1)
   Y(1)    = Ymin
   Ph(1,:) = 0
   Do i=2,NU
      Y(i)    = Ymin + (i-1)*DY
      Do IC=1,NC
         Call Linear &
           (YGO,         & ! <-- Argument grid
            AGO(:,IC),   & ! <-- Gridded function
            Y(i),        & ! <-- Interpolation point
            A(i,IC),     & ! --> Interpolated function value
            CExt=.True.)   ! <~~ Constant/linear extrapolation
      End Do
      Call Linear &
        (YGO,     & ! <-- Argument grid
         HGO,     & ! <-- Gridded function
         Y(i),    & ! <-- Interpolation point
         HY(i))     ! --> Interpolated function value
      Ph(i,:) = Ph(i-1,:) + k(:)*DY*HY(i)
      U(i,:)  = A(i,:)*Exp(Ci*Modulo(Ph(i,:),2*Pi))
   End Do
   A(1,:) = A(2,:)
   U(1,:) = U(1,:)


   !--- 5.3. Computation of multi-ray wave field in (p,xi)-representation

   YHmin = MinVal(YGO) - YB
   YHmax = Pmax + 2*YB
   DYS   = 0.125*MinVal(Lam(:))/0.05
   NY    = Nearest_Power2(1 + Ceiling((YHmax - YHmin)/DYS))
   DYH   = (YHmax - YHmin)/Real(NY-1,Double)

   If (LVrb >= 2) then
      Write(*,'(2X,A,ES10.2)')        &
         'DYS = ', DYS
      Write(*,'(2X,2(A,F6.2),A,I7)')  &
         'YHmin = ', YHmin, '   YHmax = ', YHmax, '   NY = ', NY
   End If


   Allocate(YH(0:NY-1))
   Allocate(UH(0:NY-1,NC))
   Allocate(EH(0:NY-1,NC))
   Allocate(AH(0:NY-1,NC))
   Allocate(PHH(0:NY-1,NC))

   Do i=0,NY-1
      YH(i) = YHmin + i*DYH
   End Do

   R0 = RE + Ymin

   Do IC=1,NC
      Do i=0,NY-1
         Call Linear &
           (PGO,          & ! <-- Argument grid
            EGO,          & ! <-- Gridded function
            YH(i) + RE,   & ! <-- Interpolation point
            EH(i,IC),     & ! --> Interpolated function value
            CExt=.False.)   ! <~~ Constant/linear extrapolation
         Call Linear &
           (PGO,          & ! <-- Argument grid
            AP(:,IC),     & ! <-- Gridded function
            YH(i) + RE,   & ! <-- Interpolation point
            AH(i,IC),     & ! --> Interpolated function value
            CExt=.True.)    ! <~~ Constant/linear extrapolation
      End Do
   End Do


   !--- 5.4. Computation of complex wave field in (p,xi)-representation

   Do IC=1,NC

      PHI      = 0
      AHI      = 0
      UH(0,IC) = 0

      Do i=1,NY-1
         PHI = PHI +                                      &
               k(IC)*(Sin(EH(i,IC)) + Sin(EH(i-1,IC)))*   &
                     (YH(i) - YH(i-1))/2
         If (YH(i) > MinVal(PGO) - RE) then
            AHI = AH(i,IC)
         Else
            AHI = 0
         End If
         UH(i,IC) = AHI*Exp(Ci*Modulo(PHI, 2*Pi))
      End Do

   End Do


   !--- 5.5. Inverse canonical transform

   CTSgn = -1

   Call Apodize(YH, UH, YB)

   Call Canonical_Transform &
     (CTSgn,  & ! <-- Canonical transform direction
      k,      & ! <-- Wave vectors [km^-1]
      XB,     & ! <-- X-coordinate of BP plane [km]
      R0,     & ! <-- Shift of Y-grid [km]
      YH,     & ! <-- High resolution Y-grid [km]
      UH)       ! <-> Wave function (Y,channel)


   !--- 5.6. Computation of phase and amplitude

   Do IC=1,NC

      AH(0:NY-1,IC)  = Abs(UH(0:NY-1,IC))

      Where (UH(0:NY-1,IC) /= 0.0_Double)
         PHH(0:NY-1,IC) = ATan2(AImag(UH(0:NY-1,IC)),Real(UH(0:NY-1,IC)))
      Elsewhere
         PHH(0:NY-1,IC) = 0.0
      End Where

      Call Accumulate_Phase(PHH(:,IC))

   End Do


   !--- 5.7. Phase correction

   Do IC=1,NC
      Call Linear &
        (Y,            & ! <-- Argument grid
         Ph(:,IC),     & ! <-- Gridded function
         Pmax,         & ! <-- Interpolation point
         PHI,          & ! --> Interpolated function value
         CExt=.False.)   ! <~~ Constant/linear extrapolation
      Call Linear &
        (YH,           & ! <-- Argument grid
         PHH(:,IC),    & ! <-- Gridded function
         Pmax,         & ! <-- Interpolation point
         PHHI,         & ! --> Interpolated function value
         CExt=.False.)   ! <~~ Constant/linear extrapolation
      PHH(:,IC) = PHH(:,IC) + PHI - PHHI
   End Do


   !--- 5.8. Merging complex field

   Do IC=1,NC
      Do i=1,NU
         If (Y(i) < Pmax) then
            Call Linear &
              (YH,           & ! <-- Argument grid
               PHH(:,IC),    & ! <-- Gridded function
               Y(i),         & ! <-- Interpolation point
               Ph(i,IC),     & ! --> Interpolated function value
               CExt=.False.)   ! <~~ Constant/linear extrapolation
            Call Linear &
              (YH,           & ! <-- Argument grid
               AH(:,IC),     & ! <-- Gridded function
               Y(i),         & ! <-- Interpolation point
               A(i,IC),      & ! --> Interpolated function value
               CExt=.True.)    ! <~~ Constant/linear extrapolation
            U(i,IC)  = A(i,IC)*Exp(Ci*Modulo(Ph(i,IC),2*Pi))
         End If
      End Do
   End Do


   !--- 5.9. Vacuum propagation to last screen position

   U0(:,:) = U(:,:)

   YBmax = Ymax*(XB - XGPS(Iocc))/(Xmax - XGPS(Iocc))
   IBmax = Sum(MaxLoc(Y, Mask=Y(:) < YBmax))

   U0(IBmax+1:NU,:) = 0
   Call Apodize(Y(1:IBmax), U0(1:IBmax,:), DYA)

   Do IC=1,NC
      Call Propagate &
        (U0(1:NU,IC),  & ! <-- Complex field in source plane
         k(IC),        & ! <-- Wave vector
         DY,           & ! <-- Discretization step
         Xmax-XB,      & ! <-- Distance to observation plane
         U(1:NU,IC))     ! --> Propagated field in observation plane
   End Do

   Do IC=1,NC

      A(1:NU,IC)  = Abs(U(1:NU,IC))

      Where (U(1:NU,IC) /= 0.0_Double)
         Ph(1:NU,IC) = Atan2(AImag(U(1:NU,IC)),Real(U(1:NU,IC)))
      Elsewhere
         Ph(1:NU,IC) = 0.0
      End Where

      Call Accumulate_Phase(Ph(1:NU,IC))

   End Do

   X(:) = Xmax


   !--- 5.10. Memory deallocation

   Deallocate(YGO)
   Deallocate(HGO)
   Deallocate(HY)
   Deallocate(YH)
   Deallocate(PHH)
   Deallocate(AH)
   Deallocate(EH)
   Deallocate(UH)

End If


!----------------------------------------------------------
! 6. CT IN LAST PHASE SCREEN
!----------------------------------------------------------


If (LVrb >= 1) then
   Write(*,'(2X,A)') 'CT in last phase screen.'
End If


!--- 6.1. Memory allocation

Allocate(YH(0:NU-1))
Allocate(UH(0:NU-1,NC))
Allocate(PHH(0:NU-1,NC))
Allocate(AH(0:NU-1,NC))
Allocate(AF(0:NU-1,NC))
Allocate(EH(0:NU-1,NC))

YH(:)   = Y(1:NU)
UH(:,:) = U(1:NU,:)


!--- 6.2. Fourier integral operator

R0 = RE + Ymin

CTSgn = 1

Call Canonical_Transform &
  (CTSgn,  & ! <-- Canonical transform direction
   k,      & ! <-- Wave vectors [km^-1]
   Xmax,   & ! <-- X-coordinate of BP plane [km]
   R0,     & ! <-- Shift of Y-grid [km]
   YH,     & ! <-- High resolution Y-grid [km]
   UH)       ! <-> Wave function (Y,channel)


!--- 6.3. Computation of phase and amplitude

Do IC=1,NC
   AH(:,IC)  = Abs(UH(:,IC))
   PHH(:,IC) = ATan2(AImag(UH(:,IC)),Real(UH(:,IC)))
   Call Accumulate_Phase(PHH(:,IC))
End Do


!--- 6.4. Determination of shadow border

Do IC=1,NC

   AF(NU-1,IC) = AH(NU-1,IC)

   Do i=NU-2,0,-1
      AF(i,IC) = AF(i+1,IC) + AH(i,IC)
   End Do

   Do i=NU-1,0,-1
      AF(i,IC) = AF(i,IC)/Sqrt(Real(NU-i))
   End Do

   PminC(IC) = Sum(YH(MaxLoc(AF(:,IC))))

   If (LVrb >= 2) then
      Write(*,'(2X,A,F7.3)') 'DP0  = ', PminC(IC)
   End If

End Do

Pmin = MaxVal(PminC(:))

If (LVrb >= 2) then
   Write(*,'(2X,A,F7.3)') 'Pmin = ', Pmin
End If


!--- 6.5. Computation of refraction angle

Do IC=1,NC

   WF = DYN/DY

   UH(:,IC)  = PHH(:,IC)
   Call Fourier_Filter &
     (WF,         & ! <-- Filtering window width
      UH(:,IC))     ! <-> 1D Data to be filtered
   PHH(:,IC) = UH(:,IC)

   XG = XGPS(1)

   Do i=1,NU-2
      EH(i,IC) = -(PHH(i+1,IC) - PHH(i-1,IC))/(k(IC)*(YH(i+1) - YH(i-1))) + &
                 ASin(- (XG*YH(i) + RE*(XG + Sqrt(XG**2 - YH(i)*(2*RE + YH(i))))) /  &
                        (RE**2 + XG**2))
   End Do

   EH(0,IC)    = EH(1,IC)
   EH(NU-1,IC) = EH(NU-2,IC)

End Do


!--- 6.6. Generating output

IH  = Sum(MaxLoc(YH, YH(:) < Ymax - DYA))
IL  = Sum(MinLoc(YH, YH(:) > 0))
NCT = Size(YH(IL:IH))

Allocate(PCT(NCT))
Allocate(ECT(NCT,NC))
Allocate(ACT(NCT,NC))

PCT(:)      = YH(IL:IH)
ECT(:,1:NC) = EH(IL:IH,1:NC)
ACT(:,1:NC) = AH(IL:IH,1:NC)


!--- 6.7. Memory deallocation

Deallocate(YH)
Deallocate(UH)
Deallocate(PHH)
Deallocate(AH)
Deallocate(AF)
Deallocate(EH)


!----------------------------------------------------------
! 7. PROPAGATION TO ADDITIONAL PHASE SCREEN
!----------------------------------------------------------


!--- 7.1. Propagation to additional screen

If (opt_LS) then

   If (LVrb >= 1) then
      Write(*,'(2X,A)') 'Propagation to additional screen.'
   End If


   !--- 6.1.1. Determination of propagation area

   Ymin = MinVal((PGO(:) + XLS*Sin(EGO(:)))/Cos(EGO(:)) - RE,   &
                 Mask = (-EGO(:) < ELim))
   Ymax = MaxVal((PGO(:) + XLS*Sin(EGO(:)))/Cos(EGO(:)) - RE)

   If (LVrb >= 2) then
      Write(*,'(2X,A,F6.0)') 'YGOmin = ', Ymin
   End If

   AS   = 0.25*MaxVal(XLS*2*Pi/(k(:)*DYN))

   Ymin = Ymin - DYA - AS
   Ymax = Ymax + DYA

   If (LVrb >= 2) then
      Write(*,'(2X,A,F6.0,A,F6.0,A,F6.0)') &
         'Ymin = ', Ymin, ', Ymax = ', Ymax, ', XLS = ', XLS
   End If


   !--- 7.1.2. Setting grid dimension

   NL  = Ceiling(MaxVal((Ymax-Ymin)*k(:)*Eps*2/Pi))
   NL  = Nearest_Power2(NL)
   DYL = (Ymax - Ymin)/(NL - 1)

   If (LVrb >= 2) then
      Write(*,'(2X,A,I8)') 'NL = ', NL
   End If

   Allocate(YL(NL))       ! Y-coordinate in last phase screen
   Allocate(XL(NL))       ! X-coordinate in last phase screen
   Allocate(UL0(NL,NC))   ! Complex field in last screen (Y,channel)
   Allocate(UL(NL,NC))    ! Propagated complex field in last screen (Y,channel)
   Allocate(AL(NL,NC))    ! Amplitude of complex field in last screen (Y,channel)
   Allocate(PL(NL,NC))    ! Phase of complex field in last screen (Y,channel)


   !--- 7.1.3. Interpolation

   XL(:) = XLS

   Do i=1,NL
      YL(i) = Ymin + (i-1)*DYL
   End Do

   Call Interpolate_field &
     (Y(1:NU),     & ! <-- First regular grid
      A(1:NU,:),   & ! <-- Amplitude on first grid (Y,channel)
      Ph(1:NU,:),  & ! <-- Phase on first grid (Y,channel)
      YL(:),       & ! <-- Second regular grid
      AL(:,:),     & ! --> Amplitude on second grid (Y,channel)
      PL(:,:),     & ! --> Phase on second grid (Y,channel)
      UL0)           ! --> Interpolated complex field on second grid

   !!!DEBUG
   !Call PutXY('wop-a.dat',  Y(1:NU:NU/8192),  A(1:NU:NU/8192,1),  Stat=RStat)
   !Call PutXY('wop-al.dat', YL(1:NL:NL/8192), AL(1:NL:NL/8192,1), Stat=RStat)
   !Call PutXY('wop-ph.dat', Y(1:NU:NU/8192),  Ph(1:NU:NU/8192,1), Stat=RStat)
   !Call PutXY('wop-pl.dat', YL(1:NL:NL/8192), PL(1:NL:NL/8192,1), Stat=RStat)
   !Stop

   !--- 7.1.4. Propagation

   Call Apodize(YL(1:NL), UL0(1:NL,1:NC), DYA)

   Do IC=1,NC
      Call Propagate &
        (UL0(1:NL,IC),  & ! <-- Complex field in source plane
         k(IC),         & ! <-- Wave vector
         DYL,           & ! <-- Discretization step
         XLS - Xmax,    & ! <-- Distance to observation plane
         UL(1:NL,IC))     ! --> Propagated field in observation plane
   End Do


   !--- 7.1.5. Calculation of amplitude and phase

   AL(1:NL,:) = Abs(UL(1:NL,:))
   PL(1:NL,:) = Atan2(AImag(UL(1:NL,:)),Real(UL(1:NL,:)))

   Do IC=1,NC
      Call Accumulate_Phase(PL(1:NL,IC))
   End Do

   !!!DEBUG
   !Call PutXY('wop-sal.dat', YL(1:NL:NL/(2**15)), AL(1:NL:NL/(2**15),1), Stat=RStat)
   !Call PutXY('wop-spl.dat', YL(1:NL:NL/(2**15)), PL(1:NL:NL/(2**15),1), Stat=RStat)
   !Stop '**** DEBUG UP TO HERE ****'


!--- 7.2. Use of last screen at Xmax

Else

   NL  =  NU
   DYL =  DY
   YL  => Y(1:NU)
   XL  => X(1:NU)
   AL  => A(1:NU,1:NC)
   PL  => Ph(1:NU,1:NC)
   UL0 => U0(1:NU,1:NC)
   UL  => U(1:NU,1:NC)

End If


!----------------------------------------------------------
! 8. PROPAGATION TO LEO ORBIT
!----------------------------------------------------------


If (LVrb >= 1) then
   Write(*,'(2X,A)') 'Propagation to LEO orbit.'
End If


!--- 8.1. Determination of index grid out of apodization

DYF = 0.125*MinVal(Lam(:))*(XL(1) - 6370.0*0.05)/10.0
IH  = Sum(MaxLoc(YL, YL(:) < Ymax - DYA))
IL  = Sum(MinLoc(YL, YL(:) > Ymin + DYA)) ! + AS))
DI  = Max(1, Floor(DYF/DYL))

If (LVrb >= 2) then
   Write(*,'(2X,2(A, F7.5))') 'DYF = ', DYF, '   DYL = ', DYL
End If


!--- 8.2. Propagation from last screen to LEO orbit

!!!DEBUG
!Call PutXY('test-2-wal.dat', YL(:), AL(:,1), AL(:,2), Stat=Stat)
!Call PutXY('test-2-wpl.dat', YL(:), PL(:,1), PL(:,2), Stat=Stat)


If (opt_FIO) then
   Call Forward_Propagation_FIO &
     (XL(:),            & ! <-- X coordinates of source curve
      YL(:),            & ! <-- Y coordinates of source curve
      AL(:,:),          & ! <-- Amplitude of source field (Y,channel)
      PL(:,:),          & ! <-- Accumulated phase of source field (Y,channel)
      k,                & ! <-- Wave vectors (channel)
      TR(:),            & ! <-- Time for observation curve
      XLEO,             & ! <-- X coordinates of LEO
      YLEO,             & ! <-- Y coordinates of LEO
      XImin,            & ! <-- GO estimate of minimum spatial frequency
      XImax,            & ! <-- GO estimate of maximum spatial frequency
      SR,               & ! <-- Sampling rate [Hz]
      opt_TS,           & ! <-- Time scaling
      WT,               & ! --> WP data time
      WX,               & ! --> X coordinates of observation curve
      WY,               & ! --> Y coordinates of observation curve
      WA,               & ! --> Amplitudes of observed field (time,channel)
      WP,               & ! --> Accumulated phase of observed field (time,channel)
      LVrb)               ! <~~ Verbosity mode
Else
   Call Forward_Propagation_Fresnel &
     (XL(IL:IH:DI),     & ! <-- X coordinates of source curve
      YL(IL:IH:DI),     & ! <-- Y coordinates of source curve
      AL(IL:IH:DI,:),   & ! <-- Amplitude of source field (Y,channel)
      PL(IL:IH:DI,:),   & ! <-- Accumulated phase of source field (Y,channel)
      k,                & ! <-- Wave vectors (channel)
      TR,               & ! <-- Time for observation curve
      XLEO,             & ! <-- X coordinates of LEO
      YLEO,             & ! <-- Y coordinates of LEO
      FZ,               & ! <-- Initial last Fresnal zone size [rad]
      SR,               & ! <-- Sampling rate [Hz]
      opt_TS,           & ! <-- Time scaling
      WT,               & ! --> WP data time
      WX,               & ! --> X coordinates of observation curve
      WY,               & ! --> Y coordinates of observation curve
      WA,               & ! --> Amplitudes of observed field (time,channel)
      WP,               & ! --> Accumulated phase of observed field (time,channel)
      LVrb)               ! <~~ Verbosity mode
End If

!!!DEBUG
!Call PutXY('test-2-wa.dat', WT(:), WA(:,1), WA(:,2), Stat=Stat)

!--- 8.3. Data array allocation

NW = Size(WT)

Allocate(PRLEO(1:NW))
Allocate(PRGPS(1:NW))
Allocate(WRLEO(1:NW))
Allocate(WRGPS(1:NW))
Allocate(WVLEO(1:NW))
Allocate(WVGPS(1:NW))
Allocate(WS(1:NW,NC))


!--- 8.4. Calculation of LEO and GPS coordinates and velocities

If (opt_ECEF) then
   PhiN = 0.0
Else
   PhiN =  GAST &
     (Year,    & ! <-- Year of occultation
      Month,   & ! <-- Month
      Day,     & ! <-- Day
      Hour,    & ! <-- Hour
      Minute,  & ! <-- Minute
      Second,  & ! <-- Second of beginning of occultation
      TR(Iocc))  ! <-- Second of occultation duration
End If

Do i=1,NW
   If (opt_ECEF) then
      PhiT = 0.0
   Else
      PhiT =  GAST &
        (Year,    & ! <-- Year of occultation
         Month,   & ! <-- Month
         Day,     & ! <-- Day
         Hour,    & ! <-- Hour
         Minute,  & ! <-- Minute
         Second,  & ! <-- Second of beginning of occultation
         WT(i))     ! <-- Second of occultation duration
   End If
   PRLEO(i) = RLC + AX*WX(i)      + AY*(RE + WY(i))
   PRGPS(i) = RLC + AX*XGPS(Iocc) + AY*(RE + YGPS(Iocc))
   PRLEO(i) = Rotate(PRLEO(i), PA, PhiT-PhiN)
   PRGPS(i) = Rotate(PRGPS(i), PA, PhiT-PhiN)
End Do

Call Satellite_Velocities &
  (WT,       & ! <-- Relative time of samples [sec]
   PRLEO,    & ! <-- LEO coordinates
   PRGPS,    & ! <-- GPS coordinates
   WRLEO,    & ! --> LEO coordinates from regression
   WVLEO,    & ! --> LEO velocities from regression
   WRGPS,    & ! --> GPS coordinates from regression
   WVGPS)      ! --> GPS velocities from regression


!--- 8.5. Calculation of phase excess

Do IC=1,NC
   WS(1,IC) = 1e3*(WP(1,IC)/k(IC) - Vector_Norm(WRLEO(1)-WRGPS(1)))
   Do i=NW,1,-1
      WS(i,IC) = 1e3*(WP(i,IC)/k(IC) - Vector_Norm(WRLEO(i)-WRGPS(i))) - WS(1,IC)
   End Do
End Do


!----------------------------------------------------------
! 9. MEMORY DEALLOCATION
!----------------------------------------------------------

Deallocate(ERLEO)
Deallocate(ERGPS)
Deallocate(PRLEO)
Deallocate(PRGPS)
Deallocate(XLEO)
Deallocate(YLEO)
Deallocate(XGPS)
Deallocate(YGPS)
Deallocate(WX)
Deallocate(WY)
Deallocate(WP)
Deallocate(U0)
Deallocate(U)
Deallocate(Ph)
Deallocate(A)
Deallocate(Y)
Deallocate(X)
Deallocate(RN)
Deallocate(RNI)
Deallocate(PGO)
Deallocate(EGO)
Deallocate(AGO)
Deallocate(AP)
Deallocate(k)
Deallocate(Ak)
Deallocate(Lam)
Deallocate(At)
Deallocate(PminC)

If (opt_LS) then
   Deallocate(YL)
   Deallocate(XL)
   Deallocate(UL0)
   Deallocate(UL)
   Deallocate(AL)
   Deallocate(PL)
End If


Return


!----------------------------------------------------------


Contains


!----------------------------------------------------------
Subroutine Make_Refractivity(X, Y, RN, RNI, YE)
!
! Calculation of refractivity profile in phase screen:
! 1. Calculation of refractivity profile with a lower
!    resolution;
! 2. Linear interpolation to high resolution grid.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
Real(Double), Intent(In)    :: &
   X         ! X position of phase screen
Real(Double), Intent(In)    :: &
   Y(:)      ! Y grid in phase screen
Real(Double), Intent(Out)   :: &
   RN(:)     ! Refractivity profile in phase screen
Real(Double), Intent(Out)   :: &
   RNI(:,:)  ! Imaginary part of refractivity in phase screen (Y,channel)
Real(Double), Intent(Out)   :: &
   YE        ! Y coordinate of Earth border
!----------------------------------------------------------
! Local Scalars:
!
Integer         :: i           ! Array index
Integer         :: IC          ! Channel number
Type(Cartesian) :: XY          ! Earth frame coordinates of point
Type(Cartesian) :: NGradN      ! N*Grad(N)
Integer         :: Stat        ! ECHAM_NGradN error status
!
! Local Arrays:
!
Real(Double)    :: YL(NN)      ! Low resolution grid
Real(Double)    :: RNL(NN)     ! Low resolution refractivity profile
Real(Double)    :: RNIL(NN,NC) ! Low resolution imaginary refractivity (Y,channel)
!----------------------------------------------------------
! Global variables used:
!
! NN   - number of points in lowres refractivity profiles.
! ERLC - curvature center (ECEF).
! EX   - X basis vector of occultation plane.
! EY   - Y basis vector of occultation plane.
! Ymin - lower Y coordinate of propagation area.
! Ymax - upper Y coordinate of propagation area.
!----------------------------------------------------------


YE = Ymin

!! JDH Debug 031521
!!********************************
!print*,'NN=',NN

!!if(nn.gt.100) then
!!if(nn.gt.2) then
!if(nn.ge.1) then
!!   print*,'nn greater than 100'
!!   print*,'nn greater than 2'
!   print*,'nn greater than 1'
!   stop
!end if
!!*******************************

Do i=1,NN
   YL(i)  = Ymin + Real(i-1)*(Ymax-Ymin)/Real(NN-1)
   XY     = ERLC + X*EX + (RE+YL(i))*EY
   Call Atmosphere_NGradN &
     (XY,         & ! <-- Cartesian coordinates of point
      NGradN,     & ! --> Interpolated (1 + N)*Grad(N)
      RNL(i),     & ! --> Interpolated N
      Stat,       & ! --> Error status
      RNIL(i,:))    ! ~~> Interpolated Im(N)
   If (Stat /= 0) then
      YE = Max(YE, YL(i))
   End If
End Do

Do i=1,Size(Y)
   Call Linear &
     (YL,      & ! <-- Argument grid
      RNL,     & ! <-- Gridded function
      Y(i),    & ! <-- Interpolation point
      RN(i))     ! --> Interpolated function value
   Do IC=1,NC
      Call Linear &
        (YL,           & ! <-- Argument grid
         RNIL(:,IC),   & ! <-- Gridded function
         Y(i),         & ! <-- Interpolation point
         RNI(i,IC))      ! --> Interpolated function value
   End Do
End Do

End Subroutine Make_Refractivity



!----------------------------------------------------------
Subroutine Double_Resolution(Y, Ph, U, RN, DY, NU)
!
! Doubling resolution of complex field and refractivity
! profile in phase screen:
! 1. Linear interpolation for accumulated phase and
!    refractivity between points.
! 2. Using formula U_1/2 = U0*Sqrt(U1/U0) for complex
!    field.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
Real(Double), Intent(InOut)    :: &
   Y(:)      ! Y grid in phase screen
Real(Double), Intent(InOut)    :: &
   Ph(:,:)   ! Accumulated phase of complex field (Y,channel)
Real(Double), Intent(InOut)    :: &
   RN(:)     ! Refractivity profile in phase screen
Complex(Double), Intent(InOut) :: &
   U(:,:)    ! Complex field in phase screen (Y,channel)
Real(Double), Intent(InOut)    :: &
   DY        ! Resolution of Y grid
Integer, Intent(InOut)         :: &
   NU        ! Number of points in phase screen
!----------------------------------------------------------
! Local Scalars:
!
Integer :: i   ! Array index
!----------------------------------------------------------


Do i=NU,1,-1
   Y(2*i-1)    = Y(i)
   Ph(2*i-1,:) = Ph(i,:)
   U(2*i-1,:)  = U(i,:)
   RN(2*i-1)   = RN(i)
End Do

Y(2*NU)    = Y(2*NU-1) + DY/2
Ph(2*NU,:) = Ph(2*NU-1,:)
U(2*NU,:)  = U(2*NU-1,:)
RN(2*NU)   = RN(2*NU-1)

Do i=1,NU-1
   Y(2*i)    = (Y(2*i-1) + Y(2*i+1))/2
   Ph(2*i,:) = (Ph(2*i-1,:) + Ph(2*i+1,:))/2
   Where (U(2*i+1,:) == 0)
      U(2*i,:) = 0
   Elsewhere
      U(2*i,:) = Sqrt(U(2*i-1,:)/U(2*i+1,:))*U(2*i+1,:)
   End Where
   RN(2*i)   = (RN(2*i-1) + RN(2*i+1))/2
End Do

NU = 2*NU
DY = DY/2

Do i=1,NU
   Y(i)    = Ymin + (i-1)*DY
End Do

End Subroutine Double_Resolution



!----------------------------------------------------------
Subroutine Save_Ckpt(DT)
!
! Saving checkpoint for a future restart.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
Integer, Intent(In) :: &
   DT       ! Interval from previous checkpoint (min)
!----------------------------------------------------------
! Local Scalars:
!
Logical, Save :: FirstCall = .True. ! Indicator of first call
Integer       :: ETM                ! Elapsed time in minutes
Integer       :: IU                 ! Unit number
Integer       :: Stat               ! Error status
!
! Local Arrays:
!
Integer, Save :: &
   TS(8),     & ! Time of last save
   TC(8),     & ! Current time
   ET(4)        ! Elapsed time
!----------------------------------------------------------

If (FirstCall) then
   Call Date_and_Time(Values = TS(:))
   FirstCall = .False.
End If

Call Date_and_Time(Values = TC(:))
Call Elapsed_Time &
  (TS,        & ! <-- Start time
   TC,        & ! <-- End time
   ET(1),     & ! --> Days
   ET(2),     & ! --> Hours
   ET(3),     & ! --> Minutes
   ET(4))       ! --> Seconds

ETM = ET(4)/60 + ET(3) + 60*(ET(2) + 12*ET(1))

If (ETM >= DT) then

If (LVrb >= 3) then
   Write (Line,'(2X,A30,A1)') '** Saving field **', Char(0)
   Call CPrintf(Line)
End If

Call GetFreeUnit  &
  (IU,       & ! --> Unit number
   Stat)       ! --> Error code

Open(Unit=IU,                   &
     File    = CkptName,        &
     Form    = 'UNFORMATTED',   &
     IOStat  = Stat)

If (Stat /= 0) then
   Write (Line,'(2X,A30,A1)') '** Checkpoint not taken **', Char(0)
   Call CPrintf(Line)
   Return
End If

Write(IU) NC, NU
Write(IU) X(1:NU)
Write(IU) Y(1:NU)
Write(IU) U(1:NU,1:NC)

Close(Unit=IU)

If (LVrb >= 3) then
   Write (Line,'(2X,A30,A1)') '** Field saved **', Char(0)
   Call CPrintf(Line)
End If

Call Date_and_Time(Values = TS(:))

End If

End Subroutine Save_Ckpt


!----------------------------------------------------------
Subroutine Restore_Ckpt
!
! Saving checkpoint for a future restart.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Local Scalars:
!
Integer       :: IU                 ! Unit number
Integer       :: Stat               ! Error status
!----------------------------------------------------------

If (LVrb >= 3) then
   Write (Line,'(2X,A30,A1)') '** Reading field **', Char(0)
   Call CPrintf(Line)
End If

Call GetFreeUnit  &
  (IU,       & ! --> Unit number
   Stat)       ! --> Error code

Open(Unit=IU,                 &
     File   = CkptName,       &
     Status ='OLD',           &
     Form   ='UNFORMATTED',   &
     IOStat = Stat)
!Convert='BIG_ENDIAN')

If (Stat /= 0) then
   Write (Line,'(2X,A30,A1)') '** Not restarted **', Char(0)
   Call CPrintf(Line)
   Return
End If

Read(IU) NC, NU
Read(IU) X(1:NU)
Read(IU) Y(1:NU)
Read(IU) U(1:NU,1:NC)

Close(Unit=IU)

DY         = Y(2) - Y(1)
Ymin       = Y(1)
Ymax       = Y(NU)
U0(1:NU,:) = U(1:NU,:)

A(1:NU,:)  = Abs(U(1:NU,:))
Ph(1:NU,:) = Atan2(AImag(U(1:NU,:)),Real(U(1:NU,:)))

Do IC=1,NC
   Call Accumulate_Phase(Ph(1:NU,IC))
End Do

If (LVrb >= 3) then
   Write (Line,'(2X,A30,A1)') '** Field read **', Char(0)
   Call CPrintf(Line)
End If

End Subroutine Restore_Ckpt


End Subroutine Wave_Propagate



!==========================================================
Subroutine Interpolate_field &
  (Y1,        & ! <-- First regular grid
   A1,        & ! <-- Amplitude on first grid
   P1,        & ! <-- Phase on first grid
   Y2,        & ! <-- Second regular grid
   A2,        & ! --> Amplitude on second grid
   P2,        & ! --> Phase on second grid
   U2)          ! --> Interpolated complex field on second grid
!
! Linear interpolation of amplitude and phase
! from a regular grid to another regular grid.
!----------------------------------------------------------
! (C) Copyright 2002, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 05 Dec 2002 | Original version.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In) :: &
   Y1(:)        ! First regular grid
!
Real(Double), Intent(In) :: &
   A1(:,:)      ! Amplitude on first grid
!
Real(Double), Intent(In) :: &
   P1(:,:)      ! Phase on first grid
!
Real(Double), Intent(In) :: &
   Y2(:)        ! Second grid
!
! Output arguments:
!
Real(Double), Intent(Out)    :: &
   A2(:,:)      ! Amplitude on first grid
!
Real(Double), Intent(Out)    :: &
   P2(:,:)      ! Phase on first grid
!
Complex(Double), Intent(Out) :: &
   U2(:,:)      ! Interpolated complex on second grid
!----------------------------------------------------------
! Local Parameters:
!
Complex(Double), Parameter ::  &
   Ci = (0.0_Double, 1.0_Double)  ! I = Sqrt(-1)
!
! Local Scalars:
!
Integer      :: N1    ! Dimension of first grid
Integer      :: N2    ! Dimension of second grid
Integer      :: i1    ! Index of first grid
Integer      :: i2    ! Index of second grid
Real(Double) :: DY1   ! Step of first grid
Real(Double) :: DY2   ! Step of second grid
!----------------------------------------------------------


!----------------------------------------------------------
! 1. INITIALIZATION
!----------------------------------------------------------

N1 = Size(A1,2)
N2 = Size(U2,2)

DY1 = Y1(2) - Y1(1)
DY2 = Y2(2) - Y2(1)


!----------------------------------------------------------
! 2. INTERPOLATION
!----------------------------------------------------------

Do i2=1,N2
   i1  = 1 + Floor((Y2(i2) - Y1(1))/DY1)
   If ((i1 >= 1) .and. (i1 <= N1-1)) then
      A2(i2,:) = (A1(i1,:)*(Y1(i1+1) - Y2(i2)) + A1(i1+1,:)*(Y2(i2) - Y1(i1)))/DY1
      P2(i2,:) = (P1(i1,:)*(Y1(i1+1) - Y2(i2)) + P1(i1+1,:)*(Y2(i2) - Y1(i1)))/DY1
      U2(i2,:) = A2(i2,:)*Exp(Ci*Modulo(P2(i2,:),2*Pi))
   Else
      A2(i2,:) = 0
      P2(i2,:) = 0
      U2(i2,:) = 0
   End If
End Do


End Subroutine Interpolate_field



!==========================================================
Subroutine Boundary_Condition &
  (Y,         & ! <-- Grid
   YE,        & ! <-- Earth's boundary
   X,         & ! <-- X-coordinate of phase screen
   DYA,       & ! <-- Apodization border width
   U,         & ! <-> Complex field
   YBC,       & ! --> Boundary condition location
   IR)          ! --> Reduced grid index
!
! Wave field continuation beneath the
! Earth's surface according to
! the boundary condition.
!----------------------------------------------------------
! Method:
!   Multiplication of field below surface by 1d-4,
!   reduction of grid below surface.
!----------------------------------------------------------
! (C) Copyright 2002, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0.a | 05 Dec 2002 | Original version.
!   1.1   | 03 Sen 2004 | Suppressing field in dead zone.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In) :: &
   Y(:)         ! Grid
!
Real(Double), Intent(In) :: &
   YE           ! Earth's boundary
!
Real(Double), Intent(In) :: &
   X            ! X-coordinate of phase screen
!
Real(Double), Intent(In) :: &
   DYA          ! Apodization border width
!
! Inout arguments:
!
Complex(Double), Intent(InOut) :: &
   U(:,:)       ! Complex field
!
! Output arguments:
!
Real(Double), Intent(Out) :: &
   YBC          ! Boundary condition location
!
Integer, Intent(Out) :: &
   IR           ! Reduced grid index
!----------------------------------------------------------
! Local Parameters:
!
Real(Double), Parameter :: &
   YB = 0.5     ! Boundary zone size
!
Integer, Parameter :: &
   Frac(2,3) = Reshape((/ 3,  4,      &
                          3,  8,      &
                          3, 16 /),   &
                       (/ 2, 3 /))
!
! Local Scalars:
!
Integer      :: NU  ! Dimension of field
Integer      :: i   ! Grid index
Real(Double) :: YR  ! Reduced grid limit
Integer      :: IR0 ! Initial estimation of IR
Integer      :: IRC ! Currently tested IR
!----------------------------------------------------------


!----------------------------------------------------------
! 1. INITIALIZATION
!----------------------------------------------------------

NU = Size(Y)


!----------------------------------------------------------
! 2. BOUNDARY CONDITION
!----------------------------------------------------------

If (X < 0.0) then
   YBC = 0.0
Else
   YBC = YE
End If

If (X < 0.0) then
   Do i=1,NU
      If (Y(i) < 0.0) then
         U(i,:) = U(i,:)*Exp(-(Y(i)/YB)**2)
      End If
   End Do
End If

Do i=1,NU
   If (Y(i) < YE) then
      U(i,:) = U(i,:)*1d-4
   End If
End Do


!----------------------------------------------------------
! 3. GRID REDUCTION
!----------------------------------------------------------

If (X < 0.0) then
   YR = -DYA
Else
   YR = YE - DYA
End If

IR0 = Sum(MaxLoc(Y(:), Mask = (Y(:) < YR)))

IR  = 1

Do i=1,3
   IRC = NU - Frac(1,i)*NU/Frac(2,i) + 1
   If (IRC > IR0) then
      Exit
   Else
      IR = IRC
   End If
End Do

U(1:IR-1,:) = 0


End Subroutine Boundary_Condition



End Module Wave_Propagator


