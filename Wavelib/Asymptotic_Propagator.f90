!
Module Asymptotic_Propagator
!
! Asymptotic model radio occultations
!----------------------------------------------------------
! (C) Copyright 2003, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 25 Apr 2003 | Original version.
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
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
!
Contains


!==========================================================
Subroutine Asymptotic_Propagate &
  (Year,       & ! <-- Occultation year
   Month,      & ! <-- Occultation month
   Day,        & ! <-- Occultation day
   Hour,       & ! <-- Occultation hour
   Minute,     & ! <-- Occultation begin minute
   Second,     & ! <-- Occultation begin second
   SDim,       & ! <-- Spatial dimension
   TR,         & ! <-- Relative time of samples [sec]
   Freq,       & ! <-- Frequencies [Hz]
   AFreq,      & ! <-- Frequencies for absorption [Hz]
   RLEO,       & ! <-- LEO coordinates (J2000/ECEF)
   RGPS,       & ! <-- GPS coordinates (J2000/ECEF)
   opt_ECEF,   & ! <-- Coordinates in ECEF frame
   Hmax,       & ! <-- Maximum height
   DYN,        & ! <-- Minimum vertical scale of N
   NGO,        & ! <-- Upper number of GO rays
   SR,         & ! <-- Sampling rate [Hz]
   opt_AB,     & ! <-- Modeling of absorption
   FWP,        & ! <-- Filter width for computing DPGPS/DP [km]
   EP,         & ! --> Refraction angles
   P,          & ! --> Impact parameters
   LAT,        & ! --> Logarithmic attenuation (p,channel)
   Z,          & ! --> Ray heights
   GP,         & ! --> Occultation point
   ERLC,       & ! --> Local curvature center (ECEF)
   RE,         & ! --> Local curvature radius
   WT,         & ! --> WP data time
   WA,         & ! --> WP amplitudes (time,channel)
   WS,         & ! --> WP phase excess (time,channel) [m]
   WRLEO,      & ! --> WP LEO coordinates (J2000/ECEF)
   WVLEO,      & ! --> WP LEO velocity (J2000/ECEF)
   WRGPS,      & ! --> WP GPS coordinates (J2000/ECEF)
   WVGPS,      & ! --> WP GPS velocity (J2000/ECEF)
   Stat,       & ! ~~> Error status
   Vrb)          ! <~~ Verbosity level
!
! Asymptotic modeling of radio occultation.
!----------------------------------------------------------
! Method:
!   Maslov operator.
!----------------------------------------------------------
! (C) Copyright 2003-2009, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 22 May 2003 | Original version.
!   1.1   | 21 Jul 2003 | Spline interpolation of GO functions.
!   2.0   | 20 Aug 2003 | Verbosity level.
!   2.1   | 05 Sep 2003 | 1) Imag --> AImag;
!         |             | 2) Bugs corrected;
!         |             | 3) Parameter DT.
!   3.0   | 30 Oct 2003 | DYN.
!   3.1   | 17 Dec 2003 | LAT.
!   3.2   | 25 May 2004 | AFreq.
!   3.3   | 29 Nov 2005 | Corrected processing of status.
!   3.4   | 11 Dec 2005 | Excluded variable Tdir.
!   3.5   | 24 Feb 2007 | FWP.
!   3.6   | 23 Sep 2009 | Index order changed to (time,channel).
!   3.7   | 29 Sep 2009 | Output phase defined to be positive.
!   3.8   | 08 Oct 2009 | Updated for new version of
!         |             | Interpolate_Trajectory.
!----------------------------------------------------------
! Modules used:
!
Use Coordinates, only: &
! Imported Type Definitions:
    Cartesian,      &
! Imported Routines:
    Rotate,         &
! Imported Operators:
    Operator(+),    &
    Operator(-),    &
    Operator(*)
!
Use Earth, only: &
! Imported Type Definitions:
    Geodetic,              &
! Imported Routines:
    GAST
!
Use Occ_Coordinates, only: &
! Imported Routines:
    Occ_Geometry,           &
    Occ_Geometry_ECEF,      &
    Satellite_Velocities,   &
    Interpolate_Trajectory
!
Use Ray_Problem, only: &
! Imported Routines:
    Solve_Boundary_Problem
!
Use Interpolation, only: &
! Imported Routines:
    Linear,           &
    Init_Spline,      &
    Spline,           &
    Nearest_Power2
!
Use FFTW, only: &
! Imported Routines:
    FFT1
!
Use Occ_Diffraction, only: &
! Imported Routines:
    Accumulate_Phase
!
Use IO, only: &
! Imported Routines:
    PutXY
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
Integer, Intent(In)     :: &
   SDim       ! Spatial dimension
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
Logical, Intent(In) :: &
   opt_ECEF   ! .True.  - Coordinates in ECEF frame
              ! .False. - Coordinates in absolute frame
!
Real(Double), Intent(In) :: &
   Hmax       ! Maximum height
!
Real(Double), Intent(In) :: &
   DYN        ! Minimum vertical scale of N
!
Integer, Intent(In) :: &
   NGO        ! Upper number of GO rays
!
Real(Double), Intent(In) :: &
   SR         ! Sampling rate [Hz]
!
Logical, Intent(In)      :: &
   opt_AB     ! Modeling of absorption
!
Real(Double), Intent(In) :: &
   FWP        ! Filter width for computing DPGPS/DP [km]
!
! Output arguments:
!
Real(Double), Pointer     :: &
   EP(:)      ! Refraction angles
!
Real(Double), Pointer     :: &
   P(:)       ! Impact parameters
!
Real(Double), Pointer     :: &
   LAT(:,:)   ! Logarithmic attenuation (p,channel)
!
Real(Double), Pointer     :: &
   Z(:)       ! Ray heights
!
Type(Geodetic), Intent(Out)  :: &
   GP       ! Occultation point
!
Type(Cartesian), Intent(Out)  :: &
   ERLC     ! Local curvature center (ECEF)
!
Real(Double), Intent(Out)  :: &
   RE       ! Local curvature radius
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
Integer, Intent(Out)  :: &
   Stat       ! Error status
!
! Input optional arguments:
!
Integer, Optional, Intent(In) :: &
   Vrb        ! Verbosity level
              ! 0 by default.
!----------------------------------------------------------
! Local Parameters:
!
Integer, Parameter      :: &
   NV = 5                         ! Polynomial degree for calculation of velocity
!
Complex(Double), Parameter ::  &
   Ci = (0.0_Double, 1.0_Double)  ! I = Sqrt(-1)
!
Real(Double), Parameter :: &
   DHB = 6.0_Double,           &  ! Safety border width
   DT  = 5.0_Double               ! Safety time border
!
Type(Cartesian), Parameter :: &
   PA =  Cartesian((/0,0,1/))     ! Polar axis
!
! Local Scalars:
!
! --- Geometry of occultation
!
Integer         :: N        ! Number of data
Type(Cartesian) :: RLC      ! Local curvature center (J2000)
Real(Double)    :: Phi      ! Rotation angle from ECEF to J2000
!
! --- Geometric optics
!
Real(Double)    :: Pmin     ! Minimum impact parameter
Real(Double)    :: Pmax     ! Maximum impact parameter
Real(Double)    :: Ymax     ! Maximum momentum
Real(Double)    :: PZmin    ! Minimum impact parameter for Z-grid
Real(Double)    :: PZmax    ! Maximum impact parameter for Z-grid
!
! --- Wave optics
!
Integer         :: NC       ! Number of channels
Integer         :: IC       ! Channel index
Integer         :: NH       ! Number of hi-res grid points
Real(Double)    :: YPI      ! Interpolated YP
Integer         :: FSign    ! FFT direction
Real(Double)    :: FTI      ! Interpolated phase model
Real(Double)    :: S0I      ! Interpolated satellite-to-satellite distance
Real(Double)    :: PHB      ! Base phase excess
Real(Double)    :: PB       ! Border width for increase of spectral resolution
Integer         :: IGmin    ! Index of minimum output time
Integer         :: IGmax    ! Index of maximum output time
Real(Double)    :: TGmin    ! Minimum TG
Real(Double)    :: TGmax    ! Maximum TG
Integer         :: NP       ! Number of GO solution samples
Integer         :: NWT      ! Number of wave field samples
Real(Double)    :: YTI      ! Interpolated Y
!
! --- Work variables
!
Integer         :: LVrb     ! Verbosity level
Integer         :: i        ! Array index
!
! Local Arrays:
!
!
! --- Geometry of occultation
!
Type(Cartesian), Allocatable :: &
   ERLEO(:),     & ! LEO coordinates (ECEF)
   ERGPS(:)        ! GPS coordinates (ECEF)
Type(Cartesian), Allocatable :: &
   XLEO(:),      & ! LEO coordinates from regression (J2000)
   VLEO(:),      & ! LEO velocities from regression (J2000)
   XGPS(:),      & ! GPS coordinates from regression (J2000)
   VGPS(:)         ! GPS velocities from regression (J2000)
Type(Cartesian), Allocatable :: &
   WERLEO(:),    & ! LEO coordinates (ECEF) for wave optics
   WEVLEO(:),    & ! LEO velocities (ECEF) for wave optics
   WERGPS(:),    & ! GPS coordinates (ECEF) for wave optics
   WEVGPS(:)       ! GPS velocities (ECEF) for wave optics
Real(Double)    :: &
   BG(0:NV,3),   & ! Regression coefficients for VGPS
   BL(0:NV,3)      ! Regression coefficients for VLEO
!
! --- Geometric optics
!
Real(Double), Pointer     :: &
   TP(:),        & ! Time as function of P
   YP(:),        & ! Y-coordinate as function of P
   ARP(:),       & ! Refractive amplitude as function of P
   ATP(:,:)        ! Absorptive attenuation as function of P (p,channel)
Real(Double), Allocatable :: &
   YP2(:),       & ! Spline coefficients for YP
   ARP2(:),      & ! Spline coefficients for ARP
   ATP2(:,:)       ! Spline coefficients for ATP (p,channel)
Real(Double), Pointer     :: &
   TG(:),        & ! Time for GO solution
   YT(:),        & ! Y as function of TG
   S0(:),        & ! Satellite-to-satellite distnce
   FT(:)           ! Model function F(TG)
Real(Double), Allocatable :: &
   YT2(:),       & ! Spline coefficients for YT
   S02(:),       & ! Spline coefficients for S0
   FT2(:)          ! Spline coefficients for FT
!
! --- Wave optics
!
Real(Double), Allocatable  :: &
   k(:),         & ! Wavenumbers
   Ak(:),        & ! Wavenumbers for absorption
   ZH(:),        & ! Hi-res grid of ray heights (P-PZmin)
   ARH(:),       & ! Refractive amplitude in p-representation
   ATH(:),       & ! Absorptive amplitude in p-representation
   PH(:),        & ! Eikonal of wave field
   YH(:),        & ! Hi-res grid of Y-coordinate
   AH(:)           ! Amplitude of wave field
Complex(Double), Allocatable :: &
   UH(:)           ! Complex wave field
!----------------------------------------------------------


!----------------------------------------------------------
! 0. INITIALIZATION
!----------------------------------------------------------


!--- 0.1. Determination of number of input data and channels

N  = Size(TR)
NC = Size(Freq)


!--- 0.2. Setting verbosity level

If (Present(Vrb)) then
   LVrb = Vrb
Else
   LVrb = 0
End If


!--- 0.3. Wavenumbers

Allocate(k(NC))
Allocate(Ak(NC))

k(:)  = 2*Pi*Freq(:)/C_Light
Ak(:) = 2*Pi*AFreq(:)/C_Light


!----------------------------------------------------------
! 1. COORDINATE TRANSFORMS
!----------------------------------------------------------


!--- 1.1. Calculation of curvature center and
!---      coordinates (ECEF)

Allocate(ERLEO(N))
Allocate(ERGPS(N))

If (opt_ECEF) then
   ERLEO(:) = RLEO(:)
   ERGPS(:) = RGPS(:)
   Call Occ_Geometry_ECEF &
     (ERLEO,     & ! <-- LEO coordinates (Earth frame)
      ERGPS,     & ! <-- GPS coordinates (Earth frame)
      GP,        & ! --> Geodetic coordinates of occultation point
      ERLC,      & ! --> Curvature center (Earth frame)
      RE)          ! --> Local curvature radius
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
   If (Stat /= 0) then
      Return
   End If
End If


!--- 1.2. Regression coefficients for satellite trajectories

Allocate(XLEO(N))
Allocate(VLEO(N))
Allocate(XGPS(N))
Allocate(VGPS(N))

Call Satellite_Velocities &
  (TR,       & ! <-- Relative time of samples [sec]
   ERLEO,    & ! <-- LEO coordinates (J2000)
   ERGPS,    & ! <-- GPS coordinates (J2000)
   XLEO,     & ! --> LEO coordinates from regression (J2000)
   VLEO,     & ! --> LEO velocities from regression (J2000)
   XGPS,     & ! --> GPS coordinates from regression (J2000)
   VGPS,     & ! --> GPS velocities from regression (J2000)
   Stat,     & ! ~~> Error status
   BL,       & ! ~~> LEO regression coefficients
   BG)         ! ~~> GPS regression coefficients

If (Stat /= 0) then
   Return
End If


!----------------------------------------------------------
! 2. GEOMETRIC OPTICAL PROPAGATION
!----------------------------------------------------------

If (LVrb >= 1) then
   Write(*,'(2X,A)') 'Geometric optical propagation'
End If


!--- 2.1. Solving ray boundary problem

Call Solve_Boundary_Problem &
  (SDim,      & ! <-- Spatial dimension
   TR,        & ! <-- Relative time of samples [sec]
   BL,        & ! <-- LEO regression coefficients
   BG,        & ! <-- GPS regression coefficients
   Hmax+DHB,  & ! <-- Maximum height
   ERLC,      & ! <-- Curvature center (ECEF)
   RE,        & ! <-- Local curvature radius
   DYN,       & ! <-- Minimum vertical scale of N
   NGO,       & ! <-- Upper number of GO rays
   opt_AB,    & ! <-- Modeling of absorption
   NC,        & ! <-- Number of channels for absorption
   FWP,       & ! <-- Filter width for computing DPGPS/DP [km]
   P,         & ! --> Impact parameters
   EP,        & ! --> Refraction angles
   TP,        & ! --> Time as function of P
   YP,        & ! --> Y-coordinate as function of P
   ARP,       & ! --> Refractive amplitude as function of P
   ATP,       & ! --> Absorptive attenuation as function of P (p,channel)
   TG,        & ! --> Time for GO solution
   S0,        & ! --> Satellite-to-satellite distnce
   YT,        & ! --> Y-coordinate as function of TG
   FT,        & ! --> Model function F(TG)
   Vrb)         ! <~~ Verbosity level

   !!!DEBUG
   !i = NH/4096
   !Call PutXY('test-ayh.dat', YH(1::i), AH(1::i), Stat=Stat)
   !Call PutXY('test-pyh.dat', YH(1::i), PH(1::i), Stat=Stat)
   !Call PutXY('test-FTY.dat', YT(:), FT(:), Stat=Stat)
   !Stop 'DEBUG to HERE...'

NP = Size(P)

Allocate(Z(NP))
Allocate(LAT(NP,NC))

Z(:) = P(:) - RE

Do IC=1,NC
   LAT(:,IC) = Ak(IC)*ATP(:,IC)
End Do


!--- 2.2. Computation of time for output

IGmin = Sum(MaxLoc(P(:), P(:)-RE <= Hmax))
TGMin = TP(IGmin)

IGmax = Sum(MinLoc(P(:)))
TGmax = TP(IGmax)
TGmax = TGmax + Sign(DT, TGmax - TGmin)

NWT   = 1 + Ceiling(SR*Abs(TGmax - TGmin))

Allocate(WT(NWT))          ! --> WP data time
Allocate(WA(NWT,NC))       ! --> WP amplitudes (time,channel)
Allocate(WS(NWT,NC))       ! --> WP phase excess (time,channel) [m]
Allocate(WRLEO(NWT))       ! --> WP LEO coordinates (J2000)
Allocate(WVLEO(NWT))       ! --> WP LEO velocity (J2000)
Allocate(WRGPS(NWT))       ! --> WP GPS coordinates (J2000)
Allocate(WVGPS(NWT))       ! --> WP GPS velocity (J2000)

Do i=1,NWT
   WT(i) = Min(TGmax, TGmin) + Real(i-1)/SR
End Do


!--- 2.3. Computation of trajectories in ECEF frame

Allocate(WERLEO(NWT))      ! LEO coordinates (ECEF) for wave optics
Allocate(WEVLEO(NWT))      ! LEO velocities (ECEF) for wave optics
Allocate(WERGPS(NWT))      ! GPS coordinates (ECEF) for wave optics
Allocate(WEVGPS(NWT))      ! GPS velocities (ECEF) for wave optics

Do i=1,NWT
   Call Interpolate_Trajectory  &
     (TR,          & ! <-- Relative time of samples [sec]
      BL,          & ! <-- LEO regression coefficients
      BG,          & ! <-- GPS regression coefficients
      WT(i),       & ! <-- Time
      WERLEO(i),   & ! --> LEO coordinates
      WEVLEO(i),   & ! --> LEO velocity
      WERGPS(i),   & ! --> GPS coordinates
      WEVGPS(i))     ! --> GPS velocity
End Do


!--- 2.4. Computation of trajectories in J2000 frame

Do i=1,NWT
   If (opt_ECEF) then
      Phi = 0.0
   Else
      Phi =  GAST &
        (Year,   & ! <-- Year of occultation
         Month,  & ! <-- Month
         Day,    & ! <-- Day
         Hour,   & ! <-- Hour
         Minute, & ! <-- Minute
         Second, & ! <-- Second of beginning of occultation
         WT(i))    ! <-- Second of occultation duration
   End If
   WRLEO(i) = Rotate(WERLEO(i), PA, Phi)
   WRGPS(i) = Rotate(WERGPS(i), PA, Phi)
End Do

Call Satellite_Velocities &
  (WT,       & ! <-- Relative time of samples [sec]
   WRLEO,    & ! <-- LEO coordinates (J2000)
   WRGPS,    & ! <-- GPS coordinates (J2000)
   WRLEO,    & ! --> LEO coordinates from regression (J2000)
   WVLEO,    & ! --> LEO velocities from regression (J2000)
   WRGPS,    & ! --> GPS coordinates from regression (J2000)
   WVGPS)      ! --> GPS velocities from regression (J2000)

!----------------------------------------------------------
! 3. ASYMPTOTIC WAVE SOLUTION
!----------------------------------------------------------

If (LVrb >= 1) then
   Write(*,'(2X,A)') 'Asymptotic wave solution'
End If


!--- 3.0. Initialization of spline interpolation

Allocate(YP2(NP))
Allocate(ARP2(NP))
Allocate(ATP2(NP,NC))
Allocate(YT2(NP))
Allocate(S02(NP))
Allocate(FT2(NP))

Call Init_Spline &
  (P(:),        & ! <-- Argument grid
   YP(:),       & ! <-- Gridded function
   YP2(:))        ! --> 2nd derivative of spline

Call Init_Spline &
  (P(:),        & ! <-- Argument grid
   ARP(:),      & ! <-- Gridded function
   ARP2(:))       ! --> 2nd derivative of spline

Do IC=1,NC
   Call Init_Spline &
     (P(:),        & ! <-- Argument grid
      ATP(:,IC),   & ! <-- Gridded function
      ATP2(:,IC))    ! --> 2nd derivative of spline
End Do

Call Init_Spline &
  (YT(:),       & ! <-- Argument grid
   FT(:),       & ! <-- Gridded function
   FT2(:))        ! --> 2nd derivative of spline

Call Init_Spline &
  (YT(:),       & ! <-- Argument grid
   S0(:),       & ! <-- Gridded function
   S02(:))        ! --> 2nd derivative of spline

Call Init_Spline &
  (TG(:),       & ! <-- Argument grid
   YT(:),       & ! <-- Gridded function
   YT2(:))        ! --> 2nd derivative of spline


!DEBUG
!Call PutXY('test-atp1.dat', P(:) - RE, Exp(-k(1)*ATP(:,1)), Stat=Stat)
!Call PutXY('test-atp2.dat', P(:) - RE, Exp(-k(2)*ATP(:,2)), Stat=Stat)
!Call PutXY('test-atp3.dat', P(:) - RE, Exp(-k(3)*ATP(:,3)), Stat=Stat)

!--- 3.1. Lagrange manifold dimension

Pmin  = MinVal(P)
Pmax  = MaxVal(P)
Ymax  = MaxVal(YP)

PB    = (Pmax - Pmin)/4

PZmin = Pmin - PB
PZmax = Pmax + PB


!--- 3.2. Channel processing

Channels: Do IC = 1,NC


   !--- 3.2.1. Determination of resolution

   NH = k(IC)*Ymax*(PZmax - PZmin)/Pi
   NH = Nearest_Power2(NH)

   If (LVrb >= 1) then
      Write(*,'(2X,A,I2,A,F8.4,A,I8)')       &
         'Channel ', IC,                     &
         '  Freq [GHz] = ', 1d-9*Freq(IC),   &
         '  NH = ', NH
   End If


   !--- 3.2.2. Hi-res grid

   Allocate(ZH(NH))
   Allocate(PH(NH))
   Allocate(ARH(NH))
   Allocate(ATH(NH))
   Allocate(UH(NH))
   Allocate(YH(NH))
   Allocate(AH(NH))

   Do i=1,NH
      ZH(i) = (PZmax - PZmin)*(i-1)/Real(NH-1)
   End Do


   !--- 3.2.3. Computation of solution in p-representation

   PH(1)  = 0.0
   ARH(1) = 0.0
   ATH(1) = 0.0
   UH(1)  = 0.0

   Do i=2,NH
      If (ZH(i) > Pmin - PZmin .and. &
          ZH(i) < Pmax - PZmin) then
         Call Spline &
           (P(:),            & ! <-- Argument grid
            YP(:),           & ! <-- Gridded function
            YP2(:),          & ! <-- 2nd derivative of spline
            PZmin + ZH(i),   & ! <-- Interpolation point
            YPI)               ! --> Interpolated function value
         PH(i) = PH(i-1) - YPI*(ZH(i) - ZH(i-1))
         Call Spline &
           (P(:),            & ! <-- Argument grid
            ARP(:),          & ! <-- Gridded function
            ARP2(:),         & ! <-- 2nd derivative of spline
            PZmin + ZH(i),   & ! <-- Interpolation point
            ARH(i))            ! --> Interpolated function value
         Call Spline &
           (P(:),            & ! <-- Argument grid
            ATP(:,IC),       & ! <-- Gridded function
            ATP2(:,IC),      & ! <-- 2nd derivative of spline
            PZmin + ZH(i),   & ! <-- Interpolation point
            ATH(i))            ! --> Interpolated function value
      Else
         PH(i)  = 0.0
         ARH(i) = 0.0
         ATH(i) = 0.0
      End If
      UH(i) = ARH(i)*Exp(Ci*(Modulo(k(IC)*PH(i), 2*Pi)))*Exp(-k(IC)*ATH(i))
   End Do

   !!!DEBUG
   !i = NH/4096
   !Call PutXY('test-azh.dat', ZH(1::i), ARH(1::i), Stat=Stat)
   !Call PutXY('test-pzh.dat', ZH(1::i), PH(1::i), Stat=Stat)

   !--- 3.2.4. Upper-border apodization

   Where (PZmin - RE + ZH(:) > Hmax)
      UH(:) = UH(:)*Exp(-(PZmin - RE + ZH(:) - Hmax)**2/(0.3*DHB)**2)
   End Where


   !--- 3.2.5. FIO

   FSign = 1

   Call FFT1 &
     (FSign,    & ! <-- Transform direction.
      UH(:))      ! <-> Data to be transformed.

   UH(:) = UH(:)/Sqrt(Real(NH))


   !--- 3.2.6. Determination of hi-res grid of Y

   Do i=1,NH
      YH(i) = 2*Pi*(i-1)/(k(IC)*(PZmax - PZmin))
   End Do


   !--- 3.2.7. Determination of amplitude

   AH(:) = Abs(UH(:))


   !--- 3.2.8. Computation of phase excess

   Where (UH(:) /= 0.0_Double)
      PH(:) = ATan2(AImag(UH(:)),Real(UH(:)))
   Elsewhere
      PH(:) = 0.0
   End Where

   Call Accumulate_Phase &
     (PH(:),   & ! <-> Array of (accumulated) phase
      1)         ! <~~ Phase change direction

   PH(:) = PH(:)/k(IC)

   !!!DEBUG
   !i = NH/4096
   !Call PutXY('test-py0h.dat', YH(1::i), PH(1::i), Stat=Stat)

   Do i=1,NH
      Call Spline &
        (YT(:),    & ! <-- Argument grid
         FT(:),    & ! <-- Gridded function
         FT2(:),   & ! <-- 2nd derivative of spline
         YH(i),    & ! <-- Interpolation point
         FTI)        ! --> Interpolated function value
      Call Spline &
        (YT(:),    & ! <-- Argument grid
         S0(:),    & ! <-- Gridded function
         S02(:),   & ! <-- 2nd derivative of spline
         YH(i),    & ! <-- Interpolation point
         S0I)        ! --> Interpolated function value
      PH(i) = PH(i) - FTI + PZmin*YH(i) - (S0I - S0(1))
   End Do

   PH(:) = 1d3*PH(:)

   !!!DEBUG
   !i = NH/4096
   !Call PutXY('test-ayh.dat', YH(1::i), AH(1::i), Stat=Stat)
   !Call PutXY('test-pyh.dat', YH(1::i), PH(1::i), Stat=Stat)
   !Call PutXY('test-FTY.dat', YT(:), 1d3*FT(:), Stat=Stat)
   !Stop 'DEBUG to HERE...'


   !--- 3.2.9. Interpolation of amplitude and phase excess

   Do i=1,NWT
      Call Spline &
        (TG(:),    & ! <-- Argument grid
         YT(:),    & ! <-- Gridded function
         YT2(:),   & ! <-- 2nd derivative of spline
         WT(i),    & ! <-- Interpolation point
         YTI)        ! --> Interpolated function value
      Call Linear &
        (YH(:),    & ! <-- Argument grid
         AH(:),    & ! <-- Gridded function
         YTI,      & ! <-- Interpolation point
         WA(i,IC))   ! --> Interpolated function value
      Call Linear &
        (YH(:),    & ! <-- Argument grid
         PH(:),    & ! <-- Gridded function
         YTI,      & ! <-- Interpolation point
         WS(i,IC))   ! --> Interpolated function value
   End Do

   PHB      = MinVal(WS(:,IC))
   WS(:,IC) = WS(:,IC) - PHB


   !--- 3.2.10. Memory deallocation

   Deallocate(ZH)
   Deallocate(PH)
   Deallocate(ARH)
   Deallocate(ATH)
   Deallocate(UH)
   Deallocate(YH)
   Deallocate(AH)


End Do Channels

!!!DEBUG
!Call PutXY('test-wa1.dat', WT(:), WA(:,1), Stat = Stat)
!Call PutXY('test-wa2.dat', WT(:), WA(:,2), Stat = Stat)
!Call PutXY('test-wa3.dat', WT(:), WA(:,3), Stat = Stat)
!Call PutXY('test-ws1.dat', WT(:), WS(:,1), Stat = Stat)
!Call PutXY('test-ws2.dat', WT(:), WS(:,2), Stat = Stat)
!Call PutXY('test-ws3.dat', WT(:), WS(:,3), Stat = Stat)
!Call PutXY('test-ep.dat',  EP(:), P(:)-RE, Stat = Stat)


!----------------------------------------------------------
! 4. MEMORY DEALLOCATION
!----------------------------------------------------------


Deallocate(k)
Deallocate(Ak)

Deallocate(WERLEO)
Deallocate(WEVLEO)
Deallocate(WERGPS)
Deallocate(WEVGPS)

Deallocate(YP2)
Deallocate(ARP2)
Deallocate(ATP2)
Deallocate(YT2)
Deallocate(S02)
Deallocate(FT2)


Stat = 0


End Subroutine Asymptotic_Propagate



End Module Asymptotic_Propagator


