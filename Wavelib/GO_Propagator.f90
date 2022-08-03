!
Module GO_Propagator
!
! Geometric optical propagator.
!----------------------------------------------------------
! (C) Copyright 1999-2009, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 27 Feb 1999 | Original version.
!   2.0   | 21 May 1999 | GO_Caustics.
!   3.0   | 17 Sep 1999 | GO_Caustics in Occ_Refraction.
!   4.0   | 06 Oct 1999 | Computation of ray heights.
!   5.0   | 23 May 2003 | GO_Invert.
!   5.1   | 04 Mar 2009 | Local_Profiles.
!   5.2   | 14 Mar 2009 | Local_Profiles_Vertical and
!         |             | Local_Profiles_Skew.
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
Use Atmosphere, only: &
! Imported Routines:
    Atmosphere_NGradN,        &
    Atmosphere_Constituents
!
Use Atmosphere_Rays, only: &
! Imported Routines:
    Ray_Trace
!
Use IO, only: &
! Imported Routines:
    PutXY
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Interfaces:
!
Interface Local_Profiles
   Module Procedure Local_Profiles_Vertical
   Module Procedure Local_Profiles_Skew
End Interface
!----------------------------------------------------------
!
Contains


!==========================================================
Subroutine GO_Propagate &
  (Year,      & ! <-- Occultation year
   Month,     & ! <-- Occultation month
   Day,       & ! <-- Occultation day
   Hour,      & ! <-- Occultation hour
   Minute,    & ! <-- Occultation begin minute
   Second,    & ! <-- Occultation begin second
   TR,        & ! <-- Relative time of samples [sec]
   RLEO,      & ! <-- LEO coordinates (J2000/ECEF)
   RGPS,      & ! <-- GPS coordinates (J2000/ECEF)
   Hmax,      & ! <-- Maximum height
   S,         & ! <-- Ray integration step parameter
   AFreq,     & ! <-- Absorption frequencies [Hz]
   opt_ECEF,  & ! <-- Coordinates in ECEF
   EGO,       & ! --> GO refraction angle
   PGO,       & ! --> GO impact parameter
   ZGO,       & ! --> GO ray heights
   LAT,       & ! --> Logarithmic attenuation (p,channel)
   EGR,       & ! --> GO reflected refraction angle
   PGR,       & ! --> GO reflected impact parameter
   ZGR,       & ! --> GO reflected ray heights
   GP,        & ! --> Occultation point
   ERLC,      & ! --> Local curvature center (ECEF)
   RE,        & ! --> Local curvature radius
   Vrb)         ! <~~ Verbosity level
!
! Calculation of ray propagation in atmosphere.
!----------------------------------------------------------
! Method:
!   Geometric optical ray equation.
!----------------------------------------------------------
! (C) Copyright 1999-2009, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 27 Feb 1999 | Original version.
!   2.0   | 14 Mar 1999 | Inversion.
!   2.1.  | 15 Mar 1999 | Elimination of ambiguities.
!   3.0.  | 23 Mar 1999 | Interpolated temperature.
!   4.0   | 08 Jul 1999 | Verbosity parameter.
!   4.1   | 10 Dec 1999 | Inhomogeneous grid of impact parameter.
!   5.0   | 04 Dec 2000 | QE, RNE, and TDE.
!   5.1   | 18 May 2001 | RIE.
!   6.0   | 10 Aug 2002 | Reflected rays.
!   7.0   | 23 May 2003 | Inversion in a separate subroutine.
!   8.0   | 20 Aug 2003 | Verbosity level.
!   8.1   | 24 Apr 2004 | Modified for new version of Ray_Trace.
!   9.0   | 21 May 2004 | AFreq, LAT.
!   9.1   | 14 Sep 2009 | opt_ECEF.
!   9.2   | 23 Sep 2009 | Index order changed to (p,channel).
!----------------------------------------------------------
! Modules used:
!
Use Coordinates, only: &
! Imported Type Definitions:
    Cartesian,      &
! Imported Routines:
    Vector_Normed,  &
! Imported Operators:
    Operator(+),    &
    Operator(-),    &
    Operator(*),    &
    Operator(.xx.),  &
    Assignment(=)
!
Use Earth, only: &
! Imported Type Definitions:
    Geodetic
!
Use Occ_Coordinates, only: &
! Imported Routines:
    Occ_Geometry_ECEF,     &
    Occ_Geometry,          &
    Plane_Coordinates,     &
    Plane_Basis
!
Use Occ_Refraction, only: &
! Imported Routines:
    Doppler_to_Refraction
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Integer, Intent(In)          :: &
   Year        ! Occultation year
!
Integer, Intent(In)          :: &
   Month       ! Occultation month
!
Integer, Intent(In)          :: &
   Day         ! Occultation day
!
Integer, Intent(In)          :: &
   Hour        ! Occultation hour
!
Integer, Intent(In)          :: &
   Minute      ! Occultation begin minute
!
Real(Double), Intent(In)     :: &
   Second      ! Occultation begin second
!
Real(Double), Intent(In)     :: &
   TR(:)       ! Relative time of samples [sec]
!
Type(Cartesian), Intent(In) :: &
   RLEO(:)     ! LEO coordinates (J2000)
!
Type(Cartesian), Intent(In) :: &
   RGPS(:)     ! GPS coordinates (J2000)
!
Real(Double), Intent(In) :: &
   Hmax        ! Maximum height
!
Real(Double), Intent(In) :: &
   S           ! Ray integration step parameter
!
Real(Double), Intent(In) :: &
   AFreq(:)    ! Absorption frequencies [Hz]
!
Logical, Intent(In)      :: &
   opt_ECEF    ! Coordinates in ECEF
!
! Output arguments:
!
Real(Double), Intent(Out) :: &
   EGO(:)      ! GO refraction angle
!
Real(Double), Intent(Out) :: &
   PGO(:)      ! GO impact parameter
!
Real(Double), Intent(Out) :: &
   ZGO(:)      ! GO ray heights
!
Real(Double), Intent(Out) :: &
   LAT(:,:)    ! Logarithmic attenuation (p,channel)
!
Real(Double), Intent(Out) :: &
   EGR(:)      ! GO reflected refraction angle
!
Real(Double), Intent(Out) :: &
   PGR(:)      ! GO reflected impact parameter
!
Real(Double), Intent(Out) :: &
   ZGR(:)      ! GO reflected ray heights
!
Type(Geodetic), Intent(Out) :: &
   GP          ! Occultation point
!
Type(Cartesian), Intent(Out)  :: &
   ERLC        ! Local curvature center (ECEF)
!
Real(Double), Intent(Out) :: &
   RE          ! Local curvature radius
!
! Input optional arguments:
!
Integer, Optional, Intent(In) :: &
   Vrb         ! Verbosity level
               ! 0 by default.
!----------------------------------------------------------
! Local Parameters:
!
Real(Double), Parameter :: &
   A = 0.9                  ! Inhomogeneity parameter of grid of p
!
! Local Scalars:
!
! --- Input data
!
Integer           :: N        ! Number of data read
!
! --- Geometric parameters
!
Type(Cartesian)   :: RLC      ! Local curvature center (J2000)
Type(Cartesian)   :: AX, AY   ! Occultation plane XY-basis
Type(Cartesian)   :: EX, EY   ! Occultation plane basis (ECEF)
!
! --- Variables for ray-tracing
!
Integer           :: NGO      ! Number of direct refraction angles
Integer           :: NGR      ! Number of reflected refraction angles
Type(Cartesian)   :: XT       ! Transmitter position
Type(Cartesian)   :: XN       ! Final ray point
Type(Cartesian)   :: UT       ! Ray direction at transmitter
Type(Cartesian)   :: UN       ! Ray direction at XN
Real(Double)      :: x        ! Parameter for scaling grid of p
Real(Double)      :: YP       ! Ray leveling height
Real(Double)      :: YPmin    ! Minimum height of direct ray
Real(Double)      :: YRmax    ! Maximum height of reflected ray
Real(Double)      :: Y1, Y2   ! Perigee limits for dichotomy
Real(Double)      :: Eps      ! Refraction angle
Type(Cartesian)   :: VT       ! Transmitter velocity
Type(Cartesian)   :: VN       ! Receiver velocity
Real(Double)      :: d        ! Doppler frequency shift
!
! --- Variables for absorption
!
Integer           :: NC       ! Number of channels
!
! --- Work variables
!
Integer           :: LVrb     ! Verbosity level
Integer           :: i        ! Array index
Integer           :: Stat     ! Error status
Character(Len=80) :: Line     ! Terminal line
!
!
! Local Arrays:
!
! --- Geometric parameters
!
Type(Cartesian), Allocatable :: &
   ERLEO(:),  & ! LEO coordinates (ECEF)
   ERGPS(:)     ! GPS coordinates (ECEF)
Real(Double), Allocatable :: &
   XLEO(:),   & ! X coordinates of LEO
   YLEO(:),   & ! Y coordinates of LEO
   XGPS(:),   & ! X coordinates of GPS
   YGPS(:)      ! Y coordinates of GPS
!
! --- Absorption
!
Real(Double), Allocatable  :: &
   ATP(:)       ! Absorptive attenuation
!----------------------------------------------------------


!----------------------------------------------------------
! 0. INITIALIZATION
!----------------------------------------------------------


!--- 0.1. Determination of number of data

N  = Size(TR)
NC = Size(AFreq)


!--- 0.2. Array allocation

Allocate(ERLEO(1:N))
Allocate(ERGPS(1:N))
Allocate(XLEO(1:N))
Allocate(YLEO(1:N))
Allocate(XGPS(1:N))
Allocate(YGPS(1:N))
Allocate(ATP(NC))


!--- 0.3. Verbosity status definition

If (Present(Vrb)) then
   LVrb = Vrb
Else
   LVrb = 0
End If


!----------------------------------------------------------
! 1. COORDINATE TRANSFORMS
!----------------------------------------------------------


!--- 1.1. Calculation of curvature center and
!---      coordinates (ECEF)


If (opt_ECEF) then

   ERLEO(:) = RLEO(:)
   ERGPS(:) = RGPS(:)

   Call Occ_Geometry_ECEF &
     (ERLEO,      & ! <-- LEO coordinates (ECEF)
      ERGPS,      & ! <-- GPS coordinates (ECEF)
      GP,         & ! --> Geodetic coordinates of occultation point
      ERLC,       & ! --> Curvature center (ECEF)
      RE)           ! --> Local curvature radius

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
      RE)          ! --> Local curvature radius

End If


!--- 1.2. Calculation of coordinates in occultation plane

Call Plane_Coordinates &
  (RLEO,      & ! <-- LEO coordinates (J2000)
   RGPS,      & ! <-- GPS coordinates (J2000)
   RLC,       & ! <-- Curvature center (J2000)
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
  (ERLEO,     & ! <-- LEO coordinates (J2000)
   ERGPS,     & ! <-- GPS coordinates (J2000)
   ERLC,      & ! <-- Curvature center (J2000)
   RE,        & ! <-- Local curvature radius [km]
   EX,        & ! --> Occultation plane X basis vector
   EY)          ! --> Occultation plane Y basis vector


!----------------------------------------------------------
! 2. DETERMINATION OF CRITICAL POINTS
!----------------------------------------------------------


!--- 2.1. Determination of shadow border

If (LVrb >= 1) then
   Write(*,'(2X,A)') 'Determination of shadow border.'
End If

Y1  = -20.0_Double
Y2  =  20.0_Double
Eps = 0

XT = ERLC + EX*XGPS(N) + EY*(RE+YGPS(N))

Dichotomy_Shadow: Do
   YP = (Y1 + Y2)/2
   UT = Vector_Normed(ERLC + EY*(RE+YP) - XT)
   Call Ray_Trace &
     (XT,        & ! <-- Transmitter position
      UT,        & ! <-- Transmitter ray direction
      DS = S,    & ! <-- Integration step parameter
      XN = XN,   & ! --> Ray point nearest to receiver
      UN = UN,   & ! --> Ray direction at XN
      Stat=Stat)   ! --> Error status
   If (Stat == 0) then
      Y2 = YP
   Else
      Y1 = YP
   End If
   If (LVrb >= 3) then
      Write (Line,'(2X,A,F9.4,A,F9.4,A,I1,A1)')  &
         'Y1 = ', Y1,  ', Y2 = ', Y2,  ', Stat = ', Stat, &
         Char(0)
      Call CPrintf(Line)
   End If
   If (Abs(Y2-Y1) < 0.001) then
      YPmin = Y2
      Exit Dichotomy_Shadow
   End If
End Do Dichotomy_Shadow


If (LVrb >= 3) then
   Write(*,'()')
End If


!--- 2.2. Determination of non-trapped reflection

If (LVrb >= 1) then
   Write(*,'(2X,A)') 'Determination of non-trapped reflection.'
End If

Y1  = -20.0_Double
Y2  =  YPmin
Eps = 0

XT = ERLC + EX*XGPS(N) + EY*(RE+YGPS(N))

Dichotomy_Reflection: Do
   YP = (Y1 + Y2)/2
   UT = Vector_Normed(ERLC + EY*(RE+YP) - XT)
   Call Ray_Trace &
     (XT,        & ! <-- Transmitter position
      UT,        & ! <-- Transmitter ray direction
      DS = S,    & ! <-- Integration step parameter
      XN = XN,   & ! --> Ray point nearest to receiver
      UN = UN,   & ! --> Ray direction at XN
      Stat=Stat)   ! --> Error status
   If (Stat == 1) then
      Y1 = YP
   Else
      Y2 = YP
   End If
   If (LVrb >= 3) then
      Write (Line,'(2X,A,F9.4,A,F9.4,A,I1,A1)')  &
         'Y1 = ', Y1,  ', Y2 = ', Y2,  ', Stat = ', Stat, &
         Char(0)
      Call CPrintf(Line)
   End If
   If (Abs(Y2-Y1) < 0.001) then
      YRmax = Y1
      Exit Dichotomy_Reflection
   End If
End Do Dichotomy_Reflection


If (LVrb >= 3) then
   Write(*,'()')
End If


!----------------------------------------------------------
! 3. GEOMETRIC OPTICAL PROPAGATION
!----------------------------------------------------------


!--- 3.1. Direct rays

If (LVrb >= 1) then
   Write(*,'(2X,A)') 'Direct rays.'
End If

NGO = Size(EGO)
VT  = 0
VN  = 8*EY

Do i=1,NGO
!!!   YP = Hmax + (i-1)*(YPmin - Hmax)/Real(NGO-1)
   x  = Real(NGO - i)/Real(NGO - 1)
   YP = YPmin + ((1-A)*x + A*x*x)*Hmax
   UT = Vector_Normed(ERLC + EY*(RE+YP) - XT)
   Call Ray_Trace &
     (XT,           & ! <-- Transmitter position
      UT,           & ! <-- Transmitter ray direction
      DS = S,       & ! <-- Integration step parameter
      XN = XN,      & ! --> Ray point nearest to receiver
      UN = UN,      & ! --> Ray direction at XN
      At = ATP(:),  & ! ~~> Attenuation [Int(Imag(N))]
      Stat=Stat)      ! --> Error status
   LAT(i,:) = ATP(:)*2*Pi*AFreq(:)/C_Light
   d = ((VT*UT) - (VN*UN))/(C_Light - (VT*UT))
   Call Doppler_to_Refraction &
     (XT-ERLC,  & ! <-- Transmitter position
      VT,       & ! <-- Transmitter velocity
      XN-ERLC,  & ! <-- Receiver position
      VN,       & ! <-- Receiver velocity
      d,        & ! <-- Relative Doppler frequency shift
      PGO(i),   & ! --> Impact parameter
      EGO(i))     ! --> Refraction angle
   If (LVrb >= 3) then
      Write (Line,'(2X,A,I5,A,F9.3,A,ES9.2,A1)')  &
         'i = ', i,  ', P = ', PGO(i) - RE,  ', Eps = ', EGO(i), &
         Char(0)
      Call CPrintf(Line)
   End If
End Do

If (LVrb >= 3) then
   Write(*,'()')
End If


ZGO(:) = PGO(:) - RE


!--- 3.2. Reflected rays

If (LVrb >= 1) then
   Write(*,'(2X,A)') 'Reflected rays.'
End If

NGR = Size(EGR)
VT  = 0
VN  = 8*EY

Do i=1,NGR
   x  = Real(NGR - i)/Real(NGR - 1)
   YP = YRmax - 0.2 + x*0.2
   UT = Vector_Normed(ERLC + EY*(RE+YP) - XT)
   Call Ray_Trace &
     (XT,        & ! <-- Transmitter position
      UT,        & ! <-- Transmitter ray direction
      DS = S,    & ! <-- Integration step parameter
      XN = XN,   & ! --> Ray point nearest to receiver
      UN = UN,   & ! --> Ray direction at XN
      Stat=Stat)   ! --> Error status
   d = ((VT*UT) - (VN*UN))/(C_Light - (VT*UT))
   Call Doppler_to_Refraction &
     (XT-ERLC,  & ! <-- Transmitter position
      VT,       & ! <-- Transmitter velocity
      XN-ERLC,  & ! <-- Receiver position
      VN,       & ! <-- Receiver velocity
      d,        & ! <-- Relative Doppler frequency shift
      PGR(i),   & ! --> Impact parameter
      EGR(i))     ! --> Refraction angle
   If (LVrb >= 3) then
      Write (Line,'(2X,A,I5,A,F9.3,A,ES9.2,A1)')  &
         'i = ', i,  ', P = ', PGR(i)-RE,  ', Eps = ', EGR(i), &
         Char(0)
      Call CPrintf(Line)
   End If
End Do

If (LVrb >= 3) then
   Write(*,'()')
End If


ZGR(:) = PGR(:) - RE


!----------------------------------------------------------
! 4. MEMORY DEALLOCATION
!----------------------------------------------------------


Deallocate(ERLEO)
Deallocate(ERGPS)
Deallocate(XLEO)
Deallocate(YLEO)
Deallocate(XGPS)
Deallocate(YGPS)
Deallocate(ATP)


End Subroutine GO_Propagate




!==========================================================
Subroutine GO_Invert &
  (EGO,       & ! <-- GO refraction angle
   PGO,       & ! <-- GO impact parameter
   GP,        & ! <-- Occultation point
   RE,        & ! <-- Local curvature radius
   AFreq,     & ! <-- Absorption frequency [Hz]
   Z,         & ! --> Altitudes [km]
   RN,        & ! --> Refractivities [absolute]
   TD,        & ! --> Dry temperatures [K]
   TE,        & ! --> Interpolated ECHAM temperatures [K]
   QE,        & ! --> Interpolated ECHAM humidity [kg/kg]
   RNE,       & ! --> Interpolated ECHAM refractivity [n-1]
   RAE,       & ! --> Interpolated ECHAM specific absorption (z,channel) [km^-1]
   TDE,       & ! --> Dry temperature from ECHAM refractivity [K]
   Vrb)         ! <~~ Verbosity level
!
! Inversion of refraction angle profiles
!----------------------------------------------------------
! Method:
!   Dry atmosphere inversion.
!----------------------------------------------------------
! (C) Copyright 2003-2009, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 23 May 2003 | Original version.
!   2.0   | 23 May 2003 | Verbosity level.
!   2.1   | 02 Sep 2009 | Improved notation.
!   2.2   | 23 Sep 2009 | Index order changed to (z,channel)
!----------------------------------------------------------
! Modules used:
!
Use Earth, only: &
! Imported Type Definitions:
    Geodetic
!
Use Occ_Inversion, only: &
! Imported Routines:
    Invert_Refraction
!
Use Occ_Meteoprofiles, only: &
! Imported Array Variables:
    NQ_to_TP
!
Use Interpolation, only: &
! Imported Array Variables:
    Linear
!
Use Signal, only: &
! Imported Routines:
    Monotonize
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In) :: &
   EGO(:)     ! GO refraction angle
!
Real(Double), Intent(In) :: &
   PGO(:)     ! GO impact parameter
!
Type(Geodetic), Intent(In) :: &
   GP         ! Occultation point
!
Real(Double), Intent(In) :: &
   RE         ! Local curvature radius
!
Real(Double), Intent(In) :: &
   AFreq(:)   ! Absorption frequency [Hz]
!
! Output arguments:
!
Real(Double), Intent(Out) :: &
   Z(:)       ! Altitudes [km]
!
Real(Double), Intent(Out) :: &
   RN(:)      ! Refractivities [absolute]
!
Real(Double), Intent(Out) :: &
   TD(:)      ! Dry temperatures [K]
!
Real(Double), Intent(Out) :: &
   TE(:)      ! Interpolated ECHAM temperatures [K]
!
Real(Double), Intent(Out) :: &
   QE(:)      ! Interpolated ECHAM humidity [kg/kg]
!
Real(Double), Intent(Out) :: &
   RNE(:)     ! Interpolated ECHAM refractivity [n-1]
!
Real(Double), Intent(Out) :: &
   RAE(:,:)   ! Interpolated ECHAM specific absorption (z,channel) [km^-1]
!
Real(Double), Intent(Out) :: &
   TDE(:)     ! Dry temperature from ECHAM refractivity [K]
!
! Input optional arguments:
!
Integer, Optional, Intent(In) :: &
   Vrb        ! Verbosity mode.
              ! True by default.
!----------------------------------------------------------
! Local Parameters:
!
Real(Double), Parameter :: &
   DZH   = 0.01_Double,     & ! High-resoltuion grid step
   ZHmax = 130.0_Double       ! Upper limit of high-resolution grid
!
! Local Scalars:
!
Integer           :: LVrb     ! Verbosity level
Integer           :: NGO      ! Number of direct refraction angles
Integer           :: NH       ! High-resolution grid dimension
Integer           :: i        ! Profile index
Real(Double)      :: Zmin     ! Minimum model altitude
Real(Double)      :: Zmax     ! Maximum model altitude
Integer           :: NC       ! Number of frequency channels
Real(Double)      :: ZmaxI    ! Upper integration height [km]
Real(Double)      :: DZI      ! Integration step [km]
Real(Double)      :: Tinit    ! Temperature at upper boundary [K]
Real(Double)      :: Qmin     ! Minimum specific humidity [kg/kg]
!
! Local Arrays:
!
! --- Refraction
!
Real(Double), Allocatable :: &
   PU(:)        ! Monotonized impact parameters
!
! --- Standard resolution profiles
!
Real(Double), Allocatable :: &
   RIE(:,:)     ! Profile of Im(N) (z,channel)
!
! --- High resolution profiles
!
Real(Double), Allocatable :: &
   ZH(:),     & ! High-resolution altitude grid
   RNH(:),    & ! High-resolution refractivity profile
   QH(:),     & ! High-resolution humidity profile
   TDH(:)       ! High-resolution temperature profile
!
! --- Wave vector for absorption
!
Real(Double), Allocatable :: &
   Ak(:)        ! Wave vectors
!----------------------------------------------------------


!----------------------------------------------------------
! 1. INVERSION
!----------------------------------------------------------


!--- 1.0. Setting verbosity mode

If (Present(Vrb)) then
   LVrb = Vrb
Else
   LVrb = 0
End If


!--- 1.1. Determination of number of data

NGO = Size(EGO)
NC  = Size(AFreq)


!--- 1.2. Elimination of ambiguities in eps(p)

If (LVrb >= 1) then
   Write(*,'(2X,A)') 'Elimination of ambiguities'
End If

Allocate(PU(NGO))
PU(:) = PGO(:)

!!!DEBUG
!Call PutXY('test-pgo.dat', EGO(:), PGO(:), Stat=i)

Call Monotonize(-1, 0.001_Double, PU(:), LVrb)


!--- 1.3. Inversion

If (LVrb >= 1) then
   Write(*, '(2X,A)') 'Inversion'
End If

Call Invert_Refraction &
  (PU,            & ! <-- Impact parameters [km]
   EGO,           & ! <-- Refraction angles [rad]
   GP,            & ! <-- Occultation point (geodetic)
   RE,            & ! <-- Local curvature radius [km]
   LVrb,          & ! <-- Verbosity level
   Z,             & ! --> Altitudes [km]
   RN,            & ! --> Refractivities [absolute]
   TD)              ! --> Dry temperatures [K]


!--- 1.3. Interpolation

Allocate(Ak(NC))
Allocate(RIE(NGO,NC))

Ak(:) = 2*Pi*AFreq(:)/C_Light

Do i=1,NGO
   Call Atmosphere_Constituents &
     (GP%Lambda,    & ! <-- Longiude of point [deg]
      GP%Phi,       & ! <-- Latitude of point [deg]
      Z(i),         & ! <-- Altitude of point [km]
      Zmin,         & ! --> Minimum model Z for this lon/lat
      Zmax,         & ! --> Maximum model Z for this lon/lat
      QE(i),        & ! --> Interpolated constituent
      TE(i),        & ! --> Interpolated constituent
      RNE(i),       & ! --> Interpolated Re(N)
      RIE(i,:))       ! ~~> Interpolated Im(N)
   RAE(i,:) = RIE(i,:)*Ak(:)
End Do


!--- 1.4. Computation of dry temperature

NH = 1 + Ceiling((ZHmax - Zmin)/DZH)

Allocate(ZH(NH))
Allocate(RNH(NH))
Allocate(QH(NH))
Allocate(TDH(NH))

Do i=1,NH
   ZH(i) = Zmin + Real(i-1)*(ZHmax - Zmin)/Real(NH-1)
   Call Atmosphere_Constituents &
     (GP%Lambda,    & ! <-- Longiude of point [deg]
      GP%Phi,       & ! <-- Latitude of point [deg]
      ZH(i),        & ! <-- Altitude of point [km]
      Zmin,         & ! --> Minimum model Z for this lon/lat
      Zmax,         & ! --> Maximum model Z for this lon/lat
      QH(i),        & ! --> Interpolated humidity [kg/kg]
      TDH(i),       & ! --> Interpolated temperature [K]
      RNH(i))         ! --> Interpolated Re(N)
   QH(i) = 0
End Do

ZmaxI = 120.0
DZI   = 0.015
Tinit = 0.0
Qmin  = 1d-7

Call NQ_to_TP &
  (GP%Phi,  & ! <-- Geodetic latitude [deg]
   ZH,      & ! <-- Altitude above reference ellipsoid [km]
   RNH,     & ! <-- Profile of refractivity
   QH,      & ! <-- Profile of specific humidity
   ZmaxI,   & ! <-- Upper integration height [km]
   DZI,     & ! <-- Integration step [km]
   Tinit,   & ! <-- Temperature at upper boundary [K]
   Qmin,    & ! <-- Minimum specific humidity [kg/kg]
   TDH)       ! --> Profile of temperature

Do i=1,NGO
   Call Linear &
     (ZH,      & ! <-- Argument grid
      TDH,     & ! <-- Gridded function
      Z(i),    & ! <-- Interpolation point
      TDE(i))    ! --> Interpolated function value
End Do


!----------------------------------------------------------
! 2. MEMORY DEALLOCATION
!----------------------------------------------------------

Deallocate(Ak)
Deallocate(RIE)
Deallocate(PU)
Deallocate(ZH)
Deallocate(RNH)
Deallocate(QH)
Deallocate(TDH)


End Subroutine GO_Invert



!==========================================================
Subroutine Local_Profiles_Vertical &
  (Z,         & ! <-- Altitudes [km]
   GP,        & ! <-- Occultation point
   AFreq,     & ! <-- Absorption frequency [Hz]
   TE,        & ! --> Interpolated ECHAM temperatures [K]
   QE,        & ! --> Interpolated ECHAM humidity [kg/kg]
   RNE,       & ! --> Interpolated ECHAM refractivity [n-1]
   RAE,       & ! --> Interpolated ECHAM specific absorption (z,channel) [km^-1]
   TDE,       & ! --> Dry temperature from ECHAM refractivity [K]
   Vrb)         ! <~~ Verbosity level
!
! Getting local atmospheric profiles.
!----------------------------------------------------------
! (C) Copyright 2009, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 04 Mar 2009 | Original version.
!   1.1   | 02 Sep 2009 | Improved notation.
!   1.2   | 23 Sep 2009 | Index order changed to (z,channel).
!----------------------------------------------------------
! Modules used:
!
Use Earth, only: &
! Imported Type Definitions:
    Geodetic
!
Use Occ_Inversion, only: &
! Imported Routines:
    Invert_Refraction
!
Use Occ_Meteoprofiles, only: &
! Imported Array Variables:
    NQ_to_TP
!
Use Interpolation, only: &
! Imported Array Variables:
    Linear
!
Use Signal, only: &
! Imported Routines:
    Monotonize
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In) :: &
   Z(:)       ! Altitudes [km]
!
Type(Geodetic), Intent(In) :: &
   GP         ! Occultation point
!
Real(Double), Intent(In) :: &
   AFreq(:)   ! Absorption frequency [Hz]
!
! Output arguments:
!
Real(Double), Intent(Out) :: &
   TE(:)      ! Interpolated ECHAM temperatures [K]
!
Real(Double), Intent(Out) :: &
   QE(:)      ! Interpolated ECHAM humidity [kg/kg]
!
Real(Double), Intent(Out) :: &
   RNE(:)     ! Interpolated ECHAM refractivity [n-1]
!
Real(Double), Intent(Out) :: &
   RAE(:,:)   ! Interpolated ECHAM specific absorption (z,channel) [km^-1]
!
Real(Double), Intent(Out) :: &
   TDE(:)     ! Dry temperature from ECHAM refractivity [K]
!
! Input optional arguments:
!
Integer, Optional, Intent(In) :: &
   Vrb        ! Verbosity mode.
              ! True by default.
!----------------------------------------------------------
! Local Parameters:
!
Real(Double), Parameter :: &
   DZH   = 0.01_Double,     & ! High-resoltuion grid step
   ZHmax = 130.0_Double       ! Upper limit of high-resolution grid
!
! Local Scalars:
!
Integer           :: LVrb     ! Verbosity level
Integer           :: NGO      ! Number of direct refraction angles
Integer           :: NH       ! High-resolution grid dimension
Integer           :: i        ! Profile index
Real(Double)      :: Zmin     ! Minimum model altitude
Real(Double)      :: Zmax     ! Maximum model altitude
Integer           :: NC       ! Number of frequency channels
Real(Double)      :: ZmaxI    ! Upper integration height [km]
Real(Double)      :: DZI      ! Integration step [km]
Real(Double)      :: Tinit    ! Temperature at upper boundary [K]
Real(Double)      :: Qmin     ! Minimum specific humidity [kg/kg]
!
! Local Arrays:
!
! --- Refraction
!
! --- Standard resolution profiles
!
Real(Double), Allocatable :: &
   RIE(:,:)     ! Profile of Im(N) (z,channel)
!
! --- High resolution profiles
!
Real(Double), Allocatable :: &
   ZH(:),     & ! High-resolution altitude grid
   RNH(:),    & ! High-resolution refractivity profile
   QH(:),     & ! High-resolution humidity profile
   TDH(:)       ! High-resolution temperature profile
!
! --- Wave vector for absorption
!
Real(Double), Allocatable :: &
   Ak(:)        ! Wave vectors
!----------------------------------------------------------


!----------------------------------------------------------
! 1. GETTING LOCAL PROFILES
!----------------------------------------------------------


!--- 1.0. Setting verbosity mode

If (Present(Vrb)) then
   LVrb = Vrb
Else
   LVrb = 0
End If


!--- 1.1. Determination of number of data

NGO = Size(Z)
NC  = Size(AFreq)


!--- 1.2. Interpolation

Allocate(Ak(NC))
Allocate(RIE(NGO,NC))

Ak(:) = 2*Pi*AFreq(:)/C_Light

Do i=1,NGO
   Call Atmosphere_Constituents &
     (GP%Lambda,    & ! <-- Longiude of point [deg]
      GP%Phi,       & ! <-- Latitude of point [deg]
      Z(i),         & ! <-- Altitude of point [km]
      Zmin,         & ! --> Minimum model Z for this lon/lat
      Zmax,         & ! --> Maximum model Z for this lon/lat
      QE(i),        & ! --> Interpolated constituent
      TE(i),        & ! --> Interpolated constituent
      RNE(i),       & ! --> Interpolated N
      RIE(i,:))       ! ~~> Interpolated Im(N)
   RAE(i,:) = RIE(i,:)*Ak(:)
End Do


!--- 1.4. Computation of dry temperature

NH = 1 + Ceiling((ZHmax - Zmin)/DZH)

Allocate(ZH(NH))
Allocate(RNH(NH))
Allocate(QH(NH))
Allocate(TDH(NH))

Do i=1,NH
   ZH(i) = Zmin + Real(i-1)*(ZHmax - Zmin)/Real(NH-1)
   Call Atmosphere_Constituents &
     (GP%Lambda,    & ! <-- Longiude of point [deg]
      GP%Phi,       & ! <-- Latitude of point [deg]
      ZH(i),        & ! <-- Altitude of point [km]
      Zmin,         & ! --> Minimum model Z for this lon/lat
      Zmax,         & ! --> Maximum model Z for this lon/lat
      QH(i),        & ! --> Interpolated humidity [kg/kg]
      TDH(i),       & ! --> Interpolated temperature [K]
      RNH(i))         ! --> Interpolated Re(N)
End Do

ZmaxI = 120.0
DZI   = 0.015
Tinit = 0.0
Qmin  = 1d-7

QH(:) = 0.0

Call NQ_to_TP &
  (GP%Phi,  & ! <-- Geodetic latitude [deg]
   ZH,      & ! <-- Altitude above reference ellipsoid [km]
   RNH,     & ! <-- Profile of refractivity
   QH,      & ! <-- Profile of specific humidity
   ZmaxI,   & ! <-- Upper integration height [km]
   DZI,     & ! <-- Integration step [km]
   Tinit,   & ! <-- Temperature at upper boundary [K]
   Qmin,    & ! <-- Minimum specific humidity [kg/kg]
   TDH)       ! --> Profile of temperature

Do i=1,NGO
   Call Linear &
     (ZH,      & ! <-- Argument grid
      TDH,     & ! <-- Gridded function
      Z(i),    & ! <-- Interpolation point
      TDE(i))    ! --> Interpolated function value
End Do


!----------------------------------------------------------
! 2. MEMORY DEALLOCATION
!----------------------------------------------------------

Deallocate(Ak)
Deallocate(RIE)
Deallocate(ZH)
Deallocate(RNH)
Deallocate(QH)
Deallocate(TDH)


End Subroutine Local_Profiles_Vertical



!==========================================================
Subroutine Local_Profiles_Skew &
  (GPZ,       & ! <-- Perigee locations [geodetic]
   GP,        & ! <-- Occultation point
   AFreq,     & ! <-- Absorption frequency [Hz]
   TE,        & ! --> Interpolated ECHAM temperatures [K]
   QE,        & ! --> Interpolated ECHAM humidity [kg/kg]
   RNE,       & ! --> Interpolated ECHAM refractivity [n-1]
   RAE,       & ! --> Interpolated ECHAM specific absorption (z,channel) [km^-1]
   TDE,       & ! --> Dry temperature from ECHAM refractivity [K]
   Vrb)         ! <~~ Verbosity level
!
! Getting local atmospheric profiles.
!----------------------------------------------------------
! (C) Copyright 2009, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 04 Mar 2009 | Original version.
!   1.1   | 02 Sep 2009 | Improved notation.
!   1.2   | 23 Sep 2009 | Index order changed to (z,channel).
!----------------------------------------------------------
! Modules used:
!
Use Earth, only: &
! Imported Type Definitions:
    Geodetic
!
Use Occ_Inversion, only: &
! Imported Routines:
    Invert_Refraction
!
Use Occ_Meteoprofiles, only: &
! Imported Array Variables:
    NQ_to_TP
!
Use Interpolation, only: &
! Imported Array Variables:
    Linear
!
Use Signal, only: &
! Imported Routines:
    Monotonize
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type(Geodetic), Intent(In) :: &
   GPZ(:)     ! Perigee locations [geodetic]
!
Type(Geodetic), Intent(In) :: &
   GP         ! Occultation point
!
Real(Double), Intent(In) :: &
   AFreq(:)   ! Absorption frequency [Hz]
!
! Output arguments:
!
Real(Double), Intent(Out) :: &
   TE(:)      ! Interpolated ECHAM temperatures [K]
!
Real(Double), Intent(Out) :: &
   QE(:)      ! Interpolated ECHAM humidity [kg/kg]
!
Real(Double), Intent(Out) :: &
   RNE(:)     ! Interpolated ECHAM refractivity [n-1]
!
Real(Double), Intent(Out) :: &
   RAE(:,:)   ! Interpolated ECHAM specific absorption (z,channel) [km^-1]
!
Real(Double), Intent(Out) :: &
   TDE(:)     ! Dry temperature from ECHAM refractivity [K]
!
! Input optional arguments:
!
Integer, Optional, Intent(In) :: &
   Vrb        ! Verbosity mode.
              ! True by default.
!----------------------------------------------------------
! Local Parameters:
!
Real(Double), Parameter :: &
   DZH   = 0.01_Double,   & ! High-resoltuion grid step
   ZHmax = 130.0_Double     ! Upper limit of high-resolution grid
!
! Local Scalars:
!
Integer           :: LVrb     ! Verbosity level
Integer           :: NGO      ! Number of direct refraction angles
Integer           :: NH       ! High-resolution grid dimension
Integer           :: i        ! Profile index
Real(Double)      :: Zmin     ! Minimum model altitude
Real(Double)      :: Zmax     ! Maximum model altitude
Integer           :: NC       ! Number of frequency channels
Real(Double)      :: ZmaxI    ! Upper integration height [km]
Real(Double)      :: DZI      ! Integration step [km]
Real(Double)      :: Tinit    ! Temperature at upper boundary [K]
Real(Double)      :: Qmin     ! Minimum specific humidity [kg/kg]
!
! Local Arrays:
!
! --- Refraction
!
! --- Standard resolution profiles
!
Real(Double), Allocatable :: &
   RIE(:,:)     ! Profile of Im(N) (z,channel)
!
! --- High resolution profiles
!
Type(Geodetic), Allocatable :: &
   GPZH(:)      ! High-resolution grid of perigees
Real(Double), Allocatable :: &
   RNH(:),    & ! High-resolution refractivity profile
   QH(:),     & ! High-resolution humidity profile
   TDH(:)       ! High-resolution temperature profile
!
! --- Wave vector for absorption
!
Real(Double), Allocatable :: &
   Ak(:)        ! Wave vectors
!----------------------------------------------------------


!----------------------------------------------------------
! 1. GETTING LOCAL PROFILES
!----------------------------------------------------------


!--- 1.0. Setting verbosity mode

If (Present(Vrb)) then
   LVrb = Vrb
Else
   LVrb = 0
End If


!--- 1.1. Determination of number of data

NGO = Size(GPZ)
NC  = Size(AFreq)


!--- 1.2. Interpolation

Allocate(Ak(NC))
Allocate(RIE(NGO,NC))

Ak(:) = 2*Pi*AFreq(:)/C_Light

Do i=1,NGO
   Call Atmosphere_Constituents &
     (GPZ(i)%Lambda,    & ! <-- Longiude of point [deg]
      GPZ(i)%Phi,       & ! <-- Latitude of point [deg]
      GPZ(i)%H,         & ! <-- Altitude of point [km]
      Zmin,             & ! --> Minimum model Z for this lon/lat
      Zmax,             & ! --> Maximum model Z for this lon/lat
      QE(i),            & ! --> Interpolated constituent
      TE(i),            & ! --> Interpolated constituent
      RNE(i),           & ! --> Interpolated Re(N)
      RIE(i,:))           ! ~~> Interpolated Im(N)
   RAE(i,:) = RIE(i,:)*Ak(:)
End Do


!--- 1.4. Computation of dry temperature

NH = 1 + Ceiling((ZHmax - Zmin)/DZH)

Allocate(GPZH(NH))
Allocate(RNH(NH))
Allocate(QH(NH))
Allocate(TDH(NH))

Do i=1,NH
   GPZH(i)%H = Zmin + Real(i-1)*(ZHmax - Zmin)/Real(NH-1)
   Call Linear &
     (GPZ(:)%H,         & ! <-- Argument grid
      GPZ(:)%Lambda,    & ! <-- Gridded function
      GPZH(i)%H,        & ! <-- Interpolation point
      GPZH(i)%Lambda,   & ! --> Interpolated function value
      CExt=.True.)        ! <~~ Constant/linear extrapolation
   Call Linear &
     (GPZ(:)%H,         & ! <-- Argument grid
      GPZ(:)%Phi,       & ! <-- Gridded function
      GPZH(i)%H,        & ! <-- Interpolation point
      GPZH(i)%Phi,      & ! --> Interpolated function value
      CExt=.True.)        ! <~~ Constant/linear extrapolation
   Call Atmosphere_Constituents &
     (GPZH(i)%Lambda,   & ! <-- Longiude of point [deg]
      GPZH(i)%Phi,      & ! <-- Latitude of point [deg]
      GPZH(i)%H,        & ! <-- Altitude of point [km]
      Zmin,             & ! --> Minimum model Z for this lon/lat
      Zmax,             & ! --> Maximum model Z for this lon/lat
      QH(i),            & ! --> Interpolated humidity [kg/kg]
      TDH(i),           & ! --> Interpolated temperature [K]
      RNH(i))             ! --> Interpolated Re(N)
End Do

ZmaxI = 120.0
DZI   = 0.015
Tinit = 0.0
Qmin  = 1d-7

QH(:) = 0.0

Call NQ_to_TP &
  (GP%Phi,        & ! <-- Geodetic latitude [deg]
   GPZH(:)%H,     & ! <-- Altitude above reference ellipsoid [km]
   RNH,           & ! <-- Profile of refractivity
   QH,            & ! <-- Profile of specific humidity
   ZmaxI,         & ! <-- Upper integration height [km]
   DZI,           & ! <-- Integration step [km]
   Tinit,         & ! <-- Temperature at upper boundary [K]
   Qmin,          & ! <-- Minimum specific humidity [kg/kg]
   TDH)             ! --> Profile of temperature

Do i=1,NGO
   Call Linear &
     (GPZH(:)%H,  & ! <-- Argument grid
      TDH,        & ! <-- Gridded function
      GPZ(i)%H,   & ! <-- Interpolation point
      TDE(i))       ! --> Interpolated function value
End Do


!----------------------------------------------------------
! 2. MEMORY DEALLOCATION
!----------------------------------------------------------

Deallocate(Ak)
Deallocate(RIE)
Deallocate(GPZH)
Deallocate(RNH)
Deallocate(QH)
Deallocate(TDH)


End Subroutine Local_Profiles_Skew



End Module GO_Propagator



