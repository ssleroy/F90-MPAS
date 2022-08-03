!
Module Ray_Problem
!
! Solving boundary problem for geometric
! optical rays.
!----------------------------------------------------------
! (C) Copyright 2003, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 08 May 2003 | Original version.
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
! Public Parameters:
!
Integer, Parameter :: &
   dim_2D = 2,     & ! 2D simulation
   dim_3D = 3        ! 3D simulation
!
!----------------------------------------------------------

Contains


!==========================================================
Subroutine Solve_Boundary_Problem &
  (SDim,       & ! <-- Spatial dimension
   TR,         & ! <-- Relative time of samples [sec]
   BL,         & ! <-- LEO regression coefficients
   BG,         & ! <-- GPS regression coefficients
   Hmax,       & ! <-- Maximum height
   ERLC,       & ! <-- Curvature center (ECEF)
   RE,         & ! <-- Local curvature radius
   DYN,        & ! <-- Minimum vertical scale of N
   NP,         & ! <-- Upper number of GO rays
   opt_AB,     & ! <-- Modeling of absorption
   NC,         & ! <-- Number of channels for absorption
   FWP,        & ! <-- Filter width for computing DPGPS/DP [km]
   P,          & ! --> Impact parameters
   EP,         & ! --> Refraction angles
   TP,         & ! --> Time as function of P
   YP,         & ! --> Y-coordinate as function of P
   ARP,        & ! --> Refractive amplitude (p)
   ATP,        & ! --> Absorptive attenuation (p,channel)
   TG,         & ! --> Time for GO solution
   S0,         & ! --> Satellite-to-satellite distnce
   YT,         & ! --> Y-coordinate as function of TG
   FT,         & ! --> Model function F(TG)
   Vrb)          ! <~~ Verbosity level
!
! Computation of rays for given GPS and LEO trajectories
!----------------------------------------------------------
! Method:
!   Itertive solution for time as function of
!   impact parameter.
!----------------------------------------------------------
! (C) Copyright 2003-2007, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 22 May 2003 | Original version.
!   2.0   | 20 Aug 2003 | Verbosity level.
!   2.1   | 04 Sep 2003 | Corrected bugs.
!   3.0   | 30 Oct 2003 | DYN.
!   3.1   | 04 Dec 2003 | Locate_Time moved to Occ_Coordinates.
!   3.2   | 01 Mar 2004 | Parameter Use_Doppler.
!   3.3   | 05 Mar 2004 | Monotonization of impact parameters.
!   3.4   | 07 Dec 2005 | Y-coordinate adjustment and
!         |             | time monotinization for variable
!         |             | occultation direction.
!   3.5   | 05 May 2006 | Filter width limitation.
!   3.6   | 24 Feb 2007 | FWP.
!   3.7   | 23 Sep 2009 | Index order changed to (p,channel).
!   3.8   | 08 Oct 2009 | Updated for new version of
!         |             | Interpolate_Trajectory.
!----------------------------------------------------------
! Modules used:
!
Use Coordinates, only: &
! Imported Type Definitions:
    Cartesian,      &
! Imported Routines:
    Vector_Norm,    &
    Vector_Angle,   &
! Imported Operators:
    Operator(+),    &
    Operator(-),    &
    Operator(*),    &
    Operator(.xx.)
!
Use Interpolation, only: &
! Imported Routines:
    Linear
!
Use Atmosphere_rays_adj, only: &
! Imported Routines:
    Find_Ray_P
!
Use Atmosphere_Rays, only: &
! Imported Routines:
    Ray_Trace
!
Use Occ_Coordinates, only: &
! Imported Routines:
    Impact_Parameter,        &
    Interpolate_Trajectory,  &
    Locate_Time
!
Use Occ_Refraction, only: &
! Imported Routines:
    Doppler_to_Refraction
!
Use Occ_Refraction_adj, only: &
! Imported Routines:
    Doppler_to_Refraction_adj
!
Use Signal, only: &
! Imported Routines:
    Monotonize,           &
    Sliding_Polynomial
!
Use IO, only: &
! Imported Routines:
    PutXY
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Integer, Intent(In)     :: &
   SDim        ! Spatial dimension
!
Real(Double), Intent(In)     :: &
   TR(:)       ! Relative time of samples [sec]
!
Real(Double), Intent(In)    :: &
   BG(0:,:)    ! Regression coefficients for VGPS
!
Real(Double), Intent(In)    :: &
   BL(0:,:)    ! Regression coefficients for VLEO
!
Real(Double), Intent(In)    :: &
   Hmax        ! Maximum height
!
Type(Cartesian), Intent(In)  :: &
   ERLC        ! Curvature center (ECEF)
!
Real(Double), Intent(In)     :: &
   RE          ! Local curvature radius
!
Real(Double), Intent(In)     :: &
   DYN         ! Minimum vertical scale of N
!
Integer, Intent(In)          :: &
   NP          ! Upper number of GO rays
!
Logical, Intent(In)          :: &
   opt_AB      ! Modeling of absorption
!
Integer, Intent(In)          :: &
   NC          ! Number of channels for absorption
!
Real(Double), Intent(In)     :: &
   FWP         ! Filter width for computing DPGPS/DP [km]
!
! Output arguments:
!
Real(Double), Pointer     :: &
   P(:)        ! Impact parameters
!
Real(Double), Pointer     :: &
   EP(:)       ! Refraction angles
!
Real(Double), Pointer     :: &
   TP(:)       ! Time as function of P
!
Real(Double), Pointer     :: &
   YP(:)       ! Y-coordinate as function of P
!
Real(Double), Pointer     :: &
   ARP(:)      ! Refractive amplitude as function of P
!
Real(Double), Pointer     :: &
   ATP(:,:)    ! Absorptive attenuation as function of P (p,channel)
!
Real(Double), Pointer     :: &
   TG(:)       ! Time for GO solution
!
Real(Double), Pointer     :: &
   S0(:)       ! Satellite-to-satellite distnce
!
Real(Double), Pointer     :: &
   YT(:)       ! Y-coordinate as function of TG
!
Real(Double), Pointer     :: &
   FT(:)       ! Model function F(TG)
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
   RHmin = 1.5_Double       ! Lower estimate of ray height
!
Real(Double), Parameter :: &
   A = 0.9                  ! Inhomogeneity parameter of grid of p
!
Integer, Parameter :: &
   NIT = 10                 ! Limit of iteration number
!
Logical, Parameter :: &
   Use_Doppler = .True.     ! Computation of P and EP from Doppler
!
Integer, Parameter :: &
   NPD = 3                  ! Polynomial degree for smoothing Doppler
!
! Local Scalars:
!
Integer              :: LVrb        ! Verbosity mode
Character(Len=255)   :: Line        ! Terminal line
Integer              :: N           ! Number of orbit data
Integer              :: NPR         ! Reduced number of impact parameters
Real(Double)         :: P1          ! Straight-line height for TR(N)
Real(Double)         :: PN          ! Straight-line height for TR(N)
Real(Double)         :: Pmin        ! Minimum impact parameter
Real(Double)         :: Pmax        ! Maximum impact parameter
Real(Double)         :: E_d         ! d(Eps)/d(d)
Real(Double)         :: x           ! Grid parameter
Real(Double)         :: DP          ! Grid volume
Real(Double)         :: DPM         ! Monotonization accuracy
Type(Cartesian)      :: XGPS        ! GPS coordinates
Type(Cartesian)      :: VGPS        ! GPS velocity
Type(Cartesian)      :: XLEO        ! LEO coordinates
Type(Cartesian)      :: VLEO        ! LEO velocity
Real(Double)         :: TI1         ! Time for iterations
Real(Double)         :: TI2         ! Time for iterations
Integer              :: i           ! Array index
Integer              :: j           ! Iteration index
Type(Cartesian)      :: UT          ! Transmitter ray direction
Type(Cartesian)      :: XN          ! Ray point nearest to receiver
Type(Cartesian)      :: UN          ! Ray direction at XN
Real(Double)         :: Hper        ! Ray perigee height
Real(Double)         :: DSV         ! Variable integration step
Integer              :: Stat        ! Error status
Real(Double)         :: DL          ! Spatial uncertainty
Real(Double)         :: Tmin        ! Minimum time
Real(Double)         :: Tmax        ! Maximum time
Real(Double)         :: DTM         ! Scale for monotonizing time
Integer              :: Tdir        ! Time direction flag
Real(Double)         :: Ymin        ! Y adjustment
Integer              :: WD          ! Filter width for smoothing Doppler
Real(Double)         :: W0          ! Normalizing constant for amplitude
Integer              :: WP          ! Filter width for smoothing PGPS
!
! Local Arrays:
!
Real(Double), Allocatable    :: &
   RGPS(:),               & ! GPS radius from ERLC
   RLEO(:),               & ! LEO radius from ERLC
   PGPS(:),               & ! Impact parameter at GPS
   FPGPS(:),              & ! Filtered PGPS
   DPGPS(:),              & ! D PGPS / DP
   PLEO(:),               & ! LEO radius from ERLC
   Theta(:)                 ! Satellite-to-satellite angle
!
Real(Double), Allocatable  :: &
   TMP(:),                & ! Monotonized time
   DFP(:),                & ! Doppler as function of P
   DFT(:),                & ! Doppler as function of TG
   DF0(:),                & ! Smoothed Doppler model as function of TG
   P0(:),                 & ! Impact parameter model
   EP0(:),                & ! Refraction angle model
   P0_d(:)                  ! d(P)/d(d)
!
Real(Double)    :: &
   P_XRT(6),              & ! d(P)/d(XR,XT)
   E_XRT(6)                 ! d(Eps)/d(XR,XT)
!----------------------------------------------------------


!----------------------------------------------------------
! 1. DETERMINATION OF IMPACT PARAMETER GRID
!----------------------------------------------------------


!--- 1.0. Setting verbosity level

If (Present(Vrb)) then
   LVrb = Vrb
Else
   LVrb = 0
End If


!--- 1.1. Dimension determination

N = Size(TR)


!--- 1.2. Determination of occultation direction

Call Interpolate_Trajectory  &
  (TR,      & ! <-- Relative time of samples [sec]
   BL,      & ! <-- LEO regression coefficients
   BG,      & ! <-- GPS regression coefficients
   TR(1),   & ! <-- Time
   XLEO,    & ! --> LEO coordinates
   VLEO,    & ! --> LEO velocity
   XGPS,    & ! --> GPS coordinates
   VGPS)      ! --> GPS velocity

P1 = Impact_Parameter &
  (XGPS - ERLC,    & ! <-- Ray beginning
   XLEO - ERLC,    & ! <-- Ray end
   0.0_Double)       ! <-- Refraction angle


Call Interpolate_Trajectory  &
  (TR,      & ! <-- Relative time of samples [sec]
   BL,      & ! <-- LEO regression coefficients
   BG,      & ! <-- GPS regression coefficients
   TR(N),   & ! <-- Time
   XLEO,    & ! --> LEO coordinates
   VLEO,    & ! --> LEO velocity
   XGPS,    & ! --> GPS coordinates
   VGPS)      ! --> GPS velocity

PN = Impact_Parameter &
  (XGPS - ERLC,    & ! <-- Ray beginning
   XLEO - ERLC,    & ! <-- Ray end
   0.0_Double)       ! <-- Refraction angle

If (P1 > PN) then
   TI1  = TR(1)
Else
   TI1  = TR(N)
End If

Pmin = RE + RHmin
Pmax = RE + Hmax


!--- 1.3. Determination of impact parameter grid

!NP = 1 + Ceiling((Pmax - Pmin)/DP)

Allocate(P(NP))
Allocate(EP(NP))
Allocate(TP(NP))
Allocate(YP(NP))
Allocate(ARP(NP))
Allocate(ATP(NP,NC))

Allocate(RGPS(NP))
Allocate(RLEO(NP))
Allocate(PGPS(NP))
Allocate(FPGPS(NP))
Allocate(DPGPS(NP))
Allocate(PLEO(NP))
Allocate(Theta(NP))

Allocate(DFP(NP))
Allocate(DFT(NP))
Allocate(DF0(NP))

Do i=1,NP
!   P(i) = (Pmax*(NP-i) + Pmin*(i-1))/(NP-1)
   x    = Real(NP - i)/Real(NP - 1)
   P(i) = Pmin + ((1-A)*x + A*x*x)*(Pmax-Pmin)
End Do


!--- 1.4. Setting the initial condition for integration step

!DSV = DS
DSV = DYN


!----------------------------------------------------------
! 2. GEOMETRIC OPTICAL SOLUTION
!----------------------------------------------------------


!--- 2.1. Location satellite points for
!---      given impact parameters


Profile: Do i=1,NP

!!!DEBUG
!Write(*,*) 'i = ', i

   !--- 2.1.1. Solving ray boundary problem

   Find_Ray: Do j=1,NIT

      Call Interpolate_Trajectory  &
        (TR,      & ! <-- Relative time of samples [sec]
         BL,      & ! <-- LEO regression coefficients
         BG,      & ! <-- GPS regression coefficients
         TI1,     & ! <-- Time
         XLEO,    & ! --> LEO coordinates
         VLEO,    & ! --> LEO velocity
         XGPS,    & ! --> GPS coordinates
         VGPS,    & ! --> GPS velocity
         ERLC,    & ! <~~ Curvature center (ECEF)
         Theta(i))  ! ~~> Satellite-to-satellite angle

      Call Find_Ray_P &
        (XGPS,   & ! <-- Transmitter position
         VGPS,   & ! <-- Transmitter velocity
         XLEO,   & ! <-- Receiver position
         VLEO,   & ! <-- Receiver velocity
         DSV,    & ! <-- Integration step
         ERLC,   & ! <-- Local curvature center
         P(i),   & ! <-- Impact parameter [km]
         UT,     & ! --> Transmitter ray direction
         XN,     & ! --> Ray point nearest to receiver
         UN,     & ! --> Ray direction at XN
         Hper,   & ! --> Ray perigee height
         Stat,   & ! --> Error status
         .False.)  ! <~~ Strict criterion for iterations


      If (Stat /= 0) then
         Exit Profile
      End If

      DFP(i) = ((VGPS*UT) - (VLEO*UN))/(C_Light - (VGPS*UT))

      If (Use_Doppler) then

         Call Doppler_to_Refraction &
           (XGPS-ERLC,  & ! <-- Transmitter position
            VGPS,       & ! <-- Transmitter velocity
            XN-ERLC,    & ! <-- Receiver position
            VLEO,       & ! <-- Receiver velocity
            DFP(i),     & ! <-- Relative Doppler frequency shift
            P(i),       & ! --> Impact parameter
            EP(i))        ! --> Refraction angle

      Else

         EP(i) = Vector_Angle(XGPS-ERLC,XN-ERLC)      &
                 - ACos(P(i)/Vector_Norm(XGPS-ERLC))  &
                 - ACos(P(i)/Vector_Norm(XN-ERLC))

      End If


      Call Locate_Time  &
        (P(i),    & ! <-- Impact parameter
         EP(i),   & ! <-- Refraction angle
         TR,      & ! <-- Relative time of samples [sec]
         BL,      & ! <-- LEO regression coefficients
         BG,      & ! <-- GPS regression coefficients
         ERLC,    & ! <-- Curvature center
         TI1,     & ! <-- Initial approximation
         TI2)       ! --> Located time

      DL = Abs(TI2 - TI1)*Max(Vector_Norm(VGPS), Vector_Norm(VLEO))
      TP(i) = TI2

      If (DL < 1d-4) then
         Exit Find_Ray
      Else
         TI1 = TI2
      End If

   End Do Find_Ray

   NPR = i

   If (LVrb >= 3) then
      Write (Line,'(2X,A,I5,A,F7.3,A,F6.2,A,ES14.6,A,F10.2,A1)')  &
         'i = ', i, '  P = ', P(i) - RE,  ' DSV = ', DSV,         &
         '  d = ', DFP(i), '  T = ', TP(i),                       &
         Char(0)
      Call CPrintf(Line)
   End If


   !--- 2.1.2. Computation of amplitude 

   RGPS(i) = Vector_Norm(XGPS - ERLC)
   RLEO(i) = Vector_Norm(XLEO - ERLC)
   PGPS(i) = Vector_Norm((XGPS - ERLC) .xx. (UT))
   PLEO(i) = Vector_Norm((XN   - ERLC) .xx. (UN))


   !--- 2.1.3. Computation of attenuation

   If (Opt_AB) then
      Call Ray_Trace &
        (XGPS,      & ! <-- Transmitter position
         UT,        & ! <-- Transmitter ray direction
         XLEO,      & ! <-- Receiver position
         DSV,       & ! <-- Integration step
         XN,        & ! --> Ray point nearest to receiver
         UN,        & ! --> Ray direction at XN
         ATP(i,:),  & ! ~~> Attenuation [Int(Imag(N))]
         Stat = Stat) ! --> Error status
   Else
      ATP(i,:) = 0.0
   End If


   !--- 2.1.4. Estimation of integration step

!   DSV = Min(DS, DYN/Min(0.05, Max(0.0001,Abs(EP(i)))))

End Do Profile

If (LVrb >= 3) then
   Write(*,'()')
End If

!!!DEBUG
!Call PutXY('test-dfp.dat',  DFP(1,1:NPR), P(1:NPR)-RE, Stat=i)
!Call PutXY('test-ep.dat',    EP(1:NPR),   P(1:NPR)-RE, Stat=i)
!Call PutXY('test-tp.dat',    TP(1:NPR),   P(1:NPR)-RE, Stat=i)


!--- 2.2. Shrinking arrays

P   =>   P(1:NPR)
EP  =>  EP(1:NPR)
TP  =>  TP(1:NPR)
YP  =>  YP(1:NPR)
ARP => ARP(1:NPR)
ATP => ATP(1:NPR,:)


!--- 2.3. Computation of amplitude

DP = Abs(P(NPR) - P(1))/(NPR-1)
WP = Max(Ceiling(FWP/Abs(DP)), NPD + 2)

If (LVrb >= 2) then
   Write(*,'(2X,A,I4)') 'WP = ', WP
End If

DPM = 0.010

Call Monotonize &
  (-1,        & ! <-- Array direction
   DPM,       & ! <-- Accuracy
   P,         & ! <-> Sequence to be transformed to monotonic
   LVrb)        ! <~~ Verbosity level

If (WP >= NPD) then
   Call Sliding_Polynomial &
     (P(1:NPR),      & ! <-- Time
      PGPS(1:NPR),   & ! <-- Signal samples
      WP,            & ! <-- Window width [samples]
      NPD,           & ! <-- Polynomial degree
      FPGPS(1:NPR),  & ! --> Filtered signal
      DPGPS(1:NPR))    ! ~~> Signal derivative
Else
   DPGPS(:) = 1.0
End If

!!!DEBUG
!Call PutXY('test-pg.dat',  P(1:NPR)-RE,  FPGPS(1,1:NPR)-P(1:NPR), Stat=i)
!Call PutXY('test-dpg.dat', P(1:NPR)-RE,  DPGPS(1,1:NPR),         Stat=i)


Select Case(SDim)
   Case (dim_2D)
      W0 = 2*Pi*5000.0**2*RGPS(1)
   Case (dim_3D)
      W0 = 4*Pi*5000.0**2*RGPS(1)**2
End Select


Do i=1,NPR
   ARP(i) = Sqrt(W0/(Sqrt(RGPS(i)**2-PGPS(i)**2)*    &
                     Sqrt(RLEO(i)**2-PLEO(i)**2)))*  &
            DPGPS(i)
   If (SDim == dim_3D) then
      ARP(i) = ARP(i)*Sqrt(PGPS(i)/(RGPS(i)*RLEO(i)*Sin(Theta(i))))
   End If
End Do

!!!DEBUG
!Call PutXY('test-arp.dat',  ARP(1:NPR),   P(1:NPR)-RE, Stat=i)
!Call PutXY('test-pg.dat',  PGPS(1,1:NPR) - P(1:NPR),   P(1:NPR)-RE, Stat=i)
!Call PutXY('test-pl.dat',  PLEO(1,1:NPR) - P(1:NPR),   P(1:NPR)-RE, Stat=i)
!Call PutXY('test-dpg.dat', DPGPS(1,1:NPR),   P(1:NPR)-RE, Stat=i)


!----------------------------------------------------------
! 3. DETERMINATION OF TG, Y(TG), S0(TG), F(TG)
!----------------------------------------------------------


!--- 3.1. Array allocation

Allocate(TG(NPR))
Allocate(S0(NPR))
Allocate(YT(NPR))
Allocate(FT(NPR))
Allocate(TMP(NPR))


!--- 3.2. Determination of time for GO solution

Tmin = MinVal(TP)
Tmax = MaxVal(TP)

Do i=1,NPR
   TG(i) = (Tmin*(NPR-i) + Tmax*(i-1))/(NPR-1)
End Do


!--- 3.2. Determination of smooth Doppler model

TMP(:) = TP(:)

DTM  = 10.0/Max(Vector_Norm(VGPS), Vector_Norm(VLEO))
Tdir = Nint(Sign(1.0_Double, TMP(NPR) - TMP(1)))

Call Monotonize &
  (Tdir,      & ! <-- Array direction
   DTM,       & ! <-- Accuracy
   TMP,       & ! <-> Sequence to be transformed to monotonic
   LVrb)        ! <~~ Verbosity mode

!!!DEBUG
!Call PutXY('test-dftm.dat', TMP(1:NPR), DFP(1,1:NPR), Stat=i)

Do i=1,NPR
   Call Linear &
     (TMP(1:NPR),    & ! <-- Argument grid
      DFP(1:NPR),    & ! <-- Gridded function
      TG(i),         & ! <-- Interpolation point
      DFT(i))          ! --> Interpolated function value
End Do

DTM = 20.0/Max(Vector_Norm(VGPS), Vector_Norm(VLEO))
WD  = DTM*(NPR-1)/(Tmax - Tmin)

Call Sliding_Polynomial &
  (TG(1:NPR),     & ! <-- Time
   DFT(1:NPR),    & ! <-- Signal samples
   WD,            & ! <-- Window width [samples]
   NPD,           & ! <-- Polynomial degree
   DF0(1:NPR))      ! --> Filtered signal

!!!DEBUG
!Call PutXY('test-df0.dat', TG(1:NPR), DF0(1,1:NPR), Stat=i)


!--- 3.3. Determination of smooth impact parameter model

Allocate(P0(NPR))
Allocate(EP0(NPR))
Allocate(P0_D(NPR))
Allocate(S0(NPR))

Do i=1,NPR
   Call Interpolate_Trajectory  &
     (TR,         & ! <-- Relative time of samples [sec]
      BL,         & ! <-- LEO regression coefficients
      BG,         & ! <-- GPS regression coefficients
      TG(i),      & ! <-- Time
      XLEO,       & ! --> LEO coordinates
      VLEO,       & ! --> LEO velocity
      XGPS,       & ! --> GPS coordinates
      VGPS,       & ! --> GPS velocity
      SS = S0(i))   ! ~~> Satellite-to-satellite distance
   Call Doppler_to_Refraction_adj &
     (XGPS-ERLC,  & ! <-- Transmitter position
      VGPS,       & ! <-- Transmitter velocity
      XLEO-ERLC,  & ! <-- Receiver position
      VLEO,       & ! <-- Receiver velocity
      DF0(i),     & ! <-- Relative Doppler frequency shift
      P0(i),      & ! --> Impact parameter
      EP0(i),     & ! --> Refraction angle
      P0_d(i),    & ! --> d(P)/d(d)
      P_XRT,      & ! --> d(P)/d(XR,XT)
      E_d,        & ! --> d(Eps)/d(d)
      E_XRT)        ! --> d(Eps)/d(XR,XT)
End Do

!!!DEBUG
!Call PutXY('test-ep0.dat',  EP0(1:NPR),  P0(1:NPR)-RE,   Stat=i)
!Call PutXY('test-ss0.dat',  TG(1:NPR),   S0(1:NPR),   Stat=i)


!--- 3.4. Determination of Y-coordinate

YT(1) = 0.0
Do i=2,NPR
   YT(i) = YT(i-1) - C_Light*(1/P0_d(i) + 1/P0_d(i-1))*(TG(i) - TG(i-1))/2
End Do

Ymin  = MinVal(YT)
YT(:) = YT(:) - Ymin

!!!DEBUG
!Call PutXY('test-yt.dat',  TG(1:NPR),  YT(1:NPR),   Stat=i)


!--- 3.5. Determination of Lagrange manifold

Do i=1,NPR
   Call Linear &
     (TG(1:NPR),    & ! <-- Argument grid
      YT(1:NPR),    & ! <-- Gridded function
      TP(i),        & ! <-- Interpolation point
      YP(i))          ! --> Interpolated function value
End Do

!!!DEBUG
!Call PutXY('test-yp.dat',   YP(1:NPR),   P(1:NPR)-RE, Stat=i)


!--- 3.5. Determination of model function

FT(1) = 0.0
Do i=2,NPR
   FT(i) = FT(i-1) + &
           (P0(i)   - P0_d(i)*DF0(i) +       &
            P0(i-1) - P0_d(i-1)*DF0(i-1))*   &
           (YT(i) - YT(i-1))/2
End Do

!!!DEBUG
!Call PutXY('test-fty.dat',  YT(1:NPR),  FT(1:NPR),   Stat=i)
!Call PutXY('test-df0y.dat', YT(1:NPR), DF0(1,1:NPR), Stat=i)
!Call PutXY('test-p0d.dat',  YT(1:NPR), P0_d(1:NPR), Stat=i)


!----------------------------------------------------------
! 4. MEMORY DEALLOCATION
!----------------------------------------------------------

Deallocate(RGPS)
Deallocate(RLEO)
Deallocate(PGPS)
Deallocate(FPGPS)
Deallocate(DPGPS)
Deallocate(PLEO)
Deallocate(Theta)
Deallocate(DFP)
Deallocate(DFT)
Deallocate(DF0)
Deallocate(P0)
Deallocate(EP0)
Deallocate(P0_D)


End Subroutine Solve_Boundary_Problem



End Module Ray_Problem


