!
Module Vacuum_Propagator
!
! Wave propagation in vacuum.
!----------------------------------------------------------
! (C) Copyright 2004-2006, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 27 Apr 2004 | Original version.
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Parameters:
    Double,        &
    Pi
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
Subroutine Forward_Propagation_Fresnel &
  (XS,     & ! <-- X coordinates of source curve
   YS,     & ! <-- Y coordinates of source curve
   AS,     & ! <-- Amplitude of source field (Y,channel)
   PS,     & ! <-- Accumulated phase of source field (Y,channel)
   k,      & ! <-- Wave vectors [channel]
   TR,     & ! <-- Time for observation curve
   XLEO,   & ! <-- X coordinates of LEO
   YLEO,   & ! <-- Y coordinates of LEO
   FZ,     & ! <-- Initial last Fresnal zone size [Pi rad]
   SR,     & ! <-- Sampling rate [Hz]
   opt_TS, & ! <-- Time scaling
   WT,     & ! --> WP data time
   WX,     & ! --> X coordinates of observation curve
   WY,     & ! --> Y coordinates of observation curve
   WA,     & ! --> Amplitudes of observed field (time,channel)
   WP,     & ! --> Accumulated phase of observed field (time,channel)
   Vrb)      ! <~~ Verbosity level
!
! Calculation of forward propagation of wave field in vacuum
! from given source curve to given observation curve.
!----------------------------------------------------------
! Method:
!   Forward propagation using diffractive integral,
!   calculation of observation curve from LEO orbit using
!   polynomial regression.
!----------------------------------------------------------
! (C) Copyright 1999-2009, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 26 Feb 1999 | Original version.
!   2.0   | 14 Mar 1999 | Dynamical determination of
!         |             | mumber of output data.
!   3.0   | 08 Jul 1999 | Verbosity parameter.
!   3.1   | 11 Feb 2000 | Local arrays allocatable.
!   3.2   | 17 Dec 2000 | FZ as parameter.
!   4.0   | 02 Mar 2001 | SR.
!   4.1   | 14 Jun 2002 | Correct handling of sample rate.
!   4.2   | 01 Feb 2003 | Dynamic number of channels.
!   5.0   | 20 Aug 2003 | Verbosity level.
!   5.1   | 05 Sep 2003 | Bugs corrected.
!   5.2   | 08 Nov 2005 | If SR <= 0, sampling rate unchanged.
!   5.5   | 23 Feb 2007 | opt_TS.
!   5.6   | 16 Jun 2008 | Support of both rising and setting ROs.
!   5.7   | 23 Sep 2009 | Index order changed to (time,channel).
!----------------------------------------------------------
! Modules used:
!
Use Externals, only: &
! Imported Routines:
    CPrintf
!
Use Matrix, only: &
! Imported Routines:
    Regression,         &
    Basic_Polynomials
!
Use Interpolation, only: &
! Imported Routines:
    Polynomial
!
Use Occ_Diffraction, only: &
! Imported Routines:
    Accumulate_Phase,      &
    Diffractive_Integral
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In) :: &
   XS(:)      ! X coordinates of source curve
!
Real(Double), Intent(In) :: &
   YS(:)      ! Y coordinates of source curve
!
Real(Double), Intent(In) :: &
   AS(:,:)    ! Amplitude of source field (Y,channel)
!
Real(Double), Intent(In) :: &
   PS(:,:)    ! Accumulated phase of source field (Y,channel)
!
Real(Double), Intent(In) :: &
   k(:)       ! Wave vectors [channel]
!
Real(Double), Intent(In) :: &
   TR(:)      ! Time for observation curve
!
Real(Double), Intent(In) :: &
   XLEO(:)    ! X coordinates of LEO
!
Real(Double), Intent(In) :: &
   YLEO(:)    ! Y coordinates of LEO
!
Real(Double), Intent(In) :: &
   FZ         ! Initial last Fresnal zone size [Pi rad]
!
Real(Double), Intent(In) :: &
   SR         ! Sampling rate [Hz]
!
Logical, Intent(In)      :: &
   opt_TS     ! Time scaling
              !   True  - scale time to conserve trajectory start/finish time
              !   False - no time scaling
!
! Output arguments:
!
Real(Double), Pointer :: &
   WT(:)      ! WP data time
!
Real(Double), Pointer :: &
   WX(:)      ! X coordinates of observation curve
!
Real(Double), Pointer :: &
   WY(:)      ! Y coordinates of observation curve
!
Real(Double), Pointer :: &
   WA(:,:)    ! Amplitudes of observed field (time,channel)
!
Real(Double), Pointer :: &
   WP(:,:)    ! Accumulated phase of observed field (time,channel)
!
! Input optional arguments:
!
Integer, Optional, Intent(In) :: &
   Vrb        ! Verbosity level.
              ! 0 by default.
!----------------------------------------------------------
! Local Parameters:
!
Integer, Parameter      :: &
   NP = 5        ! Polynomial degree for orbit regression
!
Real(Double), Parameter :: &
   Emax  = 0.005 ! Max parallax of interfering rays
!
! Local Scalars:
!
! --- Variables for dichotomy
!
Real(Double) :: T1   ! Upper bound of dichotomy interval
Real(Double) :: T2   ! Lower bound of dichotomy interval
Real(Double) :: TL   ! Lower time limit
Real(Double) :: TU   ! Upper time limit
Real(Double) :: XL   ! X(TL)
Real(Double) :: YL   ! Y(TL)
Real(Double) :: XU   ! X(TU)
Real(Double) :: YU   ! Y(TU)
!
! --- Wave fields
!
Real(Double) :: Pst    ! Stationary phase
Real(Double) :: DPR    ! Phase change rate
Integer      :: IC     ! Channel number
Real(Double) :: FZmax  ! Maximum Fresnel zone
!
! --- Geometry parameters
!
Integer      :: N      ! Number of observation points
Integer      :: NC     ! Number of channels
Integer      :: OCD    ! Occultation direction
Integer      :: SRG    ! Sampling rate of RO data
Integer      :: NW     ! Number of output data
Integer      :: IH     ! Index of highest point
Integer      :: IL     ! Index of lowest point
Integer      :: DI     ! Index step in observation curve
Integer      :: NSmax  ! Maximum number of subintervals
Real(Double) :: YPmin  ! Minimum Y coordinate of observation point
Real(Double) :: YPmax  ! Maximum Y coordinate of observation point
Real(Double) :: DYP    ! Resolution of observation curve
Integer      :: NS     ! Number of subintervals
Integer      :: j      ! Subinterval index
Real(Double) :: TP     ! Time of observation point
Real(Double) :: XP     ! X coordinate of observation point
Real(Double) :: YP     ! Y coordinate of observation point
!
! --- Work variables
!
Integer           :: LVrb     ! Verbosity mode
Integer           :: i        ! Sample number
Integer           :: Stat     ! Error code
Character(Len=80) :: Line     ! Terminal line
!
!
! Local Arrays:
!
! --- Variables for regression
!
Real(Double), Allocatable    :: &
   TRN(:)               ! Normed TR
Real(Double), Allocatable    :: &
   KR(:,:)              ! Regression matrix for coordinates
Real(Double)    :: &
   CX(0:NP),   &        ! Regression coefficients for X
   CY(0:NP)             ! Regression coefficients for Y
Real(Double), Allocatable  :: &
   WTN(:)               ! Normed WT
!
! --- Wave fields
!
Real(Double), Allocatable :: &
   P0(:,:),    & ! Stationary phase (time,channel)
   DP(:,:)       ! Phase addition (time,channel)
Real(Double), Allocatable :: &
   DPS(:)        ! Accumulated phase addition at TS
Real(Double), Allocatable :: &
   FZL(:)        ! Fresnel zone for current point (channel)
!----------------------------------------------------------



!----------------------------------------------------------
! 0. INITIALIZATION
!----------------------------------------------------------


!--- 0.1. Determination of verbosity mode

If (Present(Vrb)) then
   LVrb = Vrb
Else
   LVrb = 0
End If


!--- 0.2. Determination of maximum Fresnel zones

FZMax = 10*FZ


!--- 0.3. Determination of array size

N  = Size(TR)
NC = Size(k)


!--- 0.4. Determination of occultation direction

OCD = Nint(Sign(1.0_Double, YLEO(N)-YLEO(1)))


!----------------------------------------------------------
! 1. POLYNOMIAL APPROXIMATION OF OBSERVATION CURVE
!----------------------------------------------------------


!--- 1.1. Memory allocation

Allocate(TRN(N))
Allocate(KR(N,0:NP))
Allocate(FZL(NC))


!--- 1.2. Renorming time

If (OCD < 0) then
   TRN(:) = (TR(:) - TR(1))/(TR(N) - TR(1))
Else
   TRN(:) = (TR(:) - TR(N))/(TR(1) - TR(N))
End If




!--- 1.3. Calculation of regression coefficients

Call Basic_Polynomials(TRN, KR, Stat)
Call Regression(KR, XLEO, CX, Stat)
Call Regression(KR, YLEO, CY, Stat)


!----------------------------------------------------------
! 2. DETERMINATION OF OBSERVATION CURVE
!----------------------------------------------------------


!--- 2.1. Search for upper point

T1 = -1.0
T2 =  2.0

FindUpper: Do
   TU = (T1 + T2)/2
   Call Polynomial(CX, TU, XU)
   Call Polynomial(CY, TU, YU)
   Call Diffractive_Integral &
     (XS, YS,                               &
      AS(:,1), AS(:,1), PS(:,1), PS(:,1),   &
      k(1), XU, YU, FZ, FZ, FZL(1),         &
      Pst,  Stat=Stat)
   If (Stat == 2) then
      T1  = TU
   Else
      T2  = TU
   End If
   If (Abs(T1-T2) < 1d-4) then
      Exit FindUpper
   End If
End Do FindUpper


!--- 2.2. Search for lower limit

T1 = -1.0
T2 =  2.0

FindLower: Do
   TL = (T1 + T2)/2
   Call Polynomial(CX, TL, XL)
   Call Polynomial(CY, TL, YL)
   Call Diffractive_Integral &
     (XS, YS,                              &
      AS(:,1), AS(:,1), PS(:,1), PS(:,1),  &
      k(1), XL, YL, FZ, FZ, FZL(1),        &
      Pst,  Stat=Stat)
   If (Stat == 1) then
      T2  = TL
   Else
      T1  = TL
   End If
   If (Abs(T1-T2) < 1d-4) then
      Exit FindLower
   End If
End Do FindLower

If (LVrb >= 2) then
   Write(*,'(2X,A,F8.3,1X,F8.3)')     &
      'TU/TL = ', TU, TL
End If


!--- 2.3. Output data allocation

If (SR > 0.0) then
   SRG = (N-1)/Abs(TR(N) - TR(1))
   NW  = 1 + Ceiling((N-1)*(TL-TU)*(SR/SRG))
Else
   NW  = N
End If

If (LVrb >= 2) then
   Write(*,'(2X,A,I5)') 'Number of output data: ', NW
End If

Allocate(WT(1:NW))       ! --> WP data time
Allocate(WX(1:NW))       ! --> X coordinates of observation curve
Allocate(WY(1:NW))       ! --> Y coordinates of observation curve
Allocate(WA(1:NW,NC))    ! --> Amplitudes of observed field (time,channel)
Allocate(WP(1:NW,NC))    ! --> Accumulated phase of observed field (time,channel)


!--- 2.4. Determination of observation curve

Allocate(WTN(1:NW))
Allocate(P0(1:NW,NC))
Allocate(DP(1:NW,NC))

Do i=1,NW
   If (OCD < 0) then
      WTN(i) = TU    + Real(i-1)*(TL - TU)/Real(NW-1)
      If (opt_TS) then
         WT(i) = TR(1) + Real(i-1)*(TR(N)-TR(1))/Real(NW-1)
      Else
         WT(i) = WTN(i)*(TR(N) - TR(1)) + TR(1)
      End If
   Else
      WTN(i) = TL    + Real(i-1)*(TU - TL)/Real(NW-1)
      If (opt_TS) then
         WT(i) = TR(1) + Real(i-1)*(TR(N)-TR(1))/Real(NW-1)
      Else
         WT(i) = WTN(i)*(TR(1) - TR(N)) + TR(N)
      End If
   End If
   Call Polynomial(CX, WTN(i), WX(i))
   Call Polynomial(CY, WTN(i), WY(i))
End Do



!--- 2.5. Determination of maximum resolution

Call Polynomial(CY, WTN(1),  YPmin)
Call Polynomial(CY, WTN(NW), YPmax)

DYP   = Abs((YPmax - YPmin)/Real(NW-1))
NSmax = Ceiling(MaxVal(Emax*4*DYP*k(:)/(2*Pi)))
!!!TEST
NSmax = 1

If (LVrb >= 2) then
   Write(*,'(2X,A,I4)') 'NSmax = ', NSmax
End If

Allocate(DPS(0:NSmax))


!----------------------------------------------------------
! 3. FORWARD PROPAGATION
!----------------------------------------------------------


!--- 3.1. Determination of processing order

IH = Sum(MaxLoc(WY))
IL = Sum(MinLoc(WY))
DI = Sign(1, IL-IH)


!--- 3.2. Processing upper point

NS = 1

Do IC=1,NC
   Call Diffractive_Integral &
     (XS(:),      & ! <-- X-coordinates of source curve
      YS(:),      & ! <-- Y-coordinates of source curve
      AS(:,IC),   & ! <-- Amplitude of source field
      AS(:,IC),   & ! <-- Amplitude of source field
      PS(:,IC),   & ! <-- Accumulated phase of source field
      PS(:,IC),   & ! <-- Accumulated phase of source field
      k(IC),      & ! <-- Wave vector (2*Pi/Wavelength)
      WX(IH),     & ! <-- X-coordinate of observation point
      WY(IH),     & ! <-- Y-coordinate of observation point
      FZ,         & ! <-- Minimum Fresnel zone size
      FZmax,      & ! <-- Maximum Fresnel zone size
      FZL(IC),    & ! --> Last Fresnel zone size
      P0(IH,IC),  & ! --> Stationary phase in observation point
      WA(IH,IC),  & ! ~~> Amplitude in observation point
      DP(IH,IC),  & ! ~~> Addition to stationary phase
      Stat=Stat)    ! --> Error code
End Do


!--- 3.3. Phase accumulation from up to down

Do i=IH+DI,IL,DI
   Do IC=1,NC
      DPS(0) = DP(i-DI,IC)
      Do j=1,NS
         TP = WTN(i-DI) + Real(j)*(WTN(i)-WTN(i-DI))/Real(NS)
         Call Polynomial(CX, TP, XP)
         Call Polynomial(CY, TP, YP)
         Call Diffractive_Integral &
           (XS(:),      & ! <-- X-coordinates of source curve
            YS(:),      & ! <-- Y-coordinates of source curve
            AS(:,IC),   & ! <-- Amplitude of source field
            AS(:,IC),   & ! <-- Amplitude of source field
            PS(:,IC),   & ! <-- Accumulated phase of source field
            PS(:,IC),   & ! <-- Accumulated phase of source field
            k(IC),      & ! <-- Wave vector (2*Pi/Wavelength)
            XP,         & ! <-- X-coordinate of observation point
            YP,         & ! <-- Y-coordinate of observation point
            FZ,         & ! <-- Minimum Fresnel zone size
            FZmax,      & ! <-- Maximum Fresnel zone size
            FZL(IC),    & ! --> Last Fresnel zone size
            P0(i,IC),   & ! --> Stationary phase in observation point
            WA(i,IC),   & ! ~~> Amplitude in observation point
            DPS(j),     & ! ~~> Addition to stationary phase
            Stat=Stat)    ! --> Error code
      End Do
      Call Accumulate_Phase(DPS(0:NS))
      DP(i,IC) = DPS(NS)
      DPR = MaxVal(Abs(DPS(1:NS)-DPS(0:NS-1)))
      NS  = Min(NSmax, Max(NS, Ceiling(DPR/(Pi/8))))
   End Do
   If (LVrb >= 3 .and. (Modulo(i-1,Ceiling(20.0/NS)) == 0)) then
      Write (Line,'(2X,A,I6,A,F8.2,A,F8.2,A,I6,A,F5.1,1X,F5.1,A1)')  &
         'i = ', i, '   Y = ', YP, '   X = ', XP, '   NS = ', NS,    &
         '   FZL = ', MinVal(FZL(:)/Pi), MaxVal(FZL(:)/Pi),          &
         Char(0)
      Call CPrintf(Line)
   End If
End Do
If (LVrb >= 3) then
   Write (*,'()')
End If


!--- 3.4. Calculation of accumulated phase

WP(:,:) = P0(:,:) + DP(:,:)


!--- 3.5. Memory deallocation

Deallocate(TRN)
Deallocate(KR)
Deallocate(WTN)
Deallocate(P0)
Deallocate(DP)
Deallocate(DPS)
Deallocate(FZL)


End Subroutine Forward_Propagation_Fresnel



!==========================================================
Subroutine Forward_Propagation_FIO &
  (XS,     & ! <-- X coordinates of source curve
   YS,     & ! <-- Y coordinates of source curve
   AS,     & ! <-- Amplitude of source field [channel, Y]
   PS,     & ! <-- Accumulated phase of source field [channel, Y]
   k,      & ! <-- Wave vectors [channel]
   TR,     & ! <-- Time for observation curve
   XLEO,   & ! <-- X coordinates of LEO
   YLEO,   & ! <-- Y coordinates of LEO
   HGmin,  & ! <-- GO estimate of minimum spatial frequency
   HGmax,  & ! <-- GO estimate of maximum spatial frequency
   SR,     & ! <-- Sampling rate [Hz]
   opt_TS, & ! <-- Time scaling
   WT,     & ! --> WP data time
   WX,     & ! --> X coordinates of observation curve
   WY,     & ! --> Y coordinates of observation curve
   WA,     & ! --> Amplitudes of observed field (time,channel)
   WP,     & ! --> Accumulated phase of observed field (time,channel)
   Vrb)      ! <~~ Verbosity level
!
! Calculation of forward propagation of wave field in vacuum
! from vertical source line to given observation curve.
!----------------------------------------------------------
! Method:
!   Forward propagation using FIO associated with linearized
!   canonical transform from (y,eta) to (t,sigma).
!----------------------------------------------------------
! (C) Copyright 2004-2009, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 10 May 2004 | Original version.
!   1.1   | 05 Sen 2004 | Increased border of Y grid.
!   2.0   | 17 Sen 2004 | HGmin, HGmax;
!         |             | initialization of Hmax.
!   2.1   | 08 Nov 2005 | If SR <= 0, sampling rate unchanged.
!   2.2   | 06 Dec 2005 | Processing both rising and
!         |             | setting occultations.
!   2.3   | 04 May 2006 | Adjusted DYH.
!   2.4   | 09 May 2006 | Increased time boundary for shadow,
!         |             | GO shadow border estimate,
!         |             | IR/IC min/max index restriction.
!   2.5   | 07 Feb 2007 | IHL index-out-of-range check.
!   2.6   | 23 Feb 2007 | opt_TS.
!   2.7   | 23 Sep 2009 | Index order changed to (time,channel).
!----------------------------------------------------------
! Modules used:
!
Use Externals, only: &
! Imported Routines:
    CPrintf
!
Use Matrix, only: &
! Imported Routines:
    Regression,         &
    Basic_Polynomials
!
Use Interpolation, only: &
! Imported Routines:
    Polynomial,        &
    Nearest_Power2,    &
    Linear
!
Use Occ_Diffraction, only: &
! Imported Routines:
    Accumulate_Phase,      &
    Diffractive_Integral
!
Use Signal, only: &
! Imported Routines:
    Sliding_Polynomial
!
Use FFTW, only: &
! Imported Routines:
    FFT1
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In) :: &
   XS(:)      ! X coordinates of source curve
!
Real(Double), Intent(In) :: &
   YS(:)      ! Y coordinates of source curve
!
Real(Double), Intent(In) :: &
   AS(:,:)    ! Amplitude of source field [channel, Y]
!
Real(Double), Intent(In) :: &
   PS(:,:)    ! Accumulated phase of source field [channel, Y]
!
Real(Double), Intent(In) :: &
   k(:)       ! Wave vectors [channel]
!
Real(Double), Intent(In) :: &
   TR(:)      ! Time for observation curve
!
Real(Double), Intent(In) :: &
   XLEO(:)    ! X coordinates of LEO
!
Real(Double), Intent(In) :: &
   YLEO(:)    ! Y coordinates of LEO
!
Real(Double), Intent(In) :: &
   HGmin      ! GO estimate of minimum spatial frequency
!
Real(Double), Intent(In) :: &
   HGmax      ! GO estimate of maximum spatial frequency
!
Real(Double), Intent(In) :: &
   SR         ! Sampling rate [Hz]
!
Logical, Intent(In)      :: &
   opt_TS     ! Time scaling
              !   True  - scale time to conserve trajectory start/finish time
              !   False - no time scaling
!
! Output arguments:
!
Real(Double), Pointer :: &
   WT(:)      ! WP data time
!
Real(Double), Pointer :: &
   WX(:)      ! X coordinates of observation curve
!
Real(Double), Pointer :: &
   WY(:)      ! Y coordinates of observation curve
!
Real(Double), Pointer :: &
   WA(:,:)    ! Amplitudes of observed field (time,channel)
!
Real(Double), Pointer :: &
   WP(:,:)    ! Accumulated phase of observed field (time,channel)
!
! Input optional arguments:
!
Integer, Optional, Intent(In) :: &
   Vrb        ! Verbosity level.
              ! 0 by default.
!----------------------------------------------------------
! Local Parameters:
!
Integer, Parameter      :: &
   NP = 5        ! Polynomial degree for regression
!
Integer, Parameter      :: &
   NFP = 3       ! Polynomial degree for filtering phase
!
Integer, Parameter      :: &
   NPH = 2       ! Polynomial degree for phase regression
!
Real(Double), Parameter :: &
   DYR  = 0.02   ! Step of rare grid
!
Real(Double), Parameter :: &
   NPM  = 1024   ! Number of points for computing phase model
!
Real(Double), Parameter :: &
   DYFP = 5.0    ! Vertical scale for filtering phase
!
Real(Double), Parameter :: &
   Hacc = 1d-5   ! Accuracy of monotonizing frequency
!
Complex(Double), Parameter ::  &
   Ci = (0.0_Double, 1.0_Double)  ! I = Sqrt(-1)
!
! Local Scalars:
!
! --- Variables for dichotomy
!
Real(Double) :: T1   ! Upper bound of dichotomy interval
Real(Double) :: T2   ! Lower bound of dichotomy interval
Real(Double) :: TL   ! Lower time limit
Real(Double) :: TU   ! Upper time limit
Real(Double) :: XL   ! X(TL)
Real(Double) :: YL   ! Y(TL)
Real(Double) :: XU   ! X(TU)
Real(Double) :: YU   ! Y(TU)
!
! --- Wave fields
!
Real(Double) :: Pst    ! Stationary phase
Integer      :: IC     ! Channel number
Real(Double) :: FZ     ! Minimum Fresnel zones
!
! --- Geometry parameters
!
Integer      :: NS     ! Number of source field data
Integer      :: N      ! Number of observation points
Integer      :: NC     ! Number of channels
Integer      :: SRG    ! Sampling rate of RO data
Integer      :: NW     ! Number of output data
Integer      :: OCD    ! Occultation direction
Real(Double) :: YPmin  ! Minimum Y coordinate of observation point
Real(Double) :: YPmax  ! Maximum Y coordinate of observation point
Real(Double) :: YSmin  ! Minimum Y of source field
Real(Double) :: YSmax  ! Maximum Y of source field
Real(Double) :: DYS    ! Safety border width
Real(Double) :: DYP    ! Resolution of observation curve
Real(Double) :: XS0    ! X-coordinate of shifted origin
!
! --- Filtering parameters
!
Integer      :: NR     ! Input grid skip
Integer      :: NWR    ! Output grid skip
Integer      :: WFP    ! Window width for filterin phase
Integer      :: Hdir   ! Direction of frequency array
Real(Double) :: DHD    ! Eta-interval for differentiation
Integer      :: DID    ! Index-interval for differentiation
Real(Double) :: Hgc    ! Gravity center of spectrum
Integer      :: Igc    ! Index of gravity center
Integer      :: IHmax  ! Index of Hmax
Real(Double) :: Hmax   ! Maximum spatial frequency
Real(Double) :: Hmin   ! Minimum spatial frequency
Real(Double) :: HL     ! Eta where phase is lost
Integer      :: IHL    ! Index of HL
Real(Double) :: DHL    ! Interval for regression
Integer      :: IRmin  ! Beginning of regression interval
Integer      :: IRmax  ! End of regression interval
Integer      :: ICmin  ! Beginning of correction interval
Integer      :: ICmax  ! End of correction interval
Real(Double) :: YSLmin ! Minimum YS in light zone
Real(Double) :: YSLmax ! Maximum YS in light zone
Integer      :: ISLmin ! Index of YSLmin
Integer      :: ISLmax ! Index of YSLmax
Real(Double) :: HGO    ! Geometric optical H
Real(Double) :: TLGO   ! Geometric optical estimate of TL
Integer      :: IH     ! Upper border index
Integer      :: IL     ! Lower border index
Integer      :: IM     ! Middle index
Integer      :: DI     ! Index step
Real(Double) :: DHT    ! Limit increment of HT
!
! --- Work variables
!
Integer           :: LVrb     ! Verbosity mode
Integer           :: i        ! Sample number
Integer           :: Stat     ! Error code
Character(Len=80) :: Line     ! Terminal line
Real(Double)      :: SW       ! Swap variable
!
! --- Hi-res grid parameters
!
Integer           :: NH       ! Dimension of high-resolution grids
Real(Double)      :: DYH      ! Hi-res Y-grid step
Real(Double)      :: DHH      ! Hi-res eta-grid step
Real(Double)      :: XImin    ! Minimum of xi hi-res grid
Real(Double)      :: XImax    ! Maximum of xi hi-res grid
Real(Double)      :: DXIH     ! Hi-res xi-grid step
Real(Double)      :: DXIm     ! XImin/max adjustment
Integer           :: FSign    ! FFT direction
Real(Double)      :: FII      ! Interpolated 1st model
Real(Double)      :: NUI      ! Interpolated measure density
Real(Double)      :: DTH      ! Time grid step
Real(Double)      :: DTA      ! Time apodization border
!
! Local Arrays:
!
! --- Variables for regression
!
Real(Double), Allocatable    :: &
   TRN(:)               ! Normed TR
Real(Double), Allocatable    :: &
   KR(:,:)              ! Regression matrix for coordinates
Real(Double)    :: &
   CX(0:NP),   &        ! Regression coefficients for X
   CY(0:NP)             ! Regression coefficients for Y
Real(Double), Allocatable  :: &
   WTN(:)               ! Normed WT
Real(Double), Allocatable    :: &
   KPH(:,:)             ! Regression matrix for phase
Real(Double)    :: &
   CPH(0:NPH)           ! Regression coefficients for phase
!
! --- Orbit data
!
Real(Double), Allocatable  :: &
   DWX(:),            & ! X-component of velocity
   DWY(:)               ! Y-component of velocity
!
! --- Wave fields
!
Real(Double), Pointer :: &
   PAS(:),     & ! Pointer for AS
   PPS(:),     & ! Pointer for PS
   PXS(:),     & ! Pointer for XS
   PYS(:)        ! Pointer for YS
Target :: AS, PS, XS, YS
!
Real(Double), Allocatable :: &
   AFS(:),     & ! Corrected amplitude
   PFS(:),     & ! Corrected phase
   P0(:),      & ! Stationary phase
   PF(:),      & ! Filtered phase
   DPF(:),     & ! Derivative of phase
   HT(:)         ! Frequency model
!
! --- Phase models
!
Real(Double), Allocatable :: &
   FXY(:),     & ! Auxilliary function
   F(:),       & ! Derivative of first phase model (eta)
   FI(:),      & ! First phase model (eta)
   XI(:),      & ! Transformed frequency
   G(:),       & ! Derivative of second phase model (t)
   GI(:)         ! Second phase model (t)
!
! --- High-resolution wave field
!
Real(Double), Allocatable :: &
   YH(:),      & ! Y-grid for source field
   HH(:),      & ! Eta-grid for Fourier transformed field
   XIH(:),     & ! Xi-grid for Fourier transformed field
   TH(:),      & ! T-grid for propagated field
   HXIH(:),    & ! H(XIH) eta-grid for xi-grid
   AH(:),      & ! Amplitude
   PH(:),      & ! Phase
   DPH(:),     & ! Phase derivative
   DDPH(:),    & ! Monotonized phase derivative
   AFH(:),     & ! Filtered amplitude
   PFH(:),     & ! Filtered phase
   AIH(:),     & ! Interpolated amplitude
   PIH(:)        ! Interpolated phase
!
Complex(Double), Allocatable :: &
   UH(:)         ! Source field
!
! --- Debug variables
!
!Real(Double), Allocatable :: &
!   W(:)
!----------------------------------------------------------



!----------------------------------------------------------
! 0. INITIALIZATION
!----------------------------------------------------------


!--- 0.1. Determination of verbosity mode

If (Present(Vrb)) then
   LVrb = Vrb
Else
   LVrb = 0
End If


!--- 0.2. Array size determination

NS = Size(YS)
N  = Size(TR)
NC = Size(k)


!--- 0.3. Determination of occultation direction

OCD = Nint(Sign(1.0_Double, YLEO(N)-YLEO(1)))


!--- 0.3. Determination of Fresnel zones

FZ    = 60*Pi


!--- 0.4. Computation of Y-grid limits

DYS   = 0.35*Abs(YS(NS) - YS(1))
YSmin = MinVal(YS) - DYS
YSmax = MaxVal(YS) + DYS


!--- 0.5. Computation of shifted origin

XS0   = MinVal(XS)


!----------------------------------------------------------
! 1. POLYNOMIAL APPROXIMATION OF OBSERVATION CURVE
!----------------------------------------------------------


!--- 1.1. Memory allocation

Allocate(TRN(N))
Allocate(KR(N,0:NP))


!--- 1.2. Renorming time

If (OCD < 0) then
   TRN(:) = (TR(:) - TR(1))/(TR(N) - TR(1))
Else
   TRN(:) = (TR(:) - TR(N))/(TR(1) - TR(N))
End If


!--- 1.3. Calculation of regression coefficients

Call Basic_Polynomials(TRN, KR, Stat)
Call Regression(KR, XLEO, CX, Stat)
Call Regression(KR, YLEO, CY, Stat)


!----------------------------------------------------------

Channels: Do IC=1,NC

!----------------------------------------------------------
! 2. PREPROCESSING
!----------------------------------------------------------


!--- 2.1. Computation of hi-res grid

DYH  = (4.*Pi/4.)/(k(IC)*Max(Abs(HGmax), Abs(HGmin)))
NH   = Nearest_Power2(Ceiling((YSmax - YSmin)/DYH))
DYH  = (YSmax - YSmin)/NH

If (LVrb >= 1) then
   Write(*,'(2X,A,I1,2X,A,I7)') &
      'Channel ', IC, 'NH = ', NH
End If


!--- 2.2. Interpolation of complex field

Allocate(YH(NH))
Allocate(HH(NH))
Allocate(UH(NH))
Allocate(AH(NH))
Allocate(PH(NH))

Do i=1,NH
   YH(i) = YSmin + (i-1)*DYH
End Do

Call Interpolate_field &
  (YS(:),     & ! <-- First regular grid
   AS(:,IC),  & ! <-- Amplitude on first grid
   PS(:,IC),  & ! <-- Phase on first grid
   YH,        & ! <-- Second grid
   AH,        & ! --> Amplitude on second grid
   PH)          ! --> Phase on second grid

UH(:) = AH(:)*Exp(Ci*Modulo(PH(:),2*Pi))


!--- 2.3. Fourier transform

FSign = -1

Call FFT1 &
  (FSign,    & ! <-- Transform direction.
   UH(:))      ! <-> Data to be transformed.

HH(1) = 0
Do i=2,NH/2+1
   HH(i)      = (i-1)*2*Pi/(k(IC)*NH*DYH)
   HH(NH-i+2) = -HH(i)
End Do

DHH   = 2*Pi/(k(IC)*NH*DYH)

UH(:) = UH(:)*Sqrt(k(IC)/(2*Pi))*DYH

HH(:) = CShift(HH, NH/2)
UH(:) = CShift(UH, NH/2)


!--- 2.4. Computation of phase and amplitude

AH(:) = Abs(UH(:))

Where (UH(:) /= 0.0_Double)
   PH(:) = ATan2(AImag(UH(:)),Real(UH(:)))
Elsewhere
   PH(:) = 0.0
End Where

Call Accumulate_Phase &
  (PH(:),       & ! <-> Array of (accumulated) phase
   -1)            ! <~~ Phase change direction


!--- 2.5. Differentiation of phase

Allocate(DPH(NH))

DHD = 5.0d-5
DID = Ceiling(DHD/DHH)

Do i=1+DID,NH-DID
   DPH(i) = YSmin - (PH(i+DID) - PH(i-DID))/      &
                    (k(IC)*(HH(i+DID)-HH(i-DID)))
End Do

DPH(1:DID)       = DPH(1+DID)
DPH(NH-DID+1:NH) = DPH(NH-DID)

!!!DEBUG
!Call PutXY('test-2.5-ah.dat', HH(1::16),  AH(1::16), Stat=i)
!Call PutXY('test-2.5-yh.dat', HH(1::16), DPH(1::16), Stat=i)


!--- 2.6. Determination of low amplitude area

Allocate(DDPH(NH))

Hgc   = Sum(HH(:)*AH(:)**2)/Sum(AH(:)**2)
Igc   = Ceiling((Hgc - HH(1))/DHH)

DDPH(Igc) = DPH(Igc)

Do i=Igc-1,1,-1
   DDPH(i) = Min(DPH(i), DDPH(i+1))
End Do
DDPH(1:Igc) = DPH(1:Igc) - DDPH(1:Igc)

!!!DEBUG
!Call PutXY('test-2.6-dyh.dat', HH(1:Igc:16),DDPH(1:Igc:16), Stat=i)


IHL = Sum(MaxLoc(HH(1:Igc), Mask = (DDPH(:) > 10.0)))
If ((IHL < 1) .or. (IHL > Igc)) then
   IHL = 1
End If
HL  = HH(IHL) + 0.006
DHL = 0.005

IHmax = Sum(MinLoc(HH(:), HH(:) > HH(Igc) .and. AH(:) < 0.3*AH(Igc)))
Hmax  = HH(IHmax)

Hmin  = HGmin

If (LVrb >= 2) then
   Write(*,'(2X,A,F8.5,1X,F8.5,2X,A,F8.5)')  &
      'Hmin/max = ', Hmin, Hmax,             &
      'HL = ', HL
End If

Deallocate(DDPH)


!--- 2.7. Regression

IRmin = Max(1,  Ceiling((HL - HH(1))/DHH))
IRmax = Min(NH, Ceiling((HL + DHL - HH(1))/DHH))

Allocate(KPH(IRmin:IRmax,0:NPH))

Call Basic_Polynomials(HH(IRmin:IRmax), KPH, Stat)
Call Regression(KPH, PH(IRmin:IRmax), CPH, Stat)

Deallocate(KPH)


!--- 2.8. Redefinition of phase and amplitude

Allocate(AFH(NH))
Allocate(PFH(NH))

AFH(:) = AH(:)
PFH(:) = PH(:)

ICmin = Max(1,  Ceiling((Hmin - 0.001 - HH(1))/DHH))
ICmax = Min(NH, Ceiling((HL   - HH(1))/DHH))

Do i=ICmin,ICmax
   AFH(i) = Sum(AH(IRmin:IRmax))/Size(AH(IRmin:IRmax))
   Call Polynomial(CPH, HH(i), PFH(i))
End Do

Where (Hmin - 0.001 <= HH(:) .and. HH(:) < Hmax)
   AFH(:) = 1.0
Elsewhere
   AFH(:) = 0.0
End Where


!--- 2.9. Determination of light zone

Call Polynomial(CPH, HH(ICmin), PFH(ICmin), YSLmin)
YSLmin = YSmin - YSLmin/k(IC)

YSLmax = DPH(IHmax)

ISLmin = Sum(MinLoc(Abs(YS(:) - YSLmin)))
ISLmax = Sum(MinLoc(Abs(YS(:) - YSLmax)))

If (LVrb >= 2) then
   Write(*,'(2X,A,F8.3,1X,F8.3)')     &
      'YSLmin/max = ', YSLmin, YSLmax
End If


Deallocate(DPH)


!--- 2.10. Inverse Fourier transform

UH(:) = AFH(:)*Exp(Ci*Modulo(PFH(:),2*Pi))
UH(:) = CShift(UH, -NH/2)

FSign = 1

Call FFT1 &
  (FSign,    & ! <-- Transform direction.
   UH(:))      ! <-> Data to be transformed.

AFH(:) = Abs(UH(:))

Where (UH(:) /= 0.0_Double)
   PFH(:) = ATan2(AImag(UH(:)),Real(UH(:)))
Elsewhere
   PFH(:) = 0.0
End Where

Call Accumulate_Phase &
  (PFH(:))        ! <-> Array of (accumulated) phase


!--- 2.11. Interpolation of corrected phase

Allocate(AFS(NS))
Allocate(PFS(NS))

Call Interpolate_field &
  (YH,        & ! <-- First regular grid
   AFH,       & ! <-- Amplitude on first grid
   PFH,       & ! <-- Phase on first grid
   YS,        & ! <-- Second grid
   AFS,       & ! --> Amplitude on second grid
   PFS)         ! --> Phase on second grid


Deallocate(YH)
Deallocate(AFH)
Deallocate(PFH)


!----------------------------------------------------------
! 3. DETERMINATION OF OBSERVATION CURVE
!----------------------------------------------------------

If (IC == 1) then


!--- 3.1. Setting light zone

PXS => XS(ISLmin:ISLmax)
PYS => YS(ISLmin:ISLmax)
PAS => AS(ISLmin:ISLmax,IC)
PPS => PS(ISLmin:ISLmax,IC)


!--- 3.2. Search for upper point

T1 = -1.0
T2 =  2.0

FindUpper: Do
   TU = (T1 + T2)/2
   Call Polynomial(CX, TU, XU)
   Call Polynomial(CY, TU, YU)
   Call Diffractive_Integral &
     (PXS, PYS,                  &
      PAS, PAS, PPS, PPS,        &
      k(1), XU, YU, FZ, FZ, FZ,  &
      Pst,  Stat=Stat)
   If (Stat == 2) then
      T1  = TU
   Else
      T2  = TU
   End If
   If (Abs(T1-T2) < 1d-4) then
      Exit FindUpper
   End If
End Do FindUpper


!--- 3.3. Search for lower limit

T1 = -1.0
T2 =  2.0

FindLower: Do
   TL = (T1 + T2)/2
   Call Polynomial(CX, TL, XL)
   Call Polynomial(CY, TL, YL)
   Call Diffractive_Integral &
     (PXS, PYS,                  &
      PAS, PAS, PPS, PPS,        &
      k(1), XL, YL, FZ, FZ, FZ,  &
      Pst,  Stat=Stat)
   If (Stat == 1) then
      T2  = TL
   Else
      T1  = TL
   End If
   If (Abs(T1-T2) < 1d-4) then
      Exit FindLower
   End If
End Do FindLower

TL = TL + 0.1


!--- 3.4. Search for GO shadow

T1 = -1.0
T2 =  2.0

FindGOshadow: Do
   TLGO = (T1 + T2)/2
   Call Polynomial(CX, TLGO, XL)
   Call Polynomial(CY, TLGO, YL)
   HGO = (YL - YS(ISLmin))/(XL - XS(ISLmin))
   If (HGO < HGmin - 0.005) then
      T2  = TLGO
   Else
      T1  = TLGO
   End If
   If (Abs(T1-T2) < 1d-4) then
      Exit FindGOshadow
   End If
End Do FindGOshadow

TL = TLGO

If (LVrb >= 2) then
   Write(*,'(2X,A,F8.3,1X,F8.3)')     &
      'TU/TL = ', TU, TL
End If


!--- 3.5. Output data allocation

If (SR > 0.0) then
   SRG = (N-1)/Abs(TR(N) - TR(1))
   NW  = 1 + Ceiling((N-1)*(TL-TU)*(SR/SRG))
Else
   NW  = N
End If

If (LVrb >= 2) then
   Write(*,'(2X,A,I6)') 'Number of output data: ', NW
End If

Allocate(WT(1:NW))
Allocate(WX(1:NW))
Allocate(WY(1:NW))
Allocate(WA(1:NW,NC))
Allocate(WP(1:NW,NC))


!--- 3.6. Determination of observation curve

Allocate(WTN(1:NW))
Allocate(DWX(1:NW))
Allocate(DWY(1:NW))

Do i=1,NW
   If (OCD < 0) then
      WTN(i) = TU    + Real(i-1)*(TL - TU)/Real(NW-1)
      If (opt_TS) then
         WT(i) = TR(1) + Real(i-1)*(TR(N)-TR(1))/Real(NW-1)
      Else
         WT(i) = WTN(i)*(TR(N) - TR(1)) + TR(1)
      End If
   Else
      WTN(i) = TL    + Real(i-1)*(TU - TL)/Real(NW-1)
      If (opt_TS) then
         WT(i) = TR(1) + Real(i-1)*(TR(N)-TR(1))/Real(NW-1)
      Else
         WT(i) = WTN(i)*(TR(1) - TR(N)) + TR(N)
      End If
   End If
   Call Polynomial(CX, WTN(i), WX(i), DWX(i))
   Call Polynomial(CY, WTN(i), WY(i), DWY(i))
End Do


!--- 3.7. Determination of maximum resolution

Call Polynomial(CY, WTN(1),  YPmin)
Call Polynomial(CY, WTN(NW), YPmax)

DYP   = Abs((YPmax - YPmin)/Real(NW-1))


End If


!----------------------------------------------------------
! 4. COMPUTATION OF MODELS
!----------------------------------------------------------


!--- 4.1. Determination of grid skip

NR  = Ceiling(DYR/Abs(YS(2) - YS(1)))
NWR = Ceiling(NW/NPM)

If (LVrb >= 2) then
   Write(*,'(2X,A,I5,1X,I5)') 'Grid skips = ', NR, NWR
End If


!--- 4.1. Computing stationary optical path

Allocate(P0(1:NW))

Do i=1,NW,NWR
   Call Diffractive_Integral &
     (XS(1::NR),     & ! <-- X-coordinates of source curve
      YS(1::NR),     & ! <-- Y-coordinates of source curve
      AFS(1::NR),    & ! <-- Amplitude of source field
      AFS(1::NR),    & ! <-- Amplitude of source field
      PFS(1::NR),    & ! <-- Accumulated phase of source field
      PFS(1::NR),    & ! <-- Accumulated phase of source field
      k(IC),         & ! <-- Wave vector (2*Pi/Wavelength)
      WX(i),         & ! <-- X-coordinate of observation point
      WY(i),         & ! <-- Y-coordinate of observation point
      FZ,            & ! <-- Minimum Fresnel zone size
      FZ,            & ! <-- Maximum Fresnel zone size
      FZ,            & ! --> Last Fresnel zone size
      P0(i),         & ! --> Stationary phase in observation point
      Stat=Stat)       ! --> Error code
   If (LVrb >= 3 .and. (Modulo(i-1,20*NWR) == 0)) then
      Write (Line,'(2X,A,I6,A,F8.2,A,F8.2,A,A1)')        &
         'i = ', i, '   Y = ', WY(i), '   X = ', WX(i),  &
         Char(0)
      Call CPrintf(Line)
   End If
End Do
If (LVrb >= 3) then
   Write (*,'()')
End If

Do i=1,NW
   Call Linear &
     (WTN(1::NWR),  & ! <-- Argument grid
      P0(1::NWR),   & ! <-- Gridded function
      WTN(i),       & ! <-- Interpolation point
      P0(i))          ! --> Interpolated function value
End Do

Deallocate(AFS)
Deallocate(PFS)


!--- 4.2. Differentiating optical path

Allocate(PF(1:NW))
Allocate(DPF(1:NW))

WFP = Ceiling(NW*DYFP/Abs(WY(NW) - WY(1)))

Call Sliding_Polynomial &
  (WTN,       & ! <-- Time
   P0(:),     & ! <-- Signal samples
   WFP,       & ! <-- Window width [samples]
   NFP,       & ! <-- Polynomial degree
   PF(:),     & ! --> Filtered signal
   DPF(:))      ! ~~> Signal derivative

DPF(:) = DPF(:)/k(IC)


Deallocate(P0)
Deallocate(PF)


!--- 4.3. Computation of frequency model

Allocate(HT(1:NW))

HT(:) = (DWY(:)*DPF(:) - Sign(1.0_Double, DWY(:))* &
         DWX(:)*Sqrt(DWX(:)**2 + DWY(:)**2 - DPF(:)**2))/  &
        (DWX(:)**2 + DWY(:)**2)


!!!DEBUG
!Allocate(W(NW))
!W(:) = WY(:) - (WX(:)-XS0)*HT(:)/Sqrt(1 - HT(:)**2)
!Call PutXY('test-4.3-yh.dat', HT(:), W(:), Stat=Stat)
!Deallocate(W)


Deallocate(DPF)


!--- 4.4. Monotonization frequency model

Hdir = NInt(Sign(1.0_Double, WY(NW) - WY(1)))

If (Hdir == 1) then
   IH = NW
   IL = 1
   DI = -1
Else
   IH = 1
   IL = NW
   DI = 1
End If

IM = (IH + IL)/2

DHT = 0.01*Abs(Hmax - Hmin)/(NW-1)

Do i=IM+DI,IL,DI
   HT(i) = Min(HT(i-DI) - DHT, HT(i))
End Do

Do i=IM+DI,IH,-DI
   HT(i) = Max(HT(i+DI) + DHT, HT(i))
End Do

!!!DEBUG
!Allocate(W(NW))
!W(:) = WY(:) - (WX(:)-XS0)*HT(:)/Sqrt(1 - HT(:)**2)
!Call PutXY('test-4.4-yh.dat', HT(:), W(:), Stat=Stat)
!Deallocate(W)
!Stop 'DEBUG TO HERE'


!--- 4.5. Computation of transformed frequency

Allocate(FXY(1:NW))
Allocate(XI(1:NW))

FXY(:) = DWY(:) - DWX(:)*HT(:)/Sqrt(1 - HT(:)**2)

XI(1) = 0
Do i=2,NW
   XI(i) = XI(i-1) + &
           (FXY(i)+FXY(i-1))*(HT(i)-HT(i-1))/2
End Do


!--- 4.6. Computation of first phase model

Allocate(F(1:NW))
Allocate(FI(1:NW))

F(:)   = WTN(:) - &
         ((WY(:)-YSmin) - (WX(:)-XS0)*HT(:)/Sqrt(1 - HT(:)**2))/FXY(:)
FI(1)  = 0
Do i=2,NW
   FI(i) = FI(i-1) + &
           (F(i) + F(i-1))*(XI(i) - XI(i-1))/2
End Do


Deallocate(F)


!--- 4.7. Computation of the second model

Allocate(G(1:NW))
Allocate(GI(1:NW))

G(:)  = DWX(:)*Sqrt(1 - HT(:)**2) + DWY(:)*HT(:) - XI(:)

GI(1) = 0
Do i=2,NW
   GI(i) = GI(i-1) + &
           (G(i) + G(i-1))*(WTN(i) - WTN(i-1))/2
End Do


Deallocate(G)


!----------------------------------------------------------
! 5. FOURIER INTEGRAL OPERATOR
!----------------------------------------------------------


!--- 5.1. Interpolation to hi-res grid of XI

Allocate(HXIH(NH))
Allocate(XIH(NH))
Allocate(AIH(NH))
Allocate(PIH(NH))

Call Linear &
  (HT(:),   & ! <-- Argument grid
   XI(:),   & ! <-- Gridded function
   HH(1),   & ! <-- Interpolation point
   XImin)     ! --> Interpolated function value

Call Linear &
  (HT(:),   & ! <-- Argument grid
   XI(:),   & ! <-- Gridded function
   HH(NH),  & ! <-- Interpolation point
   XImax)     ! --> Interpolated function value

If (XImax < XImin) then
   SW    = XImax
   XImax = XImin
   XImin = SW
End If

DXIH  = Min((XImax-XImin)/(NH-1), 2*Pi/(k(IC)*(TL-TU+0.02)))

DXIm  = ((XImax-XImin) - (NH-1)*DXIH)/2
XImax = XImax - DXIm
XImin = XImin + DXIm


Do i=1,NH
   XIH(i) = (XImin*(NH-i) + XImax*(i-1))/(NH-1)
   Call Linear &
     (XI(:),      & ! <-- Argument grid
      HT(:),      & ! <-- Gridded function
      XIH(i),     & ! <-- Interpolation point
      HXIH(i))      ! --> Interpolated function value
End Do

Call Interpolate_field &
  (HH,        & ! <-- First regular grid
   AH,        & ! <-- Amplitude on first grid
   PH,        & ! <-- Phase on first grid
   HXIH,      & ! <-- Second grid
   AIH,       & ! --> Amplitude on second grid
   PIH)         ! --> Phase on second grid


Deallocate(AH)
Deallocate(HH)
Deallocate(PH)
Deallocate(HXIH)
Deallocate(HT)


!--- 5.2. Time grid computation and shift

Allocate(TH(NH))

Do i=1,NH
   TH(i) = (i-1)*2*Pi/(k(IC)*NH*DXIH)
End Do

DTH   = 2*Pi/(k(IC)*NH*DXIH)

DTA   = -TU + (TH(NH) - (TL - TU))/2
TH(:) = TH(:) - DTA


!--- 5.3. Superimposing 1st phase model

Do i=1,NH
   Call Linear &
     (XI(:),   & ! <-- Argument grid
      FI(:),   & ! <-- Gridded function
      XIH(i),     & ! <-- Interpolation point
      FII)          ! --> Interpolated function value
   PIH(i) = PIH(i) - k(IC)*(FII + DTA*XIH(i))
End Do


Deallocate(FI)


!--- 5.4. Multiplication with measure

Do i=1,NH
   Call Linear &
     (XI(:),        & ! <-- Argument grid
      FXY(:),       & ! <-- Gridded function
      XIH(i),       & ! <-- Interpolation point
      NUI,          & ! --> Interpolated function value
      CExt = .True.)  ! <~~ Constant/linear extrapolation
   AIH(i) = AIH(i)/NUI
End Do


Deallocate(XIH)
Deallocate(FXY)
Deallocate(XI)


!--- 5.5. Inverse Fourier transform

UH(:) = AIH(:)*Exp(Ci*Modulo(PIH(:),2*Pi))

FSign = 1

Call FFT1 &
  (FSign,    & ! <-- Transform direction.
   UH(:))      ! <-> Data to be transformed.

UH(:) = UH(:)*Sqrt(k(IC)/(2*Pi))*DXIH


!--- 5.6. Computation of phase and amplitude

AIH(:) = Abs(UH(:))

Where (UH(:) /= 0.0_Double)
   PIH(:) = ATan2(AImag(UH(:)),Real(UH(:)))
Elsewhere
   PIH(:) = 0.0
End Where

Call Accumulate_Phase &
  (PIH(:),      & ! <-> Array of (accumulated) phase
   1)             ! <~~ Phase change direction


Deallocate(UH)


!--- 5.7. Interpolation to time grid

Call Interpolate_field &
  (TH,        & ! <-- First regular grid
   AIH,       & ! <-- Amplitude on first grid
   PIH,       & ! <-- Phase on first grid
   WTN,       & ! <-- Second grid
   WA(:,IC),  & ! --> Amplitude on second grid
   WP(:,IC))    ! --> Phase on second grid


Deallocate(AIH)
Deallocate(PIH)
Deallocate(TH)


!--- 5.8. Superimposing 2nd phase model

WP(:,IC) = WP(:,IC) + k(IC)*(GI(:) + XImin*WTN(:))

Deallocate(GI)


End Do Channels


!----------------------------------------------------------
! 6. MEMORY DEALLOCATION
!----------------------------------------------------------


Deallocate(TRN)
Deallocate(KR)
Deallocate(WTN)
Deallocate(DWX)
Deallocate(DWY)



End Subroutine Forward_Propagation_FIO



!==========================================================
Subroutine Interpolate_field &
  (Y1,        & ! <-- First regular grid
   A1,        & ! <-- Amplitude on first grid
   P1,        & ! <-- Phase on first grid
   Y2,        & ! <-- Second grid
   A2,        & ! --> Amplitude on second grid
   P2)          ! --> Phase on second grid
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
   A1(:)        ! Amplitude on first grid
!
Real(Double), Intent(In) :: &
   P1(:)        ! Phase on first grid
!
Real(Double), Intent(In) :: &
   Y2(:)        ! Second grid
!
! Output arguments:
!
Real(Double), Intent(Out)    :: &
   A2(:)        ! Amplitude on first grid
!
Real(Double), Intent(Out)    :: &
   P2(:)        ! Phase on first grid
!----------------------------------------------------------
! Local Scalars:
!
Integer      :: N1    ! Dimension of first grid
Integer      :: N2    ! Dimension of second grid
Integer      :: i1    ! Index of first grid
Integer      :: i2    ! Index of second grid
Real(Double) :: DY1   ! Step of first grid
!----------------------------------------------------------


!----------------------------------------------------------
! 1. INITIALIZATION
!----------------------------------------------------------

N1 = Size(Y1)
N2 = Size(Y2)

DY1 = Y1(2) - Y1(1)


!----------------------------------------------------------
! 2. INTERPOLATION
!----------------------------------------------------------

Do i2=1,N2
   i1  = 1 + Floor((Y2(i2) - Y1(1))/DY1)
   If ((i1 >= 1) .and. (i1 <= N1-1)) then
      A2(i2) = (A1(i1)*(Y1(i1+1) - Y2(i2)) + A1(i1+1)*(Y2(i2) - Y1(i1)))/DY1
      P2(i2) = (P1(i1)*(Y1(i1+1) - Y2(i2)) + P1(i1+1)*(Y2(i2) - Y1(i1)))/DY1
   Else
      A2(i2) = 0
      P2(i2) = 0
   End If
End Do


End Subroutine Interpolate_field



End Module Vacuum_Propagator


