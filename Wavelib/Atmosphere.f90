!
Module Atmosphere
!
! Providing a model of atmosphere using grib files,
! or radio sonde vertical profile, or an analytical
! model.
!----------------------------------------------------------
! (C) Copyright 2002, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 06 Aug 2002 | Original version.
!   2.0   | 20 Sep 2009 | Horizontal gradients mode.
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Parameters:
    Double, Pi, dtr
!
Use Errors, only: &
! Imported Type Definitions:
    Error_Status,  &
! Imported Routines:
    Enter_Callee,  &
    Exit_Callee,   &
    Error
!
Use Earth, only: &
! Imported Type Definitions:
    Geodetic,        &
! Imported Routines:
    Cart_from_Geod,  &
    Geod_from_Cart
!
Use Coordinates, only: &
! Imported Type Definitions:
    Cartesian,     &
! Imported Operators:
    Operator(-),   &
    Operator(+),   &
    Operator(*),   &
    Operator(/),   &
    Assignment(=)
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Public Parameters:
!
! --- Error codes
!
Integer, Parameter :: &
   err_Wrong_Type    = 1302001
!
! --- Atmosphere model types
!
Integer, Parameter :: &
   atm_Undefined  = 0,     & ! Undefined mode
   atm_Phantom    = 1,     & ! Phantom mode
   atm_ECHAM      = 2,     & ! ECHAM mode
   atm_NCEP       = 3,     & ! NCEP mode
   atm_MPAS       = 4        ! MPAS mode
!
Character(Len=6), Parameter :: &
   Atm_Name(0:4) = &
      (/ 'undef ',      &
         'phan  ',      &
         'echam ',      &
         'ncep  ',      &
         'mpas  '  /)
!
! --- Horizontal gradients modes
!
Integer, Parameter :: &
   hg_Undefined   = 0,     & ! Undefined mode
   hg_1d          = 1,     & ! No horizontal gradients
   hg_3d          = 2        ! Full horizontal gradients
!
Character(Len=6), Parameter :: &
   HG_Name(0:2) = &
      (/ 'undef ',      &
         '1d    ',      &
         '3d    ' /)
!
! Public Scalars:
!
Integer  :: &
   Atmosphere_Type = atm_Undefined,    & ! Atmosphere model type
   HorizGrad_Mode  = hg_3d               ! Horizontal 
!
Type(Geodetic) :: &
   GP                ! Coordinates of center for HG_Mode=1d
!----------------------------------------------------------
!
Contains


!==========================================================
Subroutine Atmosphere_Init &
  (Atm_Type,      & ! <-- Type of atmosphere model
   Atm_Name,      & ! <-- Pathname(s) of atmosphere model file(s)
   HG_Mode,       & ! <-- Horizontal gradients mode
   Vrb,           & ! <-- Verbosity level
   ErrStat,       & ! <-> Error status
   XGP,           & ! <~~ Coordinates of center for HG_Mode=1d
   XFreq)           ! <~~ Frequency channels [Hz]
!
! Intialization of atmosphere model.
!----------------------------------------------------------
! Method:
!   Invoke of ECHAM_Init, NCEP_Init, or Phantom_Init.
!----------------------------------------------------------
! (C) Copyright 2002-2009, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 06 Aug 2002 | Original version.
!   2.0   | 28 Oct 2003 | Vrb is mandatory.
!   2.1   | 29 Jun 2007 | NCEP_Names.
!   2.2   | 05 May 2008 | Vrb passed to ECHAM/NCEP_Init.
!   3.0   | 20 Sep 2009 | Atm_Type.
!----------------------------------------------------------
! Modules used:
!
Use ECHAM_fields, only: &
! Imported Routines:
    ECHAM_Init
!
Use NCEP_fields, only: &
! Imported Routines:
    NCEP_Init
!
Use MPAS, only: &
! Imported Routines:
    MPAS_Init
!
Use Phantom, only: &
! Imported Routines:
    Phantom_Init
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Integer, Intent(In) :: &
   Atm_Type         ! Type of atmosphere model
   
Character(Len=*), Intent(In) :: &
   Atm_Name(:)       ! Pathname(s) of atmosphere model file(s)
!
Integer, Intent(In) :: &
   HG_Mode           ! Horizontal gradients mode
!
Integer, Intent(In)  :: &
   Vrb               ! Verbosity level
!
! InOut arguments:
!
Type(Error_Status), Pointer   :: &
   ErrStat           ! Error status
!
! Optional input arguments:
!
Type(Geodetic), Intent(In), Optional :: &
   XGP               ! Coordinates of center for HG_Mode=1d
!
Real(Double), Intent(In), Optional :: &
   XFreq(:)          ! Frequency channels [Hz]
!----------------------------------------------------------
! Global variables used:
!
!   Atmosphere_Type
!----------------------------------------------------------


!----------------------------------------------------------
! 0. INITIALIZATION
!----------------------------------------------------------


!--- 0.1. Initialization of error status

Call Enter_Callee &
  ('Atmosphere_Init',  & ! <-- User routine
   ErrStat)              ! <-> Pointer to callee status


!--- 0.2. Storing horizontal gradient mode parameters

If (Present(XGP)) then
   HorizGrad_Mode = HG_Mode
   GP             = XGP
Else
   HorizGrad_Mode = hg_3d
   GP             = Geodetic(0.0_Double, 0.0_Double, 0.0_Double)
End If


!----------------------------------------------------------
! 1. SETTING ATMOSPHERE MODEL
!----------------------------------------------------------


Select Case (Atm_Type)

   Case (atm_ECHAM)

      Atmosphere_Type = atm_ECHAM

      Call ECHAM_Init &
        (Atm_Name,      & ! <-- Pathnames of data files
         Vrb,           & ! <-- Verbosity level
         ErrStat,       & ! <-> Error status
         XFreq)           ! <~~ Frequency channels [Hz]

   Case (atm_NCEP)

      Atmosphere_Type = atm_NCEP

      Call NCEP_Init &
        (Atm_Name,      & ! <-- Pathnames of data files
         Vrb,           & ! <-- Verbosity level
         ErrStat,       & ! <-> Error status
         XFreq)           ! <~~ Frequency channels [Hz]

   Case (atm_MPAS)

      Atmosphere_Type = atm_MPAS

      Call MPAS_Init &
        (Atm_Name,      & ! <-- Pathnames of data files
         Vrb,           & ! <-- Verbosity level
         ErrStat,       & ! <-> Error status
         XFreq)           ! <~~ Frequency channels [Hz]

   Case (atm_Phantom)

      Atmosphere_Type = atm_Phantom

      Call Phantom_Init &
        (Atm_Name(1),   & ! <-- Phantom parameter file
         Vrb,           & ! <-- Verbosity level
         ErrStat,       & ! <-> Error status
         XFreq)           ! <~~ Frequency channels [Hz]

   Case Default

      Atmosphere_Type = atm_Undefined

      Call Exit_Callee &
        (ErrStat,           & ! <-> Pointer to callee status
         err_Wrong_Type,    & ! <~~ User error code
         0,                 & ! <~~ System error code
         'Wrong atm_Type')    ! <~~ Error description
      Return

End Select


!----------------------------------------------------------
! 2. EXIT
!----------------------------------------------------------


Call Exit_Callee(ErrStat)


End Subroutine Atmosphere_Init



!==========================================================
Subroutine Atmosphere_Constituents &
  (PLon,     & ! <-- Longiude of point [deg]
   PLat,     & ! <-- Latitude of point [deg]
   ZP,       & ! <-- Altitude of point [km]
   ZPmin,    & ! --> Minimum model Z for this lon/lat
   ZPmax,    & ! --> Maximum model Z for this lon/lat
   QP,       & ! --> Interpolated humidity [kg/kg]
   TP,       & ! --> Interpolated temperature [K]
   RNR,      & ! --> Interpolated Re(N)
   RNI)        ! ~~> Interpolated Im(N)
!
! Computation of interpolated humidity, temperature,
! and refractivity for given point.
!----------------------------------------------------------
! Method:
!   Invoke of ECHAM_Constituents or Phantom_Constituents.
!----------------------------------------------------------
! (C) Copyright 2002-2007, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 06 Aug 2002 | Original version.
!   1.1   | 29 Jun 2007 | Added NCEP.
!----------------------------------------------------------
! Modules used:
!
Use ECHAM_fields, only: &
! Imported Routines:
    ECHAM_Constituents
!
Use NCEP_fields, only: &
! Imported Routines:
    NCEP_Constituents

Use MPAS, only: &
! Imported Routines:
    MPAS_Constituents
!
Use Phantom, only: &
! Imported Routines:
    Phantom_Constituents
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In) :: &
   PLon          ! Latitude [deg]
!
Real(Double), Intent(In) :: &
   PLat          ! Longitude [deg]
!
Real(Double), Intent(In) :: &
   ZP            ! Altitude [km]
!
! Output arguments:
!
Real(Double), Intent(Out) :: &
   ZPmin         ! Minimum model Z for this lon/lat
!
Real(Double), Intent(Out) :: &
   ZPmax         ! Maximum model Z for this lon/lat
!
Real(Double), Intent(Out) :: &
   QP            ! Interpolated humidity [kg/kg]
!
Real(Double), Intent(Out) :: &
   TP            ! Interpolated temperature [K]
!
Real(Double), Intent(Out) :: &
   RNR           ! Interpolated Re(N)
!
Real(Double), Optional, Intent(Out) :: &
   RNI(:)        ! Interpolated Im(N)
!----------------------------------------------------------
! Global variables used:
!
!   Atmosphere_Type
!----------------------------------------------------------
! Local Scalars:
!
! --- Coordinates
!
Real(Double)      :: PGLon    ! Longitude [deg]
Real(Double)      :: PGLat    ! Latitude [deg]
!----------------------------------------------------------


!----------------------------------------------------------
! 1. COMPUTATION OF INTERPOLATION POINT COORDINATES
!----------------------------------------------------------

Select Case (HorizGrad_Mode)
   Case (hg_1d)
      PGLon = GP%Lambda
      PGLat = GP%Phi
   Case Default
      PGLon = PLon
      PGLat = PLat
End Select


!----------------------------------------------------------
! 2. COMPUTATION OF CONSTITUENTS
!----------------------------------------------------------


Select Case (Atmosphere_Type)

   Case (atm_ECHAM)

      Call ECHAM_Constituents &
        (PGLon,    & ! <-- Longiude of point [deg]
         PGLat,    & ! <-- Latitude of point [deg]
         ZP,       & ! <-- Altitude of point [km]
         ZPmin,    & ! --> Minimum model Z for this lon/lat
         ZPmax,    & ! --> Maximum model Z for this lon/lat
         QP,       & ! --> Interpolated humidity [kg/kg]
         TP,       & ! --> Interpolated temperature [K]
         RNR,      & ! --> Interpolated Re(N)
         RNI)        ! ~~> Interpolated Im(N)

   Case (atm_NCEP)

      Call NCEP_Constituents &
        (PGLon,    & ! <-- Longiude of point [deg]
         PGLat,    & ! <-- Latitude of point [deg]
         ZP,       & ! <-- Altitude of point [km]
         ZPmin,    & ! --> Minimum model Z for this lon/lat
         ZPmax,    & ! --> Maximum model Z for this lon/lat
         QP,       & ! --> Interpolated humidity [kg/kg]
         TP,       & ! --> Interpolated temperature [K]
         RNR,      & ! --> Interpolated Re(N)
         RNI)        ! ~~> Interpolated Im(N)

   Case (atm_MPAS)

      Call MPAS_Constituents &
        (PGLon,    & ! <-- Longiude of point [deg]
         PGLat,    & ! <-- Latitude of point [deg]
         ZP,       & ! <-- Altitude of point [km]
         ZPmin,    & ! --> Minimum model Z for this lon/lat
         ZPmax,    & ! --> Maximum model Z for this lon/lat
         QP,       & ! --> Interpolated humidity [kg/kg]
         TP,       & ! --> Interpolated temperature [K]
         RNR,      & ! --> Interpolated Re(N)
         RNI)        ! ~~> Interpolated Im(N)

   Case (atm_Phantom)

      Call Phantom_Constituents &
        (PGLon,    & ! <-- Longiude of point [deg]
         PGLat,    & ! <-- Latitude of point [deg]
         ZP,       & ! <-- Altitude of point [km]
         ZPmin,    & ! --> Minimum model Z for this lon/lat
         ZPmax,    & ! --> Maximum model Z for this lon/lat
         QP,       & ! --> Interpolated humidity [kg/kg]
         TP,       & ! --> Interpolated temperature [K]
         RNR,      & ! --> Interpolated Re(N)
         RNI)        ! ~~> Interpolated Im(N)

   Case Default

      Write(*,'(A,I10)') 'Wrong Atmosphere_Type in Atmosphere_Constituents:', Atmosphere_Type
      Stop

End Select


End Subroutine Atmosphere_Constituents



!==========================================================
Subroutine Atmosphere_NGradN &
  (X,        & ! <-- Cartesian coordinates of point
   NGradN,   & ! --> Interpolated (1 + N)*Grad(N)
   NP,       & ! --> Interpolated N
   Stat,     & ! --> Error status
   NIP)        ! ~~> Interpolated Im(N)
!
! Calculation of (1 + N)*Grad(N) for atmosphere model.
!----------------------------------------------------------
! Method:
!   Invoke of ECHAM_NGrad or NPhantom_NGradN
!----------------------------------------------------------
! (C) Copyright 2002-2007, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 06 Aug 2002 | Original version.
!   2.0   | 30 Apr 2004 | Added turbulence.
!   2.1   | 29 Jun 2007 | Added NCEP.
!----------------------------------------------------------
! Modules used:
!
Use ECHAM_fields, only: &
! Imported Routines:
    ECHAM_NGradN
!
Use NCEP_fields, only: &
! Imported Routines:
    NCEP_NGradN
!
!Use MPAS_fields, only: &
Use MPAS, only: &
! Imported Routines:
    MPAS_NGradN
!
Use Phantom, only: &
! Imported Routines:
    Phantom_NGradN
!
Use Turbulence, only: &
! Imported Routines:
    Turbulence_FGradF
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type(Cartesian), Intent(In) :: &
   X        ! Cartesian coordinates of point
!
! Output arguments:
!
Type(Cartesian), Intent(Out) :: &
   NGradN   ! Interpolated (1 + N)*Grad(N)
!
Real(Double), Intent(Out) :: &
   NP       ! Interpolated N
!
Integer, Intent(Out)         :: &
   Stat     ! Error status:
            !   0 - point above surface
            !   1 - point under surface
!
Real(Double), Optional, Intent(Out) :: &
   NIP(:)   ! Interpolated Im(N)
!----------------------------------------------------------
! Global variables used:
!
!   Atmosphere_Type
!----------------------------------------------------------
! Local Scalars:
!
! --- Coordinates
!
Type(Geodetic)    :: G        ! Geodetic coordintes
Type(Cartesian)   :: XG       ! Cartesian coordinates
!
! --- Regular refractivity
!
Real(Double)      :: NPR      ! Interpolated N
Type(Cartesian)   :: NGradNR  ! Interpolated (1 + N)*Grad(N)
!
! --- Turbulence
!
Real(Double)      :: FP       ! Interpolated F
Type(Cartesian)   :: GradF    ! Interpolated Grad(F)
!----------------------------------------------------------


!----------------------------------------------------------
! 1. COMPUTATION OF INTERPOLATION POINT COORDINATES
!----------------------------------------------------------

Select Case (HorizGrad_Mode)
   Case (hg_1d)
      G        = Geod_from_Cart(X)
      G%Phi    = GP%Phi
      G%Lambda = GP%Lambda
      XG       = Cart_from_Geod(G)
   Case Default
      XG       = X
End Select


!----------------------------------------------------------
! 2. COMPUTATION OF REGULAR REFRACTIVITY AND ITS GRADIENT
!----------------------------------------------------------

Select Case (Atmosphere_Type)

   Case (atm_ECHAM)

      Call ECHAM_NGradN &
        (XG,       & ! <-- Cartesian coordinates of point
         NGradNR,  & ! --> Interpolated (1 + N)*Grad(N)
         NPR,      & ! --> Interpolated N
         Stat,     & ! --> Error status
         NIP)        ! ~~> Interpolated Im(N)

   Case (atm_NCEP)

      Call NCEP_NGradN &
        (XG,       & ! <-- Cartesian coordinates of point
         NGradNR,  & ! --> Interpolated (1 + N)*Grad(N)
         NPR,      & ! --> Interpolated N
         Stat,     & ! --> Error status
         NIP)        ! ~~> Interpolated Im(N)

   Case (atm_MPAS)
!      print*,'Interpolating from MPAS - Need MPAS_NGradN stopping'
!      stop

      Call MPAS_NGradN &
        (XG,       & ! <-- Cartesian coordinates of point
         NGradNR,  & ! --> Interpolated (1 + N)*Grad(N)
         NPR,      & ! --> Interpolated N
         Stat,     & ! --> Error status
         NIP)        ! ~~> Interpolated Im(N)

!      print*,'MPAS - After MPAS_NGradN stopping'
!      print*,'MPAS - After MPAS_NGradN continuing'
!      stop

   Case (atm_Phantom)

      Call Phantom_NGradN &
        (XG,       & ! <-- Cartesian coordinates of point
         NGradNR,  & ! --> Interpolated (1 + N)*Grad(N)
         NPR,      & ! --> Interpolated N
         Stat,     & ! --> Error status
         NIP)        ! ~~> Interpolated Im(N)

   Case Default

      Write(*,'(A,I10)') 'Wrong Atmosphere_Type in Atmosphere_NGradN:', Atmosphere_Type
      Stop

End Select


!----------------------------------------------------------
! 3. COMPUTATION OF TURBULENCE PERTURBATION
!----------------------------------------------------------

Call Turbulence_FGradF &
  (XG,       & ! <-- Cartesian coordinates of point
   GradF,    & ! --> Interpolated Grad(F)
   FP)         ! --> Interpolated F


!----------------------------------------------------------
! 4. COMPUTATION OF COMBINED REFRACTIVITY AND GRADIENT
!----------------------------------------------------------

NP = NPR*(1 + FP)

NGradN = NGradNR*(1 + FP)   + &
         (1 + NP)*NPR*GradF + &
         (NPR*NGradNR/(1+NPR))*FP*(1 + FP)


End Subroutine Atmosphere_NGradN



!==========================================================
Subroutine Atmosphere_NGHN &
  (X,        & ! <-- Cartesian coordinates of point
   NGN,      & ! --> Interpolated (1 + N)*Grad(N)
   NP,       & ! --> Interpolated N
   NHN,      & ! --> Interpolated Grad x (1+N)Grad(N)
   Stat)       ! --> Error status
!
! Calculation of interpolated (1 + N)*Grad(N) for
! atmosphere model: Adjoint version.
!----------------------------------------------------------
! Method:
!   Invoke of ECHAM_NGHN or NPhantom_NGHN
!----------------------------------------------------------
! (C) Copyright 1999-2007, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 07 May 2003 | Original version.
!   1.1   | 29 Jun 2007 | Added NCEP.
!----------------------------------------------------------
! Modules used:
!
Use ECHAM_fields, only: &
! Imported Routines:
    ECHAM_NGHN
!
Use NCEP_fields, only: &
! Imported Routines:
    NCEP_NGHN
!
Use Phantom, only: &
! Imported Routines:
    Phantom_NGHN
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type(Cartesian), Intent(In)  :: &
   X               ! Cartesian coordinates of point
!
! Output arguments:
!
Type(Cartesian), Intent(Out) :: &
   NGN             ! Interpolated (1 + N)*Grad(N)
!
Real(Double), Intent(Out)    :: &
   NP              ! Interpolated N
!
Real(Double), Intent(Out)    :: &
   NHN(3,3)        ! Interpolated Grad x (1+N)Grad(N)
!
Integer, Intent(Out)         :: &
   Stat            ! Error status:
                   !   0 - point above surface
                   !   1 - point under surface
!----------------------------------------------------------
! Global variables used:
!
!   Atmosphere_Type
!----------------------------------------------------------
! Local Scalars:
!
! --- Coordinates
!
Type(Geodetic)    :: G        ! Geodetic coordintes
Type(Cartesian)   :: XG       ! Cartesian coordinates
!----------------------------------------------------------


!----------------------------------------------------------
! 1. COMPUTATION OF INTERPOLATION POINT COORDINATES
!----------------------------------------------------------

Select Case (HorizGrad_Mode)
   Case (hg_1d)
      G        = Geod_from_Cart(X)
      G%Phi    = GP%Phi
      G%Lambda = GP%Lambda
      XG       = Cart_from_Geod(G)
   Case Default
      XG       = X
End Select


!----------------------------------------------------------
! 2. COMPUTATION OF REFRACTIVITY AND ITS GRADIENT
!----------------------------------------------------------


Select Case (Atmosphere_Type)

   Case (atm_ECHAM)

      Call ECHAM_NGHN &
        (XG,       & ! <-- Cartesian coordinates of point
         NGN,      & ! --> Interpolated (1 + N)*Grad(N)
         NP,       & ! --> Interpolated N
         NHN,      & ! --> Interpolated Grad x (1+N)Grad(N)
         Stat)       ! --> Error status

   Case (atm_NCEP)

      Call NCEP_NGHN &
        (XG,       & ! <-- Cartesian coordinates of point
         NGN,      & ! --> Interpolated (1 + N)*Grad(N)
         NP,       & ! --> Interpolated N
         NHN,      & ! --> Interpolated Grad x (1+N)Grad(N)
         Stat)       ! --> Error status

   Case (atm_Phantom)

      Call Phantom_NGHN &
        (XG,       & ! <-- Cartesian coordinates of point
         NGN,      & ! --> Interpolated (1 + N)*Grad(N)
         NP,       & ! --> Interpolated N
         NHN,      & ! --> Interpolated Grad x (1+N)Grad(N)
         Stat)       ! --> Error status

   Case Default

      Write(*,'(A,I10)') 'Wrong Atmosphere_Type in Atmosphere_NGHN:', Atmosphere_Type
      Stop

End Select


End Subroutine Atmosphere_NGHN



End Module Atmosphere


