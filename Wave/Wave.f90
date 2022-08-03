!
Program Wave
!
! Wave optics propagation simulator
!----------------------------------------------------------
! Method:
!   1. Multiple phase screens in 2D + Diffractive integrals or
!      Fourier Integral Operator
!   2. Asymptotic solutions.
!----------------------------------------------------------
! Command line parameters:
!   -------------------------------------------------
!   Parameters             | Descriptions
!   -------------------------------------------------
!   -atm=type,[hg],ss,...  | Atmospheric model specification:
!                          | type:
!                          |    echam   - ECHAM girb file
!                          |    ncep    - NCEP grib file
!                          |    phantom - analytical model
!                          | hg - horizontal gradients:
!                          |    1d      - no horizontal gradients
!                          |    3d      - full horizontal gradients (default)
!                          | ss,ss,...  - pathname(s)
!   -turb=sssss            | Turbulence parameter file
!   -earth=sss[,ppp]       | sss specifies Earth shape:
!                          |    sphere - spherical
!                          |    ellips - reference elliposid WGS 1984
!                          |    geoid  - geoid
!                          |       ppp - optional path to coefficient file
!   -gps=ssss              | RO data file patname
!   -orb=ssss              | Orbit data file pathname
!   -msis=ssss             | Path to MSIS coefficient files
!   -out=ssss              | Template for output names
!                          |    by default GPS_Name or ORB_Name
!   -freq=ffff             | Frequency-channels specification:
!         gps              |    GPS frequncies
!         gal              |    GALILEO frequencies
!         glo              |    GLONASS frequencies
!         acestd           |    ACE+ standard frequencies
!         aceopt           |    ACE+ optional frequencies
!         acemax           |    ACE+ maximum frequencies
!         xxx,xxx...       |    user-specified frequencies [GHz]
!   -afreq=ffff            | Frequency-channels for absorption
!   -hmax=xxxx             | maximum height for WO (default 100)
!   -hmin=xxxx             | minimum height for WO (default -999)
!   -hgo=xxxx              | maximum height for GO (default 120)
!   -dyn=xxxx              | vertical resolution of refractivity field
!   -dx=xxxx               | step between phase screens (default 100)
!   -xls=xxxx              | last phase screen position
!   -fz=xxxx               | Fresnel zone size [Pi rad] (default 10)
!   -sr=xxxx               | sampling rate [Hz] (default 50)
!   -ts=y/n                | time scaling (default - deactivated)
!   -ecef=y/n              | all coordinates in ECEF (default J2000)
!   -wo=y/n                | wave optics (default yes)
!   -go=y/n                | geometrical optics (default yes)
!   -rh=y/n                | output of ray heights (default yes)
!   -ab=y/n                | modeling of absorption (default no)
!   -as=y/n                | asymptotic solution (default no)
!   -fio=y/n               | use FIO propagation to LEO orbit (default yes)
!   -sd=xxxx               | spatial dimension (default 2)
!   -ngo=nnnn              | number of GO rays (default 1500)
!   -fwp=xxxx              | filter width for computing DPGPS/DP [km] (default -1)
!   -v=xxxx                | verbosity level:
!                          |   -1 - silent mode
!                          |    0 - version and error info
!                          |    1 - progress messages
!                          |    2 - static messages
!                          |    3 - dynamic messages
!   -b                     | batch mode (equivalent to -v=1)
!   -dbg=xxxx              | debug mode
!   -o=                    | output files extensions:
!    XX[-],...             |    XX     - turnes specific extension on
!                          |    XX-    - turnes specific extension off
!                          |    none   - no output
!                          |    all    - output of all extensions
!                          |    min    - minimum ouput data set
!                          |    script - output data set for use
!                          |             with genvda script
!                          |    Extensions:
!                          |    WOP WCA WCE CAU GOP GRP XBP GOZ
!                          |    GON GOT MOT MOQ MON MOA MTR MTD
!                          |    NC
!   -----------------------------------
!
! Input files:
!   1. RO data/orbit description file
!   2. ECHAM/NCEP global fields in GRIB format
!
! Output files: <OutName>.ext
!   -----------------------------------
!   Extensions | Data file descriptions
!   -----------+----------------------- 
!   <wop>      | Simulated occultation data file
!   wca<c>     | CT amplitude for channes <c>
!   wce<c>     | CT BA for channel <c>
!   cau        | Caustic structure [km]
!   gop        | GO BA profile [rad, km]
!   grp        | GO reflected BA profile [rad, km]
!   xbp        | Back propagation plane position [km]
!   goz        | GO BA vs perigee height [rad, km]
!   gon        | Refractivity profile from GO BA [N-units, km]
!   got        | Dry temperature profile from GO BA [K, km]
!   mot        | Local profile of temperature [K, km]
!   moq        | Local profile of humidity [g/kg, km]
!   mon        | Local profile of refractivity [N-units]
!   moa<c>     | Local profile of specific absorption for channel <c>
!   moa<c1c2>  | Differential specific absorption for channels <c1>,<c2>
!   mtr<c>     | Transmission for channel <c>
!   mtr<c1c2>  | Differential transmission for channels <c1>,<c2>
!   mtd        | Local dry temperature [K]
!   nc         | Occultation data in atmPhs-NetCDF format
!----------------------------------------------------------
! (C) Copyright 1999-2009, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 27 Feb 1999 | Original version.
!   1.1   | 07 Mar 1999 | Hmax, better apodization.
!   1.2   | 14 Mar 1999 | Dynamic determination of
!         |             | number of output data.
!   2.0   | 20 Mar 1999 | Error check.
!   2.1   | 23 Mar 1999 | Output of interpolated temperature.
!   3.0   | 21 May 1999 | Output of caustics.
!   3.1   | 08 Jul 1999 | -b.
!   3.2   | 16 Jul 1999 | Determination of XBP.
!   3.3   | 23 Jul 1999 | Data check.
!   3.4   | 03 Oct 1999 | Output in compact format.
!   3.5   | 06 Oct 1999 | -rh.
!   3.6   | 25 Nov 2000 | Output of GO refractivity.
!   3.7   | 04 Dec 2000 | Output of interpolated humidity
!         |             | refractivity, and dry temperature.
!   3.8   | 17 Dec 2000 | -fz.
!   3.9   | 21 Mar 2001 | -nw, S=DX/2, -sr.
!   4.0   | 15 May 2001 | Absorption.
!   4.1   | 16 May 2001 | Asymptotic solution.
!   4.2   | 18 May 2001 | Output of specific absorption.
!   4.3   | 13 Oct 2001 | -mps.
!   4.4   | 17 Nov 2001 | Output of eps(z).
!   4.5   | 12 Dec 2001 | Occultation plane rotated with Earth.
!   4.6   | 27 Jan 2002 | -ngo.
!   5.0   | 09 Aug 2002 | Merged with WaveL.
!   5.1   | 14 Aug 2002 | Additional phase screen, -xls.
!   5.2   | 21 Aug 2002 | CT in last phase screen.
!   5.3   | 11 Nov 2002 | -orb.
!   6.0   | 23 May 2003 | Improved asymptotic solution.
!   6.1   | 06 Jul 2003 | -sd.
!   6.2   | 20 Aug 2003 | Verbosity level -v.
!   6.3   | 30 Oct 2003 | Dynamic integration step for -as -nw.
!   6.4   | 17 Dec 2003 | Output of log attenuation.
!   6.5   | 24 Apr 2004 | FIO propagation to LEO orbit
!         |             | in wave optics propagator.
!   6.6   | 28 Apr 2004 | Modeling turbulence.
!   6.7   | 20 May 2004 | -afreq.
!   6.8   | 25 May 2004 | Output of differential transmission;
!         |             | -dbg.
!   6.9   | 16 Sep 2004 | Pointer nullifying.
!   7.0   | 07 Dec 2005 | Modeling rising and setting occultations.
!   7.1   | 27 Feb 2006 | Corrected output of goz.
!   7.2   | 02 Apr 2006 | Longer output format for Y.
!   7.3   | 05 May 2006 | -mpsf, -ecef, -ow.
!   7.4   | 27 Nov 2006 | -out.
!   7.5   | 23 Feb 2007 | -ts.
!   7.6   | 24 Feb 2007 | -fwp.
!   7.7   | 29 Jun 2007 | -ncep.
!   7.8   | 01 Jul 2007 | Reading generic RO file format.
!   8.0   | 25 May 2008 | -o.
!   8.1   | 28 May 2008 | -v=-1 total silence.
!   8.2   | 16 Jun 2008 | -key=y/n, -msis.
!   9.0   | 24 Nov 2008 | -geoid.
!   9.1   | 27 Nov 2008 | GO inversion heights above MSL.
!   9.2   | 04 Mar 2009 | Local profile output if no simulation.
!   9.3   | 09 Mar 2009 | Occultation point definition from
!         |             | straigh-line ray touching Earth.
!   9.4   | 19 Mar 2009 | Support of GRAS NetCDF data.
!   9.5   | 26 Mar 2009 | Skew local profiles.
!   9.6   | 17 Jun 2009 | Bugs corrected.
!   9.7   | 22 Jun 2009 | Memory deallocation corrected.
!   9.8   | 02 Sep 2009 | Corrected bug: Freq used instead of AFred
!         |             | in invoke of GO_Invert and Local_Profiles.
!   9.9   | 07 Sep 2009 | Output of differential specific absorption. 
!  10.0   | 09 Sep 2009 | -freq=optmax.
!  10.1   | 14 Sep 2009 | Corrected processing -ecef=y
!         |             | for orbit data from occultation file.
!  10.2   | 15 Sep 2009 | -hmin=xxxx.
!  10.3   | 20 Sep 2009 | -atm=type,[hg],ss,...
!  10.4   | 23 Sep 2009 | Optimized code.
!  10.5   | 01 Oct 2009 | -o=NC.
!  10.6   | 20 Oct 2009 | Modified for naming GRAS data.
!  10.6   | 30 Oct 2009 | Invoke of Make_WOP_Name.
!  10.7   | 05 Nov 2009 | -hgo.
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Parameters:
    Double, Pi
!
Use Time, only: &
! Imported Routines:
    Elapsed_Time,    &
    Day_of_Year
!
Use Wave_Comline, only: &
! Imported Types:
    t_Opts,          &
! Imported Routines:
    Get_Comline,     &
    Help
!
Use Strings, only: &
! Imported Routines:
    SL => String_from_Log
!
Use Occ_GNSS, only: &
! Imported Parameters:
    Freq_GPS
!
Use Wave_Version, only: &
! Imported Parameters:
    Display_Version
!
Use Earth, only: &
! Imported Type Definitions:
    Geodetic,           &
! Imported Parameters:
    Shape_Name,         &
    shape_Sphere,       &
    shape_Ellips,       &
    shape_Geoid,        &
! Imported Routines:
    Earth_Init
!
Use Coordinates, only: &
! Imported Type Definitions:
    Cartesian
!
Use Occ_IO, only: &
! Imported Parameters:
    dt_RINEX,              &
    dt_GPSMET,             &
    dt_LEOLEO,             &
    dt_ATMPHS,             &
    dt_Name,               &
! Imported Types:
    t_MetaData,            &
! Imported Routines:
    Read_RO,               &
    Write_RO,              &
    Check_RO,              &
    Write_AtmPhs
!
Use Occ_Orbits, only: &
! Imported Routines:
    Generate_Orbit
!
Use Occ_Coordinates, only: &
! Imported Routines:
    Occ_Geometry_ECEF,     &
    Occ_Geometry,          &
    Plane_Basis
!
Use Geoid, only: &
! Imported Parameters:
    geoid_asc,        &
! Imported Routines:
    Set_Geoid_Path,   &
    Geoid_Init
!
Use MSIS, only: &
! Imported Routines:
    Set_MSIS_Path
!
Use Occ_Strings, only: &
! Imported Routines:
    Make_WOP_Name
!
Use Atmosphere, only: &
! Imported Parameters:
    atm_Undefined,   &
    atm_Phantom,     &
    atm_ECHAM,       &
    atm_NCEP,        &
    atm_MPAS,        &
    Atm_Name,        &
    hg_Undefined,    &
    hg_1d,           &
    hg_3d,           &
    HG_Name,         &
! Imported Routines:
    Atmosphere_Init
!
Use Turbulence, only: &
! Imported Routines:
    Turbulence_Init
!
Use Wave_Propagator, only: &
! Imported Routines:
    Wave_Propagate
!
Use Asymptotic_Propagator, only: &
! Imported Routines:
    Asymptotic_Propagate
!
Use Occ_Refraction, only: &
! Imported Routines:
    GO_Caustics
!
Use Interpolation, only: &
! Imported Array Variables:
    Linear
!
Use GO_Propagator, only: &
! Imported Routines:
    GO_Propagate,   &
    GO_Invert,      &
    Local_Profiles
!
Use Occ_PhaseModel, only: &
! Imported Routines:
    MSIS_Phase_Excess
!
Use IO, only: &
! Imported Routines:
    FileOpen,  &
    PutXY
!
Use Errors, only: &
! Imported Type Definitions:
    Error_Status,  &
! Imported Routines:
    Enter_Callee,  &
    Error,         &
    Display_Status
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Local Parameters:
!
Integer, Parameter :: &
   NGR = 20       ! Number of GO reflected rays
!
! Local Scalars:
!
! --- Error status
!
Type(Error_Status), Pointer :: &
                        ErrStat     ! Error status
!
! --- Command line parameters
!
Type(t_Opts)         :: Opts        ! Command line options
Logical              :: Err         ! Command line error status
!
! --- Input RO data
!
Type(t_MetaData)     :: MetaData    ! Occultation metadata
Integer              :: ID          ! Occultation ID
Integer              :: Year        ! Occultation year
Integer              :: Month       ! Occultation month
Integer              :: Day         ! Occultation day
Integer              :: Hour        ! Occultation hour
Integer              :: Minute      ! Occultation begin minute
Real(Double)         :: Second      ! Occultation begin second
Integer              :: N           ! Number of data read
Integer              :: Stat        ! Error code:
                                    !    0 - no error
                                    !   >0 - Open, Read error
                                    !  -10 - Arrays length not enough for data
Integer              :: DataType    ! RO data type
!
! --- GO refraction
!
Real(Double)         :: XBP         ! BP plane position
!
! --- Occultation geometry
!
Type(Geodetic)       :: GP          ! Occultation point
Type(Cartesian)      :: ERLC        ! Local curvature center (ECEF)
Type(Cartesian)      :: RLC         ! Curvature center (J2000)
Real(Double)         :: RE          ! Local curvature radius
Type(Cartesian)      :: AX          ! X-vector of plane basis
Type(Cartesian)      :: AY          ! Y-vector of plane basis
Integer              :: I1, I2      ! Phase model limits
!
! --- Work variables
!
Integer              :: NC          ! Number of channels
Integer              :: IC          ! Channel index
Integer              :: IC1, IC2    ! Channel indices for differential transmission
Integer              :: i           ! Work index
Integer              :: IU          ! Unit number
Integer              :: IUO, NNN, KKN  ! Unit number
Integer              :: IOStat
Logical              :: GPS_Freq    ! Frequencies = GPS frequencies
Integer              :: Data_Type   ! Output data type
Character(Len=255)   :: WOP_Name    ! Simulated GPS data file pathname
Character(Len=255)   :: WCA_Name    ! WOP CT amplitude file pathname
Character(Len=255)   :: WCE_Name    ! WOP CT refraction angle data file pathname
Character(Len=255)   :: GOP_Name    ! GO refraction file pathname
Character(Len=255)   :: GRP_Name    ! GO reflection file pathname
Character(Len=255)   :: GOZ_Name    ! GO refraction (z) file pathname
Character(Len=255)   :: GON_Name    ! Refractivity file pathname
Character(Len=255)   :: GOT_Name    ! Dry temperature file pathname
Character(Len=255)   :: MOT_Name    ! Model temperature file pathname
Character(Len=255)   :: MOQ_Name    ! Model humidity file pathname
Character(Len=255)   :: MON_Name    ! Model refractivity file pathname
Character(Len=255)   :: MTD_Name    ! Model dry temperature file pathname
Character(Len=255)   :: MOA_Name    ! Model specific absorption file pathname
Character(Len=255)   :: LAT_Name    ! Model log attenuation file pathname
Character(Len=255)   :: CAU_Name    ! Caustic file pathname
Character(Len=255)   :: XBP_Name    ! XBP file pathname
Character(Len=255)   :: mpaspath    ! MPAS file pathname
!
! Local Arrays:
!
! --- Input RO data
!
Real(Double), Pointer :: &
   TR(:),     & ! Relative time of samples [sec]
   A(:,:),    & ! Amplitudes (time,channel)
   S(:,:)       ! Phase excess (time,channel) [m]
Type(Cartesian), Pointer :: &
   RLEO(:),   & ! LEO coordinates (J2000)
   VLEO(:),   & ! LEO velocity (J2000)
   RGPS(:),   & ! GPS coordinates (J2000)
   VGPS(:)      ! GPS velocity (J2000)
Integer, Pointer      :: &
   LCF(:)       ! Lost carrier flag
Real(Double), Pointer ::&
   ROFreq(:)    ! Frequencies [Hz]
!
! --- Orbit data in ECEF
!
Type(Cartesian), Allocatable :: &
   ERLEO(:),     & ! LEO coordinates (ECEF)
   ERGPS(:)        ! GPS coordinates (ECEF)
!
! --- Simulated data
!
Real(Double), Pointer :: &
   WT(:)        ! Time
!
Type(Cartesian), Pointer :: &
   WRLEO(:),  & ! WOP LEO coordinates (J2000)
   WVLEO(:),  & ! WOP LEO velocity (J2000)
   WRGPS(:),  & ! WOP GPS coordinates (J2000)
   WVGPS(:)     ! WOP GPS velocity (J2000)
Real(Double), Pointer :: &
   WA(:,:),   & ! WOP amplitudes (time,channel)
   WS(:,:)      ! WOP phase excess (time,channel) [m]
!
Real(Double), Pointer :: &
   PCT(:),    & ! CT impact parameter 
   ECT(:,:),  & ! CT refraction angle (p,channel)
   ACT(:,:)     ! CT amplitude (p,channel)
!
! --- GO refraction
!
Real(Double), Pointer :: &
   EGO(:),    & ! GO refraction angles
   PGO(:),    & ! GO impact parameters
   ZGO(:),    & ! GO ray heights
   LAT(:,:)     ! Logarithmic attenuation (p,channel)
Real(Double), Allocatable :: &
   EGR(:),    & ! GO reflected refraction angles
   PGR(:),    & ! GO reflected impact parameters
   ZGR(:)       ! GO reflected ray heights
!
! --- Phase model
!
Real(Double), Pointer :: &
   SM(:),     & ! Phase excess model [m]
   PSM(:)       ! Model impact parameter [m]
Type(Geodetic), Allocatable :: &
   GSM(:)       ! Ray perigee locations [geo
!
! --- Local profiles
!
Type(Geodetic), Allocatable :: &
   GPZE(:)      ! Altitudes above ellipsoid/geoid [km]
Real(Double), Allocatable :: &
   Z(:),      & ! Altitudes above ellipsoid/geoid [km]
   RN(:),     & ! Refractivities [absolute]
   TD(:),     & ! Dry temperatures [K]
   TE(:),     & ! Interpolated ECHAM temperatures [K]
   QE(:),     & ! Interpolated ECHAM humidity [kg/kg]
   RNE(:),    & ! Interpolated ECHAM refractivity [n-1]
   RAE(:,:),  & ! Interpolated ECHAM specific absorption (z,channel) [km^-1]
   TDE(:)       ! Dry temperature from ECHAM refractivity [K]
!
! --- Chronometric variables
!
Integer :: TPS(8)   ! Program start time
Integer :: TPE(8)   ! Program end time
Integer :: TWS(8)   ! Wave optics start time
Integer :: TWE(8)   ! Wave optics end time
Integer :: TGS(8)   ! GO start time
Integer :: TGE(8)   ! GO end time
Integer :: ET(4)    ! Elapsed time (days, hours, minutes, seconds)
!----------------------------------------------------------


!----------------------------------------------------------
! 0. INITIALIZATION
!----------------------------------------------------------
! MPAS filename hardwired for test
!mpaspath ='MPAS/output.2018-07-18_01.nc'



! Output unit number for optional text files

   iuo=77


!--- 0.0. Error status initialization

Nullify(ErrStat)
Call Enter_Callee &
  ('Wave',        & ! <-- User routine
   ErrStat)         ! <-> Pointer to callee status


!--- 0.2. Time initialization

Call Date_and_Time(Values = TPS(:))


!--- 0.3. Getting command line


print*,' entering Get_Comline'

Call Get_Comline &
  (Opts,      & ! --> Command line options
   Err)           ! --> Command line error status

NC = Size(Opts%Freq)


!--- 0.1. Displaying program info

If (Opts%Vrb >= 0) then
   Call Display_Version()
End If


!--- 0.4. Error exit

If (Opts%GPS_Name == '' .and. Opts%ORB_Name == '') then
   Write (*,'(A/)')  'NO ORBIT DATA FILE SPECIFIED.'
   Call Help()
   Stop
End If

If (Opts%N_Files == 0) then
   Write (*,'(A/)')  'NO GRIB OR PHANTOM DATA FILE SPECIFIED.'
   Call Help()
   Stop
End If

If ((Opts%Atm_Type == atm_Phantom) .and. &
     (Opts%N_Files /= 1)) then
   Write (*,'(A/)')  'TOO MANY PHANTOM DATA FILES SPECIFIED.'
   Call Help()
   Stop
End If

If (Err) then
   Call Help()
   Stop
End If


!--- 0.5. Displaying parameters

If (Opts%Vrb >= 2) then
   Write(*,'(A,A6,   3X,A,A6,   3X,A,F6.1, 3X,A,F6.1/&
            &A,A6,   3X,A,F6.3, 3X,A,F6.1, 3X,A,F6.1/&
            &A,F6.1, 3X,A,F6.1, 3X,A,I6,   3X,A,F6.1/&
            &A,I6,   3X,A,F6.1, 3X,A,I6)')    &
      'Atm     = ',  &
         Atm_Name(Opts%Atm_Type),            & ! Atmosphere type
      'HG      = ',  &
         HG_Name(Opts%HG_Mode),              & ! Horizontal gradients mode
      'Hmax    = ', Opts%Hmax,               & ! Maximum height
      'Hmin    = ', Opts%Hmin,               & ! Minimum height
!
      'Earth   = ',  &
         Shape_Name(Opts%Earth),             & ! Earth's shape
      'DYN     = ', Opts%DYN,                & ! Minimum vertical scale of N
      'DX      = ', Opts%DX,                 & ! Step between phase screens
      'XLS     = ', Opts%XLS,                & ! Additional screen position
!
      'FZ      = ', Opts%FZ/Pi,              & ! Fresnel zone size [rad]
      'SR      = ', Opts%SR,                 & ! Sampling rate [Hz]
      'SDim    = ', Opts%SDim,               & ! Spatial dimension
      'HGO     = ', Opts%HGO,                & ! Maximum height for GO
!
      'NGO     = ', Opts%NGO,                & ! Number of GO rays
      'FWP     = ', Opts%FWP,                & ! Filter width for computing DPGPS/DP [km]
      'Vrb     = ', Opts%Vrb                   ! Option: verbosity
   Write(*,'((5(A,L1:2X)))')     &
      'opt_TS     = ', Opts%opt_TS,          & ! Time scaling
      'opt_LS     = ', Opts%opt_LS,          & ! Propagate to additional screen
      'opt_ECEF   = ', Opts%opt_ECEF,        & ! All coordinates in ECEF
      'opt_WO     = ', Opts%opt_WO,          & ! Wave optics
      'opt_GO     = ', Opts%opt_GO,          & ! Geometrical optics
      'opt_RH     = ', Opts%opt_RH,          & ! Option: output of ray heights
      'opt_AB     = ', Opts%opt_AB,          & ! Option: modeling of absorption
      'opt_AS     = ', Opts%opt_AS,          & ! Option: asymptotic solution
      'opt_FIO    = ', Opts%opt_FIO            ! FIO forward propagation
   Write(*,'(A,(99(:F8.4)))') &
      'Freq       =', 1d-9*Opts%Freq(:)        ! Frequencies
   Write(*,'(A,(99(:F8.4)))') &
      'AFreq      =', 1d-9*Opts%AFreq(:)       ! Absorption frequencies
   Write(*,'(2A)') &
      'MSIS       = ', Trim(Opts%MSIS_Path)    ! Path to MSIS coefficient files
   Write(*,'(2A)') &
      'Geoid_Path = ', Trim(Opts%Geoid_Path)   ! Path to geoid coefficient file
   Write(*,'(17(A4:1X))')  &
      'WOP', 'WCA', 'WCE', 'CAU', &
      'GOP', 'GRP', 'XBP', 'GOZ', &
      'GON', 'GOT', 'MOT', 'MOQ', &
      'MON', 'MOA', 'MTR', 'MTD', &
      ' NC'
   Write(*,'(17(A4:1X))')  &
      SL(Opts%Out%WOP), SL(Opts%Out%WCA), SL(Opts%Out%WCE), SL(Opts%Out%CAU),  &
      SL(Opts%Out%GOP), SL(Opts%Out%GRP), SL(Opts%Out%XBP), SL(Opts%Out%GOZ),  &
      SL(Opts%Out%GON), SL(Opts%Out%GOT), SL(Opts%Out%MOT), SL(Opts%Out%MOQ),  &
      SL(Opts%Out%MON), SL(Opts%Out%MOA), SL(Opts%Out%MTR), SL(Opts%Out%MTD),  &
      SL(Opts%Out%NC)
   Write(*,'()')
End If


!--- 0.6. Generating output name

If (Opts%OUT_Name /= '') then
   Continue
Else If (Opts%GPS_Name /= '') then
   Opts%OUT_Name = Opts%GPS_Name
Else If (Opts%ORB_Name /= '') then
   Opts%OUT_Name = Opts%ORB_Name
End If

Call Make_WOP_Name &
  (Opts%OUT_Name,    & ! <-- Occultation data file name
   WOP_Name)           ! --> Simulated data file name


!--- 0.7. Nullifying pointers

Nullify(TR)       ! Relative time of samples [sec]
Nullify(A)        ! Amplitudes (time,channel)
Nullify(S)        ! Phase excess (time,channel) [m]
Nullify(RLEO)     ! LEO coordinates (J2000)
Nullify(VLEO)     ! LEO velocity (J2000)
Nullify(RGPS)     ! GPS coordinates (J2000)
Nullify(VGPS)     ! GPS velocity (J2000)
Nullify(WT)       ! Time
Nullify(WRLEO)    ! LEO coordinates (J2000)
Nullify(WVLEO)    ! LEO velocity (J2000)
Nullify(WRGPS)    ! GPS coordinates (J2000)
Nullify(WVGPS)    ! GPS velocity (J2000)
Nullify(WA)       ! Amplitudes (time,channel)
Nullify(WS)       ! Phase excess (time,channel) [m]
Nullify(PCT)      ! CT impact parameter 
Nullify(ECT)      ! CT refraction angle
Nullify(ACT)      ! CT amplitude
Nullify(EGO)      ! GO refraction angles
Nullify(PGO)      ! GO impact parameters
Nullify(ZGO)      ! GO ray heights
Nullify(LAT)      ! Logarithmic attenuation



!----------------------------------------------------------
! 1. READING ORBIT DATA
!----------------------------------------------------------


!--- 1.1. Use of RO data file

If (Opts%GPS_Name /= '') then

   !--- 1.1.1. Reading RO data file

   If (Opts%Vrb >= 1) then
      Write (*,'(A/,T3,A)') 'Reading occultation data from: ', Trim(Opts%GPS_Name)
   End If

   Call Read_RO &
     (Opts%GPS_Name,  & ! <-- RO data file name
      '',             & ! <-- GPS bits file name
      MetaData,       & ! --> Occultation metadata
      ID,             & ! --> Occultation ID
      Year,           & ! --> Occultation year
      Month,          & ! --> Occultation month
      Day,            & ! --> Occultation day
      Hour,           & ! --> Occultation hour
      Minute,         & ! --> Occultation begin minute
      Second,         & ! --> Occultation begin second
      TR,             & ! --> Relative time of samples [sec]
      A,              & ! --> Amplitudes (time,channel)
      S,              & ! --> Phase excess (time,channel) [m]
      RLEO,           & ! --> LEO coordinates (J2000)
      VLEO,           & ! --> LEO velocity (J2000)
      RGPS,           & ! --> GPS coordinates (J2000)
      VGPS,           & ! --> GPS velocity (J2000)
      N,              & ! --> Number of data in file
      Stat,           & ! --> Error code
      LCF,            & ! --> Lost carrier flag
      ROFreq,         & ! --> Frequencies [Hz]
      DataType,       & ! --> Data type
      Opts%Vrb)         ! <~~ Verbosity level

   If (Stat > 0) then
      Write(*,'(A,I6)') 'DATA READ ERROR: ', Stat
      Stop
   End If


   !--- 1.1.2. Checking data

   Call Check_RO &
     (TR,        & ! <-- Relative time of samples [sec]
      A,         & ! <-- Amplitudes (time,channel)
      S,         & ! <-- Phase excess (time,channel) [m]
      RLEO,      & ! <-- LEO coordinates (J2000)
      VLEO,      & ! <-- LEO velocity (J2000)
      RGPS,      & ! <-- GPS coordinates (J2000)
      VGPS,      & ! <-- GPS velocity (J2000)
      Stat)        ! --> Error code

   If (Stat /= 0) then
      Write(*,'(A,I6)')  &
          'DATA CHECK ERROR. CODE: ', Stat
      Stop
   End If


!--- 1.2. Generating orbit data

Else If (Opts%ORB_Name /= '') then

   If (Opts%Vrb >= 1) then
      Write (*,'(A/,T3,A)') 'Reading orbit data from: ', Trim(Opts%ORB_Name)
   End If

   Call Generate_Orbit &
     (Opts%ORB_Name,  & ! <-- Name of orbit data file
      Opts%opt_ECEF,  & ! <-- Coordinates in ECEF
      ID,             & ! --> Occultation ID
      Year,           & ! --> Occultation year
      Month,          & ! --> Occultation month
      Day,            & ! --> Occultation day
      Hour,           & ! --> Occultation hour
      Minute,         & ! --> Occultation begin minute
      Second,         & ! --> Occultation begin second
      TR,             & ! --> Relative time of samples [sec]
      RLEO,           & ! --> LEO coordinates (J2000)
      VLEO,           & ! --> LEO velocity (J2000)
      RGPS,           & ! --> GPS coordinates (J2000)
      VGPS,           & ! --> GPS velocity (J2000)
      ErrStat)          ! <-> Error status

   If (Error(ErrStat)) then
      Write(*,'(A)') 'CANNOT GENERATE ORBIT DATA'
      Call Display_Status(ErrStat)
      Stop
   End If

   N = Size(TR)

   If (Opts%Vrb >= 2) then
      Write (*,'(I5,A)') N, ' data generated.'
   End If

   Metadata = t_Metadata &
     ('WOP simulation',    & ! Mission
      'N/A',               & ! Transmitting satellite
      'N/A',               & ! Receiving satellite
      'MEG WOP')             ! Calibration source

End If


!DEBUG: Reversion of occultation order
!A    => A(N:1:-1,1:2)
!S    => S(N:1:-1,1:2)
!RLEO => RLEO(N:1:-1)
!VLEO => VLEO(N:1:-1)
!RGPS => RGPS(N:1:-1)
!VGPS => VGPS(N:1:-1)


!--- 1.3. Determination of occultation geometry

Allocate(ERLEO(N))
Allocate(ERGPS(N))

If (Opts%opt_ECEF) then

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

Call Plane_Basis &
  (ERLEO,     & ! <-- LEO coordinates (ECEF)
   ERGPS,     & ! <-- GPS coordinates (ECEF)
   ERLC,      & ! <-- Curvature center (ECEF)
   RE,        & ! <-- Local curvature radius [km]
   AX,        & ! --> Occultation plane X basis vector
   AY)          ! --> Occultation plane Y basis vector


!----------------------------------------------------------
! 2. INITIALIZING ATMOSPHERIC MODEL
!----------------------------------------------------------


!--- 2.1. Initializing Geoid if required

If (Opts%Earth == shape_Geoid) then

   Call Set_Geoid_Path &
     (Opts%Geoid_Path) ! <-- Geoid coefficients search path

   Call Geoid_Init  &
     (geoid_asc,  & ! <-- Initialization mode
      ErrStat)      ! --> Error status

   If (Error(ErrStat)) then
      Write(*,'(A)') 'CANNOT INITIALIZE GEOID'
      Call Display_Status(ErrStat)
      Stop
   End If

End If


!--- 2.2. Initializing Earth

Call Earth_Init &
  (Opts%Earth,    & ! <-- Earth's shape
   Stat)            ! --> Status

If (Stat /= 0) then
   Write(*,'(A)') 'CANNOT INITIALIZE EARTH'
   Stop
End If


!--- 2.3. Setting MSIS path

Call Set_MSIS_Path &
  (Opts%MSIS_Path)  ! <-- MSIS coefficients search path


!--- 2.4. Initializing regular atmosphere model

Select Case (Opts%Atm_Type)
   Case (atm_ECHAM)
      If (Opts%Vrb >= 1) then
         Write (*,'(A/,(T3,A))') 'Reading ECHAM fields from: ', &
                             (Trim(Opts%Atm_Name(i)),i=1,Opts%N_Files)
      End If
   Case (atm_NCEP)
      If (Opts%Vrb >= 1) then
         Write (*,'(A/,(T3,A))') 'Reading NCEP fields from: ', &
                             (Trim(Opts%Atm_Name(i)),i=1,Opts%N_Files)
      End If

   Case (atm_MPAS)
      If (Opts%Vrb >= 1) then
         Write (*,'(A/,(T3,A))') 'Reading MPAS fields from: ', &
                             (Trim(Opts%Atm_Name(i)),i=1,Opts%N_Files)
      End If
   Case (atm_Phantom)
      If (Opts%Vrb >= 1) then
         Write (*,'(A/,(T3,A))') 'Reading phantom parameters from: ', &
                             Trim(Opts%Atm_Name(1))
      End If
End Select
!------------- JDH note 111020 ----------------------------------------
! Start skipping over atmopshere initialization if the case is MPAS.
! We don't want to read in the entire MPAS grid.  WE only want to read
! in the necessary subgrid in the GO_Invert step.

print*,'before Atmosphere_Init block'
print*,'N_files=',Opts%N_files

! JDH Debug replace netCDF file with the RO data to test open in Get_NETCDF
!Opts%Atm_Name(1)=Opts%GPS_Name
!if(Opts%Atm_type /= atm_MPAS) then
If (Opts%opt_AB) then
! Debug
!      Opts%GPS_Name,                   & ! <-- Pathnames of data files
   Call Atmosphere_Init &
     (Opts%Atm_Type,                   & ! <-- Type of atmosphere model
      Opts%Atm_Name(1:Opts%N_Files),   & ! <-- Pathnames of data files
      Opts%HG_Mode,                    & ! <-- Horizontal gradients mode
      Opts%Vrb,                        & ! <-- Verbosity level
      ErrStat,                         & ! <-> Error status
      XGP = GP,                        & ! <~~ Coordinates of center for HG_Mode=1d
      XFreq = Opts%AFreq)                ! <~~ Frequency channels [Hz]
Else
! Debug
!      Opts%GPS_Name,                   & ! <-- Pathnames of data files
   Call Atmosphere_Init &
     (Opts%Atm_Type,                   & ! <-- Type of atmosphere model
      Opts%Atm_Name(1:Opts%N_Files),   & ! <-- Pathnames of data files
      Opts%HG_Mode,                    & ! <-- Horizontal gradients mode
      Opts%Vrb,                        & ! <-- Verbosity level
      ErrStat,                         & ! <-> Error status
      XGP = GP)                          ! <~~ Coordinates of center for HG_Mode=1d
End If

!print*,'h=',h
!print*,'t=',t
!print*,'q=',q

If (Error(ErrStat)) then
   Write(*,'(A)') 'CANNOT INITIALIZE REGULAR ATMOSPHERE'
   Call Display_Status(ErrStat)
   Stop
End If

!----- JDH MPAS development/debug stop
!if(Opts%Atm_type == atm_MPAS) then
!   print*,'After initializing atm with MPAS'
!   print*,'Code to utilize MPAS data needs writing, stopping'
!   stop
!end if ! block around atmopsheric initialization for MPAS

!print*,'In Atmopshere_Init block'

!end if ! block around atmopsheric initialization for MPAS


!--- 2.5. Initializing turbulence model

If (Trim(Opts%TRB_Name) /= "") then

   If (Opts%Vrb >= 1) then
      Write (*,'(A/,(T3,A))') 'Reading turbulence parameters from: ', &
                          Trim(Opts%TRB_Name)
   End If

   Call Turbulence_Init &
     (Opts%TRB_Name,  & ! <-- Turbulence parameter file
      AX,             & ! <-- Occultation plane X basis vector
      AY,             & ! <-- Occultation plane Y basis vector
      Opts%Vrb,       & ! <-- Verbosity level
      ErrStat)          ! <-> Error status

   If (Error(ErrStat)) then
      Write(*,'(A)') 'CANNOT INITIALIZE TURBULENCE'
      Call Display_Status(ErrStat)
      Stop
   End If

End If


!----------------------------------------------------------
! 3. WAVE PROPAGATION
!----------------------------------------------------------


!--- 3.1. Wave propagation

If (Opts%opt_WO) then


   !--- 3.1.0. Time initialization

   Call Date_and_Time(Values = TWS(:))


   !--- 3.1.1. Wave propagation

   If (Opts%Vrb >= 1) then
      Write(*,'(A)') 'Wave propagation'
   End If

   Call Wave_Propagate &
     (Year,             & ! <-- Occultation year
      Month,            & ! <-- Occultation month
      Day,              & ! <-- Occultation day
      Hour,             & ! <-- Occultation hour
      Minute,           & ! <-- Occultation begin minute
      Second,           & ! <-- Occultation begin second
      TR,               & ! <-- Relative time of samples [sec]
      Opts%Freq,        & ! <-- Frequencies [Hz]
      Opts%AFreq,       & ! <-- Frequencies for absorption [Hz]
      RLEO,             & ! <-- LEO coordinates (J2000/ECEF)
      RGPS,             & ! <-- GPS coordinates (J2000/ECEF)
      Opts%Hmax,        & ! <-- Maximum height
      Opts%Hmin,        & ! <-- Minimum height
      Opts%DYN,         & ! <-- Minimum vertical scale of N
      Opts%DX,          & ! <-- Step between phase screens
      Opts%XLS,         & ! <-- Additional screen position
      Opts%opt_LS,      & ! <-- Propagate to additional screen
      Opts%FZ,          & ! <-- Initial last Fresnel zone size [Pi rad]
      Opts%SR,          & ! <-- Sampling rate [Hz]
      Opts%opt_TS,      & ! <-- Time scaling
      Opts%opt_AS,      & ! <-- Asymptotic solution
      Opts%opt_AB,      & ! <-- Modeling of absorption
      Opts%opt_ECEF,    & ! <-- Input/output coordinates in ECEF
      Opts%opt_FIO,     & ! <-- FIO forward propagation
      WT,               & ! --> WOP data time
      WA,               & ! --> WOP amplitudes (time,channel)
      WS,               & ! --> WOP phase excess (time,channel) [m]
      WRLEO,            & ! --> WOP LEO coordinates (J2000/ECEF)
      WVLEO,            & ! --> WOP LEO velocity (J2000/ECEF)
      WRGPS,            & ! --> WOP GPS coordinates (J2000/ECEF)
      WVGPS,            & ! --> WOP GPS velocity (J2000/ECEF)
      PCT,              & ! --> CT impact parameter
      ECT,              & ! --> CT refraction angle (p,channel)
      ACT,              & ! --> CT amplitude (p,channel)
      Stat,             & ! ~~> Error status
      Opts%Vrb)           ! <~~ Verbosity level

   print*,'After Wave_Propagate'

!----- JDH MPAS development/debug stop
if(Opts%Atm_type == atm_MPAS) then
   print*,'After Wave Propagate with MPAS'
!   print*,'Code to utilize MPAS data needs writing, stopping'
   print*,'Code to utilize MPAS data may need writing, errors may occur'
!   stop
end if

   If (Stat /= 0) then
      Write(*,'(A)') 'DATA SIZE MISMATCH IN MEDIA_WAVE_PROPAGATE.'
      Stop
   End If


   !--- 3.1.2. Elapsed time computation

   Call Date_and_Time(Values = TWE(:))

   Call Elapsed_Time &
     (TWS,       & ! <-- Start time
      TWE,       & ! <-- End time
      ET(1),     & ! --> Days
      ET(2),     & ! --> Hours
      ET(3),     & ! --> Minutes
      ET(4))       ! --> Seconds

   If (Opts%Vrb >= 1) then
      Write (*,'(A,I2.2,3(":",I2.2))') 'WOP elapsed time: ', ET(:)
   End If


   !--- 3.1.3. Writing output data

   If (Opts%Out%WCA) then

!      write(iuo,'(f7.2)',advance='NO',IOStat=IOstat) TD(1)
!      if(IOstat /=0) then
!         print*,'error in first write stopping'
!         close(iuo)
!         stop
!      end if

!      Call PutXY(GOT_Name, TD(:), Z(:), &
!                 XFmt='(F7.2)',   YFmt='(F11.5)', Stat=Stat)

      Do IC=1,NC
         nnn=size(pct)
         Write(WCA_Name,'(2A,I1)') Trim(WOP_Name), '.wca', IC
         open(iuo,file=WCA_Name)
         
         do kkn=1,nnn
            write(iuo,'(ES11.4,1x,f11.5)') ACT(kkn,IC),PCT(kkn)
         end do
         close(iuo)
         !         Call PutXY(WCA_Name, ACT(:,IC), PCT(:), &
         !                    XFmt='(ES11.4)', YFmt='(F11.5)', Stat=Stat)
      End Do
   End If

! This one doesn't seem to always work correctly.  Double check!!!
   If (Opts%Out%WCE) then
      Do IC=1,NC
         nnn=size(pct)
         Write(WCE_Name,'(2A,I1)') Trim(WOP_Name), '.wce', IC
         open(iuo,file=WCE_Name)
         
         do kkn=1,nnn
            write(iuo,'(ES11.4,1x,f11.5)') ECT(kkn,IC),PCT(kkn)
         end do
         close(iuo)
!         Call PutXY(WCE_Name, ECT(:,IC), PCT(:), &
!                    XFmt='(ES11.4)', YFmt='(F11.5)', Stat=Stat)
      End Do
   End If


!--- 3.2. Use of asymptotic Maslov solution

Else If (Opts%opt_AS) then


   !--- 3.2.0. Time initialization

   Call Date_and_Time(Values = TWS(:))


   !--- 3.2.1. Asymptotic wave propagation

   If (Opts%Vrb >= 1) then
      Write(*,'(A)') 'Asymptotic wave propagation'
   End If

   Call Asymptotic_Propagate &
     (Year,           & ! <-- Occultation year
      Month,          & ! <-- Occultation month
      Day,            & ! <-- Occultation day
      Hour,           & ! <-- Occultation hour
      Minute,         & ! <-- Occultation begin minute
      Second,         & ! <-- Occultation begin second
      Opts%SDim,      & ! <-- Spatial dimension
      TR,             & ! <-- Relative time of samples [sec]
      Opts%Freq,      & ! <-- Frequencies [Hz]
      Opts%AFreq,     & ! <-- Frequencies for absorption [Hz]
      RLEO,           & ! <-- LEO coordinates (J2000/ECEF)
      RGPS,           & ! <-- GPS coordinates (J2000/ECEF)
      Opts%opt_ECEF,  & ! <-- Coordinates in ECEF frame
      Opts%Hmax,      & ! <-- Maximum height
      Opts%DYN,       & ! <-- Minimum vertical scale of N
      Opts%NGO,       & ! <-- Upper number of GO rays
      Opts%SR,        & ! <-- Sampling rate [Hz]
      Opts%opt_AB,    & ! <-- Modeling of absorption
      Opts%FWP,       & ! <-- Filter width for computing DPGPS/DP [km]
      EGO,            & ! --> Refraction angles
      PGO,            & ! --> Impact parameters
      LAT,            & ! --> Logarithmic attenuation
      ZGO,            & ! --> Ray heights
      GP,             & ! --> Occultation point
      ERLC,           & ! --> Local curvature center (ECEF)
      RE,             & ! --> Local curvature radius
      WT,             & ! --> WOP data time
      WA,             & ! --> WOP amplitudes (time,channel)
      WS,             & ! --> WOP phase excess (time,channel) [m]
      WRLEO,          & ! --> WOP LEO coordinates (J2000/ECEF)
      WVLEO,          & ! --> WOP LEO velocity (J2000/ECEF)
      WRGPS,          & ! --> WOP GPS coordinates (J2000/ECEF)
      WVGPS,          & ! --> WOP GPS velocity (J2000/ECEF)
      Stat,           & ! ~~> Error status
      Opts%Vrb)         ! <~~ Verbosity level

   !--- 3.2.2. Elapsed time computation

   Call Date_and_Time(Values = TWE(:))

   Call Elapsed_Time &
     (TWS,       & ! <-- Start time
      TWE,       & ! <-- End time
      ET(1),     & ! --> Days
      ET(2),     & ! --> Hours
      ET(3),     & ! --> Minutes
      ET(4))       ! --> Seconds

   If (Opts%Vrb >= 1) then
      Write (*,'(A,I2.2,3(":",I2.2))') 'WOP elapsed time: ', ET(:)
   End If



!--- 3.3. Use of RO orbit data
!---      if no wave propagation requested

Else

   WT    => TR  (1:N)
   WRLEO => RLEO(1:N)
   WRGPS => RGPS(1:N)


End If


!--- 3.4. Writing simulated data

If (Opts%opt_WO .or. Opts%opt_AS) then

   If (Size(Opts%Freq) == 2) then
      GPS_Freq = (MaxVal(Abs(Opts%Freq(:)-Freq_GPS(:))) < 0.001)
   Else
      GPS_Freq = .False.
   End If

   If (GPS_Freq .and. Opts%Out%NC) then
      Data_Type = dt_ATMPHS
   Else
      Data_Type = dt_LEOLEO
   End If

   If (Opts%Out%WOP .or. Opts%Out%NC) then
      If (Opts%Vrb >= 1) then
         Write (*,'(3A/,T3,A)')                       &
            'Writing simulated RO data in ',          &
            Trim(dt_Name(Data_Type)),' format to: ',  &
            Trim(WOP_Name)
      End If
      Call Write_RO &
        (Data_Type,        & ! <-- Output data type
         WOP_Name,         & ! <-- RO data file name
         MetaData,         & ! <-- Occultation metadata
         ID,               & ! <-- Occultation ID
         Year,             & ! <-- Occultation year
         Month,            & ! <-- Occultation month
         Day,              & ! <-- Occultation day
         Hour,             & ! <-- Occultation hour
         Minute,           & ! <-- Occultation begin minute
         Second,           & ! <-- Occultation begin second
         WT,               & ! <-- Relative time of samples [sec]
         WA,               & ! <-- Amplitudes (time,channel)
         WS,               & ! <-- Phase excess (time,channel) [m]
         WRLEO,            & ! <-- LEO coordinates (J2000)
         WVLEO,            & ! <-- LEO velocity (J2000)
         WRGPS,            & ! <-- GPS coordinates (J2000)
         WVGPS,            & ! <-- GPS velocity (J2000)
         Freq=Opts%Freq,   & ! <~~ Frequencies [Hz]
         Stat=Stat)          ! --> Error code
      If (Stat /= 0) then
         Write(*,'(2X,A,I3)') '*** EXIT: Error in Write_RO: ', Stat
      End If
   End If

End If


!----------------------------------------------------------
! 4. GO PROPAGATION
!----------------------------------------------------------

If (Opts%opt_GO) then

   If (.not. Associated(EGO)) then


      !--- 4.0. Time initialization

      Call Date_and_Time(Values = TGS(:))


      !--- 4.1. GO propagation

      If (Opts%Vrb >= 1) then
         Write(*,'(A)') 'GO propagation'
      End If
   
      Allocate(EGO(Opts%NGO))
      Allocate(PGO(Opts%NGO))
      Allocate(ZGO(Opts%NGO))
      Allocate(LAT(Opts%NGO,NC))
      Allocate(EGR(NGR))
      Allocate(PGR(NGR))
      Allocate(ZGR(NGR))

      print*,'calling GO_Propogate from Wave'

      Call GO_Propagate &
        (Year,           & ! <-- Occultation year
         Month,          & ! <-- Occultation month
         Day,            & ! <-- Occultation day
         Hour,           & ! <-- Occultation hour
         Minute,         & ! <-- Occultation begin minute
         Second,         & ! <-- Occultation begin second
         WT,             & ! <-- Relative time of samples [sec]
         WRLEO,          & ! <-- LEO coordinates (J2000/ECEF)
         WRGPS,          & ! <-- GPS coordinates (J2000/ECEF)
         Opts%HGO,       & ! <-- Maximum height
         Opts%DYN,       & ! <-- Ray integration step parameter
         Opts%AFreq,     & ! <-- Absorption frequencies [Hz]
         Opts%opt_ECEF,  & ! <-- Coordinates in ECEF frame
         EGO,            & ! --> GO refraction angle
         PGO,            & ! --> GO impact parameter
         ZGO,            & ! --> GO ray heights
         LAT,            & ! --> Logarithmic attenuation
         EGR,            & ! --> GO reflected refraction angle
         PGR,            & ! --> GO reflected impact parameter
         ZGR,            & ! --> GO reflected ray heights
         GP,             & ! --> Occultation point
         ERLC,           & ! --> Local curvature center (ECEF)
         RE,             & ! --> Local curvature radius
         Opts%Vrb)         ! <~~ Verbosity mode


      If (Opts%Out%CAU) then

         CAU_Name = Trim(WOP_Name)//'.cau'

         Call GO_Caustics &
           (CAU_Name,  & ! <-- Output file name
            EGO,       & ! <-- GO refraction angle
            PGO,       & ! <-- GO impact parameter [km]
            XBP)         ! --> BPP position [km]

      End If


      !--- 4.2. Elapsed time computation

      Call Date_and_Time(Values = TGE(:))

      Call Elapsed_Time &
        (TGS,       & ! <-- Start time
         TGE,       & ! <-- End time
         ET(1),     & ! --> Days
         ET(2),     & ! --> Hours
         ET(3),     & ! --> Minutes
         ET(4))       ! --> Seconds

      If (Opts%Vrb >= 1) then
         Write (*,'(A,I2.2,3(":",I2.2))') 'GO elapsed time: ', ET(:)
      End If


   End If


   !--- 4.3. Writing output data

   If (Opts%Out%GOP) then
      nnn=size(EGO)
      GOP_Name = Trim(WOP_Name)//'.gop'
      open(iuo,file=GOP_Name)

      If (Opts%opt_RH) then
         do kkn=1,nnn
            write(iuo,'(ES11.4,1x,f11.5)') EGO(kkn),ZGO(kkn)
         end do
!         Call PutXY(GOP_Name, EGO, ZGO, &
!                    XFmt='(ES11.4)', YFmt='(F11.5)', Stat=Stat)
      Else
         do kkn=1,nnn
            write(iuo,'(ES11.4,1x,f11.5)') EGO(kkn),PGO(kkn)
         end do
!         Call PutXY(GOP_Name, EGO, PGO, &
!                    XFmt='(ES11.4)', YFmt='(F11.5)', Stat=Stat)
      End If
      close(iuo)
   End If

   If (Opts%Out%GRP) then
      If (Allocated(EGR)) then
         nnn=size(EGR)
         GRP_Name = Trim(WOP_Name)//'.grp'
         open(iuo,file=GRP_Name)
         If (Opts%opt_RH) then
            do kkn=1,nnn
               write(iuo,'(ES11.4,1x,f11.5)') EGR(kkn),ZGR(kkn)
            end do
            
!            Call PutXY(GRP_Name, EGR, ZGR, &
!                       XFmt='(ES11.4)', YFmt='(F11.5)', Stat=Stat)
         Else
            do kkn=1,nnn
               write(iuo,'(ES11.4,1x,f11.5)') EGR(kkn),PGR(kkn)
            end do
!            Call PutXY(GRP_Name, EGR, PGR, &
!                       XFmt='(ES11.4)', YFmt='(F11.5)', Stat=Stat)
         End If
         close(iuo)
      End If
   End If

   If (Opts%Out%XBP) then
      XBP_Name = Trim(WOP_Name)//'.xbp'
      Call FileOpen(  &
         Name   = XBP_Name,  & ! <-- File pathname 
         Unit   = IU,        & ! --> Unit number
         Stat   = Stat)        ! --> Error code
      Write(IU,'(F5.1)') XBP
      Close(Unit = IU)
   End If

End If


!----------------------------------------------------------
! 5. GO INVERSION
!----------------------------------------------------------

If (Opts%opt_GO) then

   If (Opts%Vrb >= 1) then
      Write(*,'(A)') 'GO inversion'
   End If

!   print*,'calling go inversion from Wave'

   !-- 5.1. Inversion

   Opts%NGO = Size(EGO)

   Allocate(Z(1:Opts%NGO))
   Allocate(RN(1:Opts%NGO))
   Allocate(TD(1:Opts%NGO))
   Allocate(TE(1:Opts%NGO))
   Allocate(QE(1:Opts%NGO))
   Allocate(RNE(1:Opts%NGO))
   Allocate(RAE(1:Opts%NGO,NC))
   Allocate(TDE(1:Opts%NGO))

!   print*,'opts%ngo=',Opts%NGO
!   print*,'ego=',ego
!   print*,'pgo=',pgo
!   print*,'z(1:10)=',z(1:10)
!   print*,'te(1:10)=',te(1:10)
!   print*,'rn=',rn

   Call GO_Invert &
     (EGO,        & ! <-- GO refraction angle
      PGO,        & ! <-- GO impact parameter
      GP,         & ! <-- Occultation point
      RE,         & ! <-- Local curvature radius
      Opts%AFreq, & ! <-- Absorption frequency [Hz]
      Z,          & ! --> Altitudes [km]
      RN,         & ! --> Refractivities [absolute]
      TD,         & ! --> Dry temperatures [K]
      TE,         & ! --> Interpolated ECHAM temperatures [K]
      QE,         & ! --> Interpolated ECHAM humidity [kg/kg]
      RNE,        & ! --> Interpolated ECHAM refractivity [n-1]
      RAE,        & ! --> Interpolated ECHAM specific absorption [km^-1]
      TDE,        & ! --> Dry temperature from ECHAM refractivity [K]
      Opts%Vrb)     ! <~~ Verbosity mode

   print*,'opts%ngo=',Opts%NGO
!   print*,'ego=',ego
!   print*,'pgo=',pgo
!   print*,'z(1:10)=',z(1:10)
!   print*,'te(1:10)=',te(1:10)
!   print*,'td(1:10)=',te(1:10)
!   print*,'rn=',rn

   !--- 5.2. Writing output data

   If (Opts%Out%GOZ) then
      nnn=size(EGO)
      GOZ_Name = Trim(WOP_Name)//'.goz'
      open(iuo,file=GOZ_Name)
      do kkn=1,nnn
         write(iuo,'(ES11.4,1x,f11.5)') EGO(kkn),Z(kkn)
      end do
!      Call PutXY(GOZ_Name, EGO(:), Z(:), &
!                 XFmt='(ES11.4)', YFmt='(F11.5)', Stat=Stat)
      close(iuo)
   End If

   If (Opts%Out%GON) then
      nnn=size(RN)
      GON_Name = Trim(WOP_Name)//'.gon'
      open(iuo,file=GON_Name)
      do kkn=1,nnn
         write(iuo,'(ES11.4,1x,f11.5)') 1d6*RN(kkn),Z(kkn)
      end do
!      Call PutXY(Trim(GON_Name), 1d6*RN(:), Z(:), &
!                 XFmt='(ES11.4)', YFmt='(F11.5)', Stat=Stat)
   End If

   If (Opts%Out%GOT) then
!      print*,'td=',td(:)
!      print*,'z=',z(:)
      nnn=size(td)
      GOT_Name = Trim(WOP_Name)//'.got'
!      open(iuo,file=GOT_Name,status='unknown')
      open(iuo,file=GOT_Name)
!      write(iuo,'(f7.2)',advance='NO',IOStat=IOstat) TD(1)
!      if(IOstat /=0) then
!         print*,'error in first write stopping'
!         close(iuo)
!         stop
!      end if

      do kkn=1,nnn
         write(iuo,'(F7.2,1x,f11.5)') TD(kkn),Z(kkn)
      end do
      close(iuo)
!      Call PutXY(GOT_Name, TD(:), Z(:), &
!                 XFmt='(F7.2)',   YFmt='(F11.5)', Stat=Stat)
   End If

   If (Opts%Out%MOT) then
      nnn=size(te)
      MOT_Name = Trim(WOP_Name)//'.mot'
      open(iuo,file=MOT_Name)
      do kkn=1,nnn
         write(iuo,'(F7.2,1x,f11.5)') TE(kkn),Z(kkn)
      end do
!      Call PutXY(MOT_Name, TE(:), Z(:), &
!                 XFmt='(F7.2)', YFmt='(F11.5)', Stat=Stat)
      close(iuo)
   End If

   If (Opts%Out%MOQ) then
      nnn=size(qe)
      MOQ_Name = Trim(WOP_Name)//'.moq'
      open(iuo,file=MOQ_Name)
      do kkn=1,nnn
         write(iuo,'(ES11.4,1x,f11.5)') 1d3*QE(kkn),Z(kkn)
      end do
      close(iuo)
!      Call PutXY(MOQ_Name, 1d3*QE(:), Z(:), &
!                 XFmt='(ES11.4)', YFmt='(F11.5)', Stat=Stat)
   End If

   If (Opts%Out%MON) then
      nnn=size(RNE)
      MON_Name = Trim(WOP_Name)//'.mon'
      open(iuo,file=MON_Name)
      do kkn=1,nnn
         write(iuo,'(ES11.4,1x,f11.5)') 1d6*RNE(kkn),Z(kkn)
      end do
      close(iuo)
!      Call PutXY(MON_Name, 1d6*RNE(:), Z(:), &
!                 XFmt='(ES11.4)', YFmt='(F11.5)', Stat=Stat)
   End If

   If (Opts%Out%MOA) then

      Do IC=1,NC
         Write(MOA_Name,'(2A,I1)') Trim(WOP_Name), '.moa', IC
         nnn=size(RAE(:,IC))
         open(iuo,file=MOA_Name)
         do kkn=1,nnn
            write(iuo,'(ES11.4,1x,f11.5)') RAE(kkn,ic),Z(kkn)
         end do
         close(iuo)
!         Call PutXY(MOA_Name, RAE(:,IC), Z(:), &
!                    XFmt='(ES11.4)', YFmt='(F11.5)', Stat=Stat)
      End Do

      Do IC1=1,NC
         Do IC2=IC1+1,NC
            nnn=size(RAE(:,IC1))
            Write(MOA_Name,'(2A,2I1)') Trim(WOP_Name), '.moa', IC1,IC2
            open(iuo,file=MOA_Name)
            do kkn=1,nnn
               write(iuo,'(ES11.4,1x,f11.5)') &
                    RAE(kkn,IC2)-RAE(kkn,IC1),Z(kkn)
            end do
            close(iuo)
!            Call PutXY(MOA_Name, (RAE(:,IC2)-RAE(:,IC1)), Z(:), &
!                    XFmt='(ES11.4)', YFmt='(F11.5)', Stat=Stat)
         End Do
      End Do

   End If

   If (Opts%Out%MTR) then

      Do IC=1,NC
         nnn=size(LAT(:,IC))
         Write(LAT_Name,'(2A,I1)') Trim(WOP_Name), '.mtr', IC
         open(iuo,file=LAT_Name)
         do kkn=1,nnn
            write(iuo,'(ES11.4,1x,f11.5)') &
                 LAT(kkn,IC)*20./Log(10.),Z(kkn)
         end do
         close(iuo)
!         Call PutXY(LAT_Name, LAT(:,IC)*20./Log(10.), ZGO(:), &
!              XFmt='(ES11.4)', YFmt='(F11.5)', Stat=Stat)
      End Do

      Do IC1=1,NC
         Do IC2=IC1+1,NC
            nnn=size(Lat(:,IC1))
            Write(LAT_Name,'(2A,2I1)') Trim(WOP_Name), '.mtr', IC1,IC2
            open(iuo,file=LAT_Name)
            do kkn=1,nnn
               write(iuo,'(ES11.4,1x,f11.5)') &
                    (LAT(kkn,IC2)-LAT(kkn,IC1))*20./Log(10.),Z(kkn)
            end do
            close(iuo)
!            Call PutXY(LAT_Name, (LAT(:,IC2)-LAT(:,IC1))*20./Log(10.), ZGO(:), &
!                       XFmt='(ES11.4)', YFmt='(F11.5)', Stat=Stat)
         End Do
      End Do

   End If

   If (Opts%Out%MTD) then
      nnn=size(tde)
      MTD_Name = Trim(WOP_Name)//'.mtd'
      open(iuo,file=MTD_Name)
      do kkn=1,nnn
         write(iuo,'(F7.2,1x,f11.5)') &
              TDE(kkn),Z(kkn)
      end do
      close(iuo)

!      Call PutXY(MTD_Name, TDE(:), Z(:), &
!                 XFmt='(F7.2)',   YFmt='(F11.5)', Stat=Stat)
   End If

End If


!----------------------------------------------------------
! 6. LOCAL PROFILES RETRIEVAL
!----------------------------------------------------------

If ((.not. Opts%opt_GO) .and. (.not. Opts%opt_WO)) then

   If (Opts%Vrb >= 1) then
      Write(*,'(A)') 'Local profiles retrieval'
   End If


   !--- 6.1. Phase model computation

   Allocate(SM(N))
   Allocate(PSM(N))
   Allocate(GSM(N))

   Call MSIS_Phase_Excess &
     (TR(:),     & ! <-- Relative time of samples [sec]
      ERLEO(:),  & ! <-- LEO coordinates (ECEF)
      ERGPS(:),  & ! <-- GPS coordinates (ECEF)
      GP,        & ! <-- Occultation point (geodetic)
      ERLC,      & ! <-- Local curvature center (ECEF)
      RE,        & ! <-- Local curvature radius
      Opts%Vrb,  & ! <-- Verbosity level
      I1, I2,    & ! --> Phase model limits
      SM(:),     & ! --> Phase excess model [m]
      PSM(:),    & ! --> Model impact parameter [m]
      GSM(:))      ! ~~> Ray perigee locations [geodetic]


   !--- 6.2. Local profiles retrieval

   Allocate(Z   (Opts%NGO))
   Allocate(RN  (Opts%NGO))
   Allocate(TD  (Opts%NGO))
   Allocate(TE  (Opts%NGO))
   Allocate(QE  (Opts%NGO))
   Allocate(RNE (Opts%NGO))
   Allocate(RAE (Opts%NGO,NC))
   Allocate(TDE (Opts%NGO))
   Allocate(GPZE(Opts%NGO))

   Do i=1,Opts%NGO
      Z(i)      = Opts%Hmax*Real(i-1)/Real(Opts%NGO - 1)
      GPZE(i)%H = Z(i)
      Call Linear &
        (GSM(:)%H,         & ! <-- Argument grid
         GSM(:)%Phi,       & ! <-- Gridded function
         GPZE(i)%H,        & ! <-- Interpolation point
         GPZE(i)%Phi,      & ! --> Interpolated function value
         CExt=.True.)        ! <~~ Constant/linear extrapolation
      Call Linear &
        (GSM(:)%H,         & ! <-- Argument grid
         GSM(:)%Lambda,    & ! <-- Gridded function
         GPZE(i)%H,        & ! <-- Interpolation point
         GPZE(i)%Lambda,   & ! --> Interpolated function value
         CExt=.True.)        ! <~~ Constant/linear extrapolation
   End Do

   Call Local_Profiles &
     (GPZE,       & ! <-- Perigee locations [geodetic]
      GP,         & ! <-- Occultation point
      Opts%AFreq, & ! <-- Absorption frequency [Hz]
      TE,         & ! --> Interpolated ECHAM temperatures [K]
      QE,         & ! --> Interpolated ECHAM humidity [kg/kg]
      RNE,        & ! --> Interpolated ECHAM refractivity [n-1]
      RAE,        & ! --> Interpolated ECHAM specific absorption [km^-1]
      TDE,        & ! --> Dry temperature from ECHAM refractivity [K]
      Opts%Vrb)     ! <~~ Verbosity level


   !--- 6.3. Writing output data

   If (Opts%Out%MOT) then
      nnn=size(te)
      MOT_Name = Trim(WOP_Name)//'.mot'
      open(iuo,file=MOT_Name)
      do kkn=1,nnn
         write(iuo,'(F7.2,3(1x,f11.5))') &
              te(kkn),gpze(kkn)%h,gpze(kkn)%phi,gpze(kkn)%Lambda
      end do
!      Call PutXY(MOT_Name, TE, GPZE(:)%H, GPZE(:)%Phi, GPZE(:)%Lambda, &
!                 XFmt='(F7.2)', YFmt='(F11.5)', Stat=Stat)
      close(iuo)
   End If

   If (Opts%Out%MOQ) then
      nnn=size(qe)
      MOQ_Name = Trim(WOP_Name)//'.moq'
      open(iuo,file=MOQ_Name)
      do kkn=1,nnn
         write(iuo,'(F7.2,3(1x,f11.5))') &
              1.0d3*qe(kkn),gpze(kkn)%H,gpze(kkn)%Phi,gpze(kkn)%Lambda
      end do
!      Call PutXY(MOQ_Name, 1d3*QE, GPZE(:)%H, GPZE(:)%Phi, GPZE(:)%Lambda, &
!                 XFmt='(ES11.4)', YFmt='(F11.5)', Stat=Stat)
      close(iuo)
   End If

   If (Opts%Out%MON) then
      nnn=size(rne)
      MON_Name = Trim(WOP_Name)//'.mon'
      open(iuo,file=MON_Name)
      do kkn=1,nnn
         write(iuo,'(F7.2,3(1x,f11.5))') &
              1.0d6*rne(kkn),gpze(kkn)%H,gpze(kkn)%Phi,gpze(kkn)%Lambda
      end do
!      Call PutXY(MON_Name, 1d6*RNE, GPZE(:)%H, GPZE(:)%Phi, GPZE(:)%Lambda, &
!                 XFmt='(ES11.4)', YFmt='(F11.5)', Stat=Stat)
      close(iuo)
   End If

   If (Opts%Out%MOA) then
      Do IC=1,NC
         nnn=size(RAE(:,IC))
         Write(MOA_Name,'(2A,I1)') Trim(WOP_Name), '.moa', IC
         open(iuo,file=MOA_Name)
         do kkn=1,nnn
            write(iuo,'(ES11.4,1x,f11.5)') &
                 rae(kkn,ic),z(kkn)
         end do
!         Call PutXY(MOA_Name, RAE(:,IC), Z(:), &
!                    XFmt='(ES11.4)', YFmt='(F11.5)', Stat=Stat)
      End Do

      Do IC1=1,NC
         Do IC2=IC1+1,NC
            nnn=size(RAE(:,IC))
            Write(MOA_Name,'(2A,2I1)') Trim(WOP_Name), '.moa', IC1,IC2
            open(iuo,file=MOA_Name)
            do kkn=1,nnn
               write(iuo,'(ES11.4,1x,f11.5)') &
                    (rae(kkn,ic2)-rae(kkn,ic1)),z(kkn)
            end do
!            Call PutXY(MOA_Name, (RAE(:,IC2)-RAE(:,IC1)), Z(:), &
!                    XFmt='(ES11.4)', YFmt='(F11.5)', Stat=Stat)
            close(iuo)
         End Do
      End Do

   End If

   If (Opts%Out%MTD) then
      nnn=size(tde)
      MTD_Name = Trim(WOP_Name)//'.mtd'
      open(iuo,file=MTD_Name)
      do kkn=1,nnn
         write(iuo,'(F7.2,1X,F11.5)') tde(kkn),z(kkn)
      end do
!      Call PutXY(MTD_Name, TDE(:), Z(:), &
!                 XFmt='(F7.2)',   YFmt='(F11.5)', Stat=Stat)
      close(iuo)
   End If

   !--- 6.4. Memory deallocation

   Deallocate(SM,    Stat=Stat)
   Deallocate(PSM,   Stat=Stat)
   Deallocate(GSM,   Stat=Stat)

   Deallocate(Z,     Stat=Stat)
   Deallocate(RN,    Stat=Stat)
   Deallocate(TD,    Stat=Stat)
   Deallocate(TE,    Stat=Stat)
   Deallocate(QE,    Stat=Stat)
   Deallocate(RNE,   Stat=Stat)
   Deallocate(RAE,   Stat=Stat)
   Deallocate(TDE,   Stat=Stat)
   Deallocate(GPZE,  Stat=Stat)

End If


!----------------------------------------------------------
! 7. EXIT
!----------------------------------------------------------

!--- 7.1. Memory deallocation

Deallocate(ERLEO, Stat=Stat)
Deallocate(ERGPS, Stat=Stat)

Deallocate(EGO,   Stat=Stat)
Deallocate(PGO,   Stat=Stat)
Deallocate(ZGO,   Stat=Stat)
Deallocate(LAT,   Stat=Stat)
Deallocate(EGR,   Stat=Stat)
Deallocate(PGR,   Stat=Stat)
Deallocate(ZGR,   Stat=Stat)


!--- 7.2. Comptation of elapsed time

Call Date_and_Time(Values = TPE(:))

Call Elapsed_Time &
  (TPS,       & ! <-- Start time
   TPE,       & ! <-- End time
   ET(1),     & ! --> Days
   ET(2),     & ! --> Hours
   ET(3),     & ! --> Minutes
   ET(4))       ! --> Seconds

!Opts%Atm_Name(1)=mpaspath
!print*,'Opts%Atm_Name=',Opts%Atm_Name(1)

If (Opts%Vrb >= 0) then
   Write (*,'(/A,I2.2,3(":",I2.2)/)') '*** Elapsed time: ', ET(:)
End If

Stop


End Program Wave


