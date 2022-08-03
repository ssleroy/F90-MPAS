!
Module Wave_Comline
!
! Processing of the command line for Wave.
!----------------------------------------------------------
! (C) Copyright 1999-2009, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 22 Feb 1999 | Original version.
!   2.0   | 25 May 2008 | Options structure.
!   2.1   | 16 Jun 2008 | Opts%opt_WO, Opts%opt_GO,
!         |             | Opts%MSIS_Path.
!   2.2   | 24 Nov 2008 | Opts%opt_Geoid, Opts%Geoid_Path.
!   2.3   | 17 Jun 2009 | Corrected output of MTD for -o=all.
!   2.4   | 20 Sep 2009 | Updated for -atm=type,[hg],ss,...
!   2.5   | 01 Oct 2009 | Opts%Out%NC.
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Parameters:
    Double,       &
    Pi
!
Use Atmosphere, only: &
! Imported Parameters:
    atm_Undefined,   &
    atm_Phantom,     &
    atm_ECHAM,       &
    atm_NCEP,        &
    atm_MPAS,        &
    hg_Undefined,    &
    hg_1d,           &
    hg_3d
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Public Type Definitions:
!
Type t_Out            ! Output extensions
   Logical :: &
      WOP,      & ! 01 Simulated occultation data
      WCA,      & ! 02 CT amplitude
      WCE,      & ! 03 CT BA for
      CAU,      & ! 04 Caustic structure [km]
      GOP,      & ! 05 GO BA profile [rad, km]
      GRP,      & ! 06 GO reflected BA profile [rad, km]
      XBP,      & ! 07 Back propagation plane position [km]
      GOZ,      & ! 08 GO BA vs perigee height [rad, km]
      GON,      & ! 09 Refractivity profile from GO BA [N-units, km]
      GOT,      & ! 10 Dry temperature profile from GO BA [K, km]
      MOT,      & ! 11 Local profile of temperature [K, km]
      MOQ,      & ! 12 Local profile of humidity [g/kg, km]
      MON,      & ! 13 Local profile of refractivity [N-units]
      MOA,      & ! 14 Local profile of specific absorption
      MTR,      & ! 15 Transmission
      MTD,      & ! 16 Local dry temperature [K]
      NC          ! 17 Simulated data in atmPhs-NetCDF format
End Type t_Out
!
Type t_Opts           ! Command line options
   Character(Len=255) :: &
      GPS_Name       ! GPS-level2 data file pathname
   Character(Len=255) :: &
      ORB_Name       ! Orbit data file pathname
   Integer            :: &
      Atm_Type       ! Atmosphere type
   Integer            :: &
      HG_Mode        ! Horizontal gradients mode
   Character(Len=255) :: &
      Atm_Name(10)   ! Atmosphere initialization file pathname(s)
   Integer            :: &
      N_Files        ! Number of atmosphere initialization files
   Character(Len=255) :: &
      TRB_Name       ! Turbulence parameter file pathname
   Integer            :: &
      Earth          ! Earth shape
   Character(Len=255) :: &
      Geoid_path     ! Path to geoid coefficient file
   Character(Len=255) :: &
      MSIS_path      ! Path to MSIS coefficient files
   Character(Len=255) :: &
      OUT_Name       ! Template for output names
   Real(Double), Pointer :: &
      Freq(:)        ! Frequencies [Hz]
   Real(Double), Pointer :: &
      AFreq(:)       ! Frequencies for absorption [Hz]
   Real(Double)       :: &
      Hmax           ! Maximum height for WO
   Real(Double)       :: &
      Hmin           ! Minimum height for WO
   Real(Double)       :: &
      DYN            ! Minimum vertical scale of N
   Real(Double)       :: &
      DX             ! Step between phase screens
   Real(Double)       :: &
      XLS            ! Additional screen position
   Logical            :: &
      opt_LS         ! Propagate to additional screen
   Real(Double)       :: &
      FZ             ! Initial last Fresnel zone size [Pi rad]
   Real(Double)       :: &
      SR             ! Sampling rate [Hz]
   Logical            :: &
      opt_TS         ! Time scaling
   Logical            :: &
      opt_ECEF       ! All coordinates in ECEF
   Logical            :: &
      opt_WO         ! Wave optics
   Logical            :: &
      opt_GO         ! Geometrical optics
   Logical            :: &
      opt_RH         ! Option: output of ray heights
   Logical            :: &
      opt_AB         ! Option: modeling of absorption
   Logical            :: &
      opt_AS         ! Option: asymptotic solution
   Logical            :: &
      opt_FIO        ! Option: FIO forward propagation
   Integer            :: &
      SDim           ! Spatial dimension
   Real(Double)       :: &
      HGO            ! Maximum height for GO
   Integer            :: &
      NGO            ! Number of geometric optical rays
   Real(Double)       :: &
      FWP            ! Filter width for computing DPGPS/DP [km]
   Integer            :: &
      Vrb            ! Verbosity level
   Type(t_Out)        :: &
      Out            ! Output options
End Type t_Opts
!----------------------------------------------------------
! Public Parameters:
!
Type(t_Out), Parameter :: &
   No_Output = t_Out ( &
      ! WOP     WCA     WCE     CAU     GOP     GRP     XBP     GOZ
      .False.,.False.,.False.,.False.,.False.,.False.,.False.,.False., &
      ! GON     GOT     MOT     MOQ     MON     MOA     MTR     MTD
      .False.,.False.,.False.,.False.,.False.,.False.,.False.,.False., &
      ! NC
      .False.), &
   Min_Output = t_Out ( &
      ! WOP     WCA     WCE     CAU     GOP     GRP     XBP     GOZ
      .True., .False.,.False.,.False.,.False.,.False.,.False.,.False., &
      ! GON     GOT     MOT     MOQ     MON     MOA     MTR     MTD
      .False.,.False.,.False.,.False.,.False.,.False.,.False.,.False., &
      ! NC
      .False.), &
   Script_Output = t_Out ( &
      ! WOP     WCA     WCE     CAU     GOP     GRP     XBP     GOZ
      .True., .False.,.False.,.False.,.True., .True., .False.,.False., &
      ! GON     GOT     MOT     MOQ     MON     MOA     MTR     MTD
      .True., .True., .True., .True., .True., .True., .True., .True., &
      ! NC
      .False.), &
   All_Output = t_Out ( &
      ! WOP     WCA     WCE     CAU     GOP     GRP     XBP     GOZ
      .True., .True., .True., .True., .True., .True., .True., .True.,  &
      ! GON     GOT     MOT     MOQ     MON     MOA     MTR     MTD
      .True., .True., .True., .True., .True., .True., .True., .True., &
      ! NC
      .False.)
!----------------------------------------------------------
!
Contains


!==========================================================
Subroutine Get_Comline &
  (Opts,      & ! --> Command line options
   Err)         ! --> Command line error status !
! Parsing command line for Wave.
!----------------------------------------------------------
! (C) Copyright 1998-2008, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 10 Jan 1998 | Original version.
!   2.0   | 08 Jul 1999 | opt_Vrb.
!   2.1   | 06 Oct 1999 | -rh.
!   2.2   | 17 Dec 2000 | -fz.
!   2.3   | 21 Mar 2001 | -sr.
!   2.4   | 15 May 2001 | -ab.
!   2.5   | 16 May 2001 | -as.
!   2.6   | 13 Oct 2001 | -mps.
!   2.7   | 27 Jan 2002 | -ngo.
!   3.0   | 06 Aug 2002 | -phantom.
!   3.1   | 10 Aug 2002 | -dyn.
!   3.2   | 14 Aug 2002 | -xls.
!   3.3   | 11 Nov 2002 | -orb.
!   3.4   | 21 May 2003 | -freq.
!   3.5   | 06 Jul 2003 | -sd.
!   3.6   | 20 Aug 2003 | Verbosity level.
!   3.7   | 28 Apr 2004 | -turb.
!   3.8   | 20 May 2004 | -afreq.
!   3.9   | 25 May 2004 | -dbg.
!   4.0   | 28 Sen 2004 | Initialization of Stat.
!   4.1   | 07 May 2006 | -mpsf, -ecef, -ow.
!   4.2   | 27 Nov 2006 | -out.
!   4.3   | 23 Feb 2007 | -ts.
!   4.4   | 24 Feb 2007 | -fwp.
!   4.5   | 27 Jun 2007 | -ncep.
!   5.0   | 26 May 2008 | Opts, -o.
!   5.1   | 16 Jun 2008 | -key=y/n, -msis.
!   5.2   | 24 Nov 2008 | -geoid.
!   5.3   | 09 Sep 2009 | -freq=acemax.
!   5.4   | 15 Sep 2009 | -hmin=xxxx.
!   5.5   | 01 Oct 2009 | -o=NC.
!   5.6   | 05 Nov 2009 | -hgo.
!----------------------------------------------------------
! Modules used:
!
Use Occ_GNSS, only: &
! Imported Parameters:
    Freq_GPS,      &
    Freq_GAL,      &
    Freq_GLO,      &
    Freq_ACESTD,   &
    Freq_ACEOPT,   &
    Freq_ACEMAX
!
Use Wave_Defaults, only: &
! Imported Parameters:
    dfl_Hmax,      &
    dfl_DYN,       &
    dfl_DX,        &
    dfl_FZ,        &
    dfl_SR,        &
    dfl_NGO
!
!Use Comline, only: &
!!Use f90_unix, only: &
!! Imported Routines:
!    IArgC,   &
!    GetArg
!
Use Strings, only: &
! Imported Routines:
    Split_Arg,     &
    String_to_Log
!
Use Debug, only: &
! Imported Scalars:
    IDbg
!
Use Earth, only: &
! Imported Parameters:
    shape_Sphere,    &
    shape_Ellips,    &
    shape_Geoid
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Output arguments:
!
Type(t_Opts)       :: &
   Opts      ! Command line options
!
Logical            :: &
   Err            ! Command line error status
!----------------------------------------------------------
! Local Scalars:
!
Integer            :: NArg   ! Number of command line arguments
Integer            :: IArg   ! Command line argument index
Character(Len=255) :: Arg    ! Command line argument
Integer            :: NFreq  ! Number of frequencies
Integer            :: Stat   ! Read error code
Integer            :: N_Sub  ! Number of suboptions
Integer            :: ISub   ! Suboption index
!
! Local Arrays:
!
Character(Len=255) :: &
   Sub_Opt(16)               ! Suboptions
!----------------------------------------------------------


!----------------------------------------------------------
! 1. SETTING VALUES BY DEFAULT
!----------------------------------------------------------

Err               = .False.
Opts%GPS_Name     = ''
Opts%ORB_Name     = ''
Opts%Atm_Type     = atm_Undefined
Opts%HG_Mode      = hg_3d
Opts%Atm_Name(:)  = ''
Opts%N_Files      = 0
Opts%TRB_Name     = ''
Opts%Earth        = shape_Ellips
Opts%Geoid_Path   = ''
Opts%MSIS_Path    = ''
Opts%OUT_Name     = ''
Opts%Hmax         = dfl_Hmax
Opts%Hmin         = -999.0
Opts%DYN          = dfl_DYN
Opts%DX           = dfl_DX
Opts%XLS          = 0.0_Double
Opts%opt_LS       = .False.
Opts%FZ           = dfl_FZ
Opts%SR           = dfl_SR
Opts%opt_TS       = .False.
Opts%opt_ECEF     = .False.
Opts%opt_WO       = .True.
Opts%opt_GO       = .True.
Opts%opt_RH       = .True.
Opts%opt_AB       = .False.
Opts%opt_AS       = .False.
Opts%opt_FIO      = .True.
Opts%SDim         = 2
Opts%HGO          = 120.0
Opts%NGO          = dfl_NGO
Opts%FWP          = -1.0
Opts%Vrb          = 3
Opts%Out          = Script_Output
IDbg              = 0

NFreq         = Size(Freq_GPS)
Allocate(Opts%Freq(NFreq))
Opts%Freq(:)  = Freq_GPS(:)
Nullify(Opts%AFreq)


!----------------------------------------------------------
! 2. PARSING COMMAND LINE
!----------------------------------------------------------

!--- 2.1. Getting number of arguments

NArg = IArgC()


print*,'in Get_Comline'

!--- 2.2. Parsing arguments

GetParameters: Do IArg = 1, NArg
   Stat = 0
   Call GetArg(IArg, Arg)

   !--- 2.2.1. Getting key arguments

   If (Arg(:1) == '-') then

      !--- 2.2.1.1. Option -atm

      If (Arg(1:5) == '-atm=') then
         Call Split_Arg(Arg, Sub_Opt, N_Sub)
         If (N_Sub > 12) then
            Write (*,'(A/2A/)')  &
               'TOO MANY SUBPARAMETERES: ',  &
               '   ', Trim(Arg)
            Err = .True.
         End If
         If (.not. Err) then
            Select Case (Trim(Sub_Opt(1)))
               Case ('echam')
                  Opts%Atm_Type = atm_ECHAM
               Case ('ncep')
                  Opts%Atm_Type = atm_NCEP
               Case ('phantom')
                  Opts%Atm_Type = atm_Phantom
               Case ('mpas')
                  Opts%Atm_Type = atm_MPAS
               Case Default
                  Write (*,'(3A/2A/)')  &
                    'UNSUPPORTED SUBOPTION "',           &
                    Trim(Sub_Opt(1)) ,'" IN PARAMETER: ',   &
                    '   ', Trim(Arg)
                  Err = .True.
            End Select
         End If
         If ((.not. Err) .and. (N_Sub > 1)) then
            Select Case (Trim(Sub_Opt(2)))
               Case ('1d')
                  Opts%HG_Mode = hg_1d
               Case ('', '3d')
                  Opts%HG_Mode = hg_3d
               Case Default
                  Write (*,'(3A/2A/)')  &
                    'UNSUPPORTED SUBOPTION "',           &
                    Trim(Sub_Opt(2)) ,'" IN PARAMETER: ',   &
                    '   ', Trim(Arg)
                  Err = .True.
            End Select
         End If
!         print*,'sub_opt=',sub_opt
!         print*,'sub_opt(2)=',sub_opt(2)
         If ((.not. Err) .and. (N_Sub > 2)) then
            Opts%Atm_Name(1:N_Sub-2) = Sub_Opt(3:N_Sub)
            Opts%N_Files = N_Sub-2
         End If
      

      !--- 2.2.1.2. Option -turb

      Else If (Arg(1:6) == '-turb=') then
         Read(Arg(7:),'(A)',IOStat=Stat) Opts%TRB_Name

      !--- 2.2.1.3. Option -earth

      Else If (Arg(1:7) == '-earth=') then
         Call Split_Arg(Arg, Sub_Opt, N_Sub)
         If (N_Sub < 1 .or. N_Sub > 2) then
            Write (*,'(A/2A/)')  &
               'WRONG NUMBER OF SUBOPTIONS: ',  &
               '   ', Trim(Arg)
            Err = .True.
         End If
         Select Case (Trim(Sub_Opt(1)))
            Case ('sphere')
               Opts%Earth = shape_Sphere
            Case ('ellips')
               Opts%Earth = shape_Ellips
            Case ('geoid')
               Opts%Earth = shape_Geoid
            Case Default
               Write (*,'(3A/2A/)')  &
                 'UNSUPPORTED SUBOPTION "',           &
                 Trim(Sub_Opt(1)) ,'" IN PARAMETER: ',   &
                 '   ', Trim(Arg)
               Err = .True.
         End Select
         If (N_Sub > 1) then
            Opts%Geoid_Path = Sub_Opt(2)
         End If

      !--- 2.2.1.4. Option -msis

      Else If (Arg(1:6) == '-msis=') then
         Opts%MSIS_Path = Arg(7:)

      !--- 2.2.1.5. Option -gps

      Else If (Arg(1:5) == '-gps=') then
         Read(Arg(6:),'(A)',IOStat=Stat) Opts%GPS_Name

      !--- 2.2.1.6. Option -orb

      Else If (Arg(1:5) == '-orb=') then
         Read(Arg(6:),'(A)',IOStat=Stat) Opts%ORB_Name

      !--- 2.2.1.7. Option -out

      Else If (Arg(1:5) == '-out=') then
         Read(Arg(6:),'(A)',IOStat=Stat) Opts%OUT_Name

      !--- 2.2.1.8. Option -freq

      Else If (Arg(1:6) == '-freq=') then
         Deallocate(Opts%Freq)
         Call Get_Freq(Arg(7:), Opts%Freq, Err, Stat)

      !--- 2.2.1.9. Option -afreq

      Else If (Arg(1:7) == '-afreq=') then
         Call Get_Freq(Arg(8:), Opts%AFreq, Err, Stat)

      !--- 2.2.1.10. Option -hmax

      Else If (Arg(1:6) == '-hmax=') then
         Read(Arg(7:),*,IOStat=Stat) Opts%Hmax

      !--- 2.2.1.11. Option -hmin

      Else If (Arg(1:6) == '-hmin=') then
         Read(Arg(7:),*,IOStat=Stat) Opts%Hmin

      !--- 2.2.1.12. Option -dyn

      Else If (Arg(1:5) == '-dyn=') then
         Read(Arg(6:),*,IOStat=Stat) Opts%DYN

      !--- 2.2.1.13. Option -dx

      Else If (Arg(1:4) == '-dx=') then
         Read(Arg(5:),*,IOStat=Stat) Opts%DX

      !--- 2.2.1.14. Option -xls

      Else If (Arg(1:5) == '-xls=') then
         Read(Arg(6:),*,IOStat=Stat) Opts%XLS
         Opts%opt_LS = .True.

      !--- 2.2.1.15. Option -fz

      Else If (Arg(1:4) == '-fz=') then
         Read(Arg(5:),*,IOStat=Stat) Opts%FZ
         Opts%FZ = Pi*Opts%FZ

      !--- 2.2.1.16. Option -sr

      Else If (Arg(1:4) == '-sr=') then
         Read(Arg(5:),*,IOStat=Stat) Opts%SR

      !--- 2.2.1.17. Option -ts

      Else If (Arg(1:4) == '-ts=') then
         Call String_to_Log &
           (Arg(5:),      & ! <-- String value
            Opts%opt_TS,  & ! --> Logical value
            Stat)           ! --> Error status

      !--- 2.2.1.18. Option -ecef

      Else If (Arg(1:6) == '-ecef=') then
         Call String_to_Log &
           (Arg(7:),        & ! <-- String value
            Opts%opt_ECEF,  & ! --> Logical value
            Stat)             ! --> Error status

      !--- 2.2.1.19. Option -wo

      Else If (Arg(1:4) == '-wo=') then
         Call String_to_Log &
           (Arg(5:),          & ! <-- String value
            Opts%opt_WO,  & ! --> Logical value
            Stat)           ! --> Error status

      !--- 2.2.1.20. Option -go

      Else If (Arg(1:4) == '-go=') then
         Call String_to_Log &
           (Arg(5:),      & ! <-- String value
            Opts%opt_GO,  & ! --> Logical value
            Stat)           ! --> Error status

      !--- 2.2.1.21. Option -rh

      Else If (Arg(1:4) == '-rh=') then
         Call String_to_Log &
           (Arg(5:),      & ! <-- String value
            Opts%opt_RH,  & ! --> Logical value
            Stat)           ! --> Error status

      !--- 2.2.1.22. Option -ab

      Else If (Arg(1:4) == '-ab=') then
         Call String_to_Log &
           (Arg(5:),      & ! <-- String value
            Opts%opt_AB,  & ! --> Logical value
            Stat)           ! --> Error status

      !--- 2.2.1.23. Option -as

      Else If (Arg(1:4) == '-as=') then
         Call String_to_Log &
           (Arg(5:),      & ! <-- String value
            Opts%opt_AS,  & ! --> Logical value
            Stat)           ! --> Error status

      !--- 2.2.1.24. Option -fio

      Else If (Arg(1:5) == '-fio=') then
         Call String_to_Log &
           (Arg(6:),       & ! <-- String value
            Opts%opt_FIO,  & ! --> Logical value
            Stat)            ! --> Error status

      !--- 2.2.1.25. Option -sd

      Else If (Arg(1:4) == '-sd=') then
         Read(Arg(5:),*,IOStat=Stat) Opts%SDim

      !--- 2.2.1.26. Option -hgo

      Else If (Arg(1:5) == '-hgo=') then
         Read(Arg(6:),*,IOStat=Stat) Opts%HGO

      !--- 2.2.1.27. Option -ngo

      Else If (Arg(1:5) == '-ngo=') then
         Read(Arg(6:),*,IOStat=Stat) Opts%NGO

      !--- 2.2.1.28. Option -fwp

      Else If (Arg(1:5) == '-fwp=') then
         Read(Arg(6:),*,IOStat=Stat) Opts%FWP

      !--- 2.2.1.29. Option -b

      Else If (Arg(1:7) == '-b') then
         Opts%Vrb = 1

      !--- 2.2.1.30. Option -v

      Else If (Arg(1:3) == '-v=') then
         Read(Arg(4:),*,IOStat=Stat) Opts%Vrb

      !--- 2.2.1.31. Option -dbg

      Else If (Arg(1:5) == '-dbg=') then
         Read(Arg(6:),*,IOStat=Stat) IDbg

      !--- 2.2.1.32. Option -o

      Else If (Arg(1:3) == '-o=') then
         Call Split_Arg(Arg, Sub_Opt, N_Sub)
         If (N_Sub > Size(Sub_Opt)) then
            Write (*,'(A/2A/)')  &
               'TOO MANY EXTENSIONS: ',  &
               '   ', Trim(Arg)
            Err = .True.
         End If
         Do ISub=1,N_Sub
            Select Case (Trim(Sub_Opt(ISub)))
               !
               ! --- Presets
               !
               Case ('none')
                  Opts%Out      = No_Output
               Case ('all')
                  Opts%Out      = All_Output
               Case ('min')
                  Opts%Out      = Min_Output
               Case ('script')
                  Opts%Out      = Script_Output
               !
               ! --- Turning extensions on
               !
               Case ('WOP')
                  Opts%Out%WOP  = .True.
               Case ('WCA')
                  Opts%Out%WCA  = .True.
               Case ('WCE')
                  Opts%Out%WCE  = .True.
               Case ('CAU')
                  Opts%Out%CAU  = .True.
               Case ('GOP')
                  Opts%Out%GOP  = .True.
               Case ('GRP')
                  Opts%Out%GRP  = .True.
               Case ('XBP')
                  Opts%Out%XBP  = .True.
               Case ('GOZ')
                  Opts%Out%GOZ  = .True.
               Case ('GON')
                  Opts%Out%GON  = .True.
               Case ('GOT')
                  Opts%Out%GOT  = .True.
               Case ('MOT')
                  Opts%Out%MOT  = .True.
               Case ('MOQ')
                  Opts%Out%MOQ  = .True.
               Case ('MON')
                  Opts%Out%MON  = .True.
               Case ('MOA')
                  Opts%Out%MOA  = .True.
               Case ('MTR')
                  Opts%Out%MTR  = .True.
               Case ('MTD')
                  Opts%Out%MTD  = .True.
               Case ('NC')
                  Opts%Out%NC   = .True.
               !
               ! --- Turning extensions off
               !
               Case ('WOP-')
                  Opts%Out%WOP  = .False.
               Case ('WCA-')
                  Opts%Out%WCA  = .False.
               Case ('WCE-')
                  Opts%Out%WCE  = .False.
               Case ('CAU-')
                  Opts%Out%CAU  = .False.
               Case ('GOP-')
                  Opts%Out%GOP  = .False.
               Case ('GRP-')
                  Opts%Out%GRP  = .False.
               Case ('XBP-')
                  Opts%Out%XBP  = .False.
               Case ('GOZ-')
                  Opts%Out%GOZ  = .False.
               Case ('GON-')
                  Opts%Out%GON  = .False.
               Case ('GOT-')
                  Opts%Out%GOT  = .False.
               Case ('MOT-')
                  Opts%Out%MOT  = .False.
               Case ('MOQ-')
                  Opts%Out%MOQ  = .False.
               Case ('MON-')
                  Opts%Out%MON  = .False.
               Case ('MOA-')
                  Opts%Out%MOA  = .False.
               Case ('MTR-')
                  Opts%Out%MTR  = .False.
               Case ('MTD-')
                  Opts%Out%MTD  = .False.
               Case ('NC-')
                  Opts%Out%NC   = .False.
               !
               ! --- Error detection
               !
               Case Default
                  Write (*,'(3A/2A/)')  &
                    'UNSUPPORTED SUBOPTION "',           &
                    Trim(Sub_Opt(ISub)) ,'" IN PARAMETER: ',   &
                    '   ', Trim(Arg)
                  Err = .True.
           End Select
         End Do
         print*,'WCA,GOT=',Opts%Out%WCA,Opts%Out%GOT

      !--- 2.2.1.33. Error status for unsupported option

      Else
         Write (*,'(A/2A/)')  &
            'UNSUPPORTED OPTION: ',  &
            '   ', Trim(Arg)
         Err = .True.
      End If

      !--- 2.2.1.34. Error status for wrong value in parameter

      If (Stat /= 0) then
         Write (*,'(A/2A/)')  &
            'BAD VALUE IN PARAMETER: ',  &
            '   ', Trim(Arg)
         Err = .True.
      End If

   !--- 2.2.2. Error status for illegal parameter

   Else
      Write (*,'(A/2A/)')  &
         'ILLEGAL PARAMETER: ',  &
         '   ', Trim(Arg)
      Err = .True.
   End If

End Do GetParameters


!--- 2.3. Processing frequencies

If (.not. Associated(Opts%AFreq)) then
   Allocate(Opts%AFreq(Size(Opts%Freq)))
   Opts%AFreq(:) = Opts%Freq(:)
Else If (Size(Opts%AFreq) /= Size(Opts%Freq)) then
   Write (*,'(A)')  &
      'DIFFERENT NUMBERS OF FREQUENCIES'
   Err = .True.
End If


!--- 2.4. Error staturs for missing arguments

If (NArg == 0) then
   Err = .True.
End If


Contains


!----------------------------------------------------------
Subroutine Get_Freq(Arg, Freq, Err, Stat)
!
! Parsing command line parameter for frequencies.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
Character(Len=*), Intent(In) :: &
   Arg         ! Command line argument
Real(Double), Pointer        :: &
   Freq(:)     ! Frequencies [Hz]
Logical, Intent(Out)         :: &
   Err         ! Error presence flag
Integer, Intent(Out)         :: &
   Stat        ! Parsing error status
!----------------------------------------------------------
! Local Parameters:
!
Integer, Parameter :: &
   NMaxFreq = 10             ! Maximum number of frequency channels
!
! Local Scalars:
!
Integer            :: i      ! Subargument index
!
! Local Arrays:
!
Character(Len=255) :: &
   FN(NMaxFreq)              ! Frequency subarguments
!----------------------------------------------------------


Err = .False.


Select Case (Arg)
   Case('gps')
      NFreq = Size(Freq_GPS)
      Allocate(Freq(NFreq))
      Freq(:) = Freq_GPS(:)
   Case('gal')
      NFreq = Size(Freq_GAL)
      Allocate(Freq(NFreq))
      Freq(:) = Freq_GAL(:)
   Case('glo')
      NFreq = Size(Freq_GLO)
      Allocate(Freq(NFreq))
      Freq(:) = Freq_GLO(:)
   Case('acestd')
      NFreq = Size(Freq_ACESTD)
      Allocate(Freq(NFreq))
      Freq(:) = Freq_ACESTD(:)
   Case('aceopt')
      NFreq = Size(Freq_ACEOPT)
      Allocate(Freq(NFreq))
      Freq(:) = Freq_ACEOPT(:)
   Case('acemax')
      NFreq = Size(Freq_ACEMAX)
      Allocate(Freq(NFreq))
      Freq(:) = Freq_ACEMAX(:)
   Case Default
      Call Split_Arg(Arg, FN, NFreq)
      If (NFreq > NMaxFreq) then
         Write (*,'(A,I2,A/2A/)')  &
            'TOO MANY FREQUENCIES (>', NMaxFreq,'): ',  &
            '   ', Trim(Arg)
         Err = .True.
      End If
      Allocate(Freq(NFreq))
      Do i=1,NFreq
         Read(FN(i),*,IOStat=Stat) Freq(i)
         If (Stat /= 0) then
            Err = .True.
            Exit
         End If
      End Do
      Freq(:) = 1d9*Freq(:)
End Select


End Subroutine Get_Freq





End Subroutine Get_Comline



!==========================================================
Subroutine Help()
!
! Displaying help screen for Wave.
!----------------------------------------------------------
! (C) Copyright 1998-2008, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 22 Feb 1999 | Original version.
!   1.1   | 08 Jul 1999 | -b.
!   1.2   | 06 Oct 1999 | -rh.
!   1.3   | 17 Dec 2000 | -fz.
!   1.4   | 21 Mar 2001 | -sr.
!   1.5   | 15 May 2001 | -ab.
!   1.6   | 16 May 2001 | -as.
!   1.7   | 13 Oct 2001 | -mps.
!   1.8   | 27 Jan 2002 | -ngo.
!   2.0   | 06 Aug 2002 | -phantom.
!   2.1   | 10 Aug 2002 | -dyn.
!   2.2   | 14 Aug 2002 | -xls.
!   2.3   | 11 Nov 2002 | -orb.
!   2.4   | 21 May 2003 | -freq.
!   2.5   | 06 Jul 2003 | -sd.
!   2.6   | 28 Apr 2004 | -turb.
!   3.8   | 20 May 2004 | -afreq.
!   3.9   | 25 May 2004 | -dbg.
!   4.0   | 05 May 2006 | -mpsf.
!   4.1   | 07 May 2006 | -ecef, -ow.
!   4.2   | 27 Nov 2006 | -out.
!   4.3   | 23 Feb 2007 | -ts.
!   4.4   | 24 Feb 2007 | -fwp.
!   4.5   | 27 Jun 2007 | -ncep.
!   5.0   | 25 May 2008 | -o.
!   5.2   | 24 Nov 2008 | -geoid.
!   5.3   | 09 Sep 2009 | -freq=acemax.
!   5.4   | 15 Sep 2009 | -hmin=xxxx.
!   5.5   | 01 Oct 2009 | -o=NC.
!   5.6   | 05 Nov 2009 | -hgo.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------


Write (*,'(A)')     &
   'Usage: ',       &
   'wave <Options>',                                                          &
   'Available options: ',                                                     &
   '   -atm=type,[hg],ss,... - Atmospheric model specification:',             &
   '                       type:',                                            &
   '                          echam   - ECHAM girb file',                     &
   '                          ncep    - NCEP grib file',                      &
   '                          phantom - analytical model',                    &
   '                       hg - horizontal gradients:',                       &
   '                          1d      - no horizontal gradients',             &
   '                          3d      - full horizontal gradients (default)', &
   '                       ss,ss,...  - pathname(s)',                         &
   '   -turb=sssss       - Turbulence parameter file',                        &
   '   -earth=sss[,ppp]  - sss specifies Earth shape:',                       &
   '                          sphere - spherical',                            &
   '                          ellips - reference elliposid WGS 1984',         &
   '                          geoid  - geoid',                                &
   '                             ppp - optional path to coefficient file',    &
   '   -gps=ssss         - RO data file patname',                             &
   '   -orb=ssss         - Orbit data file pathname',                         &
   '   -msis=ssss        - Path to MSIS coefficient files',                   &
   '   -out=ssss         - Template for output names',                        &
   '   -freq=ffff        - Frequency-channels specification:',                &
   '         gps         -     GPS',                                          &
   '         gal         -     GALILEO',                                      &
   '         glo         -     GLONASS',                                      &
   '         acestd      -     ACE+ standard',                                &
   '         aceopt      -     ACE+ optional',                                &
   '         acemax      -     ACE+ maximum',                                 &
   '         xxx,xxx...  -     up to 10 user-specified [GHz]',                &
   '   -afreq=ffff       - Frequency-channels for absorption',                &
   '   -hmax=xxxx        - maximum heigh (default 100)',                      &
   '   -hmin=xxxx        - minimum height (default -999)',                    &
   '   -dyn=xxxx         - vertical resolution of refractivity field',        &
   '   -dx=xxxx          - step between phase screens (default 100)',         &
   '   -xls=xxxx         - last phase screen position',                       &
   '   -fz=xxxx          - Fresnel zone size [Pi rad] (default 10)',          &
   '   -sr=xxxx          - sampling rate [Hz] (default 50)'
Write (*,'(A)')      &
   '   -ts=y/n           - time scaling (default - deactivated)',             &
   '   -ecef=y/n         - all coordinates in ECEF (default J2000)',          &
   '   -wo=y/n           - wave optics (default yes)',                        &
   '   -go=y/n           - geometrical optics (default yes)',                 &
   '   -rh=y/n           - output of ray heights (default yes)',              &
   '   -ab=y/n           - modeling of absorption (default no)',              &
   '   -as=y/n           - asymptotic solution (default no)',                 &
   '   -fio=y/n          - use FIO propagation to LEO orbit (default yes)',   &
   '   -sd=xxxx          - spatial dimension (default 2)',                    &
   '   -hgo=xxxx         - maximum height for GO (default 120)',              &
   '   -ngo=nnnn         - number of GO rays (default 1500)',                 &
   '   -fwp=xxxx         - filter width for DPGPS/DP [km] (default -1)',      &
   '   -v=xxxx           - verbosity level:',                                 &
   '                     -    0 - silent mode',                               &
   '                     -    1 - progress messages',                         &
   '                     -    2 - static messages',                           &
   '                     -    3 - dynamic messages',                          &
   '   -b                - batch mode (equivalent to -v=1)',                  &
   '   -dbg=xxxx         - debug mode',                                       &
   '   -o=               - output files extensions:',                         &
   '    XX[-],...             XX     - turnes specific extension on',         &
   '                          XX-    - turnes specific extension off',        &
   '                          none   - no output',                            &
   '                          all    - output of all extensions',             &
   '                          min    - minimum ouput data set',               &
   '                          script - output data set for use',              &
   '                                   with genvda script',                   &
   '                          Extensions:',                                   &
   '                          WOP WCA WCE CAU GOP GRP XBP GOZ',               &
   '                          GON GOT MOT MOQ MON MOA MTR MTD'
Write (*,'()')


End Subroutine Help



End Module Wave_Comline


