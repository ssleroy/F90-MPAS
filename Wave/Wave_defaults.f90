!
Module Wave_defaults
!
! Definitions of constants used by Wave.
!----------------------------------------------------------
! (C) Copyright 1999-2009, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 22 Feb 1999 | Original version.
!   2.0   | 15 May 2001 | Absorption frequencies.
!   2.1   | 09 Sep 2009 | ACEMAX.
!   2.2   | 30 Oct 2009 | GNSS parameters moved to Occ_GNSS.
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Parameters:
    Double, Pi
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Public Declarations:
!
Public
!
! Public Parameters:
!
Integer, Parameter :: &
   dfl_Nmax = 5000   ! Maximum number of data by default.
!
Integer, Parameter :: &
   dfl_Hmax = 100    ! Maximum height by default.
!
Real(Double), Parameter :: &
   dfl_DYN  = 0.025  ! Minimum vertical scale of N by default
!
Real(Double), Parameter :: &
   dfl_DX   = 30.0   ! Step between phase screens by default.
!
Real(Double), Parameter :: &
   dfl_FZ   = 10*Pi  ! Fresnel zone size [Pi rad].
!
Real(Double), Parameter :: &
   dfl_SR   = 50     ! Sampling rate [Hz].
!
Integer, Parameter      :: &
   dfl_NGO = 1500    ! Number of geometric optical rays
!----------------------------------------------------------
!

End Module Wave_defaults

