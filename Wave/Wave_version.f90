!
Module Wave_Version
!
! Version and date information for program Wave.
!----------------------------------------------------------
! (C) Copyright 2003, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 18 Dec 2003 | Original version
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Public Parameters:
!
Character(Len=*), Parameter :: &
   Version = '10.7.220',   & ! Program version
   Date    = '05 Nov 2009'   ! Program date
!----------------------------------------------------------
!
Contains


!==========================================================
Subroutine Display_Version()
!
! Displaying version number and program info.
!----------------------------------------------------------
! (C) Copyright 2003, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 18 Dec 2003 | Original version.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------


!----------------------------------------------------------
! 1. FIRST SECTION
!----------------------------------------------------------

Write (*,'(/5A/A)') &
   'Wave. Version ', Version, ' (c) M. E. Gorbunov. ', Date, '.',  &
   'Wave optics propagation simulator.'
Write (*,'()')


End Subroutine Display_Version



End Module Wave_Version


