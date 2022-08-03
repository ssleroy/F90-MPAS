!
Module Atmosphere_Rays
!
! Tracing of geometric optical rays for atmospheric
! refractivity field model.
!----------------------------------------------------------
! (C) Copyright 2003, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 01 Aug 2003 | Original version.
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Parameters:
    Double, dtr
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
!
Contains


!==========================================================
Subroutine Ray_Trace &
  (XT,        & ! <-- Transmitter position
   UT,        & ! <-- Transmitter ray direction
   XR,        & ! <~~ Receiver position
   DS,        & ! <-- Integration step parameter
   XN,        & ! --> Ray point nearest to receiver
   UN,        & ! --> Ray direction at XN
   At,        & ! ~~> Attenuation [Int(Imag(N))]
   Hper,      & ! ~~> Perigee altitude [km]
   Stat)        ! --> Error status
!
! Tracing ray with given initial conditions.
!----------------------------------------------------------
! Method:
!   Numerical intergration of the differential equation
!   of rays in cartesian coordinates with given intial
!   conditions, stopping at the point nearest to the
!   receiver.
!----------------------------------------------------------
! (C) Copyright 1999-2005, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 18 Feb 1999 | Original version.
!   2.0   | 22 Feb 1999 | XR optional.
!   3.0   | 06 Aug 2002 | Reflection and attenuation.
!   3.1   | 10 Aug 2002 | Check for ray trap.
!   3.2   | 06 Mar 2004 | Adaptive integration step.
!   3.3   | 27 Aug 2004 | Dynamic upper estimate
!         |             | of integration step.
!   4.0   | 08 Nov 2005 | Hper.
!----------------------------------------------------------
! Modules used:
!
Use Coordinates, only: &
! Imported Type Definitions:
    Cartesian,     &
! Imported Routines:
    Vector_Normed, &
    Vector_Norm,   &
! Imported Operators:
    Operator(*),   &
    Operator(+),    &
    Operator(-)
!
Use Earth, only: &
! Imported Type Definitions:
    Geodetic,  &
! Imported Parameters:
    R_Earth,   &
    H_atm,     &
! Imported Routines:
    Geod_from_Cart
!
Use Atmosphere, only: &
! Imported Routines:
    Atmosphere_NGradN
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type(Cartesian), Intent(In) :: &
   XT       ! Transmitter position
            ! Transmitter must be outside atmosphere:
            ! |XT| > R_Earth + H_atm
!
Type(Cartesian), Intent(In) :: &
   UT       ! Transmitter ray direction
!
Type(Cartesian), Optional, Intent(In) :: &
   XR       ! Receiver position
            ! Receiver must be outside atmosphere:
            ! |XR| > R_Earth + H_atm
!
Real(Double), Intent(In) :: &
   DS       ! Integration step parameter
!
! Output arguments:
!
Type(Cartesian), Intent(Out) :: &
   XN       ! Ray point nearest to receiver
!
Type(Cartesian), Intent(Out) :: &
   UN       ! Ray direction at XN
!
Real(Double), Optional, Intent(Out) :: &
   At(:)    ! Attenuation [Int(Imag(N))]
            ! Amplitude factor = exp(k*At)
!
Real(Double), Optional, Intent(Out) :: &
   Hper     ! Perigee altitude [km]
!
Integer, Intent(Out) :: &
   Stat     ! Error status:
            !    0 - no error
            !    1 - ray reflected
            !    2 - transmitter inside atmosphere
            !    3 - receiver inside atmosphere
            !    4 - ray trapped
!----------------------------------------------------------
! Local Parameters:
!
Real(Double), Parameter :: &
   R_atm = R_Earth + H_atm   ! Radius of atmosphere
!
! Local Scalars:
!
Type(Cartesian) :: X0     ! Ray point coordinate before step
Type(Cartesian) :: U0     ! Ray direction before step
Type(Cartesian) :: X      ! Ray point coordinate after step
Type(Cartesian) :: U      ! Ray direction after step
Type(Cartesian) :: V      ! Local vertical vector
Type(Geodetic)  :: G      ! Geodetic coordinates
Real(Double)    :: SMin   ! Integration step without collision
Real(Double)    :: SMax   ! Integration step with collision
Real(Double)    :: SE     ! Integration step
Real(Double)    :: p      ! Starting ray impact parameter
Real(Double)    :: Sa     ! Distance from transmitter to atmosphere
Real(Double)    :: Sr     ! Distance to reciever
Real(Double)    :: DYN    ! Vertical scale
Real(Double)    :: Sup    ! Upper estimate of integration step [km]
Real(Double)    :: S      ! Current integration step
Integer         :: NC     ! Number of frequency channels for
                          ! computation of attenuation
Integer         :: PStat  ! Status of refractivity computation
Logical         :: Refl   ! Reflection status
Integer         :: NRefl  ! Number of reflections
Real(Double)    :: Path   ! Ray path
!
! Local Arrays:
!
Type(Cartesian) :: DX(4)  ! Runge-Kutta intermediate dX/dt
Type(Cartesian) :: DU(4)  ! Runge-Kutta intermediate dU/dt
Real(Double)    :: NP(4)  ! Interpolated N
Real(Double), Allocatable :: &
   At0(:),              & ! Attenuation before step
   NIP(:,:)               ! Im(N)
!----------------------------------------------------------


!----------------------------------------------------------
! 0. INITIALIZATION
!----------------------------------------------------------


!--- 0.1. Parameter check

If (Vector_Norm(XT) < R_atm) then
   Stat = 2
   Return
End If

If (Present(XR)) then
   If (Vector_Norm(XR) < R_atm) then
      Stat = 3
      Return
   End If
End If


!--- 0.2. Status initialization

Stat = 0


!--- 0.3. Array initialization

If (Present(At)) then
   At(:) = 0
   NC    = Size(At)
   Allocate(NIP(4,NC))
   Allocate(At0(NC))
End If


!--- 0.4. Number of reflections

NRefl = 0


!----------------------------------------------------------
! 1. FIRST STEP FROM TRANSMITTER TO ATMOSPHERE
!----------------------------------------------------------

U  = Vector_Normed(UT)
p  = Sqrt((XT*XT) - (XT*U)**2)
Sa = -(XT*U) - Sqrt(Dim(R_atm**2, p**2))
X  = XT + Sa*U


!----------------------------------------------------------
! 2. RAY INTEGRATION
!----------------------------------------------------------

Path = 0

If (Present(Hper)) then
   G    = Geod_from_Cart(X)
   Hper = G%H
End If

Ray_Integrate: Do


   !--- 2.1. Checking for ray outgoing from atmosphere

   If (((X*U) > 0) .and. (Vector_Norm(X) > R_atm)) then
      Exit Ray_Integrate
   End If


   !--- 2.2. Estimate of integration step

   G = Geod_from_Cart(X)

   If (Present(Hper)) then
      Hper = Min(G%H, Hper)
   End If

   If (G%H > 35.0) then
      DYN = 0.5
   Else
      DYN = DS
   End If

   Sup = 50*DYN
   S   = DYN/(Abs(Vector_Normed(U)*Vector_Normed(X)) + DYN/Sup)


   !--- 2.2. Making a step of Runge-Kutta ray integration

   X0   = X
   U0   = U
   SE   = S
   SMin = 0
   SMax = S
   Refl = .False.

   If (Present(At)) then
      At0(:) = At(:)
   End If

   Step: Do

      !--- 2.2.1. Integration step

      If (.not. Present(At)) then
         DX(1) = U0
         Call Atmosphere_NGradN(X0,                DU(1), NP(1), PStat)
         DX(2) = U0 + DU(1)*(SE/2)
         Call Atmosphere_NGradN(X0 + DX(1)*(SE/2), DU(2), NP(2), PStat)
         DX(3) = U0 + DU(2)*(SE/2)
         Call Atmosphere_NGradN(X0 + DX(2)*(SE/2), DU(3), NP(3), PStat)
         DX(4) = U0 + DU(3)*SE
         Call Atmosphere_NGradN(X0 + DX(3)*SE,     DU(4), NP(4), PStat)
      Else
         DX(1) = U0
         Call Atmosphere_NGradN(X0,                DU(1), NP(1), PStat, NIP(1,:))
         DX(2) = U0 + DU(1)*(SE/2)
         Call Atmosphere_NGradN(X0 + DX(1)*(SE/2), DU(2), NP(2), PStat, NIP(2,:))
         DX(3) = U0 + DU(2)*(SE/2)
         Call Atmosphere_NGradN(X0 + DX(2)*(SE/2), DU(3), NP(3), PStat, NIP(3,:))
         DX(4) = U0 + DU(3)*SE
         Call Atmosphere_NGradN(X0 + DX(3)*SE,     DU(4), NP(4), PStat, NIP(4,:))
         At(:) = At0(:) +                                            &
                 ((1+NP(1))*NIP(1,:) + 2*(1+NP(2))*NIP(2,:) +       &
                  2*(1+NP(3))*NIP(3,:) + (1+NP(4))*NIP(4,:))*(SE/6)
      End If

      X = X0 + (DX(1) + 2*DX(2) + 2*DX(3) + DX(4))*(SE/6)
      U = U0 + (DU(1) + 2*DU(2) + 2*DU(3) + DU(4))*(SE/6)


      !--- 2.2.2. Renorming U

      U = Vector_Normed(U)*(1+NP(4))


      !--- 2.2.3. Check for reflection

      If (PStat /= 0) then
         Refl  = .True.
         Stat  = 1
         SMax  = SE
      Else
         SMin  = SE
      End If

      !--- 2.2.4. Determination of reflection point

      If (Abs(SMax - SMin) < 0.001) then

         !--- 2.2.4.1. Reflection

         If (Refl) then
            G      = Geod_from_Cart(X)
            V%X(1) = Cos(G%Phi*dtr)*Cos(G%Lambda*dtr)
            V%X(2) = Cos(G%Phi*dtr)*Sin(G%Lambda*dtr)
            V%X(3) = Sin(G%Phi*dtr)
            U      = U - 2*(U*V)*V
            NRefl  = NRefl + 1
         End If

         Path = Path + SE

         !--- 2.2.4.2. Check for ray trap

         If (Path > 6000) then
            Stat = 4
            Exit Ray_Integrate
         End If

         If (NRefl > 1) then
            Stat = 4
            Exit Ray_Integrate
         End If

         Exit Step

      Else

         !--- 2.2.4.3. Refining integration step

         SE = (SMax + SMin)/2

      End If

   End Do Step


End Do Ray_Integrate


!----------------------------------------------------------
! 3. LAST STEP FROM ATMOSPHERE TO RECEIVER
!----------------------------------------------------------

UN = Vector_Normed(U)
XN = X

If (Present(XR)) then
   Sr = (XR - X)*UN
   XN = XN + Sr*UN
End If


!----------------------------------------------------------
! 4. MEMORY DEALLOCATION
!----------------------------------------------------------

If (Present(At)) then
   Deallocate(NIP)
   Deallocate(At0)
End If


End Subroutine Ray_Trace



End Module Atmosphere_Rays


