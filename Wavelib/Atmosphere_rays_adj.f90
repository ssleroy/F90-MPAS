!
Module Atmosphere_rays_adj
!
! Adjoint version of ray-tracer for atmosphere model.
!----------------------------------------------------------
! (C) Copyright 2003, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 07 May 2003 | Original version.
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Parameters:
    Double, Pi
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
!
Contains


!==========================================================
Subroutine Ray_Trace_ZNZ0 &
  (XT,        & ! <-- Transmitter position
   UT,        & ! <-- Transmitter ray direction
   XR,        & ! <-- Receiver position
   DS,        & ! <-- Integration step parameter
   XN,        & ! --> Ray point nearest to receiver
   UN,        & ! --> Ray direction at XN
   ZN_Z0,     & ! --> d(XN,UN)/d(XT,UT)
   Hper,      & ! --> Ray perigee height
   Stat)        ! --> Error status
!
! Tracing ray with given initial conditions and calculation
! of derivative d(ZN)/d(Z0).
!----------------------------------------------------------
! Method:
!   Numerical intergration of the differential equation
!   of rays in cartesian coordinates with given intial
!   conditions, stopping at the point nearest to the
!   receiver.
!----------------------------------------------------------
! (C) Copyright 1999-2004, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   2.0   | 22 Feb 1999 | Basic non-adjoint version.
!   1.0   | 03 May 2000 | Reduced adjoint version.
!   2.0   | 07 May 2003 | Modification for Atmosphere.
!   2.1   | 06 Mar 2004 | Adaptive integration step.
!   3.3   | 27 Aug 2004 | Dynamic upper estimate
!         |             | of integration step.
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
    Operator(+),   &
    Operator(-),   &
    Operator(/),   &
    Operator(.xx.), &
    Assignment(=)
!
Use Matrix, only: &
! Imported Routines:
    Tensor_Product
!
Use Earth, only: &
! Imported Type Definitions:
    Geodetic,         &
! Imported Routines:
    Geod_from_Cart,   &
    R_Earth,          &
    H_atm
!
Use Atmosphere, only: &
! Imported Routines:
    Atmosphere_NGHN
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type(Cartesian), Intent(In) :: &
   XT         ! Transmitter position
              ! Transmitter must be outside atmosphere:
              ! |XT| > R_Earth + H_atm
!
Type(Cartesian), Intent(In) :: &
   UT         ! Transmitter ray direction
!
Type(Cartesian), Intent(In) :: &
   XR         ! Receiver position
              ! Receiver must be outside atmosphere:
              ! |XR| > R_Earth + H_atm
!
Real(Double), Intent(In) :: &
   DS         ! Integration step parameter
!
! Output arguments:
!
Type(Cartesian), Intent(Out) :: &
   XN         ! Ray point nearest to receiver
!
Type(Cartesian), Intent(Out) :: &
   UN         ! Ray direction at XN
!
Real(Double), Intent(Out)    :: &
   ZN_Z0(6,6) ! d(XN,UN)/d(XT,UT)
!
Real(Double), Intent(Out)    :: &
   Hper       ! Ray perigee height
!
Integer, Intent(Out) :: &
   Stat       ! Error status:
              !    0 - no error
              !    1 - ray collided with Earth
              !    2 - transmitter inside atmosphere
              !    3 - receiver inside atmosphere
              !    4 - ray captured
!----------------------------------------------------------
! Local Parameters:
!
Real(Double), Parameter :: &
   R_atm = R_Earth + H_atm    ! Radius of atmosphere
!
Real(Double), Parameter :: &
   E3(3,3) = Reshape((/1, 0, 0,  &
                       0, 1, 0,  &
                       0, 0, 1 /),  &
                     (/3,3/)),  & ! Unit 3x3 matrix
   E6(6,6) = Reshape((/1, 0, 0, 0, 0, 0,  &
                       0, 1, 0, 0, 0, 0,  &
                       0, 0, 1, 0, 0, 0,  &
                       0, 0, 0, 1, 0, 0,  &
                       0, 0, 0, 0, 1, 0,  &
                       0, 0, 0, 0, 0, 1 /),  &
                     (/6,6/))     ! Unit 6x6 matrix
!
! Local Scalars:
!
Type(Cartesian) :: X      ! Current ray point coordinates
Type(Geodetic)  :: G      ! Current ray point coordinates
Type(Cartesian) :: U      ! Current ray direction
Type(Cartesian) :: U1     ! Unnormed vector U
Real(Double)    :: p      ! Starting ray impact parameter
Real(Double)    :: Sa     ! Distance from transmitter to atmosphere
Real(Double)    :: DYN    ! Vertical scale
Real(Double)    :: Sup    ! Upper estimate of integration step [km]
Real(Double)    :: S      ! Current integration step
Real(Double)    :: NP     ! Interpolated N
Real(Double)    :: Sr     ! Distance to reciever
Real(Double)    :: St     ! Distance along ray
Integer         :: Mu     ! Integration substep index
!
! Local Arrays:
!
Type(Cartesian) :: &
    DX(4),    & ! Runge-Kutta intermediate dX/dt
    DU(4)       ! Runge-Kutta intermediate dU/dt
Real(Double)    :: &
    NHN(4,3,3)  ! Interpolated Grad x (1+N)Grad(N)
!
! --- Matrices for adjoint calculations
!
Real(Double) :: &
   R(6,6)                          ! Renorming matrix
Real(Double) :: &
   Sa_X(3),                      & ! d(Sa)/d(XT)
   Sa_U(3),                      & ! d(Sa)/d(U)
   Sr_X(3),                      & ! d(Sr)/d(X)
   Sr_U(3)                         ! d(Sr)/d(U)
Real(Double) :: &
   BM(4,6,6),                    & ! B^{mu} matrices
   B21(6,6), B32(6,6), B43(6,6), & ! Combinations of 2 BM
   B321(6,6), B432(6,6),         & ! Combinations of 3 BM
   B4321(6,6)                      ! Combinations of 4 BM
Real(Double) :: &
   B(6,6),                       & ! B_{n} matrix
   BN(6,6)                         ! Product of B_{n} matrices
!----------------------------------------------------------


!----------------------------------------------------------
! 0. INITIALIZATION
!----------------------------------------------------------

!--- 0.1. Parameter check

If (Vector_Norm(XT) < R_atm) then
   Stat = 2
   Return
End If

If (Vector_Norm(XR) < R_atm) then
   Stat = 3
   Return
End If

Stat = 0


!----------------------------------------------------------
! 1. FIRST STEP FROM TRANSMITTER TO ATMOSPHERE
!----------------------------------------------------------


!--- 1.1. First step

U    = Vector_Normed(UT)
p    = Sqrt((XT*XT) - (XT*U)**2)
Sa   = -(XT*U) - Sqrt(Dim(R_atm**2, p**2))
X    = XT + Sa*U


!--- 1.2. Calculation of height

G    = Geod_from_Cart(X)
Hper = G%H


!--- 1.3. Calculation of B_{1} matrix

R(:,:) = 0
R(1:3,1:3) = E3(:,:)
R(4:6,4:6) = &
   (E3(:,:) - Tensor_Product(UT%X(:),UT%X(:))/Sqrt(UT*UT))/ &
   Sqrt(UT*UT)

Sa_X(:) = -U%X(:) - &
          (-XT%X(:) + (XT*U)*U%X(:))/Sqrt(Dim(R_atm**2, p**2))
Sa_U(:) = -XT%X(:) - &
          ((XT*U)*XT%X(:))/Sqrt(Dim(R_atm**2, p**2))

B(1:3,1:3) = E3(:,:) + Tensor_Product(U%X(:),Sa_X(:))
B(1:3,4:6) = Sa*E3(:,:) + Tensor_Product(U%X(:),Sa_U(:))
B(4:6,4:6) = E3(:,:)
B(4:6,1:3) = 0

B(:,:) = MatMul(B(:,:),R(:,:))


!--- 1.4. Initialization of BN matrix

BN(:,:) = B(:,:)


!----------------------------------------------------------
! 2. RAY INTEGRATION
!----------------------------------------------------------


St = 0.0

Ray_Integrate: Do

   !--- 2.1. Checking for ray outgoing from atmosphere

   If (((X*U) > 0) .and. (Vector_Norm(X) > R_atm)) then
      Exit Ray_Integrate
   End If


   !--- 2.2. Estimate of integration step

   G = Geod_from_Cart(X)

   If (G%H > 35.0) then
      DYN = 0.5
   Else
      DYN = DS
   End If

   Sup = 50*DYN
   S   = DYN/(Abs(Vector_Normed(U)*Vector_Normed(X)) + DYN/Sup)


   !--- 2.3. Making a step of Runge-Kutta ray integration

   DX(1) = U
   Call Atmosphere_NGHN &
     (X,                    & ! <-- Cartesian coordinates of point
      DU(1),                & ! --> Interpolated (1 + N)*Grad(N)
      NP,                   & ! --> Interpolated N
      NHN(1,:,:),           & ! --> Interpolated Grad x (1+N)Grad(N)
      Stat)                   ! --> Error status
   DX(2) = U + DU(1)*(S/2)
   Call Atmosphere_NGHN &
     (X + DX(1)*(S/2),      & ! <-- Cartesian coordinates of point
      DU(2),                & ! --> Interpolated (1 + N)*Grad(N)
      NP,                   & ! --> Interpolated N
      NHN(2,:,:),           & ! --> Interpolated Grad x (1+N)Grad(N)
      Stat)                   ! --> Error status
   DX(3) = U + DU(2)*(S/2)
   Call Atmosphere_NGHN &
     (X + DX(2)*(S/2),      & ! <-- Cartesian coordinates of point
      DU(3),                & ! --> Interpolated (1 + N)*Grad(N)
      NP,                   & ! --> Interpolated N
      NHN(3,:,:),           & ! --> Interpolated Grad x (1+N)Grad(N)
      Stat)                   ! --> Error status
   DX(4) = U + DU(3)*S
   Call Atmosphere_NGHN &
     (X + DX(3)*S,          & ! <-- Cartesian coordinates of point
      DU(4),                & ! --> Interpolated (1 + N)*Grad(N)
      NP,                   & ! --> Interpolated N
      NHN(4,:,:),           & ! --> Interpolated Grad x (1+N)Grad(N)
      Stat)                   ! --> Error status

   X  = X + (DX(1) + 2*DX(2) + 2*DX(3) + DX(4))*(S/6)
   U  = U + (DU(1) + 2*DU(2) + 2*DU(3) + DU(4))*(S/6)
   St = St + S

!   If (St > 4000.0) then
!      Stat = 4
!   End If


   !--- 2.4. Calculation of height

   G    = Geod_from_Cart(X)
   Hper = Min(Hper, G%H)


   !--- 2.5. Renorming U

   U1 = U
   U  = Vector_Normed(U)*(1+NP)


   !--- 2.6. Checking for ray collision with Earth

   If (Stat /= 0) then
      UN = U
      XN = X
      Return
   End If


   !--- 2.7. Calculation of BM matrices

   BM(:,:,:) = 0
   Do Mu=1,4
      BM(Mu,1:3,4:6) = E3(:,:)
      BM(Mu,4:6,1:3) = NHN(Mu,:,:)
   End Do

   B21(:,:)   = MatMul(BM(2,:,:), BM(1,:,:))
   B32(:,:)   = MatMul(BM(3,:,:), BM(2,:,:))
   B43(:,:)   = MatMul(BM(4,:,:), BM(3,:,:))

   B321(:,:)  = MatMul(BM(3,:,:), B21(:,:))
   B432(:,:)  = MatMul(BM(4,:,:), B32(:,:))

   B4321(:,:) = MatMul(BM(4,:,:), B321(:,:))


   !--- 2.8. Calculation of B matrix

   B(:,:) = &
      E6(:,:) + &
      (S/6)*(BM(1,:,:) + 2*BM(2,:,:) + 2*BM(3,:,:) + BM(4,:,:)) + &
      (S**2/6)*(B21(:,:) + B32(:,:) + B43(:,:)) + &
      (S**3/12)*(B321(:,:) + B432(:,:)) + &
      (S**4/24)*B4321(:,:)


   !--- 2.9. Calculation of renorming matrix

   R(1:3,1:3) = E3(:,:)
   R(1:3,4:6) = 0
   R(4:6,1:3) = &
      Tensor_Product(U1%X(:),DU(4)%X(:))/ &
      ((1 + NP)*Sqrt(U1*U1))
   R(4:6,4:6) = &
      (1 + NP)*(E3(:,:) - Tensor_Product(U1%X(:),U1%X(:))/Sqrt(U1*U1))/ &
      Sqrt(U1*U1)


   !--- 2.10. Multiplication of B and C^{mu} with R

   B(:,:) = MatMul(R(:,:),B(:,:))


   !--- 2.11. Calculation of BN

   BN(:,:) = MatMul(B(:,:),BN(:,:))


End Do Ray_Integrate


!----------------------------------------------------------
! 3. LAST STEP FROM ATMOSPHERE TO RECEIVER
!----------------------------------------------------------


!--- 3.1. Last step

UN = Vector_Normed(U)
XN = X

Sr = (XR - X)*UN
XN = XN + Sr*UN


!--- 3.2. Calculation of B_{N} matrix

R(:,:) = 0
R(1:3,1:3) = E3(:,:)
R(4:6,4:6) = &
   (E3(:,:) - Tensor_Product(U%X(:),U%X(:))/Sqrt(U*U))/ &
   Sqrt(U*U)

Sr_X(:) = -UN%X(:)
Sr_U(:) = XR%X(:) - X%X(:)

B(1:3,1:3) = E3(:,:) + Tensor_Product(U%X(:),Sr_X(:))
B(1:3,4:6) = Sr*E3(:,:) + Tensor_Product(U%X(:),Sr_U(:))
B(4:6,4:6) = E3(:,:)
B(4:6,1:3) = 0

B(:,:) = MatMul(B(:,:),R(:,:))


!--- 3.3. Calculation of BN

BN(:,:) = MatMul(B(:,:),BN(:,:))


!--- 3.4. Setting ZN_Z0

ZN_Z0(:,:) = BN(:,:)


End Subroutine Ray_Trace_ZNZ0



!==========================================================
Subroutine Find_Ray_P &
  (XT,     & ! <-- Transmitter position
   VT,     & ! <-- Transmitter velocity
   XR,     & ! <-- Receiver position
   VR,     & ! <-- Receiver velocity
   DS,     & ! <-- Integration step parameter
   XLC,    & ! <-- Local curvature center
   P,      & ! <-- Impact parameter [km]
   UT,     & ! --> Transmitter ray direction
   XN,     & ! --> Ray point nearest to receiver
   UN,     & ! --> Ray direction at XN
   Hper,   & ! --> Ray perigee height
   Stat,   & ! --> Error status
   Strict)   ! <~~ Strict criterion for iterations
!
! Finding ray with prescribed impact parameter.
!----------------------------------------------------------
! Method:
!   Iterative Newton solution.
!----------------------------------------------------------
! (C) Copyright 2000-2003, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 04 May 2000 | Original version.
!   1.1   | 20 Dec 2001 | Accurate combination of
!         |             | Newton method and binary search.
!   2.0   | 07 May 2003 | Modification for Atmosphere.
!   3.0   | 31 Oct 2003 | Strict.
!   3.1   | 14 May 2006 | Temporary storage for comiler compatibility.
!----------------------------------------------------------
! Modules used:
!
Use Coordinates, only: &
! Imported Type Definitions:
    Cartesian,     &
! Imported Routines:
    Rotate,        &
    Vector_Normed, &
    Vector_Norm,   &
    Vector_Angle,  &
! Imported Operators:
    Operator(*),   &
    Operator(+),   &
    Operator(-),   &
    Operator(/),   &
    Operator(.xx.), &
    Assignment(=)
!
Use Occ_Refraction_adj, only: &
! Imported Routines:
    Ray_Refraction_adj,   &
    Ray_Derivatives
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type(Cartesian), Intent(In)  :: &
   XT        ! Transmitter position
!
Type(Cartesian), Intent(In)  :: &
   VT        ! Transmitter velocity
!
Type(Cartesian), Intent(In)  :: &
   XR        ! Receiver position
!
Type(Cartesian), Intent(In)  :: &
   VR        ! Receiver velocity
!
Real(Double), Intent(In)     :: &
   DS        ! Integration step parameter
!
Type(Cartesian), Intent(In)  :: &
   XLC       ! Local curvature center
!
Real(Double), Intent(In)     :: &
   P         ! Impact parameter [km]
!
! Output arguments:
!
Type(Cartesian), Intent(Out) :: &
   UT        ! Transmitter ray direction
!
Type(Cartesian), Intent(Out) :: &
   XN        ! Ray point nearest to receiver
!
Type(Cartesian), Intent(Out) :: &
   UN        ! Ray direction at XN
!
Real(Double), Intent(Out)    :: &
   Hper       ! Ray perigee height
!
Integer, Intent(Out)         :: &
   Stat      ! Error status
!
! Optional arguments:
!
Logical, Optional, Intent(In) :: &
   Strict    ! Strict criterion for iterations
!----------------------------------------------------------
! Local Parameters:
!
Real(Double), Parameter :: &
   Pacc = 1d-6         ! Accuracy of impact parameter
!
Integer, Parameter :: &
   Nit  = 20           ! Maximum number of iterations
!
! Local Scalars:
!
Logical      :: LStrict      ! Stric criterion for iterations
Real(Double) :: Alpha        ! Leveling angle
Real(Double) :: ER           ! Ray refraction angle
Real(Double) :: PR           ! Ray impact parameter
Real(Double) :: DeltaP       ! Impact parameter correction
Integer      :: i            ! Iteration index
Real(Double) :: D            ! Estimate of observation distance
Real(Double) :: AlphaMin     ! Lower estimate of Alpha
Real(Double) :: AlphaMax     ! Upper estimate of Alpha
Integer      :: RStat        ! Ray-tracing status
!
! Local Arrays:
!
! --- Geometric parameters
!
Real(Double) :: Z0(6)        ! Initial ray point (X0,U0)
Real(Double) :: ZN(6)        ! Final ray point (XN,UN)
!
! --- Ray derivatives
!
Real(Double)    :: ZN_Z0(6,6)   ! d(XN,UN)/d(XT,UT)
Real(Double)    :: E_Z0(6)      ! d(Eps)/d(Z0)
Real(Double)    :: E_ZN(6)      ! d(Eps)/d(ZN)
Real(Double)    :: P_Z0(6)      ! d(P)/d(Z0)
Real(Double)    :: P_ZN(6)      ! d(P)/d(ZN)
Real(Double)    :: FZ0_P(6)     ! D(Z0)/D(P)
Real(Double)    :: FE_ZN(6)     ! D(E)/D(ZN)
Type(Cartesian) :: TC           ! Temporary storage
!----------------------------------------------------------


!----------------------------------------------------------
! 1. INITIAL APPROXIMATION
!----------------------------------------------------------

Alpha = ASin(P/Vector_Norm(XT-XLC))

UT = Vector_Normed(Rotate(-(XT-XLC), (XR-XLC).xx.(XT-XLC), Alpha))

D = Sqrt((XT-XLC)*(XT-XLC) - P**2)

AlphaMin = Alpha - 20.0_Double/D
AlphaMax = Alpha + 20.0_Double/D


!----------------------------------------------------------
! 2. ITERATIVE SOLUTION
!----------------------------------------------------------

If (Present(Strict)) then
   LStrict = Strict
Else
   LStrict = .True.
End If

Stat = 1

Newton: Do i=1,Nit

!!!DEBUG
!Write(*,*) '  i = ', i

   Alpha = Vector_Angle(UT, -(XT-XLC))
   
   If ((Alpha < AlphaMin) .or. (Alpha > AlphaMax)) then
      UT = Vector_Normed(Rotate(-(XT-XLC),     &
                         (XR-XLC).xx.(XT-XLC),  &
                         (AlphaMin + AlphaMax)/2))
      Alpha = Vector_Angle(UT, -(XT-XLC))
   End If

   Call Ray_Trace_ZNZ0 &
     (XT,        & ! <-- Transmitter position
      UT,        & ! <-- Transmitter ray direction
      XR,        & ! <-- Receiver position
      DS,        & ! <-- Integration step parameter
      XN,        & ! --> Ray point nearest to receiver
      UN,        & ! --> Ray direction at XN
      ZN_Z0,     & ! --> d(XN,UN)/d(XT,UT)
      Hper,      & ! --> Ray perigee height
      RStat)       ! --> Error status

   If (.not. LStrict) then
      Stat = RStat
   End If

   If (RStat /= 0) then
      AlphaMin = Max(AlphaMin,Alpha)
      UT = Vector_Normed(Rotate(-(XT-XLC),     &
                         (XR-XLC).xx.(XT-XLC),  &
                         (AlphaMin + AlphaMax)/2))
      Cycle Newton
   End If

   Z0(1:3) = XT - XLC
   Z0(4:6) = UT
   ZN(1:3) = XN - XLC
   ZN(4:6) = UN

   Call Ray_Refraction_adj &
     (Z0,        & ! <-- Initial ray point (X0,U0)
      ZN,        & ! <-- Final ray point (XN,UN)
      VT,        & ! <-- Transmitter velocity [km/s]
      VR,        & ! <-- Receiver velocity [km/s]
      ER,        & ! --> Refraction angle [rad]
      PR,        & ! --> Impact parameter [km]
      E_Z0,      & ! --> d(Eps)/d(Z0)
      E_ZN,      & ! --> d(Eps)/d(ZN)
      P_Z0,      & ! --> d(P)/d(Z0)
      P_ZN)        ! --> d(P)/d(ZN)

   DeltaP = P - PR
   
   If (PR < P) then
      AlphaMin = Max(AlphaMin,Alpha)
   Else
      AlphaMax = Min(AlphaMax,Alpha)
   End If

   If (Abs(DeltaP) < Pacc) then
      If (LStrict) then
         Stat = 0
      End If
      Exit Newton
   Else
      If (LStrict) then
         Stat = 1
      End If
   End If

   Call Ray_Derivatives &
     (Z0,        & ! <-- Initial ray point (X0,U0)
      ZN,        & ! <-- Final ray point (XN,UN)
      ZN_Z0,     & ! <-- d(ZN)/d(Z0)
      P_Z0,      & ! <-- d(P)/d(Z0)
      P_ZN,      & ! <-- d(P)/d(ZN)
      E_Z0,      & ! <-- d(E)/d(Z0)
      E_ZN,      & ! <-- d(D)/d(ZN)
      FZ0_P,     & ! --> D(Z0)/D(P)
      FE_ZN)       ! --> D(E)/D(ZN)

   TC = FZ0_P(4:6)
   UT = UT + DeltaP*TC


End Do Newton



End Subroutine Find_Ray_P




End Module Atmosphere_rays_adj


