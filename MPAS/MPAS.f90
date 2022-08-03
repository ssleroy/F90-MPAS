!
Module MPAS
!
! Reading MPAS global gridded fields of H, P, T, and Q,
! and calculation of gridded field of N and interpolation
! coefficients for N and Q
!----------------------------------------------------------
! J. D. Hegarty 2020 based on
! (C) Copyright 2007-2009, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 26 Jun 2007 | Original version (from ECHAM_fields).
!   1.1   | 13 Feb 2008 | T_Profile, T_Matrix moved to GCM_grid.
!   1.2   | 06 Jul 2009 | Get_IDX.
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Parameters:
    Single, Double, rtd

Use IO, only: &
! Imported Routines:
    FileOpen
!
Use Errors, only: &
! Imported Type Definitions:
    Error_Status,  &
! Imported Routines:
    Enter_Callee,  &
    Exit_Callee,   &
    Error
!
!Use MPAS_IO, only: &
!! Imported Types:
!    T_MPAS_Fields

!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Public Type Definitions:
!
Type :: T_MPAS_Fields
   Real(Single), allocatable :: &
      H(:),    & ! Geopotential [gpm]
      T(:),    & ! Theta [K] , need temperature
      QR(:),   & ! Specific Humidity, need Relative humidity [%]
      U(:),    & ! U-velocity [m/s]
      V(:),    & ! V-velocity [m/s]
      OZ(:),   & ! Ozone mixing ration [kg/kg]
      LWC(:),  & ! Liquid water content [kg/kg]
      IWC(:),  & ! Ice water content [kg/kg]
      latcell(:), &
      loncell(:), &
      latvertex(:), &
      lonvertex(:), &
      latedge(:), &
      lonedge(:), &
      rdzw(:), &
      fzp(:), &
      fzm(:)
   Real(Single) :: Hsur ! Surface geopotential [gpm]
   Integer, allocatable :: &
        nedgesoncell(:),  &
        cellsonedge(:,:), &
        edgesoncell(:,:), &
        cellsoncell(:,:), &
        verticesoncell(:,:), &
        cellsonvertex(:,:)
   Integer :: ncid_sav
   !Integer :: ncid_sav
End Type T_MPAS_Fields
!
! Public Parameters:
!
Integer, Parameter :: &
   err_Bad_DataType    = 1201001, &
   err_LonLat_Mismatch = 1201002, &
   err_Lev_Mismatch    = 1201003, &
   err_No_NLev         = 1201004, &
   err_Bad_LevelType   = 1201005, &
   err_Bad_Pgrid       = 1201006, &
   err_Open_Grib       = 1201007
!----------------------------------------------------------
! --- Field states
!
Integer, Parameter :: &
   fst_Null        = 0,  &
   fst_Allocated   = 1,  &
   fst_Initialized = 2

Real(Double), Parameter :: &
   P_up = 1e-5_Double   ! Pressure at uppermost half level [mbar]
!
! Public Scalars:
!
Integer :: NM         ! Number of MSIS pressure levels
Integer :: NLev       ! Number of model levels
Integer :: Data_Type  ! Data representation type
!
! Public Arrays:
!
! --- Grids and gridded GCM fields
!
Type(T_MPAS_Fields) :: &
             Fields   ! Structure with fields
!
Real(Double), Allocatable, Public :: &
   Freq(:)            ! Frequency channels [Hz]
!
Integer, Allocatable :: &
   PStat(:,:,:)       ! Dynamical field status
!
!
! Private Scalars:
!
Logical, Private :: &
   Undefined_Pointers = .TRUE. ! Pointer status indicator.
!----------------------------------------------------------
!
Contains

!==========================================================
Subroutine MPAS_Init &
  (Pathnames, & ! <-- Pathnames of data files
   Vrb,       & ! <-- Verbosity level
   ErrStat,   & ! <-> Error status
   XFreq,     & ! <~~ Frequency channels [Hz]
   Mask)        ! <~~ Mask of fields to read
!
! Reading global fields of surface geopotential,
! surface pressure, temperature and humidity,
! calculation of refractive index and interpolation
! coefficients of refractivie index and humidity.
!----------------------------------------------------------
! Method:
!----------------------------------------------------------
!J. D. Hegarty 2020, based on
! (C) Copyright 2007-2009, M. E. Gorbunov.
!
Use Earth, only: &
! Imported Routines:
    GCLat_from_GDLat
!
Use MSIS, only: &
! Imported Parameters:
    mod_MSIS,             &
! Imported Routines:
    MSIS_Init,            &
    MSIS_Num_Levels,      &
    MSIS_Pressure_Levels
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Character(Len=*), Intent(In)  :: &
   Pathnames(:)   ! Pathnames of data files
!
Integer, Intent(In)  :: &
   Vrb            ! Verbosity level
!
! InOut arguments:
!
Type(Error_Status), Pointer   :: &
   ErrStat        ! Error status
!
! Optional input arguments:
!
Real(Double), Intent(In), Optional :: &
   XFreq(:)       ! Frequency channels [Hz]
!
Logical, Intent(In), Optional :: &
   Mask(:)        ! Mask of fields to read
!----------------------------------------------------------
! Local Scalars:
!
!Integer   :: Date     ! Data set date
character(len=19) :: xtime     ! Data set date
Integer   :: Mon      ! Data set month
Integer   :: Time     ! Data set time
Integer   :: IFile    ! Data file number
Integer   :: NC       ! Number of frequency channels
Integer   :: ErrCode  ! Error code

!----------------------------------------------------------
! 0. INITIALIZATION
!----------------------------------------------------------

!--- 0.1. Error status

Call Enter_Callee &
  ('MPAS_Init',  & ! <-- User routine
   ErrStat)         ! <-> Pointer to callee status

!----------------------------------------------------------
! 1. READING DATA
!----------------------------------------------------------

Read: Do IFile = 1, Size(Pathnames)

! I don't think we will want to read more than 1 file at a time
!  can set ifile=1 and remove the do loop to force it to do only 1 time
   print*,'reading netcdf file',Pathnames(Ifile)
!!   Call Get_GRIB &
   Call Get_NETCDF_MPAS &
     (PathNames(IFile),  & ! <-- Pathname of data file
      Vrb,               & ! <-- Verbosity level
      Data_Type,         & ! --> Data representation type
      Fields,            & ! <-> Structure with fields
      xtime,             & ! --> Field date
      Time,              & ! --> Field time
      ErrStat)             ! --> Error code

   If (Error(ErrStat)) then
      Call Exit_Callee(ErrStat)
      Return
   End If

End Do Read

!----------------------------------------------------------
! 7. MSIS INITIALIZATION
!----------------------------------------------------------

read(xtime(6:7),*) mon

print*,'Mon=',Mon

Call MSIS_Init &
  (Mon,          & ! <-- Month
   'MSIS??.asc', & ! <-- Generic name of MSIS coefficient files
   mod_MSIS,     & ! <-- Initialization mode
   Vrb,          & ! <-- Verbosity level
   ErrStat)        ! --> Error code

If (Error(ErrStat)) then
   Call Exit_Callee(ErrStat)
   Return
End If

!print*,'In MPAS_INIT, NLEV=', NLev
!
Call Exit_Callee(ErrStat)

End Subroutine MPAS_Init

!==========================================================
Subroutine Interpolate_Refractivity_MPAS &
  (PLon,     & ! <-- Longiude of point [deg]
   PLat,     & ! <-- Latitude of point [deg]
   ZP,       & ! <-- Altitude of point [km]
   Zmin,     & ! --> Minimum model Z for this lon/lat
   Zmax,     & ! --> Maximum model Z for this lon/lat
   NP,       & ! --> Interpolated N
   NG,       & ! ~~> Interpolated dN/d(alt,lat,lon)
   NH,       & ! ~~> Interpolated hessian matrix of N
   NIP)        ! ~~> Interpolated Im(N)
!
! Interpolation of refractivity calculated from MPAS meteorological profiles.
!----------------------------------------------------------
! Method:
!   1. Vertical spline interpolation to given altitude.
!   2. Summation of 2D Lon/Lat interpolation series:
!   f(z,Lat,Lon) = sum WG f(z,Lat ,Lon )
!                   i    i       i    i
!   Index i enumerates points of the interpolation
!   subgrid. For rectangular grid, weigting function WG
!   is calculated for polynomial interpolation with the
!   polynomial power defined by the interpolation subgrid
!   dimension NI. For icosahedral grid WG is defined as
!   symplectic weights of the vertices of surrounding
!   triangle.
!----------------------------------------------------------
! J. D. Hegarty, 2021 based on code from
! (C) Copyright 1998-2004, M. E. Gorbunov.
!
! Modules used:
!
Use NetCDF, only: &
! Imported Parameters:
    NF90_NoWrite,             &
    NF90_NoErr,               &
    NF90_BYTE,                &
    NF90_CHAR,                &
    NF90_SHORT,               &
    NF90_INT,                 &
    NF90_FLOAT,               &
    NF90_DOUBLE,              &
    NF90_GLOBAL,              &
! Do I need to specifically import all the routines that I use?
! Imported Routines:
    NF90_Open,                &
    NF90_Inquire,             &
    NF90_Inquire_Dimension,   &
    NF90_INQ_DIMID,           &
    NF90_Inquire_Variable,    &
    NF90_INQ_VARID,           &
    NF90_Inq_Attname,         &
    NF90_Inquire_Attribute,   &
    NF90_Get_Att,             &
    NF90_Get_Var,             &
    NF90_Close

Use Occ_Meteoprofiles, only: &
! Imported Routines:
    Q_from_Qrel           

!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
!----------------------------------------------------------
! Public Parameters:
!
Integer, Parameter :: &
   err_NoFile   = 0702001,  & ! File not found
   err_MissVar  = 0702002,  & ! Missing variable/attributes
   err_ReadErr  = 0702003,  & ! Data read error
   err_FmtErr   = 0702004,  & ! Data format error
   err_NoLCF    = 0702005,  & ! No LCF
   err_NoFreq   = 0702006,  & ! No frequencies
   err_WrtErr   = 0702007,  & ! Write error
   err_DataType = 0702008,  & ! Bad data type
   err_NoData   = 0702009,  & ! No data
   err_Dimen    = 0702010,  & ! Inconsistent dimensions
   err_Channel  = 0702011,  & ! Inconsistent numbers of channels
   err_Coord    = 0702012,  & ! Bad coordinates
   err_FileName = 0702013     ! Bad file name
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
   Zmin          ! Minimum model Z for this lon/lat
!
Real(Double), Intent(Out) :: &
   Zmax          ! Maximum model Z for this lon/lat
!
Real(Double), Intent(Out) :: &
   NP            ! Interpolated N
!
Real(Double), Optional, Intent(Out) :: &
   NG(3)         ! Interpolated gradient dN/d(alt,lat,lon)
!
Real(Double), Optional, Intent(Out) :: &
   NH(3,3)       ! Interpolated hessian
                 ! (d/d(alt,lat,lon))x(d/d(alt,lat,lon))N
!
Real(Double), Optional, Intent(Out) :: &
   NIP(:)        ! Interpolated Im(N)

real, parameter :: pi_const = 2.0 * asin(1.0)

!----------------------------------------------------------
! Local Variables:
!
integer :: ncid, iostat,ndims,unlimid,nd,nrel
integer :: varid, nattrs,xtype
integer :: idxm, last_idxm,j
integer :: maxedges, ncells, nvertices
Character(Len=255) :: Name         ! Parameter name
real :: platr,plonr,clatd,clond
Real(Double) :: pmb
Real(Double) :: q
Real(Double) :: TR
Real(Double) :: QRR
integer :: nearestcell
integer :: sourcecells(3)
integer :: sourcevertices(3)
Integer, dimension(13)  :: dimids, lend, lendh, lends
!integer :: sourceedges(maxedges,npts)
real :: cellweights(3)
real :: vertexweights(3)
real :: pointinterp(3)
!      real :: edgeweights (maxedges,npts)
real :: vertcoords(3,3)
real :: vec(3)
real :: mpaswach(3)
real :: t2m(3), q2(3)
! ---- JDH debug -----
real :: thetatest
integer :: thetadimid
!real, allocatable :: theta (:,:,:)
!--------------------------------
Integer        :: NGP     ! Subgrid dimension
Integer        :: IGP     ! Subgrid index
Integer        :: KLon    ! Longitude subgrid index
Integer        :: KLat    ! Latitude subgrid index
Integer        :: m1, m2  ! Interpolation triangle indices
Real(Double)   :: FP      ! Interpolated Ln(N)
Integer        :: i       ! Work index
Integer        :: NC      ! Number of frequency channels
Integer        :: IC      ! Channel number
Integer        :: Stat    ! Error status
!
!  These should all be 3, for the 3 MPAS cell points
!Real(Double), Allocatable :: &
!   ZNmin(:),          & ! Minimum model Z on subgrid
Real(Double) :: &
   ZNmin(3),          & ! Minimum model Z on subgrid
   ZNmax(3),          & ! Maximum model Z on subgrid
   FV(3),             & ! Vertically interpolated Ln(N) on subgrid
   FZ(3),             & ! Vertical gradient on subgrid
   F2Z(3)               ! Second vertial derivative on subgrid
Real(Double) :: &
   FG(3),             & ! Interpolated dln(N)/d(alt,lat,lon)
   FH(3,3)              ! Interpolated hessian
                        ! (d/d(alt,lat,lon))x(d/d(alt,lat,lon))ln(N)
Real(Double), Allocatable :: &
   FVI(:,:),          & ! Vertically interpolated Ln(Im(N)) on subgrid
   FIP(:)               ! Interpolated Ln(Im(N))

real, allocatable :: ti(:,:) ! working array for 3 points (level,3) temperature (k)
real, allocatable :: qi(:,:)  ! working arrays for 3 points (level,3) specific humidity (g/g)
real, allocatable :: hi(:,:) ! working array for 3 points (level,3) geopotential height (m) 
real, allocatable :: pi(:,:) ! working array for 3 points (level,3) pressure (Pa)
! Not sure I need the below arrays as the surface geopotential will be the first
! level of the geopotential profile 
real ::  hs(3)         ! working array for 3 points surface geopotential height (m)

! interpolated profiles from MPAS
real, allocatable :: t(:) ! temperature (k)
real, allocatable :: qr(:)  ! specific humidity (g/g)
real, allocatable :: h(:) ! geopotential height (m) 
real, allocatable :: p(:) ! pressure (mb) - converted from Pa before passed into subroutines
 
real, allocatable :: t_int(:) ! temperature (k)
real, allocatable :: q_int(:)  ! specific humidity (g/g)
real, allocatable :: p_int(:) ! pressure (Pa) 
real ::  hsur,tv,hscl     ! surface geopotential, virtual potential temp, scale height 

!----------------------------------------------------------
! 1. INITIALIZATION
!----------------------------------------------------------

!print*,'in interpolate_refractivity fields%ncid_sav=',fields%ncid_sav
ncid = fields%ncid_sav
!print*,'ncid=',ncid
!print*,'fields%loncell(10)=',fields%loncell(10)

!latcell = fields%latcell

! The above line would probably require latcell to be allocated.
! Can I just enter the fields%latcell into the subroutines? 
!latcell => fields%latcell

!print*,size(fields%latcell)
!print*,'in interpolate_refractivity_MPAS, latcell=',fields%latcell(10),fields%latcell(100)
!print*,'fields%loncell(10)=',fields%loncell(10)
!!loncell => fields%loncell
!print*,size(fields%loncell)
!print*,'in interpolate_refractivity_MPAS, loncell=',fields%loncell(10),fields%loncell(100)
!print*,size(fields%latvertex)
!print*,size(fields%lonvertex)
!print*,size(fields%latedge)
!-------------------------------------------------------------------------------------
!  2.0 Find 3 MPAS Horizontal Cells and Calculate Cell Weights
!-------------------------------------------------------------------------------------

!print*,'plon,plat,zp=',plon,plat,zp

! MPAS latitude go from -pi/2 to +pi/2
! MPAS longitudes go from 0 to 2 pi

plonr=plon*pi_const/180.0
platr=plat*pi_const/180.0

ncells=size(fields%verticesoncell,2)
maxedges=size(fields%verticesoncell,1)
nvertices=size(fields%latvertex)

!print*,'ncells,maxedges,nvertices=',ncells,maxedges,nvertices

last_idxm=1

idxm = nearest_cell(platr,plonr, last_idxm, nCells, & 
     maxEdges,fields%nEdgesOnCell,fields%cellsOnCell,fields%latCell,fields%lonCell)

!print*,'nearest cell, idxm=',idxm

last_idxm = 1
idxm = nearest_vertex(platr, plonr, last_idxm, nCells, & 
     nVertices,maxEdges, 3, fields%nEdgesOnCell, fields%verticesOnCell, & 
     fields%cellsOnVertex, fields%latCell, fields%lonCell, & 
     fields%latVertex, fields%lonVertex )

if (idxm <= 0) then
   sourceCells(:) = 1
   cellWeights(:) = 0.0
else

!! convert_lx returns a 3 element vector, pointInterp is an array of 3
!! SourceCells an array of 3 x ncells.  Not sure if they can be assigned
!! in this way.  If not loop it.
   call convert_lx(platr,plonr,6371229.0,vec)
   pointinterp(1:3)=vec(1:3)
!                if(idebug==1) then
!                   print*,'pointinterp=',pointinterp
!                   stop
!                end if

   sourceCells(:) = fields%cellsOnVertex(:,idxm)
   do j=1,3
      call convert_lx (fields%latCell(fields%cellsOnVertex(j,idxm)), &
           fields%lonCell(fields%cellsOnVertex(j,idxm)),6371229.0, vec)
      vertCoords(1:3,j) = vec(1:3)
   end do

!! This call creates cellweights for all three vertices at once in a single call
!! to the function.  I may need to re-write if this syntax doesn't work in a single line
   call mpas_wachspress_coordinates(3,vertCoords,pointInterp,mpaswach)
   cellWeights(:)=mpaswach(:) 
   last_idxm = idxm
end if

!print*,'sourcecells=', sourcecells
!print*,'cellweights=',cellweights

!----------------------------------------------------------------------------------------------
! 3.0 Read MPAS Meteorological Profiles for 3 Horizontal Cells
!---------------------------------------------------------------------------------------------- 
!
! Read theta, specific humidity, pressure and geopotential profiles from MPAS
!
!------------------ theta ----------------------
name='theta'
IOStat=NF90_INQ_VARID(NCID,name,varid)
If (IOStat /= NF90_NoErr) then
   print*,'nf90_inq_varid IOStat=',IOStat
   print*,'Interpolate_Refractivity_MPAS name=',name
   Stat = err_ReadErr
   Return
End If
!print*,'varid=',varid
! The first call to NF90_Inquire _Variable gets the
! number of dimensions, the second call gets the dimension
! IDs from there I think I can get the actual diminesion perhaps using
! Inquire dimension.  Try to follow scan_inpur.F in convert_mpas.F code.
   IOStat = NF90_Inquire_Variable(NCID,varid,name,xtype,ndims,dimids,nattrs )
   If (IOStat /= NF90_NoErr) then
      Stat = err_ReadErr
      Return
   End If

!print*,'after NF90_Inquire_Variable, ndims,nattrss,name=',ndims,nattrs,name,dimids

do nd=1,ndims
   IOStat=NF90_INQUIRE_DIMENSION(NCID,dimids(nd),name,lend(nd))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_inq_dimension IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
!   print*,'name,lend=',name,lend(nd)
end do

!allocate(theta(lend(1),lend(2),lend(3)))
!
!IOStat=NF90_GET_VAR(NCID,varid,theta)
!If (IOStat /= NF90_NoErr) then
!   print*,'nf90_get_var IOStat=',IOStat
!   Stat = err_ReadErr
!   Return
!End If
!print*,'theta(50,5000,1)=',theta(50,5000,1)

allocate(ti(lend(1),3))

do j=1,3
   IOStat=NF90_GET_VAR(NCID,varid,ti(:,j),start=(/1,sourcecells(j),1/), &
        count=(/lend(1),1,1/))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_get_var IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
end do

!----------------- pressure ---------------------------
name='pressure'
IOStat=NF90_INQ_VARID(NCID,name,varid)
If (IOStat /= NF90_NoErr) then
   print*,'nf90_inq_varid IOStat=',IOStat
   print*,'Interpolate_Refractivity_MPAS name=',name
   Stat = err_ReadErr
   Return
End If
!print*,'varid=',varid
   IOStat = NF90_Inquire_Variable(NCID,varid,name,xtype,ndims,dimids,nattrs )
   If (IOStat /= NF90_NoErr) then
      Stat = err_ReadErr
      Return
   End If

!print*,'after NF90_Inquire_Variable, ndims,nattrss,name=',ndims,nattrs,name,dimids

do nd=1,ndims
   IOStat=NF90_INQUIRE_DIMENSION(NCID,dimids(nd),name,lend(nd))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_inq_dimension IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
!   print*,'name,lend=',name,lend(nd)
end do

allocate(pi(lend(1),3))

do j=1,3
   IOStat=NF90_GET_VAR(NCID,varid,pi(:,j),start=(/1,sourcecells(j),1/), &
        count=(/lend(1),1,1/))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_get_var IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
end do

!--------------------- specific humidity ----------------------
!name='qv'
! read in relative humidity instead.
name='relhum'
IOStat=NF90_INQ_VARID(NCID,name,varid)
If (IOStat /= NF90_NoErr) then
   print*,'nf90_inq_varid IOStat=',IOStat
   print*,'Interpolate_Refractivity_MPAS name=',name
   Stat = err_ReadErr
   Return
End If
!print*,'varid=',varid
   IOStat = NF90_Inquire_Variable(NCID,varid,name,xtype,ndims,dimids,nattrs )
   If (IOStat /= NF90_NoErr) then
      Stat = err_ReadErr
      Return
   End If

!print*,'after NF90_Inquire_Variable, ndims,nattrss,name=',ndims,nattrs,name,dimids

do nd=1,ndims
   IOStat=NF90_INQUIRE_DIMENSION(NCID,dimids(nd),name,lend(nd))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_inq_dimension IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
!   print*,'name,lend=',name,lend(nd)
end do

allocate(qi(lend(1),3))

do j=1,3
   IOStat=NF90_GET_VAR(NCID,varid,qi(:,j),start=(/1,sourcecells(j),1/), &
        count=(/lend(1),1,1/))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_get_var IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
end do

!-------------------- geopotential height ---------------------------------
name='zgrid'
IOStat=NF90_INQ_VARID(NCID,name,varid)
If (IOStat /= NF90_NoErr) then
   print*,'nf90_inq_varid IOStat=',IOStat
   print*,'Interpolate_Refractivity_MPAS name=',name
   Stat = err_ReadErr
   Return
End If
!print*,'varid=',varid
   IOStat = NF90_Inquire_Variable(NCID,varid,name,xtype,ndims,dimids,nattrs )
   If (IOStat /= NF90_NoErr) then
      Stat = err_ReadErr
      Return
   End If

!print*,'after NF90_Inquire_Variable, ndims,nattrss,name=',ndims,nattrs,name,dimids

do nd=1,ndims
   IOStat=NF90_INQUIRE_DIMENSION(NCID,dimids(nd),name,lendh(nd))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_inq_dimension IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
!   print*,'name,lendh=',name,lendh(nd)
end do

! only get up to the number of half levels
 allocate(hi(lend(1),3))
!allocate(hi(lendh(1),3))

do j=1,3
   IOStat=NF90_GET_VAR(NCID,varid,hi(:,j),start=(/1,sourcecells(j),1/), &
        count=(/lend(1),1,1/))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_get_var IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
end do

!-------------------- t2m ---------------------------------
name='t2m'
IOStat=NF90_INQ_VARID(NCID,name,varid)
If (IOStat /= NF90_NoErr) then
   print*,'nf90_inq_varid IOStat=',IOStat
   print*,'Interpolate_Refractivity_MPAS name=',name
   Stat = err_ReadErr
   Return
End If
!print*,'varid=',varid
   IOStat = NF90_Inquire_Variable(NCID,varid,name,xtype,ndims,dimids,nattrs )
   If (IOStat /= NF90_NoErr) then
      Stat = err_ReadErr
      Return
   End If

!print*,'after NF90_Inquire_Variable, ndims,nattrss,name=',ndims,nattrs,name,dimids

do nd=1,ndims
   IOStat=NF90_INQUIRE_DIMENSION(NCID,dimids(nd),name,lends(nd))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_inq_dimension IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
!   print*,'name,lends=',name,lends(nd)
end do

do j=1,3
   IOStat=NF90_GET_VAR(NCID,varid,t2m(j),start=(/sourcecells(j),1/))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_get_var IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
end do

!-------------------- q2 ---------------------------------
name='q2'
IOStat=NF90_INQ_VARID(NCID,name,varid)
If (IOStat /= NF90_NoErr) then
   print*,'nf90_inq_varid IOStat=',IOStat
   print*,'Interpolate_Refractivity_MPAS name=',name
   Stat = err_ReadErr
   Return
End If
!print*,'varid=',varid
   IOStat = NF90_Inquire_Variable(NCID,varid,name,xtype,ndims,dimids,nattrs )
   If (IOStat /= NF90_NoErr) then
      Stat = err_ReadErr
      Return
   End If

!print*,'after NF90_Inquire_Variable, ndims,nattrss,name=',ndims,nattrs,name,dimids

do nd=1,ndims
   IOStat=NF90_INQUIRE_DIMENSION(NCID,dimids(nd),name,lends(nd))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_inq_dimension IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
!   print*,'name,lends=',name,lends(nd)
end do

do j=1,3
   IOStat=NF90_GET_VAR(NCID,varid,q2(j),start=(/sourcecells(j),1/))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_get_var IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
end do

!-------------  End Read 3 cells of MPAS meteorological profiles  - ------------
! 35974 is sourcecells(1) and there are 3 interpolation points

! Debug print

!do j=1,3
!   print*,'j=',j
!   do i=1, lend(1)
!!      print*,'i,theta,hi,pi,ti,qi=',i,hi(i,j),pi(i,j),theta(i,sourcecells(j),1),ti(i,j),qi(i,j)
!      print*,'i,hi,pi,ti,qi=',i,hi(i,j),pi(i,j),ti(i,j),qi(i,j)
!   end do
!end do

! Following Gorbunov Code, lets do vertical regridding and calculate refractivity for each of 
! the cell profiles, then can apply MPAS cell weights later. 


! At this point let's convert layer relative humidity to 
! specific humidity using pi, ti and qi

do j=1,3
   do i=1,lend(1)

      !   If (Qi(i) >= 0.0 .and. Q(i) <= 100.0) then
      If (Qi(i,j) >= 0.0 .and. Qi(i,j) <= 110.0) then
         tv = ti(i,j) * (pi(i,j)/1e05)**0.286
         pmb=pi(i,j)/100.0
         
         TR   = Real(tv, Double)
         QRR  = Real(Qi(i,j),Double)
         
         Q = Q_from_Qrel &
              (TR,    & ! <-- Temperature [K]
              PMB,  & ! <-- Pressure [mb]
              QRR)     ! <-- Relative humidity [%]
      Else
         Q = 0.0
      End If

! Replace relative humidity with specific humidity
!      if(i.eq.20) then
!         print*,'PMB,TR,QRR,Q=',PMB,TR,QRR,Q
!      end if
      qi(i,j)=q
      
   End Do
end do

allocate(p(lend(1)))
allocate(t(lend(1)))
allocate(qr(lend(1)))
allocate(h(lend(1)))

t(1)=0.0
qr(1)=0.0

If (Allocated(Freq) .and. Present(NIP)) then
   NC = Size(Freq)
   Allocate(FVI(3,NC))
   Allocate(FIP(NC))
End If

do j=1,3
   t(1)=t2m(j)
   qr(1)=q2(j)
   hsur=hi(1,j)
   h(1)=hi(1,j)

! 4.2 Vertical Regridding to geopotential Height Grid (MPAS zgrid)

   do i=2,lend(1)
      p(i)=fields%fzm(i)*pi(i,j) + fields%fzp(i)*pi(i-1,j)
      t(i)=fields%fzm(i)*ti(i,j) + fields%fzp(i)*ti(i-1,j)
      qr(i)=fields%fzm(i)*qi(i,j) + fields%fzp(i)*qi(i-1,j)
      h(i)=hi(i,j)
   end do

! Calculate surface pressure using the hypsometric equation
!     First get the virtual temperature
   tv = ti(1,j) * (pi(1,j)/1e05)**0.286
! Approximate virtual temperature factors are in Celcius 
   tv = (tv-273.15)*(1 + 0.61*qi(1,j)*1000.0) + 273.15
!     Scale height
   hscl=29.3 * tv
   p(1)=p(2)*exp((h(2) - h(1))/hscl)

! Convert theta to temperature 
   do i=1,lend(1)
      t(i)=t(i) * (p(i)/1e05)**0.286
   end do

! Convert pressure to mb from Pa 
   do i=1,lend(1)
      p(i)=p(i)/100.0
   end do

! Debug print
!   print*,'Vertical Regridded MPAS meteorology profile at j=1',j
!!   print*,'lat,lon=',fields%latCell(sourcecells(j)),fields%lonCell(sourcecells(j))
   clatd=fields%latCell(sourcecells(j)) * 180.0 / pi_const
   clond=fields%lonCell(sourcecells(j)) * 180.0 / pi_const
!   print*,'lat,lon=',clatd,clond,plon,plat
!   do i=1,lend(1)
!      print*,h(i),p(i),t(i),qr(i)
!   end do

   If (Allocated(Freq) .and. Present(NIP)) then
      Call Vertical_Interpolation_LnN_MPAS &
        (clatd,         & ! latitude of MPAS cell (deg)
         clond,         & ! longitude of MPAS cell (deg)
         hsur,         & ! surface geopotential
         h,            & ! geopotential
         t,            & ! temperature
         qr,           & ! specific humidity
         p,            & ! pressure [Pa]
         ZP,                 & ! <-- Altitude [km]
         ZNmin(j),         & ! --> Minimum model Z for this lon/lat
         ZNmax(j),         & ! --> Maximum model Z for this lon/lat
         FV(j),            & ! --> Interpolated Ln(N)
         FZ(j),            & ! --> Interpolated dLn(N)/dZ
         F2Z(j),           & ! --> Interpolated d2Ln(N)/dZ2
         FVI(j,:))           ! ~~> Interpolated Ln(NI)
   Else
      Call Vertical_Interpolation_LnN_MPAS &
        (clatd,         & ! latitude of MPAS cell (deg)
         clond,         & ! longitude of MPAS cell (deg)
         hsur,         & ! surface geopotential
         h,            & ! geopotential
         t,            & ! temperature
         qr,           & ! specific humidity
         p,            & ! pressure [Pa]
         ZP,                 & ! <-- Altitude [km]
         ZNmin(j),         & ! --> Minimum model Z for this lon/lat
         ZNmax(j),         & ! --> Maximum model Z for this lon/lat
         FV(j),            & ! --> Interpolated Ln(N)
         FZ(j),            & ! --> Interpolated dLn(N)/dZ
         F2Z(j)) ! --> Interpolated d2Ln(N)/dZ2
   End If

end do   ! j- loop through 3 points

! Note, in the above j loop through 3 MPAS points the refractivity interpolated to ZP 
! for each of the 3 MPAS points and is calculated and saved in FV array.  The horizonatl weights sumweights
! are calculated before this loop. 


!!-------------------------------------------------------------------------------------------------------
!! 4. Perform Horizontal Interpolation from MPAS Cells to RO Point
!!-------------------------------------------------------------------------------------------------------
!--- 4.1. Function and height limits

!
FP   = Sum(cellweights(:)*FV(:))
Zmin = Sum(cellweights(:)*ZNmin(:))
Zmax = Sum(cellweights(:)*ZNmax(:))
!
!
!!--- 4.2. Gradient
!
If (Present(NG)) then
!   FG(1) = Sum(WG(:,0,0)*FZ(:))
   FG(1) = Sum(cellweights(:)*FZ(:))
!   print*,'cell weighted vertical gradient calculated in ln fg1=',fg(1)
! longitude and latitude gradients
!   FG(2) = Sum(WG(:,0,1)*FV(:))
!   FG(3) = Sum(WG(:,1,0)*FV(:))
! In NCEP version of code NG1 is set to 2 which in the latitude and longitude
! interpolation routines makes these particular weights equal to zero.
   fg(2)=0.0
   fg(3)=0.0
End If
!
!
!!--- 4.3. Hessian
!
! The Hessian portion of the code appears to be only utilized
! in the adjoint version.  It is not clear in the Wave scripts how to run
! the adjoint.
If (Present(NH)) then
   print*,'Hessian interpolation specified'

!   FH(1,1) = Sum(WG(:,0,0)*F2Z(:))
!   FH(2,2) = Sum(WG(:,0,2)*FV(:))
!   FH(3,3) = Sum(WG(:,2,0)*FV(:))
!   FH(1,2) = Sum(WG(:,0,1)*FZ(:))
!   FH(2,1) = FH(1,2)
!   FH(1,3) = Sum(WG(:,1,0)*FZ(:))
!   FH(3,1) = FH(1,3)
!   FH(3,2) = Sum(WG(:,1,1)*FV(:))
!   FH(2,3) = FH(3,2)
End If
!
!
!--- 4.4. Imaginary part
!
If (Allocated(Freq) .and. Present(NIP)) then
   Do IC=1,NC
!      FIP(IC) = Sum(WG(:,0,0)*FVI(:,IC))
      FIP(IC) = Sum(cellweights(:)*FVI(:,IC))
   End Do
End If
!
!
!----------------------------------------------------------
! 5. TRANSFORM FROM LN(N) TO N
!----------------------------------------------------------
!
NP = Exp(FP)

!print*,'np,zmin,zmax=',np,zmin,zmax
!
! Note dlnN/dalt = dN/N/dalt, so N*dN/N/dalt = dN/dalt
If (Present(NG)) then
   NG(:) = NP*FG(:)
!   print*,'vertical gradient N*dlnN/dalt is',NG(1)
End If
!
If (Present(NH)) then
   print*,'need to figure out Hessian Interpolation'
!   Do i=1,3
!      NH(i,:) = NP*(FH(i,:) + FG(i)*FG(:))
!   End Do
End If
!
If (Present(NIP)) then
   If (Allocated(Freq)) then
      NIP(:) = Exp(FIP(:))
   Else
      NIP(:) = 0.0_Double
   End If
End If

!
!
!!----------------------------------------------------------
!! 6. MEMORY DEALLOCATION
!!----------------------------------------------------------
!
! These are not allocatable anymore as there will always be 3 points
!Deallocate(IDX,   Stat=Stat)
!Deallocate(WG,    Stat=Stat)
!Deallocate(ZNmin, Stat=Stat)
!Deallocate(ZNmax, Stat=Stat)
!Deallocate(FV,    Stat=Stat)
!Deallocate(FZ,    Stat=Stat)
!Deallocate(F2Z,   Stat=Stat)
!

deallocate (ti, stat=stat) ! working array for 3 points (level,3) temperature (k)
deallocate (qi, stat=stat)  ! working arrays for 3 points (level,3) specific humidity (g/g)
deallocate (hi, stat=stat) ! working array for 3 points (level,3) geopotential height (m) 
deallocate (pi, stat=stat) ! working array for 3 points (level,3) pressure (Pa)

! interpolated profiles from MPAS
deallocate (t, stat=stat) ! temperature (k)
deallocate (qr, stat=stat)  ! specific humidity (g/g)
deallocate (h, stat=stat) ! geopotential height (m) 
deallocate (p, stat=stat) ! pressure (mb) - converted from Pa before passed into subroutines
 
deallocate (t_int, stat=stat) ! temperature (k)
deallocate (q_int, stat=stat) ! specific humidity (g/g)
deallocate (p_int, stat=stat) ! pressure (Pa) 

If (Allocated(Freq) .and. Present(NIP)) then
   Deallocate(FVI, Stat=Stat)
   Deallocate(FIP, Stat=Stat)
End If

End Subroutine Interpolate_Refractivity_MPAS

!---------------------------------------------------------------------

integer function nearest_cell(target_lat, target_lon, start_cell, & 
     nCells, maxEdges,nEdgesOnCell, cellsOnCell, latCell, lonCell)

  implicit none

  real, intent(in) :: target_lat, target_lon
  integer, intent(in) :: start_cell
  integer, intent(in) :: nCells, maxEdges
  integer, dimension(nCells), intent(in) :: nEdgesOnCell
  integer, dimension(maxEdges,nCells), intent(in) :: cellsOnCell
  real, dimension(nCells), intent(in) :: latCell, lonCell
      
  integer :: i
  integer :: iCell
  integer :: current_cell
  real :: current_distance, d
  real :: nearest_distance
!  real :: sphere_distance

!      print*,' IN NEAREST CELL FUNCTION'

!      print*,'target_lat=',target_lat
!      print*,'target_lon=',target_lon
!      print*,'start_cell=',start_cell
!      print*,'ncells,maxedges,nedgesoncell=',
!     &     ncells,maxedges,nedgesoncell(10)

!      print*,'max / min cell lat=',maxval(latcell),minval(latcell)
!      print*,'max / min cell lon=',maxval(loncell),minval(loncell)

  nearest_cell = start_cell
  current_cell = -1
  
  do while (nearest_cell /= current_cell)
     current_cell = nearest_cell
     current_distance = sphere_distance(latCell(current_cell), & 
          lonCell(current_cell), target_lat,target_lon, 1.0)
     nearest_cell = current_cell
     nearest_distance = current_distance
     do i = 1, nEdgesOnCell(current_cell)
        iCell = cellsOnCell(i,current_cell)
        if (iCell > 0 .and. iCell <= nCells) then
           d = sphere_distance(latCell(iCell), lonCell(iCell), & 
                target_lat, target_lon, 1.0)
           if (d < nearest_distance) then
              nearest_cell = iCell
              nearest_distance = d
           end if
        else
           nearest_cell = 0
           return
        end if
     end do
  end do
      
end function nearest_cell

!------------------------------------------------------------------------------

real function sphere_distance(lat1, lon1, lat2, lon2, radius)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Compute the great-circle distance between (lat1, lon1) and (lat2, lon2) on a
  !   sphere with given radius.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  implicit none
  
  real, intent(in) :: lat1, lon1, lat2, lon2, radius

  real :: arg1

  arg1 = sqrt( sin(0.5*(lat2-lat1))**2 + &  
       cos(lat1)*cos(lat2)*sin(0.5*(lon2-lon1))**2 )
  sphere_distance = 2.*radius*asin(arg1)
      
end function sphere_distance

!----------------------------------------------------------

integer function nearest_vertex( target_lat, target_lon, & 
     start_vertex, nCells, nVertices, maxEdges, & 
     vertexDegree, nEdgesOnCell, verticesOnCell, & 
     cellsOnVertex, latCell, lonCell, &
     latVertex, lonVertex )

  implicit none

  real, intent(in) :: target_lat, target_lon
  integer, intent(in) :: start_vertex
  integer, intent(in) :: nCells, nVertices, maxEdges, vertexDegree
  integer, dimension(nCells), intent(in) :: nEdgesOnCell
  integer, dimension(maxEdges,nCells), intent(in) :: verticesOnCell
  integer, dimension(vertexDegree,nVertices), intent(in) :: cellsOnVertex
  real, dimension(nCells), intent(in) :: latCell, lonCell
  real, dimension(nVertices), intent(in) :: latVertex, lonVertex

  integer :: i, cell1, cell2, cell3, iCell
  integer :: iVtx
  integer :: current_vertex
  real :: cell1_dist, cell2_dist, cell3_dist
  real :: current_distance, d
  real :: nearest_distance
!      real :: sphere_distance
!      real :: mpas_triangle_signed_area_sphere
!      real :: mpas_arc_length
!!      real, dimension(nVertices) :: mpas_wachspress_coordinates

  nearest_vertex = start_vertex
  current_vertex = -1

!      print*,'nearest_vertex=',nearest_vertex
!      print*,'target_lat,target_lon=',target_lat,target_lon
!      print*,'cellsonvertex=',cellsonvertex(:,10000)

  do while (nearest_vertex /= current_vertex)
     current_vertex = nearest_vertex
     current_distance = sphere_distance(latVertex(current_vertex), & 
          lonVertex(current_vertex),target_lat,target_lon,1.0)
     nearest_vertex = current_vertex
     nearest_distance = current_distance
     cell1 = cellsOnVertex(1,current_vertex)
     if (cell1 <= 0) then
        nearest_vertex = 0
        return
     end if
     cell1_dist = sphere_distance(latCell(cell1), lonCell(cell1), & 
          target_lat, target_lon, 1.0)
     cell2 = cellsOnVertex(2,current_vertex)
     if (cell2 <= 0) then
        nearest_vertex = 0
        return
     end if
     cell2_dist = sphere_distance(latCell(cell2), lonCell(cell2), & 
          target_lat, target_lon, 1.0)
     if (vertexDegree == 3) then
        cell3 = cellsOnVertex(3,current_vertex)
        if (cell3 <= 0) then
           nearest_vertex = 0
           return
        end if
        cell3_dist = sphere_distance(latCell(cell3), lonCell(cell3), & 
             target_lat, target_lon, 1.0)
     end if
     if (vertexDegree == 3) then
        if (cell1_dist < cell2_dist) then
           if (cell1_dist < cell3_dist) then
              iCell = cell1
           else
              iCell = cell3
           end if
        else
           if (cell2_dist < cell3_dist) then
              iCell = cell2
           else
              iCell = cell3
           end if
        end if
     else
        if (cell1_dist < cell2_dist) then
           iCell = cell1
        else
           iCell = cell2
        end if
     end if
     do i = 1, nEdgesOnCell(iCell)
        iVtx = verticesOnCell(i,iCell)
        d = sphere_distance(latVertex(iVtx), lonVertex(iVtx), & 
             target_lat, target_lon, 1.0)
        if (d < nearest_distance) then
           nearest_vertex = iVtx
           nearest_distance = d
        end if
     end do
  end do
  
end function nearest_vertex

!----------------------------------------------------------

!***********************************************************************
!
!  function mpas_wachspress_coordinates
!
!> \brief Compute the barycentric Wachspress coordinates for a polygon
!> \author  Phillip Wolfram
!> \date    01/26/2015
!> \details
!>  Computes the barycentric Wachspress coordinates for a polygon with nVertices
!>  points in R3, vertCoords for a particular pointInterp with normalized radius.
!>  Follows Gillette, A., Rand, A., Bajaj, C., 2011.
!>  Error estimates for generalized barycentric interpolation.
!>  Advances in computational mathematics 37 (3), 417–439.
!>  Optimized version of mpas_wachspress_coordinates uses optional cached B_i areas
!------------------------------------------------------------------------
subroutine mpas_wachspress_coordinates(nVertices, vertCoords, & 
     pointInterp, mpaswach)
      
  implicit none
      
! input points
  integer, intent(in) :: nVertices
  real, dimension(3, nVertices), intent(in) :: vertCoords
  real, dimension(3), intent(in) :: pointInterp
!      real, dimension(nVertices), optional, intent(in) :: areaBin
      
! output
!      real, dimension(nVertices) :: mpas_wachspress_coordinates
  real, dimension(nVertices) :: mpaswach

! parameters
! The line below creates a real type with 12 decimal places
  integer, parameter :: RHI = selected_real_kind(12)
      
! computational intermediates
  real(kind=RHI), dimension(nVertices) :: wach ! The wachpress area-product
  real(kind=RHI) :: wach_total ! The wachpress total weight
  integer :: i, j           ! Loop indices
  integer :: im1, i0, ip1   ! im1 = (i-1), i0 = i, ip1 = (i+1)
      
! triangle areas to compute wachspress coordinate
  real, dimension(nVertices) :: areaA
  real, dimension(nVertices) :: areaB
  real :: radiusLocal
!      real :: mpas_triangle_signed_area_sphere 
      
  radiusLocal = sqrt(sum(vertCoords(:,1)**2))
      
!      if (.not. present(areaBin)) then
! compute areas
  do i = 1, nVertices
! compute first area B_i
! get vertex indices
     im1 = mod(nVertices + i - 2, nVertices) + 1
     i0  = mod(nVertices + i - 1, nVertices) + 1
     ip1 = mod(nVertices + i    , nVertices) + 1
            
! precompute B_i areas
! always the same because B_i independent of xp,yp,zp
! (COULD CACHE AND USE RESULT FROM ARRAY FOR FURTHER OPTIMIZATION)
     areaB(i) = mpas_triangle_signed_area_sphere(vertCoords(:, im1), & 
          vertCoords(:, i0), vertCoords(:, ip1), radiusLocal)
  end do
!      else
!! assign areas
!         do i = 1, nVertices
!            areaB(i) = areaBin(i)
!         end do
!      end if
      
! compute areas
  do i = 1, nVertices
!     compute first area B_i
            ! get vertex indices
     im1 = mod(nVertices + i - 2, nVertices) + 1
     i0  = mod(nVertices + i - 1, nVertices) + 1
     ip1 = mod(nVertices + i    , nVertices) + 1
         
!     compute A_ij areas
! must be computed each time
     areaA(i0) = mpas_triangle_signed_area_sphere(pointInterp, & 
          vertCoords(:, i0), vertCoords(:, ip1), radiusLocal)
         
! precomputed B_i areas, cached
  end do
      
! for each vertex compute wachpress coordinate
  do i = 1, nVertices
     wach(i) = areaB(i)
     do j = (i + 1), (i + nVertices - 2)
        i0  = mod(nVertices + j - 1, nVertices) + 1
!     accumulate products for A_ij subareas
        wach(i) = wach(i) * areaA(i0)
     end do
  end do
      
! get summed weights for normalization
  wach_total = 0
  do i = 1, nVertices
     wach_total = wach_total + wach(i)
  end do
      
! compute lambda
!      mpas_wachspress_coordinates= 0.0
  mpaswach= 0.0
  do i = 1, nVertices
!         mpas_wachspress_coordinates(i) = real(wach(i)/wach_total)
     mpaswach(i) = real(wach(i)/wach_total)
  end do

!      print*,'mpaswach=',mpaswach
      
!      end function mpas_wachspress_coordinates
end subroutine mpas_wachspress_coordinates
!----------------------------------------------------------------

!***********************************************************************
!
!  routine mpas_triangle_signed_area_sphere
!
!> \brief   Calculates area of a triangle on a sphere
!> \author  Matthew Hoffman
!> \date    13 January 2015
!> \details
!>  This routine calculates the area of a triangle on the surface of a sphere.
!>  Uses the spherical analog of Heron's formula.
!>  Copied from mesh generator.  A CCW winding angle is positive.
!-----------------------------------------------------------------------
real function mpas_triangle_signed_area_sphere(a, b, c, radius)

  implicit none
      
!-----------------------------------------------------------------
! input variables
!-----------------------------------------------------------------
  real, dimension(3), intent(in) :: a, b, c !< Input: 3d (x,y,z) points forming the triangle in which to calculate the bary weights
  real, intent(in) :: radius !< sphere radius
 
!-----------------------------------------------------------------
! local variables
!-----------------------------------------------------------------
  real :: ab, bc, ca, semiperim, tanqe
!  real :: mpas_arc_length
  real, dimension(3) :: ablen, aclen, Dlen
      
  ab = mpas_arc_length(a(1), a(2), a(3), b(1), b(2), b(3))/radius
  bc = mpas_arc_length(b(1), b(2), b(3), c(1), c(2), c(3))/radius
  ca = mpas_arc_length(c(1), c(2), c(3), a(1), a(2), a(3))/radius
  semiperim = 0.5 * (ab + bc + ca)
      
  tanqe = sqrt(max(0.0,tan(0.5 * semiperim) &
       * tan(0.5 * (semiperim - ab)) &
       * tan(0.5 * (semiperim - bc)) * tan(0.5 * (semiperim - ca))))
      
  mpas_triangle_signed_area_sphere = 4.0 * radius * & 
       radius * atan(tanqe)
      
! computing correct signs (in similar fashion to mpas_sphere_angle)
  ablen(1) = b(1) - a(1)
  ablen(2) = b(2) - a(2)
  ablen(3) = b(3) - a(3)
      
  aclen(1) = c(1) - a(1)
  aclen(2) = c(2) - a(2)
  aclen(3) = c(3) - a(3)
      
  dlen(1) =   (ablen(2) * aclen(3)) - (ablen(3) * aclen(2))
  dlen(2) = -((ablen(1) * aclen(3)) - (ablen(3) * aclen(1)))
  dlen(3) =   (ablen(1) * aclen(2)) - (ablen(2) * aclen(1))
      
  if ((Dlen(1)*a(1) + Dlen(2)*a(2) + Dlen(3)*a(3)) < 0.0) then
     mpas_triangle_signed_area_sphere =  & 
          -mpas_triangle_signed_area_sphere
  end if
      
end function mpas_triangle_signed_area_sphere

!------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FUNCTION MPAS_ARC_LENGTH
!
! Returns the length of the great circle arc from A=(ax, ay, az) to
!    B=(bx, by, bz). It is assumed that both A and B lie on the surface of the
!    same sphere centered at the origin.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function mpas_arc_length(ax, ay, az, bx, by, bz)
      
  implicit none
      
  real, intent(in) :: ax, ay, az, bx, by, bz
      
  real :: r, c
  real :: cx, cy, cz
      
  cx = bx - ax
  cy = by - ay
  cz = bz - az
      
  r = sqrt(ax*ax + ay*ay + az*az)
  c = sqrt(cx*cx + cy*cy + cz*cz)
      
  mpas_arc_length = r * 2.0 * asin(c/(2.0*r))
      
end function mpas_arc_length

!----------------------------------------------------------

subroutine convert_lx(lat, lon, radius, vec)

  implicit none
      
  real, intent(in) :: lat, lon, radius
      
  real, dimension(3) :: vec
      
  vec(1) = radius * cos(lon) * cos(lat)
  vec(2) = radius * sin(lon) * cos(lat)
  vec(3) = radius * sin(lat)
      
!      print*,'lat,lon,radius=',lat,lon,radius
      
!      print*,'vec=',vec

end subroutine convert_lx

!-----------------------------------------------------------------

function index2d(irank, idx) result(i)
      
  implicit none
      
  integer, intent(in) :: irank, idx
      
  integer :: i
      
  i = irank * (idx - 1) + 1
      
end function index2d

!==========================================================
Subroutine MPAS_NGradN &
  (X,        & ! <-- Cartesian coordinates of point
   NGradN,   & ! --> Interpolated (1 + N)*Grad(N)
   NP,       & ! --> Interpolated N
   Stat,     & ! --> Error status
   NIP)        ! ~~> Interpolated Im(N)
!
! Calculation of interpolated (1 + N)*Grad(N) for
! MPAS global fields.
!----------------------------------------------------------
! Method:
!   Calculation of gradient in geodetic coordinates and
!   transform to Cartesian coordinates.
!----------------------------------------------------------
! J. D. Hegarty 2020 based on:
! (C) Copyright 1999-2001, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 25 Jan 1999 | Original version.
!   2.0   | 18 Feb 1999 | Argument NP.
!   2.1   | 20 Feb 1999 | Relaxation to MSIS profile.
!   2.2   | 10 Mar 1999 | Stat=1 when G%H < Zmin.
!   3.0   | 22 Jun 1999 | No combination with N_MSIS.
!   4.0   | 15 May 2001 | Argument NIP.
!----------------------------------------------------------
! Modules used:
!
Use Coordinates, only: &
! Imported Type Definitions:
    Cartesian
!
Use Earth, only: &
! Imported Type Definitions:
    Geodetic,  &
! Imported Routines:
    Geod_from_Cart,  &
    Jacobian_GC
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
! Local Scalars:
!
Type(Geodetic) :: G     ! Geodetic coordinates of point
Real(Double)   :: PLon  ! Longitude of point
Real(Double)   :: PLat  ! Latitude of point
Real(Double)   :: ZP    ! Altitude of point
Real(Double)   :: Zmin  ! Minimum model Z
Real(Double)   :: Zmax  ! Maximum model Z
!
! Local Arrays:
!
Real(Double)   :: &
   NGP(3)     ! Interpolated grad(N) in geodetic coordinates
Real(Double)   :: &
   JGC(3,3)   ! Jacoubian d(Geodetic)/d(Cartesian)
!----------------------------------------------------------

!----------------------------------------------------------
! 1. INTERPOLATION OF NCEP-MSIS FIELD
!----------------------------------------------------------

G    = Geod_from_Cart(X)
PLon = G%Lambda
PLat = G%Phi
ZP   = G%H

!print*,'calling interpolate_refractivity_mpas'
!print*,'plon,plat,zp=',plon,plat,zp

Call Interpolate_Refractivity_MPAS &
  (PLon,     & ! <-- Longiude of point [deg]
   PLat,     & ! <-- Latitude of point [deg]
   ZP,       & ! <-- Altitude of point [km]
   Zmin,     & ! --> Minimum model Z for this lon/lat
   Zmax,     & ! --> Maximum model Z for this lon/lat
   NP,       & ! --> Interpolated N
   NGP,      & ! ~~> Interpolated dN/d(alt,lat,lon)
   NIP = NIP)  ! ~~> Interpolated Im(N)


!----------------------------------------------------------
! 2. TRANSFORM TO CARTESIAN COORDINATES
!----------------------------------------------------------

JGC    = Jacobian_GC(G)

NGradN%X(:) = (1+NP)*MatMul(NGP(:), JGC(:,:))


!----------------------------------------------------------
! 3. STATUS DEFINITION
!----------------------------------------------------------

If (G%H >= Zmin) then
   Stat = 0
Else
   Stat = 1
End If


End Subroutine MPAS_NGradN

!==========================================================
!Subroutine MPAS_NGHN &
!  (X,        & ! <-- Cartesian coordinates of point
!   NGN,      & ! --> Interpolated (1 + N)*Grad(N)
!   NP,       & ! --> Interpolated N
!   NHN,      & ! --> Interpolated Grad x (1+N)Grad(N)
!   Stat)       ! --> Error status
!!
!! Calculation of interpolated (1 + N)*Grad(N) for
!! MPAS global fields: Adjoint version.
!!----------------------------------------------------------
!! Method:
!!   Calculation of gradient in geodetic coordinates and
!!   transform to Cartesian coordinates.
!!----------------------------------------------------------
!! (C) Copyright 1999-2000, M. E. Gorbunov.
!! Modified for MPAS November 2020, J. Hegarty.
!!
!! Version |    Date     | Comment
!!---------+-------------+--------------------
!!   3.0   | 22 Jun 1999 | Basic non-adjoint version.
!!   1.0   | 03 May 2000 | Version with hessian.
!!----------------------------------------------------------
!! Modules used:
!!
!Use Coordinates, only: &
!! Imported Type Definitions:
!    Cartesian
!!
!Use Earth, only: &
!! Imported Type Definitions:
!    Geodetic,  &
!! Imported Routines:
!    Geod_from_Cart,  &
!    Jacobian_GC,     &
!    Hessian_GC
!!----------------------------------------------------------
!Implicit None
!!----------------------------------------------------------
!! Input arguments:
!!
!Type(Cartesian), Intent(In)  :: &
!   X               ! Cartesian coordinates of point
!!
!! Output arguments:
!!
!Type(Cartesian), Intent(Out) :: &
!   NGN             ! Interpolated (1 + N)*Grad(N)
!!
!Real(Double), Intent(Out)    :: &
!   NP              ! Interpolated N
!!
!Real(Double), Intent(Out)    :: &
!   NHN(3,3)        ! Interpolated Grad x (1+N)Grad(N)
!!
!Integer, Intent(Out)         :: &
!   Stat            ! Error status:
!                   !   0 - point above surface
!                   !   1 - point under surface
!!----------------------------------------------------------
!! Local Scalars:
!!
!Type(Geodetic) :: G     ! Geodetic coordinates of point
!Real(Double)   :: PLon  ! Longitude of point
!Real(Double)   :: PLat  ! Latitude of point
!Real(Double)   :: ZP    ! Altitude of point
!Real(Double)   :: Zmin  ! Minimum model Z
!Real(Double)   :: Zmax  ! Maximum model Z
!Integer        :: i, j  ! Dimension indexes
!!
!! Local Arrays:
!!
!Real(Double)   :: &
!   NG(3)        ! Interpolated Grad(N) in geodetic coordinates
!Real(Double)   :: &
!   NH(3,3)      ! Interpolated Hess(N) in geodetic coordinates
!Real(Double)   :: &
!   JGC(3,3),  & ! Jacobian d(Geodetic)/d(Cartesian)
!   HGC(3,3,3)   ! Hessian d2(Geodetic)/d(Cartesian)2
!!----------------------------------------------------------

!!----------------------------------------------------------
!! 1. INTERPOLATION OF NCEP-MSIS FIELD
!!----------------------------------------------------------

!!--- 1.1. Coordinate calculation
!
!G    = Geod_from_Cart(X)
!PLon = G%Lambda
!PLat = G%Phi
!ZP   = G%H
!
!
!!--- 1.2. Interpolation and calculation
!!---      of T, P, and Q derivatives
!
!Call Interpolate_Refractivity_MPAS &
!  (PLon,     & ! <-- Longiude of point [deg]
!   PLat,     & ! <-- Latitude of point [deg]
!   ZP,       & ! <-- Altitude of point [km]
!   Zmin,     & ! --> Minimum model Z for this lon/lat
!   Zmax,     & ! --> Maximum model Z for this lon/lat
!   NP,       & ! --> Interpolated N
!   NG,       & ! --> Interpolated dN/d(alt,lat,lon)
!   NH)         ! --> Interpolated hessian matrix of N

!!----------------------------------------------------------
!! 2. TRANSFORM TO CARTESIAN COORDINATES
!!----------------------------------------------------------
!
!
!!--- 2.1. Transform of gradient
!
!JGC    = Jacobian_GC(G)
!
!NGN%X(:) = (1+NP)*MatMul(NG(:), JGC(:,:))
!
!
!!--- 2.2. Transform of hessian
!
!HGC = Hessian_GC(G)
!
!Do i=1,3
!   Do j=1,3
!      NHN(i,j) = &
!         Sum(JGC(:,i)*NG(:))*Sum(JGC(:,j)*NG(:)) + &
!         (1+NP)*Sum(HGC(:,i,j)*NG(:))            + &
!         (1+NP)*Sum(JGC(:,i)*MatMul(NH(:,:),JGC(:,j)))
!   End Do
!End Do

!!----------------------------------------------------------
!! 3. STATUS DEFINITION
!!----------------------------------------------------------
!
!If (G%H >= Zmin) then
!   Stat = 0
!Else
!   Stat = 1
!End If
!
!
!End Subroutine MPAS_NGHN

!==========================================================
Subroutine Get_NETCDF_MPAS &
  (PathName,  & ! <-- Pathname of data file
   Vrb,       & ! <-- Verbosity level
   Data_Type, & ! --> Data representation type
   Fields,    & ! <-> Structure with fields
   xtime,      & ! --> Field date
   Time,      & ! --> Field time
   ErrStat)     ! --> Error code

!   XLon,      & ! --> Longitude grid [deg]
!   XLat,      & ! --> Latitude grid [deg]

!
! Reading MPAS global fields from NETCDF data file.
!----------------------------------------------------------
! J. D. Hegarty, 2020 based on (C) Copyright 2007-2009, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 26 Jun 2007 | Original version (from ECHAM_fields).
!   1.1   | 29 Mar 2008 | Using grid limits from GDS.
!   2.0   | 05 May 2008 | Vrb.
!   2.1   | 30 Mar 2009 | Corrected error code.
!----------------------------------------------------------
! Modules used:
!
!---- JDH 120120
! I don't link in the netCDF libraries developed from Gorbunov in the .mak
! files.  Perhaps these are standard routines that come from my netCDF library.
Use NetCDF, only: &
! Imported Parameters:
    NF90_NoWrite,             &
    NF90_NoErr,               &
    NF90_BYTE,                &
    NF90_CHAR,                &
    NF90_SHORT,               &
    NF90_INT,                 &
    NF90_FLOAT,               &
    NF90_DOUBLE,              &
    NF90_GLOBAL,              &
! Do I need to specifically import all the routines that I use?
! Imported Routines:
    NF90_Open,                &
    NF90_Inquire,             &
    NF90_Inquire_Dimension,   &
    NF90_INQ_DIMID,           &
    NF90_Inquire_Variable,    &
    NF90_INQ_VARID,           &
    NF90_Inq_Attname,         &
    NF90_Inquire_Attribute,   &
    NF90_Get_Att,             &
    NF90_Get_Var,             &
    NF90_Close

!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
!----------------------------------------------------------
! Public Parameters:
!
Integer, Parameter :: &
   err_NoFile   = 0702001,  & ! File not found
   err_MissVar  = 0702002,  & ! Missing variable/attributes
   err_ReadErr  = 0702003,  & ! Data read error
   err_FmtErr   = 0702004,  & ! Data format error
   err_NoLCF    = 0702005,  & ! No LCF
   err_NoFreq   = 0702006,  & ! No frequencies
   err_WrtErr   = 0702007,  & ! Write error
   err_DataType = 0702008,  & ! Bad data type
   err_NoData   = 0702009,  & ! No data
   err_Dimen    = 0702010,  & ! Inconsistent dimensions
   err_Channel  = 0702011,  & ! Inconsistent numbers of channels
   err_Coord    = 0702012,  & ! Bad coordinates
   err_FileName = 0702013     ! Bad file name
!----------------------------------------------------------
! Input arguments:
!
Character(Len=*), Intent(In)     :: &
   PathName     ! Pathname of data file
!
Integer, Intent(In)     :: &
   Vrb          ! Verbosity level
!
! Output arguments:
!
Integer, Intent(Out)  :: &
   Data_Type    ! Data representation type
!
Type(T_MPAS_Fields), Intent(InOut) :: &
   Fields       ! Structure with fields
!
!Real(Single), Pointer :: &
!   XLon(:,:,:)  ! Longitude grid [deg]
!
! Indices of XLon:
!   ILon, 1,  1  - for rectangular grid
!   j1,   j2, jd - for icosahedral
!
!Real(Single), Pointer :: &
!   XLat(:,:,:)  ! Latitude grid [deg]
!
! Indices of XLat:
!   1,  ILat, 1  - for rectangular grid
!   j1, j2,   jd - for icosahedral
!
!
!Integer       :: &
!   Date         ! Field date [YYYYMMDD]
                ! MM = 01..12 - month
                ! MM = 13..16 - season (DJF,MAM,JJA,SON)

Character(len=19) :: xtime
!
Integer       :: &
   Time         ! Field time [hhmm]
!
Type(Error_Status), Pointer :: &
   ErrStat      ! Error status
!----------------------------------------------------------
! Local Scalars:
!
! --- Work variables
!
Integer              :: IOStat      ! I/O operation status
Integer              :: Stat      ! I/O operation status
! try
!Integer, target      :: NCID
Integer              :: NCID        ! Identifier of NetCDF file
Integer              :: NDims       ! Number of dimensions
Integer              :: NVars       ! Number of variables
Integer              :: NEVars      ! Number of extracted variables
Integer              :: NAttrs      ! Number of global attributes
Integer              :: NEAttrs     ! Number of extracted attributes
Integer              :: UnlimID     ! ID of unlimited dimension
Integer              :: IDim        ! Dimension index
Integer              :: IVar        ! Variable index
Integer              :: IAttr       ! Atrtribute index
Integer              :: i           ! Work index
Integer            :: ErrCode      ! Error code
Character(Len=256) :: ZPathName    ! 0-terminated pathname
Character(Len=256) :: ErrMsg       ! Error message
Integer            :: N            ! Total grid points
Integer            :: NLon         ! Longitudinal dimension
Integer            :: NLat         ! Latitudinal dimension
Integer            :: NSP          ! Triangular spectral truncation
Integer            :: NI           ! Triangle per diamond side
Integer            :: ND           ! Nubmer of diamonds
Integer            :: NLev         ! Number of levels
Integer            :: La1          ! Start latitude
Integer            :: Lo1          ! Start longitude
Integer            :: La2          ! Final latitude
Integer            :: Lo2          ! Final longitude
Integer            :: Center       ! Center ID
Integer            :: Subcenter    ! Subcenter ID
Integer            :: Grid         ! Grid ID
Integer            :: Code         ! Parameter/units ID
Integer            :: varid        ! netcdf variable ID
Character(Len=255) :: Name         ! Parameter name
Character(Len=255) :: Comment      ! Parameter comment
Integer            :: Level        ! Model level
Integer            :: LevelType    ! Level type
Integer            :: Table        ! Parameter table ID
Integer            :: Year         ! Year
Integer            :: Month        ! Month of year
Integer            :: Day          ! Day of month
Integer            :: Hour         ! Hour of day
Integer            :: Minute       ! Minute of hour
Integer            :: NV           ! Number of vertical coordinate parameters
Logical            :: HybridLevel  ! Hybrid level indicator
Integer            :: ILev         ! Level index
Integer, dimension(13)  :: dimids, lend
!Integer, Pointer   :: DimId(:)    ! Dimension indices
!Integer,           :: DimId(:)    ! Dimension indices
Integer            :: XType       ! External type

!
! --- Parameter codes
!
Integer   :: fcd_Geopotential   ! Surface geopotential code
Integer   :: fcd_Temperature    ! Temperature code
Integer   :: fcd_Humidity       ! Humidity code
Integer   :: fcd_Pressure       ! Surface pressure code
Integer   :: fcd_LogPressure    ! Ln(Surface pressure) code
Integer   :: fcd_U              ! U-velocity code
Integer   :: fcd_V              ! V-velocity code
Integer   :: fcd_Ozone          ! Ozone mixing ratio code
Integer   :: fcd_LWC            ! Liquid water content code
Integer   :: fcd_IWC            ! Ice water content code
!
! Local Arrays:
!
Real, Allocatable :: &
   F(:,:,:)                     ! Decoded GRIB data

! MPAS arrays
!integer, allocatable :: nedgesoncell(:)
!integer, allocatable :: cellsonedge(:,:)
!integer, allocatable :: edgesoncell(:,:)
!integer, allocatable :: cellsoncell(:,:)
!integer, allocatable :: verticesoncell(:,:)
!integer, allocatable :: cellsonvertex (:,:)

!real, allocatable :: latcell(:)
!real, allocatable :: loncell(:)
!real, allocatable :: latvertex(:)
!real, allocatable :: lonvertex(:)
!real, allocatable :: latedge(:)
!real, allocatable :: lonedge(:)
!real, allocatable :: rdzw(:)
!real, allocatable :: fzp(:)
!real, allocatable :: fzm(:)

!----------------------------------------------------------
! 0. INITIALIZATION
!----------------------------------------------------------

Call Enter_Callee &
  ('Get_NETCDF_MPAS',   & ! <-- User routine
   ErrStat)        ! <-> Pointer to callee status

!----------------------------------------------------------
! 1. OPENING NETCDF File
!----------------------------------------------------------

Stat    = 0
N       = 0
NEVars  = 0
NEAttrs = 0

!--- 1.1. Opening file

print*,'In GET_NETCDF_MPAS pathname = ', pathname

IOStat = NF90_Open   &
  (PathName,             &
   NF90_NoWrite,     &
   NCID)

print*,'NF90_NoErr,IOStat=',NF90_NoErr,IOStat

If (IOStat /= NF90_NoErr) then
   Stat = err_NoFile
   Return
End If

!fields%ncid => ncid
!fields%ncid_sav => ncid
! make just a regular integer
fields%ncid_sav = ncid

!print*,'in get_netdcf ncid=',ncid,fields%ncid
print*,'in get_netdcf ncid=',ncid,fields%ncid_sav

!Write(ZPathName, '(2A)') Trim(PathName), Char(0)

!--- 1.2. Inquiring file contents
IOStat = NF90_Inquire(NCID,NDims,NVars,NAttrs,UnlimID)

If (IOStat /= NF90_NoErr) then
   Stat = err_ReadErr
   Return
End If

!If (LVrb >= 4) then
   Write(*,'(4(A,I5:/))')                 &
      'Number of dimensions:   ', NDims,  &
      'Number of variables:    ', NVars,  &
      'Number of attributes:   ', NAttrs, &
      'Unlimited dimension ID: ', Unlimid
!End If


!------------------ xtime ----------------------
name='xtime'
IOStat=NF90_INQ_VARID(NCID,name,varid)
If (IOStat /= NF90_NoErr) then
   print*,'nf90_inq_varid IOStat=',IOStat
   print*,'Get_NETCDF_MPAS name=',name
   Stat = err_ReadErr
   Return
End If
print*,'varid=',varid
   IOStat = NF90_Inquire_Variable(NCID,varid,name,xtype,ndims,dimids,nattrs )
   If (IOStat /= NF90_NoErr) then
      Stat = err_ReadErr
      Return
   End If

print*,'after NF90_Inquire_Variable, ndims,nattrss,name=',ndims,nattrs,name,dimids

do nd=1,ndims

   IOStat=NF90_INQUIRE_DIMENSION(NCID,dimids(nd),name,lend(nd))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_inq_dimension IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If

   print*,'name,lend=',name,lend(nd)

end do

!allocate(fields%latcell(lend(1)))

   IOStat=NF90_GET_VAR(NCID,varid,xtime)
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_get_var IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If

!   print*,'fields%latcell=',fields%latcell(10),fields%latcell(100)
   print*,'xtime=',xtime

!   if(ndims.eq.2) then
!      print*,'debug for xtime, stopping'
!   end if

!--------------------------------------------------------------------

!------------------ latCell ----------------------
name='latCell'
IOStat=NF90_INQ_VARID(NCID,name,varid)
If (IOStat /= NF90_NoErr) then
   print*,'nf90_inq_varid IOStat=',IOStat
   print*,'GET_NETCDF_MPAS name=', name
   Stat = err_ReadErr
   Return
End If
print*,'varid=',varid
! The first call to NF90_Inquire _Variable gets the
! number of dimensions, the second call gets the dimension
! IDs from there I think I can get the actual diminesion perhaps using
! Inquire dimension.  Try to follow scan_inpur.F in convert_mpas.F code.
   IOStat = NF90_Inquire_Variable(NCID,varid,name,xtype,ndims,dimids,nattrs )
   If (IOStat /= NF90_NoErr) then
      Stat = err_ReadErr
      Return
   End If

print*,'after NF90_Inquire_Variable, ndims,nattrss,name=',ndims,nattrs,name,dimids

do nd=1,ndims

   IOStat=NF90_INQUIRE_DIMENSION(NCID,dimids(nd),name,lend(nd))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_inq_dimension IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If

   print*,'name,lend=',name,lend(nd)

end do

!allocate(latcell(lend(1)))
allocate(fields%latcell(lend(1)))

!   IOStat=NF90_GET_VAR(NCID,varid,latcell)
   IOStat=NF90_GET_VAR(NCID,varid,fields%latcell)
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_get_var IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
!   print*,'latcell=',latcell(10),latcell(100)
!   print*,'fields%latcell=',fields%latcell(10),fields%latcell(100)

!--------------------------------------------------------------------


!------------------ lonCell ----------------------
name='lonCell'
IOStat=NF90_INQ_VARID(NCID,name,varid)
If (IOStat /= NF90_NoErr) then
   print*,'nf90_inq_varid IOStat=',IOStat
   print*,'GET_NETCDF_MPAS name=',name
   Stat = err_ReadErr
   Return
End If
print*,'varid=',varid
! The first call to NF90_Inquire _Variable gets the
! number of dimensions, the second call gets the dimension
! IDs from there I think I can get the actual diminesion perhaps using
! Inquire dimension.  Try to follow scan_inpur.F in convert_mpas.F code.
   IOStat = NF90_Inquire_Variable(NCID,varid,name,xtype,ndims,dimids,nattrs )
   If (IOStat /= NF90_NoErr) then
      Stat = err_ReadErr
      Return
   End If

print*,'after NF90_Inquire_Variable, ndims,nattrss,name=',ndims,nattrs,name,dimids

do nd=1,ndims

   IOStat=NF90_INQUIRE_DIMENSION(NCID,dimids(nd),name,lend(nd))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_inq_dimension IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If

   print*,'name,lend=',name,lend(nd)

end do

!allocate(loncell(lend(1)))
allocate(fields%loncell(lend(1)))

!   IOStat=NF90_GET_VAR(NCID,varid,loncell)
   IOStat=NF90_GET_VAR(NCID,varid,fields%loncell)
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_get_var IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
!   print*,'loncell=',loncell(10),loncell(100)
!   print*,'loncell=',fields%loncell(10),fields%loncell(100)

!--------------------------------------------------------------------

!------------------ latVertex ----------------------
name='latVertex'
IOStat=NF90_INQ_VARID(NCID,name,varid)
If (IOStat /= NF90_NoErr) then
   print*,'nf90_inq_varid IOStat=',IOStat
   print*,'GET_NETCDF_MPAS name=',name
   Stat = err_ReadErr
   Return
End If
print*,'varid=',varid
   IOStat = NF90_Inquire_Variable(NCID,varid,name,xtype,ndims,dimids,nattrs )
   If (IOStat /= NF90_NoErr) then
      Stat = err_ReadErr
      Return
   End If

print*,'after NF90_Inquire_Variable, ndims,nattrss,name=',ndims,nattrs,name,dimids

do nd=1,ndims

   IOStat=NF90_INQUIRE_DIMENSION(NCID,dimids(nd),name,lend(nd))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_inq_dimension IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If

   print*,'name,lend=',name,lend(nd)

end do

!allocate(latvertex(lend(1)))
allocate(fields%latvertex(lend(1)))

!   IOStat=NF90_GET_VAR(NCID,varid,latvertex)
   IOStat=NF90_GET_VAR(NCID,varid,fields%latvertex)
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_get_var IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
!   print*,'latvertex=',latvertex(10),latvertex(100)
!   print*,'latvertex=',fields%latvertex(10),fields%latvertex(100)

!--------------------------------------------------------------------

!------------------ lonVertex ----------------------
name='lonVertex'
IOStat=NF90_INQ_VARID(NCID,name,varid)
If (IOStat /= NF90_NoErr) then
   print*,'nf90_inq_varid IOStat=',IOStat
   print*,'GET_NETCDF_MPAS name=',name
   Stat = err_ReadErr
   Return
End If
print*,'varid=',varid
   IOStat = NF90_Inquire_Variable(NCID,varid,name,xtype,ndims,dimids,nattrs )
   If (IOStat /= NF90_NoErr) then
      Stat = err_ReadErr
      Return
   End If
print*,'after NF90_Inquire_Variable, ndims,nattrss,name=',ndims,nattrs,name,dimids
do nd=1,ndims
   IOStat=NF90_INQUIRE_DIMENSION(NCID,dimids(nd),name,lend(nd))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_inq_dimension IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
   print*,'name,lend=',name,lend(nd)
end do

!allocate(lonvertex(lend(1)))
allocate(fields%lonvertex(lend(1)))

!   IOStat=NF90_GET_VAR(NCID,varid,lonvertex)
   IOStat=NF90_GET_VAR(NCID,varid,fields%lonvertex)
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_get_var IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
!   print*,'lonvertex=',lonvertex(10),lonvertex(100)
!   print*,'lonvertex=',fields%lonvertex(10),fields%lonvertex(100)

!--------------------------------------------------------------------

!------------------ latEdge ----------------------
name='latEdge'
IOStat=NF90_INQ_VARID(NCID,name,varid)
If (IOStat /= NF90_NoErr) then
   print*,'nf90_inq_varid IOStat=',IOStat
   print*,'GET_NETCDF_MPAS name=',name
   Stat = err_ReadErr
   Return
End If
print*,'varid=',varid
   IOStat = NF90_Inquire_Variable(NCID,varid,name,xtype,ndims,dimids,nattrs )
   If (IOStat /= NF90_NoErr) then
      Stat = err_ReadErr
      Return
   End If
print*,'after NF90_Inquire_Variable, ndims,nattrss,name=',ndims,nattrs,name,dimids
do nd=1,ndims
   IOStat=NF90_INQUIRE_DIMENSION(NCID,dimids(nd),name,lend(nd))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_inq_dimension IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
   print*,'name,lend=',name,lend(nd)
end do

!allocate(latedge(lend(1)))
allocate(fields%latedge(lend(1)))

!   IOStat=NF90_GET_VAR(NCID,varid,latedge)
   IOStat=NF90_GET_VAR(NCID,varid,fields%latedge)
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_get_var IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
!   print*,'latedge=',latedge(10),latedge(100)
!   print*,'latedge=',fields%latedge(10),fields%latedge(100)

!--------------------------------------------------------------------

!------------------ lonEdge ----------------------
name='lonEdge'
IOStat=NF90_INQ_VARID(NCID,name,varid)
If (IOStat /= NF90_NoErr) then
   print*,'nf90_inq_varid IOStat=',IOStat
   print*,'GET_NETCDF_MPAS name=',name
   Stat = err_ReadErr
   Return
End If
print*,'varid=',varid
   IOStat = NF90_Inquire_Variable(NCID,varid,name,xtype,ndims,dimids,nattrs )
   If (IOStat /= NF90_NoErr) then
      Stat = err_ReadErr
      Return
   End If
print*,'after NF90_Inquire_Variable, ndims,nattrss,name=',ndims,nattrs,name,dimids
do nd=1,ndims
   IOStat=NF90_INQUIRE_DIMENSION(NCID,dimids(nd),name,lend(nd))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_inq_dimension IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
   print*,'name,lend=',name,lend(nd)
end do

!allocate(lonedge(lend(1)))
allocate(fields%lonedge(lend(1)))

!   IOStat=NF90_GET_VAR(NCID,varid,lonedge)
   IOStat=NF90_GET_VAR(NCID,varid,fields%lonedge)
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_get_var IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
!   print*,'lonedge=',lonedge(10),lonedge(100)
!   print*,'lonedge=',fields%lonedge(10),fields%lonedge(100)

!--------------------------------------------------------------------
!------------------ nEdgesOnCell ----------------------
name='nEdgesOnCell'
IOStat=NF90_INQ_VARID(NCID,name,varid)
If (IOStat /= NF90_NoErr) then
   print*,'nf90_inq_varid IOStat=',IOStat
   print*,'GET_NETCDF_MPAS name=',name
   Stat = err_ReadErr
   Return
End If
print*,'varid=',varid
   IOStat = NF90_Inquire_Variable(NCID,varid,name,xtype,ndims,dimids,nattrs )
   If (IOStat /= NF90_NoErr) then
      Stat = err_ReadErr
      Return
   End If
print*,'after NF90_Inquire_Variable, ndims,nattrss,name=',ndims,nattrs,name,dimids
do nd=1,ndims
   IOStat=NF90_INQUIRE_DIMENSION(NCID,dimids(nd),name,lend(nd))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_inq_dimension IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
   print*,'name,lend=',name,lend(nd)
end do

!allocate(nedgesoncell(lend(1)))
allocate(fields%nedgesoncell(lend(1)))

   IOStat=NF90_GET_VAR(NCID,varid,fields%nedgesoncell)
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_get_var IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If

!   print*,'nedgesoncell=',nedgesoncell(10),nedgesoncell(100)
!   print*,'nedgesoncell=',fields%nedgesoncell(10),fields%nedgesoncell(100)

!--------------------------------------------------------------------

!------------------ cellsOnEdge ----------------------
name='cellsOnEdge'
IOStat=NF90_INQ_VARID(NCID,name,varid)
If (IOStat /= NF90_NoErr) then
   print*,'nf90_inq_varid IOStat=',IOStat
   print*,'GET_NETCDF_MPAS name=',name
   Stat = err_ReadErr
   Return
End If
print*,'varid=',varid
   IOStat = NF90_Inquire_Variable(NCID,varid,name,xtype,ndims,dimids,nattrs )
   If (IOStat /= NF90_NoErr) then
      Stat = err_ReadErr
      Return
   End If
print*,'after NF90_Inquire_Variable, ndims,nattrss,name=',ndims,nattrs,name,dimids
do nd=1,ndims
   IOStat=NF90_INQUIRE_DIMENSION(NCID,dimids(nd),name,lend(nd))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_inq_dimension IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
   print*,'name,lend=',name,lend(nd)
end do

!allocate(cellsonedge(lend(1),lend(2)))
allocate(fields%cellsonedge(lend(1),lend(2)))

   IOStat=NF90_GET_VAR(NCID,varid,fields%cellsonedge)
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_get_var IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
!   print*,'cellsonedge=',cellsonedge(1,10),cellsonedge(1,100)
!   print*,'cellsonedge=',fields%cellsonedge(1,10),fields%cellsonedge(1,100)

!--------------------------------------------------------------------

!------------------ edgesOnCell ----------------------
name='edgesOnCell'
IOStat=NF90_INQ_VARID(NCID,name,varid)
If (IOStat /= NF90_NoErr) then
   print*,'nf90_inq_varid IOStat=',IOStat
   print*,'GET_NETCDF_MPAS name=',name
   Stat = err_ReadErr
   Return
End If
print*,'varid=',varid
   IOStat = NF90_Inquire_Variable(NCID,varid,name,xtype,ndims,dimids,nattrs )
   If (IOStat /= NF90_NoErr) then
      Stat = err_ReadErr
      Return
   End If
print*,'after NF90_Inquire_Variable, ndims,nattrss,name=',ndims,nattrs,name,dimids
do nd=1,ndims
   IOStat=NF90_INQUIRE_DIMENSION(NCID,dimids(nd),name,lend(nd))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_inq_dimension IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
   print*,'name,lend=',name,lend(nd)
end do

!allocate(edgesoncell(lend(1),lend(2)))
allocate(fields%edgesoncell(lend(1),lend(2)))

   IOStat=NF90_GET_VAR(NCID,varid,fields%edgesoncell)
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_get_var IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If

!   print*,'edgesoncell=',edgesoncell(1,10),edgesoncell(1,100)
!   print*,'edgesoncell=',fields%edgesoncell(1,10),fields%edgesoncell(1,100)

!--------------------------------------------------------------------

!------------------ cellsOnCell ----------------------
name='cellsOnCell'
IOStat=NF90_INQ_VARID(NCID,name,varid)
If (IOStat /= NF90_NoErr) then
   print*,'nf90_inq_varid IOStat=',IOStat
   print*,'GET_NETCDF_MPAS name=',name
   Stat = err_ReadErr
   Return
End If
print*,'varid=',varid
   IOStat = NF90_Inquire_Variable(NCID,varid,name,xtype,ndims,dimids,nattrs )
   If (IOStat /= NF90_NoErr) then
      Stat = err_ReadErr
      Return
   End If
print*,'after NF90_Inquire_Variable, ndims,nattrss,name=',ndims,nattrs,name,dimids
do nd=1,ndims
   IOStat=NF90_INQUIRE_DIMENSION(NCID,dimids(nd),name,lend(nd))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_inq_dimension IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
   print*,'name,lend=',name,lend(nd)
end do

!allocate(cellsoncell(lend(1),lend(2)))
allocate(fields%cellsoncell(lend(1),lend(2)))

!   IOStat=NF90_GET_VAR(NCID,varid,cellsoncell)
   IOStat=NF90_GET_VAR(NCID,varid,fields%cellsoncell)
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_get_var IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
!   print*,'cellsoncell=',cellsoncell(1,10),cellsoncell(1,100)
!   print*,'cellsoncell=',fields%cellsoncell(1,10),fields%cellsoncell(1,100)

!-----------------------------------------------------------------

!------------------ verticesOnCell ----------------------
name='verticesOnCell'
IOStat=NF90_INQ_VARID(NCID,name,varid)
If (IOStat /= NF90_NoErr) then
   print*,'nf90_inq_varid IOStat=',IOStat
   print*,'GET_NETCDF_MPAS name=',name
   Stat = err_ReadErr
   Return
End If
print*,'varid=',varid
   IOStat = NF90_Inquire_Variable(NCID,varid,name,xtype,ndims,dimids,nattrs )
   If (IOStat /= NF90_NoErr) then
      Stat = err_ReadErr
      Return
   End If
print*,'after NF90_Inquire_Variable, ndims,nattrss,name=',ndims,nattrs,name,dimids
do nd=1,ndims
   IOStat=NF90_INQUIRE_DIMENSION(NCID,dimids(nd),name,lend(nd))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_inq_dimension IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
   print*,'name,lend=',name,lend(nd)
end do

!allocate(verticesoncell(lend(1),lend(2)))
allocate(fields%verticesoncell(lend(1),lend(2)))

   IOStat=NF90_GET_VAR(NCID,varid,fields%verticesoncell)
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_get_var IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
!   print*,'verticesoncell=',verticesoncell(1,10),verticesoncell(1,100)
!   print*,'verticesoncell=',fields%verticesoncell(1,10),fields%verticesoncell(1,100)

!-----------------------------------------------------------------

!------------------ cellsOnVertex ----------------------
name='cellsOnVertex'
IOStat=NF90_INQ_VARID(NCID,name,varid)
If (IOStat /= NF90_NoErr) then
   print*,'nf90_inq_varid IOStat=',IOStat
   print*,'GET_NETCDF_MPAS name=',name
   Stat = err_ReadErr
   Return
End If
print*,'varid=',varid
   IOStat = NF90_Inquire_Variable(NCID,varid,name,xtype,ndims,dimids,nattrs )
   If (IOStat /= NF90_NoErr) then
      Stat = err_ReadErr
      Return
   End If
print*,'after NF90_Inquire_Variable, ndims,nattrss,name=',ndims,nattrs,name,dimids
do nd=1,ndims
   IOStat=NF90_INQUIRE_DIMENSION(NCID,dimids(nd),name,lend(nd))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_inq_dimension IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
   print*,'name,lend=',name,lend(nd)
end do

!allocate(cellsonvertex(lend(1),lend(2)))
allocate(fields%cellsonvertex(lend(1),lend(2)))

   IOStat=NF90_GET_VAR(NCID,varid,fields%cellsonvertex)
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_get_var IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If

!   print*,'cellsonvertex=',cellsonvertex(1,10),cellsonvertex(1,100)
!   print*,'cellsonvertex=',fields%cellsonvertex(1,10),fields%cellsonvertex(1,100)

!-----------------------------------------------------------------

!------------------ rdzw ----------------------
! This field is never used  - not in some history files
!name='rdzw'
!IOStat=NF90_INQ_VARID(NCID,name,varid)
!If (IOStat /= NF90_NoErr) then
!   print*,'nf90_inq_varid IOStat=',IOStat
!   print*,'GET_NETCDF_MPAS name=',name
!
!   Stat = err_ReadErr
!   Return
!End If
!print*,'varid=',varid
!   IOStat = NF90_Inquire_Variable(NCID,varid,name,xtype,ndims,dimids,nattrs )
!   If (IOStat /= NF90_NoErr) then
!      Stat = err_ReadErr
!      Return
!   End If
!print*,'after NF90_Inquire_Variable, ndims,nattrss,name=',ndims,nattrs,name,dimids
!do nd=1,ndims
!   IOStat=NF90_INQUIRE_DIMENSION(NCID,dimids(nd),name,lend(nd))
!   If (IOStat /= NF90_NoErr) then
!      print*,'nf90_inq_dimension IOStat=',IOStat
!      Stat = err_ReadErr
!      Return
!   End If
!   print*,'name,lend=',name,lend(nd)
!end do

!!allocate(rdzw(lend(1)))
!allocate(fields%rdzw(lend(1)))
!
!   IOStat=NF90_GET_VAR(NCID,varid,fields%rdzw)
!   If (IOStat /= NF90_NoErr) then
!      print*,'nf90_get_var IOStat=',IOStat
!      Stat = err_ReadErr
!      Return
!   End If
!!   print*,'rdzw=',rdzw(10),rdzw(100)
!   print*,'rdzw=',fields%rdzw(10),fields%rdzw(100)

!-----------------------------------------------------------------

!------------------ fzp ----------------------
name='fzp'
IOStat=NF90_INQ_VARID(NCID,name,varid)
If (IOStat /= NF90_NoErr) then
   print*,'nf90_inq_varid IOStat=',IOStat
   print*,'GET_NETCDF_MPAS name=',name

   Stat = err_ReadErr
   Return
End If
print*,'varid=',varid
   IOStat = NF90_Inquire_Variable(NCID,varid,name,xtype,ndims,dimids,nattrs )
   If (IOStat /= NF90_NoErr) then
      Stat = err_ReadErr
      Return
   End If
print*,'after NF90_Inquire_Variable, ndims,nattrss,name=',ndims,nattrs,name,dimids
do nd=1,ndims
   IOStat=NF90_INQUIRE_DIMENSION(NCID,dimids(nd),name,lend(nd))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_inq_dimension IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
   print*,'name,lend=',name,lend(nd)
end do

!allocate(fzp(lend(1)))
allocate(fields%fzp(lend(1)))

   IOStat=NF90_GET_VAR(NCID,varid,fields%fzp)
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_get_var IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
!   print*,'fzp=',fzp(10),fzp(100)
!   print*,'fzp=',fields%fzp(10),fields%fzp(100)

!-----------------------------------------------------------------
!------------------ fzm ----------------------
name='fzm'
IOStat=NF90_INQ_VARID(NCID,name,varid)
If (IOStat /= NF90_NoErr) then
   print*,'nf90_inq_varid IOStat=',IOStat
   print*,'GET_NETCDF_MPAS name=',name

   Stat = err_ReadErr
   Return
End If
print*,'varid=',varid
   IOStat = NF90_Inquire_Variable(NCID,varid,name,xtype,ndims,dimids,nattrs )
   If (IOStat /= NF90_NoErr) then
      Stat = err_ReadErr
      Return
   End If
print*,'after NF90_Inquire_Variable, ndims,nattrss,name=',ndims,nattrs,name,dimids
do nd=1,ndims
   IOStat=NF90_INQUIRE_DIMENSION(NCID,dimids(nd),name,lend(nd))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_inq_dimension IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
   print*,'name,lend=',name,lend(nd)
end do

!allocate(fzm(lend(1)))
allocate(fields%fzm(lend(1)))

   IOStat=NF90_GET_VAR(NCID,varid,fields%fzm)
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_get_var IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
!   print*,'fzm=',fzm(10),fzm(100)
!   print*,'fzm=',fields%fzm(10),fields%fzm(100)
!--------------------------------------------------------------------

!   print*,' get_netcdf_mpas subroutine under construction ---- Returning'

   print*,'end of get_netcdf_mpas --- fields%ncid_sav,ncid=',fields%ncid_sav,ncid

Return

End Subroutine Get_NETCDF_MPAS

!==========================================================
Subroutine Vertical_Interpolation_LnN_MPAS &
  (clat,     & ! latitude of MPAS cell (deg)
   clon,     & ! longitude of MPAS cell (deg)
   Hsur,     & ! surface geopotential[m]
   H,        & ! geopotential [m]
   T,        & ! temperature [k]
   QR,       & ! specific humidity [m]
   P,        & ! Pressure [mb]
   ZP,       & ! <-- Altitude [km]
   Zmin,     & ! --> Minimum model Z for this lon/lat
   Zmax,     & ! --> Maximum model Z for this lon/lat
   LnNP,     & ! --> Interpolated Ln(N)
   LnNZ,     & ! --> Interpolated dLn(N)/dZ
   LnNZ2,    & ! --> Interpolated d2Ln(N)/dZ2
   LnNIP)      ! ~~> Interpolated Ln(NI)
!
! Vertical interpolation of refractivity profiles from NCEP.
!----------------------------------------------------------
! Method:
!   Vertical spline interpolation to given altitude.
!----------------------------------------------------------
! (C) Copyright 1998-2004, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 12 Apr 2000 | Extracted from
!         |             | Interpolate_Refractivity.
!   2.0   | 14 May 2001 | Optional argument LnNIP.
!   3.0   | 07 Apr 2002 | Interpolation of data on
!         |             | icosahedral grid included.
!   3.1   | 16 Apr 2002 | Field status.
!   3.2   | 26 Aug 2004 | Corrected initialization of NC.
!----------------------------------------------------------
! Modules used:
!
Use MSIS, only: &
! Imported Routines:
    MSIS_Num_Levels
!    MSIS_Pressure_Levels,  &
!    MSIS_Geop,             &
!    MSIS_Refractivity
Use Interpolation, only: &
! Imported Routines:
    Init_Spline,  &
    Spline
!
Use Earth, only: &
! Imported Rourines:
    GCLat_from_GDLat, &
! Imported Type Definitions:
    Geodetic
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real, Intent(In) :: clat  ! MPAS cell latitude
Real, Intent(In) :: clon  ! MPAS cell longitude
Real, Intent(In) ::   Hsur ! surface geopotential[m]
Real, Intent(In) ::   H(:)    ! geopotential [m]
Real, Intent(In) ::   T(:)    ! temperature [k]
Real, Intent(In) ::   QR(:)   ! specific humidity [g/g]
Real, Intent(In) ::   P(:)   ! pressure [mb]
!
Real(Double), Intent(In) :: &
   ZP            ! Altitude [km]
!
! Output arguments:
!
Real(Double), Intent(Out) :: &
   Zmin          ! Minimum model Z for this lon/lat
!
Real(Double), Intent(Out) :: &
   Zmax          ! Maximum model Z for this lon/lat
!
Real(Double), Intent(Out) :: &
   LnNP          ! Interpolated Ln(N)
!
Real(Double), Intent(Out) :: &
   LnNZ          ! Interpolated dLn(N)/dZ
!
Real(Double), Intent(Out) :: &
   LnNZ2         ! Interpolated d2Ln(N)/dZ2
!
Real(Double), Optional, Intent(Out) :: &
   LnNIP(:)      ! Interpolated Ln(NI)[channel]
!----------------------------------------------------------
! Local Scalars:
!
Integer        :: NLev    ! Number of levels
Type(Geodetic) :: G       ! Geodetic coordinates of point
Real(Double)   :: PGCL    ! Geocentric latitude of point
Real(Double)   :: doub_clat ! Double precision Geodetic latitude of point
Integer        :: NC      ! Number of channels
Integer        :: IC      ! Channel number
Real(Double)   :: LnNIZ   ! Interpolated 1st derivative of LnNI
Real(Double)   :: LnNIZ2  ! Interpolated 2nd derivative of LnNI
Real(Double)   ::  Zsur   ! Surface model level for this lat/lon
Real(Double)   :: PR      ! Double preceision pressure level

!
! Local Arrays:
!
Real(Double), allocatable :: &
   ZLL(:),       & ! Cross section of Z
   FLL(:),       & ! Cross section of LnN or LnNI
   DLL(:)          ! Cross section of D2N or D2NI

Real(Double), allocatable :: &
     Z(:), &          ! Model altitudes (km)
     LnN(:), &        ! Log of refractive index
     D2N(:), &        ! Derivative of refractive index
     LnNI (:,:), &    ! Log of imaginary refractive index
     D2NI (:,:)    ! Derivative of imaginary refractive index
     
!----------------------------------------------------------
! 1. CALCULATION OF GRID SIZE
!----------------------------------------------------------

NLev = Size(T)

! This needs to be done in calling routine to allocate Z, lnN, D2nN 
PR   = Real(P(Nlev),Double)
NM = MSIS_Num_Levels(PR, P_up)

!print*,'nlev,nm=',nlev,nm

! Need to call function to determine nm, because it isn't yet set.  
! The nm is needed to allocate vertical arrays.

! The surface geopotential 

! Not sure if this will work yet
If (Allocated(Freq)) then
   NC = Size(Freq)
End If

!----------------------------------------------------------
! 2. VERTICAL INTERPOLATION
!----------------------------------------------------------


!--- 2.1. Vertical profile allocation
! Not sure what this is for.  Where are these declared?

!If (PStat(I1, I2, ID) == fst_Null) then

! The following allocation does not make sense to me
!
!   Allocate(  Z(I1, I2, ID)%P(1:NLev+NM))
!   Allocate(LnN(I1, I2, ID)%P(1:NLev+NM))
!   Allocate(D2N(I1, I2, ID)%P(1:NLev+NM))

! In the NCEP version of the code add a print*,size(z) or print*,shape(z)
! to see what it should be. 

! The allocation of Z does not work
   Allocate(Z(1:NLev+NM))
   Allocate(LnN(1:NLev+NM))
   Allocate(D2N(1:NLev+NM))

!   print*,'in vertical interpolation, size(z)=',size(z)

   If (Allocated(Freq)) then
!
!      Allocate(LnNI(I1, I2, ID)%M(1:NLev+NM, NC))
!      Allocate(D2NI(I1, I2, ID)%M(1:NLev+NM, NC))
      Allocate(LnNI(1:NLev+NM, NC))
      Allocate(D2NI(1:NLev+NM, NC))
   End If
!
!   PStat(I1, I2, ID) = fst_Allocated
!
!End If

!--- 2.2. Vertical profile intialization

!If (PStat(I1, I2, ID) == fst_Allocated) then
G    = Geodetic(0.0_Double, clat, clon)
doub_clat=clat
PGCL = GCLat_from_GDLat(doub_clat)

   If (Allocated(Freq)) then

      print*,'allocated freq calling Make_Refractivity_MPAS'
      print*,'nc=',nc
      Call Make_Refractivity_MPAS &
        (Hsur,        & ! <-- Surface geopotential [gpm]
         H,     & ! <-- Geopotential [gpm]
         T,     & ! <-- Temperature [K]
         QR,     & ! <-- Relative humidity [%]
         P,      & ! <-- Pressure [Pa]
         G,                         & ! <-- Geodetic coordinates
         PGCL,                      & ! <-- Geocentric latitude [rad]
         LnN,   & ! --> Ln of refractive index
         Zsur,        & ! --> Surface altitude [km]
         Z,   & ! --> Altitudes of model levels [km]
         LnNI)   ! ~~> Ln of dispersive imaginary part of refractive index

      Do IC=1,NC
         Call  Init_Spline   &
           (Z,            & ! <-- Argument grid
            LnNI(:,IC),   & ! <-- Gridded function
            D2NI(:,IC))     ! --> 2nd derivative of spline
      End Do

   Else

!      print*,'calling Make_Refractivity_MPAS freq not allocated'

      Call Make_Refractivity_MPAS &
        (Hsur,        & ! <-- Surface geopotential [gpm]
         H,     & ! <-- Geopotential [gpm]
         T,     & ! <-- Temperature [K]
         QR,     & ! <-- Relative humidity [%]
         P,      & ! Pressure [Pa]
         G,                         & ! <-- Geodetic coordinates
         PGCL,                      & ! <-- Geocentric latitude [rad]
         LnN,   & ! --> Ln of refractive index
         Zsur,        & ! --> Surface altitude [km]
         Z)     ! --> Altitudes of model levels [km]
!
   End If
   Call  Init_Spline   &
     (Z,        & ! <-- Argument grid
      LnN,      & ! <-- Gridded function
      D2N)        ! --> 2nd derivative of spline

!   PStat(I1, I2, ID) = fst_Initialized

!End If
!--- 2.3. Vertical spline interpolation of LnN

!!ZLL => Z  (I1, I2, ID)%P(:)
!!FLL => LnN(I1, I2, ID)%P(:)
!!DLL => D2N(I1, I2, ID)%P(:)

! The above are no longer pointers

ZLL = Z
FLL = LnN
DLL = D2N

!print*,'after fill of zll, fll and dll'

!print*,'size(z)=',size(z)
!print*,'size(zll)=',size(zll)
!print*,'size(fll)=',size(fll)
!print*,'size(dll)=',size(dll)
!print*,'zp=',zp


Call Spline  &
  (ZLL,       & ! <-- Argument grid
   FLL,       & ! <-- Gridded function
   DLL,       & ! <-- 2nd derivative of spline
   ZP,        & ! <-- Interpolation point
   LnNP,      & ! --> Interpolated function value
   LnNZ,      & ! --> Interpolated 1st derivative
   LnNZ2)       ! --> Interpolated 2nd derivative

!--- 2.4. Calculation of limits of Z
!
!!Zmin = Zsur(I1, I2, ID)
!!Zmax = Max(ZLL(1), ZLL(NLev))

Zmin = Zsur
Zmax = Max(ZLL(1), ZLL(NLev))
!print*,'zmin,zmax=',zmin,zmax

!------------------------------------------------------------------------------

!--- 2.5. Vertical spline interpolation of LnNI

If (Allocated(Freq) .and. Present(LnNIP)) then

   Do IC=1,NC

!      FLL => LnNI(I1, I2, ID)%M(:, IC)
!      DLL => D2NI(I1, I2, ID)%M(:, IC)

      FLL = LnNI(:, IC)
      DLL = D2NI(:, IC)

      Call Spline  &
        (ZLL,         & ! <-- Argument grid
         FLL,         & ! <-- Gridded function
         DLL,         & ! <-- 2nd derivative of spline
         ZP,          & ! <-- Interpolation point
         LnNIP(IC),   & ! --> Interpolated function value
         LnNIZ,       & ! --> Interpolated 1st derivative
         LnNIZ2)        ! --> Interpolated 2nd derivative

   End Do

! JDH Not sure if I will need this
!deallocate(lnNI)
!deallocate(D2NI)

End If

deallocate(Z)
deallocate(LnN)
deallocate(D2N)

deallocate(ZLL)
deallocate(FLL)
deallocate(DLL)

End Subroutine Vertical_Interpolation_LnN_MPAS

!==========================================================
Subroutine Make_Refractivity_MPAS &
  (Hsur,   & ! <-- Surface geopotential [gpm]
   H,      & ! <-- Geopotential [gpm]
   T,      & ! <-- Temperature [K]
   QR,     & ! <-- Specific humidity [%]
   P,      & ! <-- Pressure [Pa]
   G,      & ! <-- Geodetic coordinates
   GCLat,  & ! <-- Geocentric latitude [rad]
   LnN,    & ! --> Ln of refractive index
   Zsur,   & ! --> Surface altitude [km]
   Z,      & ! --> Altitudes of model levels [km]
   LnNI)     ! ~~> Ln of dispersive imaginary part of refractive index
!
! Calculation of refractivity profile and geometrical
! altitudes of isobaric levels from profiles of
! geopotential, temperature, and humidity.
!----------------------------------------------------------
! Method:
!----------------------------------------------------------
! (C) Copyright 2007, M. E. Gorbunov.
!  Adapted for MPAS 2021, by J. D. Hegarty
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 28 Jun 2007 | Original version.
!   1.1   | 10 Oct 2007 | Bug corrected.
!   2.0   | 18 Oct 2007 | MSIS level correction for matching
!         |             | dry temperature to MSIS temperature.
!   2.1   | 15 Sep 2008 | 
!----------------------------------------------------------
! Modules used:
!
Use Earth, only: &
! Imported Type Definitions:
    Geodetic,  &
! Imported Parameters:
    g_ave,     &
! Imported Routines:
    Alt_from_Geop,    &
    Gravity
!
Use Defaults, only: &
! Imported Parameters:
    Rd
!
Use Occ_Meteoprofiles, only: &
! Imported Parameters:
    eps, aq,  bq,          &
! Imported Routines:
    Q_from_Qrel,           &
    N_from_TPQ,            &
    T_from_NPQ,            &
    NQ_to_TP
!
Use MSIS, only: &
! Imported Routines:
!    MSIS_Num_Levels,       &
    MSIS_Pressure_Levels,  &
    MSIS_Geop,             &
    MSIS_Refractivity
!
Use Advanced_MPM93Model, only: &
! Imported Routines:
    AdvMPM93Model
!
!DEBUG
Use IO, only: &
! Imported Routines:
    PutXY,      &
    DeleteFile
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Single), Intent(In)   :: &
   Hsur           ! Surface geopotential [gpm]
!
Real(Single), Intent(In)   :: &
   H(:)           ! Geopotential [gpm]
!
Real(Single), Intent(In)   :: &
   T(:)           ! Temperature profile [K]
!
Real(Single), Intent(In)   :: &
   QR(:)          ! Specific humidity [%]
Real(Single), Intent(In)   :: &
   P(:)          ! Pressure [mb]
!
Type(Geodetic), Intent(In) :: &
   G              ! Geodetic coordinates
!
Real(Double), Intent(In)   :: &
   GCLat          ! Geocentric latitude [rad]
!
! Output arguments:
!
Real(Double), Intent(Out)  :: &
   LnN(:)         ! Ln of refractive index
!
Real(Double), Intent(Out)  :: &
   Zsur           ! Surface altitude [km]
!
Real(Double), Intent(Out)  :: &
   Z(:)           ! Altitudes of model levels [km]
!
! Output optional arguments:
!
Real(Double), Intent(Out), Optional :: &
   LnNI(1:,:)     ! Ln of dispersive imaginary part of refractive index
!----------------------------------------------------------
! Local Parameters:
!
Logical, Parameter :: &
   Dbg_mode = .False.      ! Debugging mode
!
! Local Scalars:
!
Integer         :: i       ! Level index
Real(Double)    :: TR      ! Temperature (Double)
Real(Double)    :: QRR     ! Relative humidity (Double)
Real(Double)    :: PR      ! Pressure Double, mb
Type(Geodetic)  :: GM      ! Point in vertical profile
Real(Double)    :: Zmax    ! Upper integration height [km]
Real(Double)    :: DZI     ! Integration step [km]
Real(Double)    :: Tinit   ! Temperature at upper boundary [K]
Real(Double)    :: Qmin    ! Minimum specific humidity [kg/kg]
Real(Double)    :: DZC     ! MSIS Z-level correction
Integer         :: IC      ! Channel number
Real(Double)    :: E       ! Water vapur pressure
Complex(Double) :: ZN      ! Dispersive part of refractivity [N-units]
!
!
! Local Arrays:
!
! --- Work arrays
!
!Real(Double) :: &
!   Q(1:NLev),         & ! Specific humidity [kg/kg]
!   GA(0:NLev)           ! Gravity acceleration [m/s**2]
!
Real(Double), allocatable :: &
! I think all the MSIS pressures are in mb
   PM(:),          & ! MSIS level pressure [Pa]
   HM(:),          & ! MSIS level geopotential [gpkm]
   ZM(:),          & ! MSIS level altitude [km]
   TM(:)             ! MSIS level temperature [K]
Real(Double), Allocatable :: &
   NNC(:),            & ! Local MPAS profile of N
   QNC(:),            & ! Local MPAS profile of Q
   TNC(:),            & ! Local MPAS profile of T
   DT (:)               ! Finite difference of T for debug.
!----------------------------------------------------------
! Global variables used:
!
!   NM           ! Number of MSIS half levels.
!----------------------------------------------------------

!----------------------------------------------------------
! 1. CALCULATION OF NCEP HALF AND FULL LEVELS
!----------------------------------------------------------

!print*,'In Make_Refractivity_MPAS'

Zsur = Alt_from_Geop(1e-3_Double*Hsur, GCLat)

!print*,'hsur,zsur=',hsur,zsur

Nlev=size(T)
!print*,'nlev=',nlev
Do i = 1,NLev
   Z(i)   = Alt_from_Geop(1d-3*H(i), GCLat)

! Already have Q in specific humidity so don't need

!   If (QR(i) >= 0.0 .and. QR(i) <= 100.0) then
!      TR   = Real(T(i), Double)
!      QRR  = Real(QR(i),Double)
!      Q(i) = Q_from_Qrel &
!                 (TR,    & ! <-- Temperature [K]
!                  P(i),  & ! <-- Pressure [mb]
!                  QRR)     ! <-- Relative humidity [%]
!   Else
!      Q(i) = 0.0
!   End If

   TR   = Real(T(i), Double)
   QRR  = Real(QR(i),Double)
   PR   = Real(P(i),Double)
! Convert pressure from Pascals to millibars which is expected in 
! function below
!   PR=PR/100.0
   LnN(i) = Log(N_from_TPQ(TR, PR, QRR ))

!   print*,'i,Z(i),PR,TR,QR,LnN=',i,Z(i),PR,TR,QRR,LnN(i)

End Do

!----------------------------------------------------------
! 2. CALCULATION OF MSIS LEVELS
!----------------------------------------------------------

! This needs to be done in calling routine to allocate Z, lnN, D2nN 
!PR   = Real(P(Nlev),Double)
!NM = MSIS_Num_Levels(PR, P_up)

!print*,'NM=',NM

allocate(pm(nm))
allocate(hm(nm))
allocate(zm(nm))
allocate(tm(nm))

!!--- 2.1. Calculation of MSIS pressure levels
!
Call MSIS_Pressure_Levels &
     (PR,     & ! <-- Basic pressure level [mb]
     PM(NM:1:-1))   ! --> MSIS pressure levels [mb]

!print*,'after MSIS_Pressure_levels'
!print*,'pm=',pm

!!--- 2.2. Calculation of MSIS level altitudes
!
Call MSIS_Geop &
  (G,           & ! <-- Geodetic coordinates
   PM(1:NM),    & ! <-- Pressure levels [mb]
   HM(1:NM),    & ! --> Geopotential heigths [gpkm]
   ZM(1:NM))      ! --> Altitudes [km]

!print*,'nlev,nm=',nlev,nm
!print*,'size(z)=',size(z)

Z(NLev+1:NLev+NM) = ZM(1:NM)

!print*,'hm=',hm
!print*,'zm=',zm

!--- 2.3. Calculation of MSIS refractivities

GM = G
!
Do i=NLev+1,NLev+NM
   GM%H = Z(i)
   Call MSIS_Refractivity &
     (GM,     & ! <-- Geodetic coordinates
      LnN(i))   ! --> Refractivity
   LnN(i) = Log(LnN(i))
!   print*,'i,LnN=',i,LnN(i)
End Do

!----------------------------------------------------------
! 3. CORRECTION OF MSIS LEVELS
!----------------------------------------------------------

!--- 3.1. Estimate of MSIS temperature at Z(NLev)

Allocate(NNC(1:NLev+NM))
Allocate(QNC(1:NLev+NM))
Allocate(TNC(1:NLev+NM))
Allocate(DT (1:NLev+NM))

TNC(1:NLev) = T(1:NLev)

Zmax  = 120.0 ! MaxVal(ZH)
DZI   = 3.0
Tinit = 0.0
Qmin  = 1d-7

NNC(:)      = Exp(LnN(:))
QNC(:)      = Qmin
QNC(1:NLev) = QR(1:NLev)

Call NQ_to_TP &
  (G%Phi, & ! <-- Geodetic latitude [deg]
   Z,     & ! <-- Altitude above reference ellipsoid [km]
   NNC,   & ! <-- Profile of refractivity
   QNC,   & ! <-- Profile of specific humidity
   Zmax,  & ! <-- Upper integration height [km]
   DZI,   & ! <-- Integration step [km]
   Tinit, & ! <-- Temperature at upper boundary [K]
   Qmin,  & ! <-- Minimum specific humidity [kg/kg]
   TNC)     ! --> Profile of temperature

!--- 3.2. Z-levels correction

DZC = 1d-3*((T(NLev) - TNC(NLev))*Rd*Exp(LnN(NLev))/  &
             Gravity(Z(NLev),GCLat))*                 &
      (LnN(NLev+1)      - LnN(NLev)) /                &
      (Exp(LnN(NLev+1)) - Exp(LnN(NLev)))

Z(NLev+1:NLev+NM) = ZM(1:NM) + DZC

Deallocate(NNC)
Deallocate(QNC)
Deallocate(TNC)
Deallocate(DT )

!----------------------------------------------------------
! 4. CALCULATION OF ABSORPTION PROFILE
!----------------------------------------------------------

If (Present(LnNI) .and. Allocated(Freq)) then

   print*,'Make_Refractivity --- allocated freq'

   !--- 4.1. Absorption computation at MSIS levels

   Do i=1,NM
      TM(i) = T_from_NPQ &
        (Exp(LnN(NLev+i)),  & ! <-- Refractivity [dimensionless]
         PM(i),             & ! <-- Pressure [mb]
         0.0_Double)          ! <-- Specific humidity [kg/kg]
      Do IC=1,Size(Freq)
         Call AdvMPM93Model &
           (Freq(IC),     & ! <-- Frequency [Hz]
            PM(i),        & ! <-- Pressure [mbar]
            TM(i),        & ! <-- Temperature [K]
            0.0_Double,   & ! <-- Specific humidity [kg/kg]
            0.0_Double,   & ! <-- Liquid water content [g/m^3]
            0.0_Double,   & ! <-- Rain rate [mm/h]
            0.0_Double,   & ! <-- Ice water content [g/m^3]
            ZN)             ! ~~> Complex refracitivity [N-units]
         LnNI(NLev+i,IC) = Log(Real(1e-6_Double*ZN,Double))
      End Do
   End Do

   !--- 4.2. Absorption computation at NCEP levels

   Do i=1,NLev
      Do IC=1,Size(Freq)
         QRR  = Real(QR(i),Double)
         PR  = Real(P(i),Double)
!         E = P(i)*Q(i)/(aq + bq*Q(i))
         E = PR*QRR/(aq + bq*QRR)
         Call AdvMPM93Model &
           (Freq(IC),            & ! <-- Frequency [Hz]
            PR,                & ! <-- Pressure [mbar]
            Real(T(i), Double),  & ! <-- Temperature [K]
            E,                   & ! <-- Water vapor pressure [mbar]
            0.0_Double,          & ! <-- Liquid water content [g/m^3]
            0.0_Double,          & ! <-- Rain rate [mm/h]
            0.0_Double,          & ! <-- Ice water content [g/m^3]
            ZN)                    ! ~~> Complex refracitivity [n-1]
         LnNI(i,IC) = Log(Real(1e-6_Double*ZN,Double))
      End Do
   End Do

! JDH test put in stop just to see if code makes it here
!   print*,'Make_Refractivity --- stopping in allocated freq if block'
!   stop

End If

deallocate(pm)
deallocate(hm)
deallocate(zm)
deallocate(tm)

End Subroutine Make_Refractivity_MPAS

!-----------------------------------------------------------------

Subroutine MPAS_Constituents &
  (PLon,     & ! <-- Longiude of point [deg]
   PLat,     & ! <-- Latitude of point [deg]
   ZP,       & ! <-- Altitude of point [km]
   ZPmin,    & ! --> Minimum model Z for this lon/lat
   ZPmax,    & ! --> Maximum model Z for this lon/lat
   QP,       & ! --> Interpolated specific humidity [kg/kg]
   TP,       & ! --> Interpolated temperature [K]
   RNR,      & ! --> Interpolated Re(N)
   RNI)        ! ~~> Interpolated Im(N)
!
! Computation of humidity, temperature,
! and complex refractivity for phantom.
!----------------------------------------------------------
! Method:
!   Linear interpolation of vertical profiles.
!----------------------------------------------------------
! (C) Copyright 2002-2007, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 04 Aug 2002 | Original version
!         |             | (from ECHAM_Constituents).
!   1.1   | 02 Jul 2007 | Conversion relative -> specific
!         |             | humidity.
!----------------------------------------------------------
! Modules used:
!
Use NetCDF, only: &
! Imported Parameters:
    NF90_NoWrite,             &
    NF90_NoErr,               &
    NF90_BYTE,                &
    NF90_CHAR,                &
    NF90_SHORT,               &
    NF90_INT,                 &
    NF90_FLOAT,               &
    NF90_DOUBLE,              &
    NF90_GLOBAL,              &
! Do I need to specifically import all the routines that I use?
! Imported Routines:
    NF90_Open,                &
    NF90_Inquire,             &
    NF90_Inquire_Dimension,   &
    NF90_INQ_DIMID,           &
    NF90_Inquire_Variable,    &
    NF90_INQ_VARID,           &
    NF90_Inq_Attname,         &
    NF90_Inquire_Attribute,   &
    NF90_Get_Att,             &
    NF90_Get_Var,             &
    NF90_Close

Use MSIS, only: &
! Imported Routines:
    MSIS_Num_Levels
!    MSIS_Pressure_Levels,  &
!    MSIS_Geop,             &
!    MSIS_Refractivity
Use Interpolation, only: &
! Imported Routines:
    Init_Spline,  &
    Linear,       &
    Spline
!
Use Earth, only: &
! Imported Routines
    GCLat_from_GDLat, &
! Imported Type Definitions:
    Geodetic
Use Occ_Meteoprofiles, only: &
! Imported Routines:
    Q_from_Qrel
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Public Parameters:
!
Integer, Parameter :: &
   err_NoFile   = 0702001,  & ! File not found
   err_MissVar  = 0702002,  & ! Missing variable/attributes
   err_ReadErr  = 0702003,  & ! Data read error
   err_FmtErr   = 0702004,  & ! Data format error
   err_NoLCF    = 0702005,  & ! No LCF
   err_NoFreq   = 0702006,  & ! No frequencies
   err_WrtErr   = 0702007,  & ! Write error
   err_DataType = 0702008,  & ! Bad data type
   err_NoData   = 0702009,  & ! No data
   err_Dimen    = 0702010,  & ! Inconsistent dimensions
   err_Channel  = 0702011,  & ! Inconsistent numbers of channels
   err_Coord    = 0702012,  & ! Bad coordinates
   err_FileName = 0702013     ! Bad file name
!----------------------------------------------------------
! Input arguments:
!
!--------------------------------------------------------

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
! Local Scalars:
!
real, parameter :: pi_const = 2.0 * asin(1.0)
!----------------------------------------------------------
! Local Variables:
!
integer :: ncid, iostat,ndims,unlimid,nd
integer :: varid, nattrs,xtype
integer :: idxm, last_idxm,j
integer :: maxedges, ncells, nvertices
Character(Len=255) :: Name         ! Parameter name
real :: platr,plonr,clatd,clond
integer :: nearestcell
integer :: sourcecells(3)
integer :: sourcevertices(3)
Integer, dimension(13)  :: dimids, lend, lendh, lends
!integer :: sourceedges(maxedges,npts)
real :: cellweights(3)
real :: vertexweights(3)
real :: pointinterp(3)
!      real :: edgeweights (maxedges,npts)
real :: vertcoords(3,3)
real :: vec(3)
real :: mpaswach(3)
real :: t2m(3), q2(3)

!----------------------------------------------------------
! Local Scalars:
!
Logical        :: SplInt  ! Spline interpolation indicator
Integer        :: NLev    ! Number of levels
Integer        :: i
Integer        :: Stat    ! Error status
Type(Geodetic) :: G       ! Geodetic coordinates of point
Real(Double)   :: pmb  ! Temporary pressure in mb for Q from relhum calculation
Real(Double)   :: q ! Double specific humidy for Q from relhum calculation
Real(Double)   :: QRR ! double precision specific humidity
Real(Double)   :: TR ! Double precision temperature
Real(Double)   :: doub_clat ! Double precision Geodetic latitude of point
Real(Double)   :: PGCL    ! Geocentric latitude of point
Real(Double)      :: QRP  ! Interpolated relative humidity [%].
Real(Double)   ::  Zsur   ! Surface model level for this lat/lon
Real(Double)   ::  PR   ! Double precision Pressure at top of model atmosphere[mbar]
! In the NCEP_Consituents PP there is no space for PP ro be passed out.
! It is calculated in 1 call to interpolate_constituents but then never used.
!Real(Double)      :: PP   ! Interpolated pressure [mbar].
!
! I think the weight values can just be dimentioned to 
! 3 since there are 3 MPAS points.
Real(Double) :: &
   ZNmin(3),          & ! Minimum model Z on subgrid
   ZNmax(3),          & ! Maximum model Z on subgrid
   FV(3),             & ! Vertically interpolated F on subgrid
   FVT(3),             & ! Vertically interpolated T on subgrid
   FVQ(3),             & ! Vertically interpolated Q on subgrid
   FZ(3),             & ! Vertically interpolated DF/DZ
   LnPV(3)             ! Vertically interpolated Ln(P[mbar])

Real(Double), Allocatable :: &
   LnPLL(:)             ! Vertical profile of Ln(P)
! The following three variables were originally pointers.
Real(Double), Allocatable :: &
   ZLL(:),            & ! Cross sections of Z
   FLL(:),            & ! Cross sections of F
   DLL(:),            & ! Cross sections of D2F
   TLL(:),            & ! Cross sections of T
   QLL(:)               ! Cross sections of Q
!----------------------------------------------------------

real, allocatable :: ti(:,:) ! working array for 3 points (level,3) temperature (k)
real, allocatable :: qi(:,:)  ! working arrays for 3 points (level,3) specific humidity (g/g)
real, allocatable :: hi(:,:) ! working array for 3 points (level,3) geopotential height (m) 
real, allocatable :: pi(:,:) ! working array for 3 points (level,3) pressure (Pa)
! Not sure I need the below arrays as the surface geopotential will be the first
! level of the geopotential profile 
real ::  hs(3)         ! working array for 3 points surface geopotential height (m)

! interpolated profiles from MPAS
real, allocatable :: t(:) ! temperature (k)
real, allocatable :: qr(:)  ! specific humidity (g/g)
real, allocatable :: h(:) ! geopotential height (m) 
real, allocatable :: p(:) ! pressure (mb) - converted from Pa before passed into subroutines
 
real, allocatable :: t_int(:) ! temperature (k)
real, allocatable :: q_int(:)  ! specific humidity (g/g)
real, allocatable :: p_int(:) ! pressure (Pa) 
real ::  hsur,tv,hscl     ! surface geopotential, virtual potential temp, scale height 

Real(Double), allocatable :: &
     Z(:), &          ! Model altitudes (km)
     LnN(:), &        ! Log of refractive index
     D2N(:)        ! Derivative of refractive index
!     LnNI (:,:), &    ! Log of imaginary refractive index
!     D2NI (:,:)    ! Derivative of imaginary refractive index

! D2F was not passed into interpolate_consituents in NCEP_Constituents, ignore for now.
!SplInt = Present(D2F)

!----------------------------------------------------------
! 2. CALCULATION OF INTERPOLATION WEIGHTS
!----------------------------------------------------------

ncid = fields%ncid_sav

plonr=plon*pi_const/180.0
platr=plat*pi_const/180.0

ncells=size(fields%verticesoncell,2)
maxedges=size(fields%verticesoncell,1)
nvertices=size(fields%latvertex)

last_idxm=1

idxm = nearest_cell(platr,plonr, last_idxm, nCells, & 
     maxEdges,fields%nEdgesOnCell,fields%cellsOnCell,fields%latCell,fields%lonCell)

last_idxm = 1
idxm = nearest_vertex(platr, plonr, last_idxm, nCells, & 
     nVertices,maxEdges, 3, fields%nEdgesOnCell, fields%verticesOnCell, & 
     fields%cellsOnVertex, fields%latCell, fields%lonCell, & 
     fields%latVertex, fields%lonVertex )

if (idxm <= 0) then
   sourceCells(:) = 1
   cellWeights(:) = 0.0
else
!! convert_lx returns a 3 element vector, pointInterp is an array of 3
!! SourceCells an array of 3 x ncells.  Not sure if they can be assigned
!! in this way.  If not loop it.
   call convert_lx(platr,plonr,6371229.0,vec)
   pointinterp(1:3)=vec(1:3)
   sourceCells(:) = fields%cellsOnVertex(:,idxm)
   do j=1,3
      call convert_lx (fields%latCell(fields%cellsOnVertex(j,idxm)), &
           fields%lonCell(fields%cellsOnVertex(j,idxm)),6371229.0, vec)
      vertCoords(1:3,j) = vec(1:3)
   end do

!! This call creates cellweights for all three vertices at once in a single call
!! to the function.  I may need to re-write if this syntax doesn't work in a single line
   call mpas_wachspress_coordinates(3,vertCoords,pointInterp,mpaswach)
   cellWeights(:)=mpaswach(:) 
   last_idxm = idxm
end if

!----------------------------------------------------------------------------------------------
! 3.0 Read MPAS Meteorological Profiles for 3 Horizontal Cells
!---------------------------------------------------------------------------------------------- 
!
! Read theta, specific humidity, pressure and geopotential profiles from MPAS
!
!------------------ theta ----------------------
name='theta'
IOStat=NF90_INQ_VARID(NCID,name,varid)
If (IOStat /= NF90_NoErr) then
   print*,'nf90_inq_varid IOStat=',IOStat
   print*,'MPAS_Constituents name=',name
   Stat = err_ReadErr
   Return
End If
!print*,'varid=',varid
! The first call to NF90_Inquire _Variable gets the
! number of dimensions, the second call gets the dimension
! IDs from there I think I can get the actual diminesion perhaps using
! Inquire dimension.  Try to follow scan_inpur.F in convert_mpas.F code.
   IOStat = NF90_Inquire_Variable(NCID,varid,name,xtype,ndims,dimids,nattrs )
   If (IOStat /= NF90_NoErr) then
      Stat = err_ReadErr
      Return
   End If

do nd=1,ndims
   IOStat=NF90_INQUIRE_DIMENSION(NCID,dimids(nd),name,lend(nd))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_inq_dimension IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
!   print*,'name,lend=',name,lend(nd)
end do

allocate(ti(lend(1),3))

do j=1,3
   IOStat=NF90_GET_VAR(NCID,varid,ti(:,j),start=(/1,sourcecells(j),1/), &
        count=(/lend(1),1,1/))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_get_var IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
end do

!----------------- pressure ---------------------------
name='pressure'
IOStat=NF90_INQ_VARID(NCID,name,varid)
If (IOStat /= NF90_NoErr) then
   print*,'nf90_inq_varid IOStat=',IOStat
   print*,'MPAS_Constituents name=',name

   Stat = err_ReadErr
   Return
End If
!print*,'varid=',varid
   IOStat = NF90_Inquire_Variable(NCID,varid,name,xtype,ndims,dimids,nattrs )
   If (IOStat /= NF90_NoErr) then
      Stat = err_ReadErr
      Return
   End If

!print*,'after NF90_Inquire_Variable, ndims,nattrss,name=',ndims,nattrs,name,dimids

do nd=1,ndims
   IOStat=NF90_INQUIRE_DIMENSION(NCID,dimids(nd),name,lend(nd))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_inq_dimension IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
!   print*,'name,lend=',name,lend(nd)
end do

allocate(pi(lend(1),3))

do j=1,3
   IOStat=NF90_GET_VAR(NCID,varid,pi(:,j),start=(/1,sourcecells(j),1/), &
        count=(/lend(1),1,1/))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_get_var IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
end do

!--------------------- specific humidity ----------------------
!name='qv'
! read in relative humidity instead
name='relhum'
IOStat=NF90_INQ_VARID(NCID,name,varid)
If (IOStat /= NF90_NoErr) then
   print*,'nf90_inq_varid IOStat=',IOStat
   print*,'MPAS_Constituents name=',name

   Stat = err_ReadErr
   Return
End If
!print*,'varid=',varid
   IOStat = NF90_Inquire_Variable(NCID,varid,name,xtype,ndims,dimids,nattrs )
   If (IOStat /= NF90_NoErr) then
      Stat = err_ReadErr
      Return
   End If

!print*,'after NF90_Inquire_Variable, ndims,nattrss,name=',ndims,nattrs,name,dimids

do nd=1,ndims
   IOStat=NF90_INQUIRE_DIMENSION(NCID,dimids(nd),name,lend(nd))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_inq_dimension IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
!   print*,'name,lend=',name,lend(nd)
end do

allocate(qi(lend(1),3))

do j=1,3
   IOStat=NF90_GET_VAR(NCID,varid,qi(:,j),start=(/1,sourcecells(j),1/), &
        count=(/lend(1),1,1/))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_get_var IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
end do

!-------------------- geopotential height ---------------------------------
name='zgrid'
IOStat=NF90_INQ_VARID(NCID,name,varid)
If (IOStat /= NF90_NoErr) then
   print*,'nf90_inq_varid IOStat=',IOStat
   print*,'MPAS_Constituents name=',name

   Stat = err_ReadErr
   Return
End If
!print*,'varid=',varid
   IOStat = NF90_Inquire_Variable(NCID,varid,name,xtype,ndims,dimids,nattrs )
   If (IOStat /= NF90_NoErr) then
      Stat = err_ReadErr
      Return
   End If

!print*,'after NF90_Inquire_Variable, ndims,nattrss,name=',ndims,nattrs,name,dimids

do nd=1,ndims
   IOStat=NF90_INQUIRE_DIMENSION(NCID,dimids(nd),name,lendh(nd))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_inq_dimension IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
!   print*,'name,lendh=',name,lendh(nd)
end do

! only get up to the number of half levels
 allocate(hi(lend(1),3))
!allocate(hi(lendh(1),3))

do j=1,3
   IOStat=NF90_GET_VAR(NCID,varid,hi(:,j),start=(/1,sourcecells(j),1/), &
        count=(/lend(1),1,1/))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_get_var IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
end do

!-------------------- t2m ---------------------------------
name='t2m'
IOStat=NF90_INQ_VARID(NCID,name,varid)
If (IOStat /= NF90_NoErr) then
   print*,'nf90_inq_varid IOStat=',IOStat
   print*,'MPAS_Constituents name=',name

   Stat = err_ReadErr
   Return
End If
!print*,'varid=',varid
   IOStat = NF90_Inquire_Variable(NCID,varid,name,xtype,ndims,dimids,nattrs )
   If (IOStat /= NF90_NoErr) then
      Stat = err_ReadErr
      Return
   End If

!print*,'after NF90_Inquire_Variable, ndims,nattrss,name=',ndims,nattrs,name,dimids

do nd=1,ndims
   IOStat=NF90_INQUIRE_DIMENSION(NCID,dimids(nd),name,lends(nd))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_inq_dimension IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
!   print*,'name,lends=',name,lends(nd)
end do

do j=1,3
   IOStat=NF90_GET_VAR(NCID,varid,t2m(j),start=(/sourcecells(j),1/))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_get_var IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
end do
!-------------------- q2 ---------------------------------
name='q2'
IOStat=NF90_INQ_VARID(NCID,name,varid)
If (IOStat /= NF90_NoErr) then
   print*,'nf90_inq_varid IOStat=',IOStat
   print*,'MPAS_Constituents name=',name

   Stat = err_ReadErr
   Return
End If
!print*,'varid=',varid
   IOStat = NF90_Inquire_Variable(NCID,varid,name,xtype,ndims,dimids,nattrs )
   If (IOStat /= NF90_NoErr) then
      Stat = err_ReadErr
      Return
   End If

!print*,'after NF90_Inquire_Variable, ndims,nattrss,name=',ndims,nattrs,name,dimids

do nd=1,ndims
   IOStat=NF90_INQUIRE_DIMENSION(NCID,dimids(nd),name,lends(nd))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_inq_dimension IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
!   print*,'name,lends=',name,lends(nd)
end do

do j=1,3
   IOStat=NF90_GET_VAR(NCID,varid,q2(j),start=(/sourcecells(j),1/))
   If (IOStat /= NF90_NoErr) then
      print*,'nf90_get_var IOStat=',IOStat
      Stat = err_ReadErr
      Return
   End If
end do

!-------------  End Read 3 cells of MPAS meteorological profiles  - ------------

! Following Gorbunov Code, lets do vertical regridding and calculate refractivity for each of 
! the cell profiles, then can apply MPAS cell weights later. 

allocate(p(lend(1)))
allocate(t(lend(1)))
allocate(qr(lend(1)))
allocate(h(lend(1)))

! Convert relative humidity to specific humidity

do j=1,3
   do i=1,lend(1)

      !   If (Qi(i) >= 0.0 .and. Q(i) <= 100.0) then
      If (Qi(i,j) >= 0.0 .and. Qi(i,j) <= 110.0) then
         tv = ti(i,j) * (pi(i,j)/1e05)**0.286
         pmb=pi(i,j)/100.0
         
         TR   = Real(tv, Double)
         QRR  = Real(Qi(i,j),Double)
         
         Q = Q_from_Qrel &
              (TR,    & ! <-- Temperature [K]
              PMB,  & ! <-- Pressure [mb]
              QRR)     ! <-- Relative humidity [%]
      Else
         Q = 0.0
      End If

! Replace relative humidity with specific humidity
! Debug print
!      if(i.eq.20) then
!         print*,'PMB,TR,QRR,Q=',PMB,TR,QRR,Q
!      end if
      qi(i,j)=q
      
   End Do
end do

t(1)=0.0
qr(1)=0.0

! Loop through 3 points and interpolate to the plon,plat location

do j=1,3

   if(j.eq.1) then
      NLev = Size(T)
      Allocate(FLL(NLev))
      Allocate(TLL(NLev))
      Allocate(QLL(NLev))
      Allocate(LnPLL(NLev))

   end if

   t(1)=t2m(j)
   qr(1)=q2(j)
   hsur=hi(1,j)
   h(1)=hi(1,j)

! 4.2 Vertical Regridding to geopotential Height Grid (MPAS zgrid)

   do i=2,lend(1)
      p(i)=fields%fzm(i)*pi(i,j) + fields%fzp(i)*pi(i-1,j)
      t(i)=fields%fzm(i)*ti(i,j) + fields%fzp(i)*ti(i-1,j)
      qr(i)=fields%fzm(i)*qi(i,j) + fields%fzp(i)*qi(i-1,j)
      h(i)=hi(i,j)
   end do

! Calculate surface pressure using the hypsometric equation
!     First get the virtual temperature
   tv = ti(1,j) * (pi(1,j)/1e05)**0.286
! Approximate virtual temperature factors are in Celcius 
   tv = (tv-273.15)*(1 + 0.61*qi(1,j)*1000.0) + 273.15
!     Scale height
   hscl=29.3 * tv
   p(1)=p(2)*exp((h(2) - h(1))/hscl)

! Convert theta to temperature 
   do i=1,lend(1)
      t(i)=t(i) * (p(i)/1e05)**0.286
   end do

! Convert pressure to mb from Pa 
   do i=1,lend(1)
      p(i)=p(i)/100.0
   end do

   clatd=fields%latCell(sourcecells(j)) * 180.0 / pi_const
   clond=fields%lonCell(sourcecells(j)) * 180.0 / pi_const

!----------------------------------------------------------
! 3. VERTICAL INTERPOLATION
!----------------------------------------------------------

G    = Geodetic(0.0_Double, clatd, clond)
doub_clat=clatd
PGCL = GCLat_from_GDLat(doub_clat)

! This needs to be done each time through loop because
! P(Nlev) may change 
PR   = Real(P(Nlev),Double)
NM = MSIS_Num_Levels(PR, P_up)

allocate(Z(Nlev+NM))
allocate(LnN(Nlev+NM))
allocate(D2N(Nlev+NM))

! Not clear what LnPLL is used for
!If (Present(PP)) then
!   LnPLL(:) = Log(P(:))
!End If

! JDH Debug 032521
!if(nm.lt.100) then 

   print*,'calling Make_Refractivity_MPAS'
!   print*,'plon,plat,clond,clatd=',plon,plat,clond,clatd
!   print*,'j,hsur,pgcl=',j,hsur,pgcl
!   print*,'t=',t
!   print*,'qr=',qr
!   print*,'p=',p
!   print*,'nlev,nm=',nlev,nm
!   print*,'size(z)=',size(z)

   if(size(z).lt.100) then
      print*,'In MPAS_Constituents stopping before Call_Make_Refractivity_MPAS'
      stop
   end if

!end if

      Call Make_Refractivity_MPAS    &
        (Hsur,       &
         H,    &
         T,    &
         QR,    &
         P,     &
         G,                       &
         PGCL,                    &
         LnN,  &
         Zsur,       &
         Z)

      Call Init_Spline (Z,LnN,D2N)
!      PStat(I1, I2, ID) = fst_Initialized
!   End If

!   ZLL    => Z(I1, I2, ID)%P(1:NLev)
!   FLL(1:NLev) = Real(F(I1, I2, ID, 1:NLev), Double)
   ZLL = Z
!   FLL(1:NLev) = Real(F(1:NLev), Double)
   TLL(1:NLev) = Real(T(1:NLev), Double)
   QLL(1:NLev) = Real(QR(1:NLev), Double)

!   If (SplInt) then

!      If (.not. Associated(D2F(I1, I2, ID)%P)) then
!         Allocate(D2F(I1, I2, ID)%P(1:NLev))
!         Call Init_Spline   &
!           (  Z(I1, I2, ID)%P(1:NLev),    &
!            FLL,                          &
!            D2F(I1, I2, ID)%P(1:NLev))
!      End If
!
!      DLL => D2F(I1, I2, ID)%P(1:NLev)
!      Call Spline(ZLL, FLL, DLL, ZP, FV(IGP), FZ(IGP))

!   Else

      Call Linear(ZLL, TLL, ZP, FVT(j), CExt =.True.)
      Call Linear(ZLL, QLL, ZP, FVQ(j), CExt =.True.)
      FZ(j) = 0

!   End If

! Not sure what lnPLL is used for, not passed out of NCEP_Constituents
!   If (Present(PP)) then
!      Call Linear(ZLL, LnPLL, ZP, LnPV(j), CExt =.True.)
!   End If

   ZNmin(j) = Zsur
   ZNmax(j) = Max(ZLL(1), ZLL(NLev))

   FLL = LnN
   DLL = D2N

!print*,'after fill of zll, fll and dll'

!print*,'size(z)=',size(z)
!print*,'size(zll)=',size(zll)
!print*,'size(fll)=',size(fll)
!print*,'size(dll)=',size(dll)
!print*,'zp=',zp

Call Spline  &
  (ZLL,       & ! <-- Argument grid
   FLL,       & ! <-- Gridded function
   DLL,       & ! <-- 2nd derivative of spline
   ZP,        & ! <-- Interpolation point
   FV(j),      & ! --> Interpolated function value
   FZ(j))        ! --> Interpolated 1st derivative

! NM may chage each time through loop so have to allocate and then
! deallocate each time through.

deallocate(Z)
deallocate(LnN)
deallocate(D2N)

End Do     ! j loop through 3 MPAS points

!----------------------------------------------------------
! 4. HORIONTAL INTERPOLATION
!----------------------------------------------------------

TP   = Sum(cellweights(:)*FVT(:))
QP   = Sum(cellweights(:)*FVQ(:))
RNR   = Sum(cellweights(:)*FV(:))
ZPmin = Sum(cellweights(:)*ZNmin(:))
ZPmax = Sum(cellweights(:)*ZNmax(:))

!TP   = Sum(WG(:)*FVT(:))
!QP   = Sum(WG(:)*FVQ(:))
!ZPmin = Sum(WG(:)*ZNmin(:))
!ZPmax = Sum(WG(:)*ZNmax(:))

! DF was not passed into Interpolate_Consituents from NCEP_Constituents so this is just a place holder
! in case we decide we need the vertical gradients here.

!If (Present(DF)) then
!!   DF   = Sum(WG(:)*FZ(:))
!   DF   = Sum(cellweights(:)*FZ(:))
!End If


! Not sure what PP was used for.  It is not passed out of NCEP_Constituents
!If (Present(PP)) then
!!   PP   = Exp(Sum(WG(:)*LnPV(:)))
!   PP   = Exp(Sum(cellweights(:)*LnPV(:)))
!End If

!----------------------------------------------------------
! 6. MEMORY DEALLOCATION
!----------------------------------------------------------

Deallocate(FLL)
Deallocate(DLL)
Deallocate(TLL)
Deallocate(QLL)
!Deallocate(LnPV)
Deallocate(LnPLL)

!deallocate(Z)
!deallocate(LnN)
!deallocate(D2N)

deallocate (ti, stat=stat) ! working array for 3 points (level,3) temperature (k)
deallocate (qi, stat=stat)  ! working arrays for 3 points (level,3) specific humidity (g/g)
deallocate (hi, stat=stat) ! working array for 3 points (level,3) geopotential height (m) 
deallocate (pi, stat=stat) ! working array for 3 points (level,3) pressure (Pa)

! interpolated profiles from MPAS
deallocate (t, stat=stat) ! temperature (k)
deallocate (qr, stat=stat)  ! specific humidity (g/g)
deallocate (h, stat=stat) ! geopotential height (m) 
deallocate (p, stat=stat) ! pressure (mb) - converted from Pa before passed into subroutines
 
deallocate (t_int, stat=stat) ! temperature (k)
deallocate (q_int, stat=stat) ! specific humidity (g/g)
deallocate (p_int, stat=stat) ! pressure (Pa) 

End Subroutine MPAS_Constituents

!----------------------------------------------------------

End Module MPAS

