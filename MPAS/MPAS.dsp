# Microsoft Developer Studio Project File - Name="MPAS" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=MPAS - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "MPAS.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "MPAS.mak" CFG="MPAS - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "MPAS - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "MPAS - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 1
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "MPAS - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /compile_only /include:"../Lib/Release" /include:"../GRIB/Release" /include:"../ICO/Release" /include:"../Models/Release" /include:"../occlib/Release" /include:"../GCM/Release" /nologo /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD BASE RSC /l 0x419 /d "NDEBUG"
# ADD RSC /l 0x419 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "MPAS - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /dbglibs /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /check:bounds /compile_only /dbglibs /debug:full /include:"../Lib/Debug" /include:"../GRIB/Debug" /include:"../ICO/Debug" /include:"../Models/Debug" /include:"../occlib/Debug" /include:"../GCM/Debug" /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD BASE RSC /l 0x419 /d "_DEBUG"
# ADD RSC /l 0x419 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "MPAS - Win32 Release"
# Name "MPAS - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\NCEP_1DVar.f90

!IF  "$(CFG)" == "NCEP - Win32 Release"

DEP_F90_NCEP_=\
	"..\GRIB\Release\GRIB.mod"\
	"..\Lib\Release\Coordinates.mod"\
	"..\Lib\Release\Defaults.mod"\
	"..\Lib\Release\Earth.mod"\
	"..\Lib\Release\Externals.mod"\
	".\Release\NCEP_fields.mod"\
	".\Release\NCEP_fields_adj.mod"\
	

!ELSEIF  "$(CFG)" == "NCEP - Win32 Debug"

DEP_F90_NCEP_=\
	"..\GRIB\Debug\GRIB.mod"\
	"..\Lib\Debug\Coordinates.mod"\
	"..\Lib\Debug\Defaults.mod"\
	"..\Lib\Debug\Earth.mod"\
	"..\Lib\Debug\Externals.mod"\
	".\Debug\NCEP_fields.mod"\
	".\Debug\NCEP_fields_adj.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\MPAS_fields.f90

!IF  "$(CFG)" == "MPAS - Win32 Release"

DEP_F90_NCEP_F=\
	"..\GCM\Release\GCM_grid.mod"\
	"..\GRIB\Release\GRIB.mod"\
	"..\ICO\Release\ICO_grid.mod"\
	"..\ICO\Release\ICO_interpolation.mod"\
	"..\Lib\Release\Coordinates.mod"\
	"..\Lib\Release\Defaults.mod"\
	"..\Lib\Release\Earth.mod"\
	"..\Lib\Release\Errors.mod"\
	"..\Lib\Release\Interpolation.mod"\
	"..\Lib\Release\IO.MOD"\
	"..\Models\Release\Advanced_MPM93Model.mod"\
	"..\Models\Release\MSIS.mod"\
	"..\occlib\Release\Occ_Meteoprofiles.mod"\
	".\Release\NCEP_IO.mod"\
	

!ELSEIF  "$(CFG)" == "NCEP - Win32 Debug"

DEP_F90_NCEP_F=\
	"..\GCM\Debug\GCM_grid.mod"\
	"..\GRIB\Debug\GRIB.mod"\
	"..\ICO\Debug\ICO_grid.mod"\
	"..\ICO\Debug\ICO_interpolation.mod"\
	"..\Lib\Debug\Coordinates.mod"\
	"..\Lib\Debug\Defaults.mod"\
	"..\Lib\Debug\Earth.mod"\
	"..\Lib\Debug\Errors.mod"\
	"..\Lib\Debug\Interpolation.mod"\
	"..\Lib\Debug\IO.MOD"\
	"..\Models\Debug\Advanced_MPM93Model.mod"\
	"..\Models\Debug\MSIS.mod"\
	"..\occlib\Debug\Occ_Meteoprofiles.mod"\
	".\Debug\NCEP_IO.MOD"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\NCEP_fields_adj.f90

!IF  "$(CFG)" == "NCEP - Win32 Release"

DEP_F90_NCEP_FI=\
	"..\GCM\Release\GCM_grid.mod"\
	"..\GRIB\Release\GRIB.mod"\
	"..\ICO\Release\ICO_grid.mod"\
	"..\ICO\Release\ICO_interpolation.mod"\
	"..\Lib\Release\Defaults.mod"\
	"..\Lib\Release\Earth.mod"\
	"..\Lib\Release\Interpolation_adj.mod"\
	"..\Models\Release\Advanced_MPM93Model.mod"\
	"..\Models\Release\MSIS.mod"\
	"..\occlib\Release\Occ_Meteoprofiles.mod"\
	"..\occlib\Release\Occ_Meteoprofiles_adj.mod"\
	".\Release\NCEP_fields.mod"\
	

!ELSEIF  "$(CFG)" == "NCEP - Win32 Debug"

DEP_F90_NCEP_FI=\
	"..\GCM\Debug\GCM_grid.mod"\
	"..\GRIB\Debug\GRIB.mod"\
	"..\ICO\Debug\ICO_grid.mod"\
	"..\ICO\Debug\ICO_interpolation.mod"\
	"..\Lib\Debug\Defaults.mod"\
	"..\Lib\Debug\Earth.mod"\
	"..\Lib\Debug\Interpolation_adj.mod"\
	"..\Models\Debug\Advanced_MPM93Model.mod"\
	"..\Models\Debug\MSIS.mod"\
	"..\occlib\Debug\Occ_Meteoprofiles.mod"\
	"..\occlib\Debug\Occ_Meteoprofiles_adj.mod"\
	".\Debug\NCEP_fields.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\NCEP_IO.f90

!IF  "$(CFG)" == "NCEP - Win32 Release"

DEP_F90_NCEP_I=\
	"..\GCM\Release\GCM_grid.mod"\
	"..\GRIB\Release\GRIB.mod"\
	"..\ICO\Release\ICO_grid.mod"\
	"..\Lib\Release\Defaults.mod"\
	"..\Lib\Release\Errors.mod"\
	"..\Lib\Release\IO.MOD"\
	

!ELSEIF  "$(CFG)" == "NCEP - Win32 Debug"

DEP_F90_NCEP_I=\
	"..\GCM\Debug\GCM_grid.mod"\
	"..\GRIB\Debug\GRIB.mod"\
	"..\ICO\Debug\ICO_grid.mod"\
	"..\Lib\Debug\Defaults.mod"\
	"..\Lib\Debug\Errors.mod"\
	"..\Lib\Debug\IO.MOD"\
	

!ENDIF 

# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# End Group
# End Target
# End Project
