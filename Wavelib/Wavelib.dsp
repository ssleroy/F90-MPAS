# Microsoft Developer Studio Project File - Name="Wavelib" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=Wavelib - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "Wavelib.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "Wavelib.mak" CFG="Wavelib - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "Wavelib - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "Wavelib - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 1
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "Wavelib - Win32 Release"

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
# ADD BASE F90 /compile_only /include:"Release/" /nologo /warn:nofileopt
# ADD F90 /compile_only /dbglibs /include:"Release/" /include:"../Lib/Release/" /include:"../ECHAM/Release" /include:"../NCEP/Release" /include:"../occlib/Release" /include:"../FFTW/Release" /include:"../Phantom/Release" /include:"../GRIB/Release" /include:"../MSIS/Release" /math_library:fast /nologo /traceback /warn:nofileopt /fast
# SUBTRACT F90 /assume:noaccuracy_sensitive /automatic
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /YX /FD /c
# ADD CPP /nologo /G6 /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /YX /FD /c
# ADD BASE RSC /l 0x419
# ADD RSC /l 0x419
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "Wavelib - Win32 Debug"

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
# ADD BASE F90 /compile_only /debug:full /include:"Debug/" /nologo /warn:nofileopt
# ADD F90 /check:bounds /check:format /check:power /check:overflow /check:underflow /compile_only /dbglibs /debug:full /include:"Debug/" /include:"../Lib/Debug" /include:"../ECHAM/Debug" /include:"../NCEP/Debug" /include:"../occlib/Debug" /include:"../FFTW/Debug" /include:"../Phantom/Debug" /include:"../GRIB/Debug" /include:"../MSIS/Debug" /nologo /traceback /warn:nofileopt /warn:unused
# ADD BASE CPP /nologo /W3 /GX /Z7 /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /Zi /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /YX /FD /c
# ADD BASE RSC /l 0x419
# ADD RSC /l 0x419
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "Wavelib - Win32 Release"
# Name "Wavelib - Win32 Debug"
# Begin Source File

SOURCE=.\Asymptotic_Propagator.f90

!IF  "$(CFG)" == "Wavelib - Win32 Release"

DEP_F90_ASYMP=\
	"..\FFTW\Release\FFTW.mod"\
	"..\Lib\Release\Coordinates.mod"\
	"..\Lib\Release\Defaults.mod"\
	"..\Lib\Release\Earth.mod"\
	"..\Lib\Release\Externals.mod"\
	"..\Lib\Release\Interpolation.mod"\
	"..\Lib\Release\IO.MOD"\
	"..\occlib\Release\Occ_Coordinates.mod"\
	"..\occlib\Release\Occ_Diffraction.mod"\
	".\Release\Ray_Problem.mod"\
	

!ELSEIF  "$(CFG)" == "Wavelib - Win32 Debug"

DEP_F90_ASYMP=\
	"..\FFTW\Debug\FFTW.mod"\
	"..\Lib\Debug\Coordinates.mod"\
	"..\Lib\Debug\Defaults.mod"\
	"..\Lib\Debug\Earth.mod"\
	"..\Lib\Debug\Externals.mod"\
	"..\Lib\Debug\Interpolation.mod"\
	"..\Lib\Debug\IO.MOD"\
	"..\occlib\Debug\Occ_Coordinates.mod"\
	"..\occlib\Debug\Occ_Diffraction.mod"\
	".\Debug\Ray_Problem.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\Atmosphere.f90

!IF  "$(CFG)" == "Wavelib - Win32 Release"

DEP_F90_ATMOS=\
	"..\ECHAM\Release\ECHAM_fields.mod"\
	"..\Lib\Release\Coordinates.mod"\
	"..\Lib\Release\Defaults.mod"\
	"..\Lib\Release\Errors.mod"\
	"..\NCEP\Release\NCEP_fields.mod"\
	"..\Phantom\Release\Phantom.mod"\
	"..\Phantom\Release\Turbulence.mod"\
	

!ELSEIF  "$(CFG)" == "Wavelib - Win32 Debug"

DEP_F90_ATMOS=\
	"..\ECHAM\Debug\ECHAM_fields.mod"\
	"..\Lib\Debug\Coordinates.mod"\
	"..\Lib\Debug\Defaults.mod"\
	"..\Lib\Debug\Errors.mod"\
	"..\NCEP\Debug\NCEP_fields.mod"\
	"..\Phantom\Debug\Phantom.mod"\
	"..\Phantom\Debug\Turbulence.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\Atmosphere_rays.f90

!IF  "$(CFG)" == "Wavelib - Win32 Release"

DEP_F90_ATMOSP=\
	"..\Lib\Release\Coordinates.mod"\
	"..\Lib\Release\Defaults.mod"\
	"..\Lib\Release\Earth.mod"\
	".\Release\Atmosphere.mod"\
	

!ELSEIF  "$(CFG)" == "Wavelib - Win32 Debug"

DEP_F90_ATMOSP=\
	"..\Lib\Debug\Coordinates.mod"\
	"..\Lib\Debug\Defaults.mod"\
	"..\Lib\Debug\Earth.mod"\
	".\Debug\Atmosphere.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\Atmosphere_rays_adj.f90

!IF  "$(CFG)" == "Wavelib - Win32 Release"

DEP_F90_ATMOSPH=\
	"..\Lib\Release\Coordinates.mod"\
	"..\Lib\Release\Defaults.mod"\
	"..\Lib\Release\Earth.mod"\
	"..\Lib\Release\Matrix.mod"\
	"..\occlib\Release\Occ_Refraction_adj.mod"\
	".\Release\Atmosphere.mod"\
	

!ELSEIF  "$(CFG)" == "Wavelib - Win32 Debug"

DEP_F90_ATMOSPH=\
	"..\Lib\Debug\Coordinates.mod"\
	"..\Lib\Debug\Defaults.mod"\
	"..\Lib\Debug\Earth.mod"\
	"..\Lib\Debug\Matrix.mod"\
	"..\occlib\Debug\Occ_refraction_adj.mod"\
	".\Debug\Atmosphere.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\GO_Propagator.f90

!IF  "$(CFG)" == "Wavelib - Win32 Release"

DEP_F90_GO_PR=\
	"..\Lib\Release\Coordinates.mod"\
	"..\Lib\Release\Defaults.mod"\
	"..\Lib\Release\Earth.mod"\
	"..\Lib\Release\Externals.mod"\
	"..\Lib\Release\Interpolation.mod"\
	"..\Lib\Release\IO.MOD"\
	"..\Lib\Release\Signal.mod"\
	"..\occlib\Release\Occ_Coordinates.mod"\
	"..\occlib\Release\Occ_Inversion.mod"\
	"..\occlib\Release\Occ_Meteoprofiles.mod"\
	"..\occlib\Release\Occ_Refraction.mod"\
	".\Release\Atmosphere.mod"\
	".\Release\Atmosphere_Rays.mod"\
	

!ELSEIF  "$(CFG)" == "Wavelib - Win32 Debug"

DEP_F90_GO_PR=\
	"..\Lib\Debug\Coordinates.mod"\
	"..\Lib\Debug\Defaults.mod"\
	"..\Lib\Debug\Earth.mod"\
	"..\Lib\Debug\Externals.mod"\
	"..\Lib\Debug\Interpolation.mod"\
	"..\Lib\Debug\IO.MOD"\
	"..\Lib\Debug\Signal.mod"\
	"..\occlib\Debug\Occ_Coordinates.mod"\
	"..\occlib\Debug\Occ_Inversion.mod"\
	"..\occlib\Debug\Occ_Meteoprofiles.mod"\
	"..\occlib\Debug\Occ_Refraction.mod"\
	".\Debug\Atmosphere.mod"\
	".\Debug\Atmosphere_Rays.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\Ray_Problem.f90

!IF  "$(CFG)" == "Wavelib - Win32 Release"

DEP_F90_RAY_P=\
	"..\Lib\Release\Coordinates.mod"\
	"..\Lib\Release\Defaults.mod"\
	"..\Lib\Release\Externals.mod"\
	"..\Lib\Release\Interpolation.mod"\
	"..\Lib\Release\IO.MOD"\
	"..\Lib\Release\Signal.mod"\
	"..\occlib\Release\Occ_Coordinates.mod"\
	"..\occlib\Release\Occ_Refraction.mod"\
	"..\occlib\Release\Occ_Refraction_adj.mod"\
	".\Release\Atmosphere_Rays.mod"\
	".\Release\Atmosphere_rays_adj.mod"\
	

!ELSEIF  "$(CFG)" == "Wavelib - Win32 Debug"

DEP_F90_RAY_P=\
	"..\Lib\Debug\Coordinates.mod"\
	"..\Lib\Debug\Defaults.mod"\
	"..\Lib\Debug\Externals.mod"\
	"..\Lib\Debug\Interpolation.mod"\
	"..\Lib\Debug\IO.MOD"\
	"..\Lib\Debug\Signal.mod"\
	"..\occlib\Debug\Occ_Coordinates.mod"\
	"..\occlib\Debug\Occ_Refraction.mod"\
	"..\occlib\Debug\Occ_refraction_adj.mod"\
	".\Debug\Atmosphere_Rays.mod"\
	".\Debug\Atmosphere_rays_adj.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\Vacuum_Propagator.f90

!IF  "$(CFG)" == "Wavelib - Win32 Release"

DEP_F90_VACUU=\
	"..\FFTW\Release\FFTW.mod"\
	"..\Lib\Release\Defaults.mod"\
	"..\Lib\Release\Externals.mod"\
	"..\Lib\Release\Interpolation.mod"\
	"..\Lib\Release\IO.MOD"\
	"..\Lib\Release\Matrix.mod"\
	"..\Lib\Release\Signal.mod"\
	"..\occlib\Release\Occ_Diffraction.mod"\
	

!ELSEIF  "$(CFG)" == "Wavelib - Win32 Debug"

DEP_F90_VACUU=\
	"..\FFTW\Debug\FFTW.mod"\
	"..\Lib\Debug\Defaults.mod"\
	"..\Lib\Debug\Externals.mod"\
	"..\Lib\Debug\Interpolation.mod"\
	"..\Lib\Debug\IO.MOD"\
	"..\Lib\Debug\Matrix.mod"\
	"..\Lib\Debug\Signal.mod"\
	"..\occlib\Debug\Occ_Diffraction.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\Wave_Propagator.f90

!IF  "$(CFG)" == "Wavelib - Win32 Release"

DEP_F90_WAVE_=\
	"..\FFTW\Release\FFTW.mod"\
	"..\Lib\Release\Coordinates.mod"\
	"..\Lib\Release\Debug.mod"\
	"..\Lib\Release\Defaults.mod"\
	"..\Lib\Release\Earth.mod"\
	"..\Lib\Release\Externals.mod"\
	"..\Lib\Release\Interpolation.mod"\
	"..\Lib\Release\IO.MOD"\
	"..\Lib\Release\Signal.mod"\
	"..\Lib\Release\Time.mod"\
	"..\occlib\Release\FIO.MOD"\
	"..\occlib\Release\Occ_Coordinates.mod"\
	"..\occlib\Release\Occ_Diffraction.mod"\
	".\Release\Atmosphere.mod"\
	".\Release\Atmosphere_Rays.mod"\
	".\Release\Vacuum_Propagator.mod"\
	

!ELSEIF  "$(CFG)" == "Wavelib - Win32 Debug"

DEP_F90_WAVE_=\
	"..\FFTW\Debug\FFTW.mod"\
	"..\Lib\Debug\Coordinates.mod"\
	"..\Lib\Debug\Debug.mod"\
	"..\Lib\Debug\Defaults.mod"\
	"..\Lib\Debug\Earth.mod"\
	"..\Lib\Debug\Externals.mod"\
	"..\Lib\Debug\Interpolation.mod"\
	"..\Lib\Debug\IO.MOD"\
	"..\Lib\Debug\Signal.mod"\
	"..\Lib\Debug\Time.mod"\
	"..\occlib\Debug\FIO.MOD"\
	"..\occlib\Debug\Occ_Coordinates.mod"\
	"..\occlib\Debug\Occ_Diffraction.mod"\
	".\Debug\Atmosphere.mod"\
	".\Debug\Atmosphere_Rays.mod"\
	".\Debug\Vacuum_Propagator.mod"\
	

!ENDIF 

# End Source File
# End Target
# End Project
