# Microsoft Developer Studio Project File - Name="Wave" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=Wave - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "Wave.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "Wave.mak" CFG="Wave - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "Wave - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "Wave - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 1
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "Wave - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /include:"Release/" /nologo /warn:nofileopt
# ADD F90 /check:noflawed_pentium /compile_only /dbglibs /include:"Release/" /include:"../Lib/Release" /include:"../occlib/Release" /include:"../ECHAM/Release" /include:"../GRIB/Release" /include:"../Models/Release" /include:"../Phantom/Release" /include:"../FFTW/Release" /include:"../Wavelib/Release" /include:"../GCM/Release" /math_library:fast /nologo /traceback /warn:nofileopt /fast
# SUBTRACT F90 /assume:noaccuracy_sensitive /automatic /check:bounds /check:format /check:output_conversion /check:overflow /check:underflow
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /G6 /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x407 /d "NDEBUG"
# ADD RSC /l 0x407 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib /nologo /stack:0x10000000 /subsystem:console /machine:I386 /nodefaultlib:"libc"

!ELSEIF  "$(CFG)" == "Wave - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /debug:full /include:"Debug/" /nologo /warn:nofileopt
# ADD F90 /compile_only /dbglibs /debug:full /include:"Debug/" /include:"../Lib/Debug" /include:"../occlib/Debug" /include:"../ECHAM/Debug" /include:"../Models/Debug" /include:"../Phantom/Debug" /include:"../FFTW/Debug" /include:"../Wavelib/Debug" /include:"../GRIB/Debug" /include:"../GCM/Debug" /nologo /traceback /warn:nofileopt /warn:unused
# SUBTRACT F90 /check:bounds /check:noflawed_pentium
# ADD BASE CPP /nologo /W3 /Gm /GX /Zi /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /Gm /GX /Zi /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x407 /d "_DEBUG"
# ADD RSC /l 0x407 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib /nologo /stack:0x1000000 /subsystem:console /incremental:no /debug /machine:I386 /nodefaultlib:"libc.lib" /pdbtype:sept

!ENDIF 

# Begin Target

# Name "Wave - Win32 Release"
# Name "Wave - Win32 Debug"
# Begin Source File

SOURCE=.\Wave.f90

!IF  "$(CFG)" == "Wave - Win32 Release"

DEP_F90_WAVE_=\
	"..\Lib\Release\Coordinates.mod"\
	"..\Lib\Release\Defaults.mod"\
	"..\Lib\Release\Earth.mod"\
	"..\Lib\Release\Errors.mod"\
	"..\Lib\Release\Geoid.mod"\
	"..\Lib\Release\IO.MOD"\
	"..\Lib\Release\Strings.mod"\
	"..\Lib\Release\Time.mod"\
	"..\Models\Release\MSIS.mod"\
	"..\occlib\Release\Occ_Coordinates.mod"\
	"..\occlib\Release\Occ_IO.mod"\
	"..\occlib\Release\Occ_Orbits.mod"\
	"..\occlib\Release\Occ_Refraction.mod"\
	"..\Phantom\Release\Turbulence.mod"\
	"..\Wavelib\Release\Asymptotic_Propagator.mod"\
	"..\Wavelib\Release\Atmosphere.mod"\
	"..\Wavelib\Release\GO_Propagator.mod"\
	"..\Wavelib\Release\Wave_Propagator.mod"\
	".\Release\Wave_Comline.mod"\
	".\Release\Wave_Defaults.mod"\
	".\Release\Wave_Version.mod"\
	

!ELSEIF  "$(CFG)" == "Wave - Win32 Debug"

DEP_F90_WAVE_=\
	"..\Lib\Debug\Coordinates.mod"\
	"..\Lib\Debug\Defaults.mod"\
	"..\Lib\Debug\Earth.mod"\
	"..\Lib\Debug\Errors.mod"\
	"..\Lib\Debug\Geoid.mod"\
	"..\Lib\Debug\IO.MOD"\
	"..\Lib\Debug\Strings.mod"\
	"..\Lib\Debug\Time.mod"\
	"..\Models\Debug\MSIS.mod"\
	"..\occlib\Debug\Occ_Coordinates.mod"\
	"..\occlib\Debug\Occ_IO.mod"\
	"..\occlib\Debug\Occ_Orbits.mod"\
	"..\occlib\Debug\Occ_Refraction.mod"\
	"..\Phantom\Debug\Turbulence.mod"\
	"..\Wavelib\Debug\Asymptotic_Propagator.mod"\
	"..\Wavelib\Debug\Atmosphere.mod"\
	"..\Wavelib\Debug\GO_Propagator.mod"\
	"..\Wavelib\Debug\Wave_Propagator.mod"\
	".\Debug\Wave_Comline.mod"\
	".\Debug\Wave_defaults.mod"\
	".\Debug\Wave_Version.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\Wave_comline.f90

!IF  "$(CFG)" == "Wave - Win32 Release"

DEP_F90_WAVE_C=\
	"..\Lib\Release\Comline.mod"\
	"..\Lib\Release\Debug.mod"\
	"..\Lib\Release\Defaults.mod"\
	"..\Lib\Release\Earth.mod"\
	"..\Lib\Release\Strings.mod"\
	".\Release\Wave_Defaults.mod"\
	

!ELSEIF  "$(CFG)" == "Wave - Win32 Debug"

DEP_F90_WAVE_C=\
	"..\Lib\Debug\Comline.mod"\
	"..\Lib\Debug\Debug.mod"\
	"..\Lib\Debug\Defaults.mod"\
	"..\Lib\Debug\Earth.mod"\
	"..\Lib\Debug\Strings.mod"\
	".\Debug\Wave_defaults.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\Wave_defaults.f90

!IF  "$(CFG)" == "Wave - Win32 Release"

DEP_F90_WAVE_D=\
	"..\Lib\Release\Defaults.mod"\
	

!ELSEIF  "$(CFG)" == "Wave - Win32 Debug"

DEP_F90_WAVE_D=\
	"..\Lib\Debug\Defaults.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\Wave_version.f90
# End Source File
# End Target
# End Project
