#!/bin/sh

occ=$1
wop=$2
opts=$3

#waveoldx=/nas/CMAQ_Scratch/p2291/gorbunov_2013/F90/Wave/Wave.x
wavenewx=/nas/CMAQ_Scratch/p2291/gorbunov_2013/F90/Wave/Wave.x

#   -atm=echam,3d,ecmfc20070715121.grb \
#   -atm=ncep,3d,NCEP_FNL1X1/fnl_20140426_00_00.grib2 \
#   -atm=ncep,3d,NCEP_FNL1X1/fnl_20140426_00_00.grb \
#   -go=n \


optsbase="\
   -atm=mpas,3d,MPAS/output.2018-07-18_01.nc \
   -earth=geoid \
   -gps=$occ \
   -out=$wop \
   -freq=gps \
   -afreq=gps \
   -hmax=60 \
   -hmin=0 \
   -dyn=0.1
   -dx=20 \
   -sr=50 \
   -ts=n \
   -ecef=n \
   -wo=y \
   -go=n \
   -rh=y \
   -ab=n \
   -as=n \
   -fio=y \
   -sd=2 \
   -ngo=300 \
   -fwp=-1 \
   -v=3 \
   -o=all"

export LD_LIBRARY_PATH="/nas/home/jhegarty/mynetcf/lib"

#$waveoldx $optsbase $opts
#mv $wop* old

$wavenewx $optsbase $opts
#mv $wop* new

#files=`cd old ; ls -1 *`

#for file in $files ; do
#echo --- $file ---
#diff -a -q old/$file new/$file
#done

#rm -f old/* new/*
