#!/bin/bash

#PBS -N wave 
#PBS -A uate0003
#PBS -l walltime=06:00:00
#PBS -q regular
#PBS -j oe
#PBS -k eod
#PBS -m abe
#PBS -M jhegarty@aer.com
#PBS -l select=1:ncpus=1:mpiprocs=1

occ=IN_DATA/atmPhs_C006.2014.115.23.51.G19_2013.3520_nc
wop=wop_output
opts="-v=5"

#waveoldx=/nas/CMAQ_Scratch/p2291/gorbunov_2013/F90/Wave/Wave.x
wavenewx=/glade/u/home/hegarty/Gorbunov_2013/F90/Wave/Wave.x

#   -atm=mpas,3d,IN_DATA/output.2018-07-18_01.nc \
#   -atm=echam,3d,ecmfc20070715121.grb \
#   -atm=ncep,3d,NCEP_FNL1X1/fnl_20140426_00_00.grib2 \
#   -atm=ncep,3d,NCEP_FNL1X1/fnl_20140426_00_00.grb \
#   -go=n \


optsbase="\
   -atm=mpas,3d,X5WP/history_200_gw.2019-10-16_22.00.00.nc \
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

# Do I need to do this on Cheyenne?
#export LD_LIBRARY_PATH="/nas/home/jhegarty/mynetcf/lib"

#$waveoldx $optsbase $opts
#mv $wop* old

echo $optsbase
echo $opts

# for larger MPAS output files
ulimit -s unlimited

./Wave.x $optsbase $opts
#mv $wop* new

#files=`cd old ; ls -1 *`

#for file in $files ; do
#echo --- $file ---
#diff -a -q old/$file new/$file
#done

#rm -f old/* new/*
