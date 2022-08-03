#!/bin/sh

occ=$1
wop=$2
opts=$3

waveoldx=/d/F90.old/Wave/Debug/Wave.exe
wavenewx=/d/F90/Wave/Debug/Wave.exe

optsbase="\
   -atm=echam,3d,ecmfc20070715121.grb \
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
   -go=y \
   -rh=y \
   -ab=n \
   -as=n \
   -fio=y \
   -sd=2 \
   -ngo=300 \
   -fwp=-1 \
   -v=3 \
   -o=all"

$waveoldx $optsbase $opts
mv $wop* old

$wavenewx $optsbase $opts
mv $wop* new

files=`cd old ; ls -1 *`

for file in $files ; do
echo --- $file ---
diff -a -q old/$file new/$file
done

rm -f old/* new/*
