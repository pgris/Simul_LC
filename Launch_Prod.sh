#!/bin/bash

fieldname=$1
fieldid=$2
thedir=$3
stretch=$4
color=$5
zmin=$6
zmax=$7
zstep=$8
zrandom=$9
T0random=${10}

#for season in {0..9}
for season in 0
do
   #echo "Welcome $season times"
   python Loop_Prod.py --fieldname $fieldname --fieldid $fieldid --zmin $zmin --zmax $zmax --stretch $stretch --color $color --season $season --nevts 250 --dirmeas $thedir --zstep ${zstep} --zrandom ${zrandom} --T0random ${T0random}
done