#!/bin/bash

fieldname=$1
fieldid=$2
thedir=$3
stretch=$4
color=$5
zmin=$6
zmax=$7
for season in {0..9}
do
   #echo "Welcome $season times"
   python Loop_Prod.py --fieldname $fieldname --fieldid $fieldid --zmin $zmin --zmax $zmax --stretch $stretch --color $color --season $season --nevts 1000 --dirmeas $thedir
done