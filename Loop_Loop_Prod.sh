#!/bin/bash

shift=0.1
thedir=Mean_Obs_newrefs
#thedir=DD

#for field in 290 744 1427 2412 2786 
for field in 120
do
#sh Launch_Prod.sh DD ${field} DD -999. -999. 0. 1.3 0.1
for z in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8
do
zmax=`echo $z + $shift | bc`
sh Launch_Prod.sh DD ${field} ${thedir} -2.0 0.2 $z $zmax 0.1 no no
done
#sh Launch_Prod.sh DD ${field} DD 2.0 -0.2 0. 1.1 0.1
done