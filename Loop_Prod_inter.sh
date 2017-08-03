#!/bin/bash

step=0.1
stretch=-2.0
color=0.2
num=$(awk 'BEGIN{for(i=0.3;i<=0.4;i+=0.1)print i}')
for n in $num
do
  sum=$(awk "BEGIN {print $n+$step; exit}")
  echo $n,$sum
  python simu_for_cadence.py --fieldid 309 --zmin $n --zmax $sum --season 1 --nevts 100 --stretch $stretch --color $color
done
#python simu_for_cadence.py --fieldid 309 --zmin 0. --zmax 0.1 --season 1 --nevts 1000