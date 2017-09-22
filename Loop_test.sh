#!/bin/bash

for fieldid in 120 121 122 123
do
python simu_for_cadence.py --zmin 0.95 --zmax 0.96 --nevts 1 --fieldname DD --fieldid ${fieldid} --season 0 --stretch 2.0 --color -0.2 --dirmeas Mean_Obs_newrefs --T0random no --zrandom no
done