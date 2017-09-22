#!/bin/bash

#for fieldid in {132..151}
for fieldid in 1427
do
for season in 0
do
python Loop_Prod.py --fieldname DD_${season}_1 --fieldid ${fieldid} --zmin 0. --zmax 1.1 --stretch -999. --color -999. --season 0 --nevts 1000 --dirmeas Mean_Obs_newrefs --T0random yes
#python Loop_Prod.py --fieldname DD_${season}_2 --fieldid ${fieldid} --zmin 0. --zmax 1.1 --stretch 0. --color 0. --season 0 --nevts 1000 --dirmeas Mean_Obs_newrefs --T0random yes
#python Loop_Prod.py --fieldname DD --fieldid ${fieldid} --zmin 0. --zmax 1.1 --stretch 2.0 --color -0.2 --season 0 --nevts 1000 --dirmeas Mean_Obs_newrefs --T0random yes
#python Loop_Prod.py --fieldname DD --fieldid ${fieldid} --zmin 0. --zmax 1.1 --stretch -2.0 --color 0.2 --season 0 --nevts 1000 --dirmeas Mean_Obs_newrefs --T0random yes
done
done