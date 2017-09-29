#!/bin/bash

for field in 744
do
sh Launch_Prod.sh DD ${field} DD -999. -999. 0. 1.3 0.1
#sh Launch_Prod.sh DD ${field} DD -2.0 0.2 0. 1.1 0.1
#sh Launch_Prod.sh DD ${field} DD 2.0 -0.2 0. 1.1 0.1
done