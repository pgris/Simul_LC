#!/bin/bash

for field in 290 1427 2412 2786
do
sh Launch_Prod.sh DD ${field} DD -999. -999. 0. 1.1 0.05
done