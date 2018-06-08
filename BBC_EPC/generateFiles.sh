#!/bin/bash

for X in `cat ../runs.dat`
do
    hadd out/run${X}.root out/out_${X}_*root
    #root -b -l -q qcent.C\($X\)
    #root -b -l -q coef.C\($X\)
done
