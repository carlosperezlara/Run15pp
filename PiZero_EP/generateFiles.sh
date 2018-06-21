#!/bin/bash

for X in `cat ../runs.emcal.bbc.dat`
do
    hadd out/run${X}.root out/out_${X}_*root
    #hadd outERT/run${X}.root outERT/out_${X}_*root
    #hadd outERT060/run${X}.root outERT060/out_${X}_*root
done
