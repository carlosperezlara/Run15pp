#!/bin/bash

for X in `cat ../runs.emcal.dat`
do
    hadd out/run${X}.root out/out_${X}_*root
done
