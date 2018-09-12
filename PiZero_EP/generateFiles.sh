#!/bin/bash

for X in `cat ../runs.bbc.dat`
do
    rm out/run${X}.root
    hadd out/run${X}.root out/out_${X}_*root
done
rm allfiles/allNOM.root
hadd allfiles/allNOM.root out/run*root

exit 1



for Y in A0 A1 D0 D1 T0 T1 FA0 FA1 FD0 FD1 FT0 FT1
do
    for X in `cat ../runs.bbc.dat`
    do
	rm out${Y}/run${X}.root
	hadd out${Y}/run${X}.root out${Y}/out_${X}_*root
    done
    rm allfiles/all${Y}.root
    hadd allfiles/all${Y}.root out${Y}/run*root
done

for Y in ERT ERTA0 ERTA1 ERTD0 ERTD1 ERTT0 ERTT1 ERTFA0 ERTFA1 ERTFD0 ERTFD1 ERTFT0 ERTFT1
do
    for X in `cat ../runs.emcal.bbc.dat`
    do
	rm out${Y}/run${X}.root
	hadd out${Y}/run${X}.root out${Y}/out_${X}_*root
    done
    rm allfiles/all${Y}.root
    hadd allfiles/all${Y}.root out${Y}/run*root
done


