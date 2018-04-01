#!/bin/tcsh

@ N=$1 + 1
set A = `head -$N segments.dat| tail -1`
echo $A
@ B=-1

echo "FIRST: BUILD <QX> <QY>"
./ep $A $B
echo "SECOND: BUILD DPSI"
./ep $A $B
echo "THIRD: ENJOY"
./ep $A $B

#./ep $A $B

exit
