#!/bin/tcsh

@ N=$1 + 1
set A = `head -$N segments.dat| tail -1`
echo $A
@ B=-1

root ploteventplane.C\(\"$A\"\)
