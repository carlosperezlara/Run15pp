#!/bin/tcsh

@ N=$1 + 1
set A = `head -$N segments.dat| tail -1`
echo $A
@ B=-1

./Run_BBC_EPC $A $B

exit