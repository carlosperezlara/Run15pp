#!/bin/tcsh

@ N=$1 + 1
set A = `head -$N segments.dat| tail -1`
echo $A
@ B=-1

#./Run_MX_EPC $A $B
#./Run_BBC_EPC $A $B
#./Run_PiZero $A $B

./Run_PiZero_EP $A $B D0
./Run_PiZero_EP $A $B D1
./Run_PiZero_EP $A $B A0
./Run_PiZero_EP $A $B A1
./Run_PiZero_EP $A $B T0
./Run_PiZero_EP $A $B T1

#./Run_Charged_QC $A $B

exit
