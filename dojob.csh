#!/bin/tcsh

@ N=$1 + 1
set A = `head -$N segments.dat| tail -1`
echo $A
@ B=-1

./Run_MX_EPC $A $B
#./Run_BBC_EPC $A $B
#./Run_PiZero $A $B
#./Run_PiZero_EP $A $B
#./Run_Charged_QC $A $B

exit
