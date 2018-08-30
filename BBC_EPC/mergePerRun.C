for X in `cat runs.bbc.dat`
do
rm out/out${X}.root
hadd out/out${X}.root out/${X}_*.root
done
