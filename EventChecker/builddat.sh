#!bash

rm table.dat
for X in `cat ../sss.dat`
do
    root -b -l -q multi.C\(\"$X\"\)
done
