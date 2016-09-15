#!/bin/bash
prog=monthmean
gfortran -o $prog.exe $prog.f90

sy=1979
ly=2011
HS="N"
#mon=( "Jan" "Feb" "Mar" "Apr" "May" "Jun" "Jul" "Aug" "Sep" "Oct" "Nov" "Dec" ) 

Din=/work/irudeva/tracking/UOM/MERRA/out/merra.
Dout=/work/irudeva/Arctic/MERRA/monthmean/
reg=( "" ".90_225" ".160_340" ".180_280" ".180_250" )
echo $reg
#reg=.BeaufortS
nreg=3
lreg=3
while [ $nreg -le $lreg ]
do

fout=${Dout}av.$sy-$ly$HS${reg[$nreg]}.dat
fssn=${Dout}cycstat.jja$sy-$ly$HS${reg[$nreg]}.dat
fmay=${Dout}cycstat.may$sy-$ly$HS${reg[$nreg]}.dat
fjun=${Dout}cycstat.jun$sy-$ly$HS${reg[$nreg]}.dat
fjul=${Dout}cycstat.jul$sy-$ly$HS${reg[$nreg]}.dat
faug=${Dout}cycstat.aug$sy-$ly$HS${reg[$nreg]}.dat
fsep=${Dout}cycstat.sep$sy-$ly$HS${reg[$nreg]}.dat



y=$sy
while [ $y -le $ly ]
do

fin=$Din$y$HS.trk
fouty=${Dout}av.$y$HS${reg[$nreg]}

#./$prog.exe<<mark
#$fin
#$y
#$fouty
#$nreg
#mark

y=$[$y+1]
done

gfortran -o av.exe av_m.f90
./av.exe<<mark
$Dout
$HS
$fout
$fssn
$fmay
$fjun
$fjul
$faug
$fsep
$sy
$ly
$nreg
mark
rm av.exe
#rm $Dout????$HS

nreg=$[$nreg+1]
done
rm $prog.exe

exit
