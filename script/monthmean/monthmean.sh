#!/bin/bash
prog=monthmean
gfortran -o $prog.exe $prog.f90

sy=1979
ly=2016
HS="N"
d="ncep"
lev="slp"
#mon=( "Jan" "Feb" "Mar" "Apr" "May" "Jun" "Jul" "Aug" "Sep" "Oct" "Nov" "Dec" )

Din=~/work/Projects/Arctic/DATA/NCEP/trk/
Dout=~/work/Projects/Arctic/output/trkmonmean/
reg=( "" ".90_225" ".160_340" ".180_280" ".180_250" )
echo $reg
#reg=.BeaufortS
nreg=0
lreg=0
while [ $nreg -le $lreg ]
do

fout=${Dout}av.${lev}.$d.$sy-$ly$HS${reg[$nreg]}.dat
fssn=${Dout}cycstat.${lev}.$d.jja$sy-$ly$HS${reg[$nreg]}.dat
fmay=${Dout}cycstat.${lev}.$d.may$sy-$ly$HS${reg[$nreg]}.dat
fjun=${Dout}cycstat.${lev}.$d.jun$sy-$ly$HS${reg[$nreg]}.dat
fjul=${Dout}cycstat.${lev}.$d.jul$sy-$ly$HS${reg[$nreg]}.dat
faug=${Dout}cycstat.${lev}.$d.aug$sy-$ly$HS${reg[$nreg]}.dat
fsep=${Dout}cycstat.${lev}.$d.sep$sy-$ly$HS${reg[$nreg]}.dat

y=$sy
while [ $y -le $ly ]
do

fin=${Din}trkdat.${lev}.$d.$y
fouty=${Dout}av.${lev}.$d.$y$HS${reg[$nreg]}

./$prog.exe<<mark
$fin
$y
$fouty
$nreg
mark

y=$[$y+1]
done

#gfortran -o av.exe av_m.f90
#./av.exe<<mark
# #$Dout
# $HS
# $fout
# $fssn
# $fmay
# $fjun
# $fjul
# $faug
# $fsep
# $sy
# $ly
# $nreg
# mark
# rm av.exe
# #rm $Dout????$HS

nreg=$[$nreg+1]
done
rm $prog.exe

exit
