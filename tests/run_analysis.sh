#!/bin/bash

file='frame_current_0.000000'

for dir1 in N52* N160* N256* 
do
    if [ -d $dir1 ]
    then
	number_pins=`echo $dir1 | cut -d'_' -f2 | cut -c 3-`
	pin_force=`echo $dir1 | cut -d'_' -f3 | cut -c3-`

	cd $dir1

	#echo $dir1' has '$number_pins' pins and strength fp= '$pin_force

	$HOME/Codes-Scripts/Vortex_2D/vortex_geometry -f $file -p $pin_force -N $number_pins -s 36.5 --quiet >> ../data.txt

	cd -
    fi
done
