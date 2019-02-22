#!/bin/bash

file='frame_current_0.000000'

DATE=$(date +"%Y%m%d%H%M")
image_dir="eps_images_$DATE"
mkdir $image_dir


for dir1 in N*
do
    if [ -d $dir1 ]
    then
	num_particles=`echo $dir1 | cut -d'_' -f1`
	number_pins=`echo $dir1 | cut -d'_' -f2 | cut -c 3-`
	pin_force=`echo $dir1 | cut -d'_' -f3 | cut -c3-`

	cd $dir1

	#run movie converter
	if [ ! -e $file ]
	then
	    #movie converter commands
	    echo "aack!"
	fi

	#add the s(k) and gplot analysis here
	#$HOME/Codes-Scripts/StructureFactor/calc_struct_factor $file < $HOME/Codes-Scripts/StructureFactor/input_file

	#gnuplot $HOME/Codes-Scripts/StructureFactor/gplot_eps_png

	#echo $dir1' has '$number_pins' pins and strength fp= '$pin_force
	
	#$HOME/Codes-Scripts/Vortex_2D/vortex_geometry -f $file -p $pin_force -N $number_pins -s 36.5 --quiet >> ../$image_dir/$num_particles'_'$number_pins'_data.txt'
	#echo $num_particles'_'$number_pins'_data.txt'
	$HOME/Codes-Scripts/Vortex_2D/tests/contour_plots/contourf3d_demo2.py
	mv histogram.eps ../$image_dir/histogram_$num_particles'_'$number_pins'_'$pin_force.eps
	cd -
    fi
#exit #just run test case
done