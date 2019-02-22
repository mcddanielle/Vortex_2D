#!/bin/bash

for dir1 in N*
do
    if [ -d $dir1 ]
    then

	python python*/hea*py $dir1

    fi
done
