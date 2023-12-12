#!/bin/bash
DIR=$(dirname "$(realpath "$0")")
# pogoBlockGreedy3d test_anisotropic.pogo-inp;
# pogoSolve3d test_anisotropic.pogo-inp <<< o;

# ************* create the scanning grid
nx=20
ny=20

pogoBlockGreedy3d woven_test_0901.pogo-inp;

# ************* create the base model
for ((step=1; step<=(nx+1)*(ny+1); step++))
# for step in "${array[@]}"
do
	export step;
 
	# pogoBlockGreedy3d woven_test_090${step}.pogo-inp;
    	
    	# Find the file
	filename=$(find . -name '*.pogo-block' -print -quit)
	# Check if file was found
	if [ ! -z "$filename" ]; then
    		echo "Found file: $filename"
    		# Rename the file
    		new_filename="woven_test_090${step}.pogo-block"  # Change new_extension to whatever you want
    		mv "$filename" "$new_filename"
    		echo "Renamed file to: $new_filename"
	else
    		echo "No .pogo-block file found."
	fi
	#
	# if (step==150); then
	# 	pogoSolve3d woven_test_090${step}.pogo-inp --setOutputFile "woven_test_${step}" <<< o;
	# else 
	
	pogoSolve3d woven_test_090${step}.pogo-inp --setFieldSaveOff --setOutputFile "woven_test_${step}" <<< o;
	# fi
    	rm woven_test_090${step}.pogo-inp;
done


# Define the full path to the directory containing the MATLAB script
matlabScriptDir="/home/xiaoyu/pogo_work/Read_results"

# Get the current directory of the Bash script
currentPath=$(pwd)

# Run the MATLAB script with the current path as an argument
# First change to the directory containing the MATLAB script, then execute the script
matlab -batch "cd('$matlabScriptDir'); fx_Cscanread('$currentPath')"


