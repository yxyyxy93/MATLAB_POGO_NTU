#!/bin/bash

# Starting and ending seeds
START_SEED=483
END_SEED=485

# Loop over the seed range
for (( SEED=$START_SEED; SEED<=$END_SEED; SEED+=4 )); do
    # Construct the TARGET_FOLDERS array
    TARGET_FOLDERS=(
        "base_model_shiftseed_$SEED"
        "base_model_shiftseed_$(($SEED+1))"
        "base_model_shiftseed_$(($SEED+2))"
        "base_model_shiftseed_$(($SEED+3))"
    )
	
	nx=16;
	ny=16;
	dx=0.25e-3;
	dy=0.25e-3;
		
    # Copy the master script to each target folder
    for folder in "${TARGET_FOLDERS[@]}"; do
		# Get the full path of the current directory
		currentPath=$(pwd)
		modelpath="${currentPath}/${folder}.mat"
		folderpath="${currentPath}/${folder}"
		# Check if the model file exists
		if [ -f "$modelpath" ]; then
			echo "Model file found: $modelpath"
		else
			echo "Model file does not exist: $modelpath"
		fi
		# Check if the folder exists
		if [ -d "$folderpath" ]; then
			echo "Folder found: $folderpath"
		else
			echo "Folder does not exist: $folderpath"
		fi
		# Run the MATLAB script with the paths as an argument
		# First change to the directory containing the MATLAB script, then execute the script
		matlab -batch "addpath(genpath('/home/xiaoyu/pogo_work/utlis_pogo_xiaoyu')); fx_SaveStruct('$folderpath', '$modelpath', $nx, $ny, $dx, $dy)"
    done
	
    wait # Wait for all background jobs to finish
    # Optional: wait or perform any other operations between runs
    sleep 1
done

