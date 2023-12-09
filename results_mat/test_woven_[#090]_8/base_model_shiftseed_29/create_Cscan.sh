#!/bin/bash
DIR=$(dirname "$(realpath "$0")")
# Now $DIR contains the path to the directory where the script is located
# Use $DIR to reference files or directories relative to the script's location

# Get the full path of the current directory
currentPath=$(pwd)

# Extract the name of the current folder
currentFolder=$(basename "$currentPath")

# Extract the last part of the folder path (parent folder of the current folder)
lastFolderPath=$(dirname "$currentPath")
modelpath="${lastFolderPath}/${currentFolder}.mat"

# Display the results
echo "Current path: $currentPath"
echo "model Path: $modelpath"

# Run the MATLAB script with the paths as an argument
# First change to the directory containing the MATLAB script, then execute the script
matlab -batch "addpath('/home/xiaoyu/pogo_work/utlis_pogo_xiaoyu/'); fx_immersion_Cscan('$currentPath', '$modelpath')"
