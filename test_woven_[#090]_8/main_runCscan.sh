#!/bin/bash

# Function to copy a script to a folder
copy_script_to_folder() {
    local script_path=$1
    local target_folder=$2

    cp -f "$script_path" "$target_folder"
}

# Path to the master script that you want to copy
MASTER_SCRIPT_PATH="run_woven_8l_test.sh"

# List of target folders
TARGET_FOLDERS=(
"base_model_shiftseed_26"
"base_model_shiftseed_27" 
"base_model_shiftseed_28"
"base_model_shiftseed_29"
"base_model_shiftseed_30"
)

# Copy the master script to each target folder
for folder in "${TARGET_FOLDERS[@]}"; do
    copy_script_to_folder "$MASTER_SCRIPT_PATH" "$folder"
done

# Run the scripts in each folder
for dir in "${TARGET_FOLDERS[@]}"; do
    (cd "$dir" && bash "$MASTER_SCRIPT_PATH") &
done

wait # Wait for all background jobs to finish
