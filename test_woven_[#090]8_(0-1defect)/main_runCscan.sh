#!/bin/bash

# Function to copy a script to a folder
copy_script_to_folder() {
    local script_path=$1
    local target_folder=$2
    # Check if the folder exists
    if [ ! -d "$target_folder" ]; then
        # Folder does not exist, so create it
        mkdir "$target_folder"
    fi
    cp -f "$script_path" "$target_folder"
}

# Path to the master script that you want to copy
MASTER_SCRIPT_PATH="run_woven_8l_test.sh"

# Starting and ending seeds
START_SEED=501
END_SEED=600

# Loop over the seed range
for (( SEED=$START_SEED; SEED<=$END_SEED; SEED+=4 )); do
    # Construct the TARGET_FOLDERS array
    TARGET_FOLDERS=(
        "base_model_shiftseed_$SEED"
        "base_model_shiftseed_$(($SEED+1))"
        "base_model_shiftseed_$(($SEED+2))"
        "base_model_shiftseed_$(($SEED+3))"
    )

    # Copy the master script to each target folder
    for folder in "${TARGET_FOLDERS[@]}"; do
        copy_script_to_folder "$MASTER_SCRIPT_PATH" "$folder"
		copy_script_to_folder "focusing_delays.txt" "$folder"
    done

    # Run the scripts in each folder
    for dir in "${TARGET_FOLDERS[@]}"; do
        (cd "$dir" && bash "$MASTER_SCRIPT_PATH") &
    done

    wait # Wait for all background jobs to finish

    # Optional: wait or perform any other operations between runs
    sleep 1
done

