# Define the full path to the directory containing the MATLAB script
matlabScriptDir="/home/xiaoyu/pogo_work/Read_results"

# Get the current directory of the Bash script
currentPath=$(pwd)

# Run the MATLAB script with the current path as an argument
# First change to the directory containing the MATLAB script, then execute the script
matlab -batch "cd('$matlabScriptDir'); addpath(genpath('/home/xiaoyu/pogo_work/utlis_pogo_xiaoyu')); fx_Cscanread('$currentPath')"

