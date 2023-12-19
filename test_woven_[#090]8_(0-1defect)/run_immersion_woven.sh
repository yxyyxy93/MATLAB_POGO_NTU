#!/bin/bash

# Define MATLAB executable path
# MATLAB_RUN="/usr/local/MATLAB/R20xxa/bin/matlab"

# Run MATLAB script with a path and all its subdirectories
matlab -nodisplay -nosplash -nodesktop -r " \
addpath(genpath('/home/xiaoyu/pogo_work/utlis_pogo_xiaoyu')); \
run('immersion_woven.m'); \
exit; \
"


