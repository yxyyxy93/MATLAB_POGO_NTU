#!/bin/bash
# These are the command lines to run the pogo files

for file in *.pogo-inp; do
	echo "Processing $file"
	# Your commands here
	pogoBlockGreedy3d "$file"
	pogoSolve3d "$file" <<< o;
done




