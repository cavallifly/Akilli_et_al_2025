#!/bin/bash

#SBATCH --job-name plotPuppy
#SBATCH -n 1                    # Number of cores. For now 56 is the number max of core available
#SBATCH -t 4-00:00              # Runtime in D-HH:MM
#SBATCH -o 02_plotTopDomDomainsPileUps_using_coolpuppy.out # File to which STDOUT will be written
#SBATCH -e 02_plotTopDomDomainsPileUps_using_coolpuppy.out # File to which STDERR will be written 


if [ "$#" -ne 2 ]; then
        echo "Please specify 1) input folder containing clpy files and 2) output folder for .svg files"
	exit 1
fi

echo "coolpuppy plotpup.py"
for FILE in $1*".clpy"; do
	filename=${FILE##*/}
	filename=${filename/".clpy"/}
	output=$2$filename".svg"

	if [ -f "$output" ];
	then
            echo "$output already exists"
            continue
        fi
	touch ${output}
	
	cmd="plotpup.py --input_pups $FILE -o $output"
        echo "$cmd"
        #touch"$ouput"
        
        $cmd
done
