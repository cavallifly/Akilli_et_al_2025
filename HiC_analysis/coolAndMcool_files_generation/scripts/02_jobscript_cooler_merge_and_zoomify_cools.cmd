#!/bin/bash
#SBATCH --job-name MergeZoomify 
#SBATCH -n 1                    # Number of cores. For now 56 is the number max of core available
#SBATCH --mem=16000             # allocated memory per CPU
#SBATCH --partition=computepart # specify queue partiton
#SBATCH -t 10-00:00              # Runtime in D-HH:MM
#SBATCH -o cooler_merge_and_zoomify_cools.out   # File to which STDOUT will be written
#SBATCH -e cooler_merge_and_zoomify_cools.out   # File to which STDERR will be written


coolFile=$1

#outFile=${coolFile%.cool}_500bp.mcool
#if [[ ! -e ${outFile} ]];
#then
    #conda run -n cooler cooler zoomify -r 1000,2000,5000,10000,15000,20000,40000,50000,100000 ${coolFile} -o ${coolFile%.cool}.mcool --balance 2> 02_cool2mcool_${coolFile%.cool}.out
#    conda run -n cooler cooler zoomify -r 500 ${coolFile} -o ${outFile} --balance 2> 02_cool2mcool_${outFile%.mcool}.out    
#fi

outFile=${coolFile%.cool}_1000bp.mcool
if [[ ! -e ${outFile} ]];
then

    ### For loop analysis
    conda run -n cooler cooler zoomify -r 1000,2000,4000,5000,8000,10000,20000,40000,80000,100000,200000 ${coolFile} -o ${coolFile%.cool}.mcool --balance 2> 02_cool2mcool_${coolFile%.cool}.out
fi

#outFile=${coolFile%.cool}_5000bp.mcool
#if [[ ! -e ${outFile} ]];
#then
    ### For plot
#    conda run -n cooler cooler zoomify -r 5000 ${coolFile} -o ${coolFile%.cool}.mcool --balance 2> 02_cool2mcool_${coolFile%.cool}.out    
#fi

