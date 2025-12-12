#!/bin/bash
#SBATCH --job-name DumpText
#SBATCH -n 1                    # Number of cores. For now 56 is the number max of core available
#SBATCH --mem=16000             # allocated memory per CPU
#SBATCH --partition=computepart # specify queue partiton
#SBATCH -t 10-00:00              # Runtime in D-HH:MM
#SBATCH -o cooler_dump_text_matrices.out   # File to which STDOUT will be written
#SBATCH -e cooler_dump_text_matrices.out   # File to which STDERR will be written

#outDir=balanced_text_matrices
outDir=observed_text_matrices
mkdir -p ${outDir} 

for resolution in 2000 4000;
do
   for coolFile in $(ls -1 *WD*merge*.mcool );
   do
       outFile=${outDir}/${coolFile%.mcool}_at_${resolution}bp.tab
       if [[ -e ${outFile} ]];
       then
	   continue
       fi
       touch ${outFile}
       conda run -n cooler cooler dump --balanced --join ${coolFile}::/resolutions/${resolution} > ${outFile}
       
   done # Close cycle over $coolFile
done # Close cycle over $resolution
