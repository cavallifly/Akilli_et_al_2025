#!/bin/bash

#SBATCH --job-name TADPileUps
##SBATCH -n 1                    # Number of cores. For now 56 is the number max of core available
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH -t 4-00:00              # Runtime in D-HH:MM
#SBATCH -o 01_TopDomDomainsPileUps_using_coolpuppy_1000bp.out # File to which STDOUT will be written
#SBATCH -e 01_TopDomDomainsPileUps_using_coolpuppy_1000bp.out # File to which STDERR will be written 

states="All PcG Active Null Het"

resolution=1000

for state in ${states}
do
    
    inDomainFile=TopDom_domains_hic_WD_NOnub_merge_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp_with_state.bed
    outDomainFile=TopDom_domains_hic_WD_NOnub_merge_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp_with_state_${state}state_for_${resolution}bp.bed    
    if [[ ! -e ${outDomainFile} ]];
    then
	grep ${state} ${inDomainFile} | awk -v r=${resolution} '{printf("%s\t%s\t%s\n",$1,int($2/r)*r,int($3/r)*r)}' > ${outDomainFile}
	if [[ $state == "All" ]];
	then
	    awk -v r=${resolution} '{printf("%s\t%s\t%s\n",$1,int($2/r)*r,int($3/r)*r)}' ${inDomainFile} > ${outDomainFile}	
	fi
    fi
    
    for condition in WD_NOnub WD_SUMOnub ;
    do
	
	for mcoolFile in $(ls -1 ../01_cool_files/hic_*${condition}*.mcool | sed "s,/, ,g" | awk '{print $NF}')
	do
	    coolFile=../01_cool_files/${mcoolFile}::resolutions/${resolution}	    
	    expFile=expectedCis_${mcoolFile%.mcool}_${resolution}bp.tsv
	    outFile=${mcoolFile%.mcool}_at_${resolution}bp_${state}state_local_rescaled.clpy
	    echo $coolFile $expFile $outFile

	    if [[ ! -e ${expFile} ]];
	    then
		touch ${expFile}		
		cooltools expected-cis ${coolFile} -o ${expFile} --nproc 8 --ignore-diags 0

		head ${expFile}    
	    fi
	    
	    if [[ ! -e ${outFile} ]];
	    then
		touch ${outFile}
		coolpup.py ${coolFile} ${outDomainFile} --rescale --local --expected ${expFile} -o ${outFile}
	    fi	    
	    ls -lrtha ${outFile}
	done
    done
done
