#!/bin/bash
echo "Consider each of the analysed replicate samples obtained from merging all the lanes"
# Sample name assayName_targetName_cellName_condition_replicateName_laneName_assembly_tag
# e.g., chipseq_H3K27me3_WD_NOnub_Rep1_L1_dm6_NA
replicates=$(ls -1 | grep cnr | sed "s,_, ,g" | awk '{for(i=1;i<=5;i++){printf("%s_",$i)};printf("%s_","all");for(i=7;i<=(NF-1);i++){printf("%s_",$i)};printf("%s\n",$NF)}' | uniq | grep Rep | grep all)

echo "I am going to remove PCR-duplicates for each of the following replicates:"
echo ${replicates}
echo "using sambamba markdup in version"
echo "$(sambamba markdup -v 2>&1 | grep -vi unrec | awk '{if(NF>0){print $0}}')"

for replicate in ${replicates}
do
    echo ${replicate}
    echo "Analysing $replicate"

    cd $replicate

    for bamFile in $(ls -1 *bam 2> /dev/null | grep -v _dedup.bam) ;
    do
	if [[ -e ${bamFile%.bam}_dedup.bam ]];
	then	    
	    continue
	fi
	
	echo "4 - Removing PCR-duplicates using sambamba"
	sambamba markdup --remove-duplicates --tmpdir /zssd/scratch/tmp --show-progress --nthreads 8 ${bamFile} ${bamFile%.bam}_dedup.bam
	# Usage: sambamba-markdup [options] <input.bam> [<input2.bam> [...]] <output.bam>
        # By default, marks the duplicates without removing them
	# Options: --remove-duplicates
        #              remove duplicates instead of just marking them
	#          --tmpdir=TMPDIR
        #              specify directory for temporary files
	#          --show-progress
        #              show progressbar in STDERR
        #          --nthreads=NTHREADS
        #              number of threads to use
	
    done # Close cycle over $fastqFile1

    cd .. # Exit $sample dir
    
done # Close cycle over $sample
