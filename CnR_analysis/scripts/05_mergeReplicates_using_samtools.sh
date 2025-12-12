echo "Consider all the analysed samples"
# Sample name assayName_targetName_cellName_condition_replicateName_laneName_assembly_tag
# e.g., cnr_H3K27me3_WD_NOnub_Rep1_all_dm6_NA
replicates=$(ls -1 | grep cnr | grep _Rep | grep _all_) # | sed "s,_, ,g" | awk '{for(i=1;i<=4;i++){printf("%s_",$i)};printf("%s_","merge");for(i=6;i<=(NF-1);i++){printf("%s_",$i)};printf("%s\n",$NF)}' | uniq)
echo ${replicates}

for replicate in ${replicates}
do
    if [[ ! -d ${replicate} ]];
    then
	echo "${replicate} is not a directory for a sample"
	continue
    fi 
    
    echo ${replicate}
    
    assayName=$(echo ${replicate} | sed "s,_, ,g" | awk '{print $1}')
    targetName=$(echo ${replicate} | sed "s,_, ,g" | awk '{print $2}')
    cellName=$(echo ${replicate} | sed "s,_, ,g" | awk '{print $3}')
    condition=$(echo ${replicate} | sed "s,_, ,g" | awk '{print $4}')
    replicateName=merge
    laneName=$(echo ${replicate} | sed "s,_, ,g" | awk '{print $6}')
    assembly=$(echo ${replicate} | sed "s,_, ,g" | awk '{print $7}')
    tag=$(echo ${replicate} | sed "s,_, ,g" | awk '{print $8}')
    outDir=${assayName}_${targetName}_${cellName}_${condition}_${replicateName}_${laneName}_${assembly}_${tag}
    
    echo "These are the files merged in this replicate"
    bamFiles=$(ls -1 ${PWD}/${assayName}_${targetName}_${cellName}_${condition}_*${laneName}_${assembly}_${tag}/*sorted_dedup.bam)
    echo $bamFiles | awk '{for(i=1;i<=NF;i++){printf("%s\n", $i)}}'
    
    if [[ -d ${outDir} ]];
    then
	echo "Replicates already merged in ${outDir}!"
	continue
    fi
    echo "Creating merged sample ${outDir}"
    
    mkdir -p ${outDir}
    cd ${outDir}

    outBam=${outDir}.sorted_dedup.bam    
    echo "1 - Creating merged ${outBam} using samtools $(samtools 2>&1 | grep Version)"
    samtools merge -@ 16 ${outBam} ${bamFiles}
    # -@ Set number of threads.
	
    echo "2 - Sorting and indexing ${outBam} using samtools $(samtools 2>&1 | grep Version)"
    samtools sort -@ 8 -o sorted_${outBam} -O bam ${outBam}
    mv -v sorted_${outBam} ${outBam}
    # -o Write the final sorted output to FILE, rather than to standard output.
    # -@ Set number of sorting and compression threads. By default, operation is single-threaded.
    # -O Write the final output as sam, bam, or cram.
    samtools index -@ 16 ${outBam}
    # -@ Set number of threads.

    cd .. # Exit ${outDir}
    echo ""
    
done # Close cycle over $replicates
