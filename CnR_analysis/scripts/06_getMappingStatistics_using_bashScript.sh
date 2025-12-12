echo "Consider all the analysed samples"
# Sample name assayName_targetName_cellName_condition_replicateName_laneName_assembly_tag
# e.g., chipseq_H3K27me3_WD_NOnub_Rep1_L1_dm6_NA
replicates=$(ls -1 | grep cnr_)
echo ${replicates}

for replicate in ${replicates}
do
    if [[ ! -d ${replicate} ]];
    then
	continue
    fi

    if [[ -e mappingStatistics.log ]];
    then
	checkIfDone=$(grep -w ${replicate} mappingStatistics.log | wc -l)
	if [[ ${checkIfDone} -eq 1 ]];
	then
	    continue
	fi
    fi
	
    cd ${replicate}
    pwd

    replicateName=$(echo ${replicate} | sed "s,_, ,g" | awk '{print $5}')
    
    if [[ ${replicateName} == "merge" ]];
    then
	totalReads=NA
	mappedReads=NA
	filteredReads=NA
	uniqueReads=$(samtools view -c *bam 2> /dev/null)
	
	echo ${replicate} ${totalReads} ${mappedReads} ${filteredReads} ${uniqueReads} >> ../mappingStatistics.log
	cd ..
	echo ""
	continue
    fi
    
    
    
    
    logBowtie2=$(ls -1 *_bowtie2.log 2> /dev/null)
    if [[ -e ${logBowtie2} ]];
    then
	cat ${logBowtie2}
	
	echo "+ Total is the Total number of reads in the input fastq(s)"
	totalReads=$(grep "reads; of these:" ${logBowtie2} | awk '{print $1}')
	
	echo "+ Mapped is the Number of mapped reads/read-pairs"
	checkSingleEnd=$(grep "Aligning single-end sample" ${logBowtie2} | wc -l | awk '{print $1}')
	if [[ ${checkSingleEnd} -eq 1 ]];
	then
	    mappedReads=$(grep aligned ${logBowtie2} | grep -v "0 times" | awk '{s+=$1}END{print s}')
	else
	    mappedReads=$(grep overall ${logBowtie2} | sed "s/\%//g" | awk -v t=${totalReads} '{print int($1*t/100)}')
	fi
    else
	echo "logBowtie2 file doesn't exist!"
	totalReads=NA
	mappedReads=NA
    fi

    rm -fr _tmp
    for file in $(ls -1 *sorted_dedup.bam *.sorted.bam 2> /dev/null);
    do
	echo $file $(samtools view -c ${file} 2> /dev/null) >> _tmp
    done
    cat _tmp
	
    echo "+ Filtered is the Number of reads/read-pairs after filtering for MAPQ>30"
    filteredReads=$(grep -v dedup _tmp | awk 'BEGIN{v=0}{v=$2}END{print v}')

    echo "+ Unique is the Number of reads/read-pairs after removing duplicates"  
    uniqueReads=$(grep dedup _tmp | awk 'BEGIN{v=0}{v=$2}END{print v}')
    
    echo ${replicate} ${totalReads} ${mappedReads} ${filteredReads} ${uniqueReads} >> ../mappingStatistics.log
    #rm -fr _tmp
    
    cd .. # Exit ${replicate}
    echo ""

done # Close cycle over $replicates
