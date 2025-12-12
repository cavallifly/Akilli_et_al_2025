echo "Consider all the analysed samples"
# Sample name assayName_targetName_cellName_condition_replicateName_laneName_assembly_tag
# e.g., chipseq_H3K27me3_WD_NOnub_Rep1_L1_dm6_NA
replicates=$(ls -1 | grep cnr)
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
    replicateName=$(echo ${replicate} | sed "s,_, ,g" | awk '{print $5}')
    laneName=$(echo ${replicate} | sed "s,_, ,g" | awk '{print $6}')
    assembly=$(echo ${replicate} | sed "s,_, ,g" | awk '{print $7}')
    tag=$(echo ${replicate} | sed "s,_, ,g" | awk '{print $8}')

    #effectiveGenomeSize=$(faCount /zssd/scratch/DBs/*/UCSC/${assembly}/Sequence/WholeGenomeFasta/genome.fa | grep total | awk '{print $2}')
    effectiveGenomeSize=$(awk -v a=${assembly} '{if($1==a){print $2}}' /home/common_pipelines/utils/list_of_effective_genomic_size.tab)
    if [[ ${effectiveGenomeSize} == "" ]];
    then
	effectiveGenomeSize=$(/home/michael.szalay/anaconda3/envs/subread/bin/faCount /zssd/scratch/DBs/*/UCSC/${assembly}/Sequence/WholeGenomeFasta/genome.fa | grep total | awk '{print $2}')
	echo ${assembly} ${effectiveGenomeSize} >> /home/common_pipelines/utils/list_of_effective_genomic_size.tab
    fi
    echo "Effective genome size ${effectiveGenomeSize} of ${assembly} from https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html"
    
    mkdir -p ${replicate}
    cd ${replicate}

    #echo $assayName $targetName $cellName $condition $replicateName $laneName $assembly $tag
    inBam=$(ls -1 *sorted_dedup.bam)
    if [[ ${replicateName} == "merge" ]];
    then
	inBam=$(ls -1 *.bam)
    fi
    
    outBw=${replicate}_normRPKM.bw
    if [[ -e ${outBw} ]];
    then
        echo "${outBw} already present in ${replicate}!"
    else
	echo "3 - Create .bigWig for ${replicate} using DeepTools"
	/home/michael.szalay/anaconda3/envs/deeptools/bin/bamCoverage --outFileFormat bigwig --ignoreDuplicates --normalizeUsing RPKM --numberOfProcessors max/2 --effectiveGenomeSize ${effectiveGenomeSize} -v --bam ${inBam} --extendReads 0 --outFileName ${outBw} 
	# DeepTools2 Version 3.5.0
	# https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html
	# --binSize             Size of the bins, in bases, for the output of the bigwig/bedgraph file. (Using default: 50)
	# --bam                 Using all BAM files in the folder    
	# --outFileFormat       Output file type. Either “bigwig” or “bedgraph”.
	### --samFlagExclude 128 second in pair (0x80)* https://broadinstitute.github.io/picard/explain-flags.html -> Ask for explanations to Giorgio, and i case use it in the analysis
	# --outFileName         Output file name.
	# --numberOfProcessors  Number of processors to use. Type “max/2” to use half the maximum number of processors or “max” to use all available processors.
	# --ignoreDuplicates    If set, reads that have the same orientation and start position will be considered only once. If reads are paired, the mate’s position also has to coincide to ignore a read.
	# --normalizeUsing      Use one of the entered methods to normalize the number of reads per bin. RPKM = Reads Per Kilobase per Million mapped reads.
	# --effectiveGenomeSize The effective genome size is the portion of the genome that is mappable. Large fractions of the genome are stretches of NNNN that should be discarded. Also, if repetitive regions were not included in the mapping of reads, the effective genome size needs to be adjusted accordingly. The value for mm10/GRCm38 was taken from here https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html.
	# -v                    Verbose mode    
    fi
    
    outBw=${replicate}_normNONE.bw
    if [[ -e ${outBw} ]];
    then
        echo "${outBw} already present in ${replicate}!"
    else
	echo "3 - Create .bigWig for ${replicate} using DeepTools"
	/home/michael.szalay/anaconda3/envs/deeptools/bin/bamCoverage --outFileFormat bigwig --ignoreDuplicates --normalizeUsing None --numberOfProcessors max/2 --effectiveGenomeSize ${effectiveGenomeSize} -v --bam ${inBam}  --extendReads 0 --outFileName ${outBw} 
	# DeepTools2 Version 3.5.0
	# https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html
	# --binSize             Size of the bins, in bases, for the output of the bigwig/bedgraph file. (Using default: 50)
	# --bam                 Using all BAM files in the folder    
	# --outFileFormat       Output file type. Either “bigwig” or “bedgraph”.
	### --samFlagExclude 128 second in pair (0x80)* https://broadinstitute.github.io/picard/explain-flags.html -> Ask for explanations to Giorgio, and i case use it in the analysis
	# --outFileName         Output file name.
	# --numberOfProcessors  Number of processors to use. Type “max/2” to use half the maximum number of processors or “max” to use all available processors.
	# --ignoreDuplicates    If set, reads that have the same orientation and start position will be considered only once. If reads are paired, the mate’s position also has to coincide to ignore a read.
	# --normalizeUsing      Use one of the entered methods to normalize the number of reads per bin. RPKM = Reads Per Kilobase per Million mapped reads.
	# --effectiveGenomeSize The effective genome size is the portion of the genome that is mappable. Large fractions of the genome are stretches of NNNN that should be discarded. Also, if repetitive regions were not included in the mapping of reads, the effective genome size needs to be adjusted accordingly. The value for mm10/GRCm38 was taken from here https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html.
	# -v                    Verbose mode    
    fi    
    
    cd .. # Exit ${replicate}
    echo
    
done # Close cycle over $replicates
