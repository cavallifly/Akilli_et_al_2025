#!/bin/bash
echo "# Remember to load ChIPseq_pipeline environment: conda activate rnaseq_pipeline"

assembly=$1

ls -1 /zssd/scratch/DBs/*/*/ | grep -v ":" | awk '{if(NF!=0){print $0}}' > /home/common_pipelines/_available_assemblies
checkAssembly=$(awk -v a=${assembly} 'BEGIN{c=0}{if(a==$1){c=1}}END{print c}' /home/common_pipelines/_available_assemblies)
echo
if [[ $checkAssembly != 1 ]];
then
    echo "The requested assembly ${assembly} is not in the database."
    
    echo "Possible errors:"
    echo "1. You didn't provide the assembly as the first in-line argument."
    echo "2. The target assembly is not yet in the database."
    echo "Please download it from"
    echo "https://support.illumina.com/sequencing/sequencing_software/igenome.html"
    echo "in /zssd/scratch/DBs/"
    ls -1 /zssd/scratch/DBs/ | grep -v "tar\|txt"
    exit
else
    echo "The requested assembly is available at: $(ls -1 /zssd/scratch/DBs/*/*/* | grep $assembly | sed "s/://g")."
fi
echo
echo "$(bowtie2 --version 2>&1 | awk '{if(NF!=0)printf("%s ",$0)}END{printf(" is used for mapping")}')"
indexDir=$(ls -1 /zssd/scratch/DBs/*/*/${assembly}/Sequence/Bowtie2Index/*.fa | sed "s,.fa,,g")
echo "Baseline of the bowtie2 index used"
echo "${indexDir}"
echo

#samples=$(ls -1 raw_data/*/ | grep _1 | grep -v _2.fq | sed -e "s,raw_data/, ,g" -e "s,.gz,,g" -e "s,_1.fq,,g" -e "s,_1.fastq,,g")
samples=$(ls -1 raw_data | grep cnr)
echo "Samples to analyse ${samples}"

for sample in ${samples}
do

    assayName=$(echo ${sample} | sed "s,_, ,g" | awk '{print $1}')
    targetName=$(echo ${sample} | sed "s,_, ,g" | awk '{print $2}')
    cellName=$(echo ${sample} | sed "s,_, ,g" | awk '{print $3}')
    condition=$(echo ${sample} | sed "s,_, ,g" | awk '{print $4}')
    replicateName=$(echo ${sample} | sed "s,_, ,g" | awk '{print $5}')
    laneName=$(echo ${sample} | sed "s,_, ,g" | awk '{print $6}')
    #assembly=$(echo ${sample} | sed "s,_, ,g" | awk '{print $7}')
    tag=$(echo ${sample} | sed "s,_, ,g" | awk '{print $8}')

    replicate=${assayName}_${targetName}_${cellName}_${condition}_${replicateName}_${laneName}_${assembly}_${tag}

    if [[ -d ${replicate} ]];
    then
	echo "${replicate} already mapped. Going to the next sample."
	continue
    fi
    echo "Analysing $replicate"
    
    mkdir -p ${replicate}
    cd $replicate

    for fastqFile1 in $(ls -1 ../raw_data/${sample}/*_1.* 2> /dev/null) ;
    do
	# Mapping	    
	fastqFileName2=$(echo ${fastqFile1} | sed -e "s/_1\./_2\./g")
        fastqFile2=$(ls -1 $fastqFileName2 2> /dev/null | awk 'BEGIN{v="NA"}{v=$1}END{print v}')
        if [[ "${fastqFile2}" == "${fastqFile1}" ]];
        then
            if [[ $fastqFile2 != "NA" ]];
            then
                echo "The .fastq files for read1 and read2 are the same! Please, check these files:"
              	echo "Fastq file for read 1 ${fastqFile1}"
                echo "Fastq file for read 2 ${fastqFile2}"
		exit
	    fi
        fi
	echo "Fastq file for read 1 ${fastqFile1}"
	echo "Fastq file for read 2 ${fastqFile2}"

	outSam=${replicate}.sam
	outBam=${replicate}.bam
	
	if [[ ${fastqFile2} == "NA" ]];
	then
	    echo "Aligning single-end sample"
	    bowtie2 -x ${indexDir} -U ${fastqFile1} -S ${outSam} -p 8 --local --very-sensitive-local --no-unal --no-discordant --phred33 -I 10 -X 700 &>> ${replicate}_bowtie2.log
	    # Version 2.2.4
	    # (default) look for multiple alignments, report best, with MAPQ
	    # https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
	    # -x Bowtie2 index
	    # -U Read 1
	    # -S Output sam file
	    # -p number of alignment threads to launch
	    # --local When this option is specified, Bowtie 2 performs local read alignment. In this mode, Bowtie 2 might "trim" or "clip" some read characters
	    #         from one or both ends of the alignment if doing so maximizes the alignment score.
	    # --very-sensitive-local Same as: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50     
	    # --no-unal Suppress SAM records for reads that failed to align.
	    # --no-discordant By default, bowtie2 looks for discordant alignments if it cannot find any concordant alignments. A discordant alignment is an
	    #                 alignment where both mates align uniquely, but that does not satisfy the paired-end constraints (--fr/--rf/--ff, -I, -X).
	    #                 This option disables that behavior.
	    # --phred33 Input qualities are ASCII chars equal to the Phred quality plus 33. This is also called the "Phred+33" encoding, which is used by the very latest Illumina pipelines
	    # -I 10 -X 700 The minimum and maximum fragment length for valid paired-end alignments. E.g. if -I 60 is specified and a paired-end alignment consists
	    #              of two 20-bp alignments in the appropriate orientation with a 20-bp gap between them, that alignment is considered valid (as long as -X is also satisfied).
	    #              A 19-bp gap would not be valid in that case. If trimming options -3 or -5 are also used, the -I constraint is applied with respect to the untrimmed mates.
	    #              The larger the difference between -I and -X, the slower Bowtie 2 will run. This is because larger differences between -I and -X require that Bowtie 2 scan a larger
	    #              window to determine if a concordant alignment exists. For typical fragment length ranges (200 to 400 nucleotides), Bowtie 2 is very efficient.
	else
	    echo "Aligning paired-end sample"	    
	    bowtie2 -x ${indexDir} -1 ${fastqFile1} -2 ${fastqFile2} -S ${outSam} -p 8 --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 &>> ${replicate}_bowtie2.log
	    # Version 2.2.4
	    # (default) look for multiple alignments, report best, with MAPQ
	    # -x Bowtie2 index
	    # -1 Read 1
	    # -2 Read 2
	    # -S Output sam file
	    # -p number of alignment threads to launch
	    # --no-mixed If Bowtie 2 cannot find a paired-end alignment for a pair, by default it will go on to look for unpaired alignments for the constituent mates.
	    #            This is called "mixed mode." To disable mixed mode, set the --no-mixed option. Bowtie 2 runs a little faster in --no-mixed mode, but will only
	    #            consider alignment status of pairs per se, not individual mates.
	    # --no-discordant : suppress discordant alignments for paired reads
	fi

	echo "2 - Filtering: Discard the reads in ${outSam} with map quality (mapq) < 30"
	samtools view -@ 8 -h -F 4 -b -q 30 ${outSam} > ${outBam}
	# Version 1.9
	# -b     Output in the bam format
	# -F     4 Remove ummaped reads see https://broadinstitute.github.io/picard/explain-flags.html
	# -h     Include the header in the output
	# -q     Skip alignments with MAPQ smaller than INT

	echo "3 - Sorting ${outBam} using sambamba"
	sambamba sort --tmpdir /zssd/scratch/tmp/ --show-progress --nthreads 12 ${outBam}
	
	
    done # Close cycle over $fastqFile1

    cd .. # Exit $replicate dir
    
done # Close cycle over $sample
