#!/bin/bash
#echo "# Remember to load subread environment: conda activate subread"

assembly=$1

ls -1 /zssd/scratch/DBs/*/*/ | grep -v ":" | awk '{if(NF!=0){print $0}}' > /home/common_pipelines/hic/_available_assemblies
checkAssembly=$(awk -v a=${assembly} 'BEGIN{c=0}{if(a==$1){c=1}}END{print c}' /home/common_pipelines/hic/_available_assemblies)
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
echo "$(conda run -n subread subread-align -v 2>&1 | awk '{if(NF!=0) print $0}' | awk '{if(NR==1)printf("%s is used for mapping",$0)}')"
indexDir=$(ls -1 /zssd/scratch/DBs/*/*/${assembly}/Sequence/subreadIndex/subread_${assembly}* | grep log | sed "s,.log,,g")
echo "Baseline of the subread index used"
echo "${indexDir}"
echo


#samples=$(ls -1 raw_data/*/ | grep _1 | grep -v _2.fq | sed -e "s,raw_data/, ,g" -e "s,.gz,,g" -e "s,_1.fq,,g")
samples=$(ls -1 raw_data | grep rnaseq | sed -e "s,/, ,g" | awk '{print $NF}')
echo $samples

for sample in ${samples}
do
    if [[ -d ${sample} ]];
    then
	echo "${sample} already mapped. Going to the next sample."
	continue
    fi
    echo "Analysing $sample"
    
    mkdir -p ${sample}
    cd $sample

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

	outFile=${sample}.bam

	if [[ ! ${fastqFile2} == "" ]];
	then
	    echo "Mapping data as paired-end experiment"
	    conda run -n subread subread-align -i ${indexDir} -r ${fastqFile1} -R ${fastqFile2} -t 0 -o ${outFile} -T 10  --sortReadsByCoordinates
	    #Version 2.0.6
	    #Usage:
	    #./subread-align [options] -i <index_name> -r <input> -t <type> -o <output>
	    ## Mandatory arguments:
	    #  -i <string>       Base name of the index.
	    #  -r <string>       Name of an input read file. If paired-end, this should be
	    #                    the first read file (typically containing "R1"in the file
	    #                    name) and the second should be provided via "-R".
	    #                    Acceptable formats include gzipped FASTQ, FASTQ, gzipped
	    #                    FASTA and FASTA.
	    #                    These formats are identified automatically.   
	    #  -t <int>          Type of input sequencing data. Its values include
	    #                      0: RNA-seq data
	    #                      1: genomic DNA-seq data.
	    # number of CPU threads
	    #  -T <int>          Number of CPU threads used, 1 by default.
	    #  -o <string>       Name of an output file. By default, the output is in BAM
	    #                    format. Omitting this option makes the output be written to
	    #                    STDOUT.
	    #  -R <string>       Name of the second read file in paired-end data (typically
	    #                    containing "R2" the file name).
	    # read order
	    #  --sortReadsByCoordinates Output location-sorted reads. This option is
	    #                    applicable for BAM output only. A BAI index file is also
	    #                    generated for each BAM file so the BAM files can be directly
	    #                    loaded into a genome browser.
	else
	    echo "Mapping data as single-end experiment"
	    conda run -n subread subread-align -i ${indexDir} -r ${fastqFile1} -t 0 -o ${outFile} -T 10  --sortReadsByCoordinates
	fi
	outBam=${sample}.bam
	echo "2 - Sorting and indexing ${outBam}"
	samtools sort -@ 8 -o sorted_${outBam} -O bam ${outBam}
	mv -v sorted_${outBam} ${outBam}
	# -o Write the final sorted output to FILE, rather than to standard output.
	# -@ Set number of sorting and compression threads. By default, operation is single-threaded.
	# -O Write the final output as sam, bam, or cram.
	samtools index -@ 8 ${outBam}
	# -@ Set number of threads.

	
    done # Close cycle over $fastqFile1

    cd .. # Exit $sample dir
    
done # Close cycle over $sample
