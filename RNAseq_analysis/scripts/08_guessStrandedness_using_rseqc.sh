

#https://rseqc.sourceforge.net/#download-rseqc
#This program is used to “guess” how RNA-seq sequencing were configured, particulary how reads were
#stranded for strand-specific RNA-seq data, through comparing the “strandness of reads” with the “standness of transcripts”.
#The “strandness of reads” is determiend from alignment, and the “standness of transcripts” is determined from annotation.
#For non strand-specific RNA-seq data, “strandness of reads” and “standness of transcripts” are independent.
#For strand-specific RNA-seq data, “strandness of reads” is largely determined by “standness of transcripts”. See below 3 examples for details.
#You don’t need to know the RNA sequencing protocol before mapping your reads to the reference genome. Mapping your RNA-seq reads as if they were non-strand specific,
#this script can “guess” how RNA-seq reads were stranded.
#For pair-end RNA-seq, there are two different ways to strand reads (such as Illumina ScriptSeq protocol):
# In the first case the read1 is always mapped concordantly to the strand of the gene: If this happens for the majority of the reads, it is a case of "forward stranded" mapping, that is option -s 1 in featurCounts
#1++,1--,2+-,2-+
#read1 mapped to ‘+’ strand indicates parental gene on ‘+’ strand
#read1 mapped to ‘-‘ strand indicates parental gene on ‘-‘ strand
#read2 mapped to ‘+’ strand indicates parental gene on ‘-‘ strand
#read2 mapped to ‘-‘ strand indicates parental gene on ‘+’ strand

# In the second case the read1 is always mapped opposite to the strand of the gene: If this happens for the majority of the reads, it is a case of "reverse stranded" mapping, that is option -s 2 in featurCounts
#1+-,1-+,2++,2--
#read1 mapped to ‘+’ strand indicates parental gene on ‘-‘ strand
#read1 mapped to ‘-‘ strand indicates parental gene on ‘+’ strand
#read2 mapped to ‘+’ strand indicates parental gene on ‘+’ strand
#read2 mapped to ‘-‘ strand indicates parental gene on ‘-‘ strand

#Options:
#--version
#show program’s version number and exit
#-h, --help
#show this help message and exit
#-i INPUT_FILE, --input-file=INPUT_FILE
#Input alignment file in SAM or BAM format
#-r REFGENE_BED, --refgene=REFGENE_BED
#Reference gene model in bed fomat.
#-s SAMPLE_SIZE, --sample-size=SAMPLE_SIZE
#Number of reads sampled from SAM/BAM file. default=200000
#-q MAP_QUAL, --mapq=MAP_QUAL
#Minimum mapping quality (phred scaled) for an alignment to be considered as “uniquely mapped”. default=30

assembly=$1

bed12File=$(ls /zssd/scratch/DBs/*/*/${assembly}/Annotation/Genes/genes.bed12)

for bamFile in $(ls -1 ./*/*bam | head -1);
do
    echo $bamFile

    /home/michael.szalay/anaconda3/envs/rseqc/bin/infer_experiment.py -r ${bed12File} -i ${bamFile}
    
done
