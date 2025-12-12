echo "Consider all the analysed samples"
# Sample name assayName_targetName_cellName_condition_replicateName_laneName_assembly_tag
# e.g., cnr_H3K27me3_WD_NOnub_Rep1_L1_dm6_NA

if [[ $1 == "" ]];
then
    echo "Please, give as unique parameter the name of the configuration table."
    echo "The configuration table has 5 columns:"
    echo "1 - Directory name of the sample"
    echo "2 - single-end/paired-end experiment flag"
    echo "3 - cutoff value of the q (false discovery rate) parameter (e.g. 0.01)"
    echo "4 - narrow or broad flag to do a narrowpeak or broadpeak calling analysis"
    echo "5 - Directory name of the control (optional)"     
    exit
fi
inTable=$1

nlines=$(wc -l ${inTable} | awk '{print $1+1}')
for nReplicate in $(seq 1 1 ${nlines});
do    
    # Get the parameters from the $inTable
    replicate=$(awk -v n=${nReplicate} '{if(NR==n){print $1}}' ${inTable})
    experimentFlag=$(awk -v n=${nReplicate} '{if(NR==n){print $2}}' ${inTable})
    qCutoff=$(awk -v n=${nReplicate} '{if(NR==n){print $3}}' ${inTable})
    peakFlag=$(awk -v n=${nReplicate} '{if(NR==n){print $4}}' ${inTable})
    controlDir=$(awk -v n=${nReplicate} '{if(NR==n){print $5}}' ${inTable})
    if [[ $replicate == "" ]];
    then
	continue
    fi
    
    assembly=$(echo ${replicate} | sed "s,_, ,g" | awk '{print $7}')

    #effectiveGenomeSize=$(faCount /zssd/scratch/DBs/*/UCSC/${assembly}/Sequence/WholeGenomeFasta/genome.fa | grep total | awk '{print $2}')
    effectiveGenomeSize=$(awk -v a=${assembly} '{if($1==a){print $2}}' /home/common_pipelines/utils/list_of_effective_genomic_size.tab)
    if [[ ${effectiveGenomeSize} == "" ]];
    then
	effectiveGenomeSize=$(/home/michael.szalay/anaconda3/envs/subread/bin/faCount /zssd/scratch/DBs/*/UCSC/${assembly}/Sequence/WholeGenomeFasta/genome.fa | grep total | awk '{print $2}')
	echo ${assembly} ${effectiveGenomeSize} >> /home/common_pipelines/utils/list_of_effective_genomic_size.tab
    fi
    echo "Effective genome size ${effectiveGenomeSize} of ${assembly} from https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html"
    
    echo ${replicate} ${experimentFlag} ${peakFlag} ${controlDir}

    cd ${replicate}
   
    inBam=$(ls -1 *sorted_dedup*.bam)
    controlBam=$(ls -1 ../${controlDir}/*.bam 2> /dev/null)      
    
    if [[ ${peakFlag} == "narrow" ]];
    then
	if [[ ! -e ${replicate}_peaks.narrowPeak ]];	   
	then

	    if [[ $experimentFlag == "single-end" ]];
	    then
		if [[ $controlDir == "" ]];
		then
		    echo "4 - Call the peaks for ${name} with macs3"
		    /home/michael.szalay/anaconda3/envs/macs/bin/macs3 callpeak -t ${inBam} -f BAM -g ${effectiveGenomeSize} -q ${qCutoff} -n ${replicate} --outdir . --tempdir /zssd/scratch/tmp/
		    # Version 3.0.0b1
		    # -t The IP data file (this is the only REQUIRED parameter for MACS). It is the only required parameter.
		    # -f Format of input file.
		    # -g Mappable genome size. For mouse mm10.
		    # -q q-value (minimum FDR) cutoff. We use default 0.01.
		    # -n The prefix string for output files.
		    # --outdir Name of the dir to store the output.
		    # --SPMR If True, MACS will SAVE signal per million reads for fragment pileup profiles. It won't interfere with computing pvalue/qvalue during peak calling, since internally MACS3 keeps
		    #  	 using the raw pileup and scaling factors between larger and smaller dataset to calculate statistics measurements. If you plan to use the signal output in bedGraph to call peaks
		    #        using bdgcmp and bdgpeakcall, you shouldn't use this option because you will end up with different results. However, this option is recommended for displaying normalized pileup
                    #        tracks across many datasets. Require -B to be set. Default: False
		    # --bdg  Whether or not to save extended fragment pileup, and local lambda tracks (two files) at every bp into a bedGraph file. DEFAULT: False
		    # --nomodel Whether or not to build the shifting model. If True, MACS will not build model. by default it means shifting size = 100, try to set extsize to change it.
		    #           It's highly recommended that while you have many datasets to process and you plan to compare different conditions, aka differential calling, 
		    #	    use both 'nomodel' and 'extsize' to make signal files from
                    #           different datasets comparable. DEFAULT: False
		    # --tempdir TEMPDIR     Optional directory to store temp files.
		    # Zhang Y, Liu T, Meyer CA, Eeckhoute J, Johnson DS, Bernstein BE, Nusbaum C, Myers RM, Brown M, Li W, Liu XS. (2008)
		    # Model-based Analysis of ChIP-Seq (MACS), Genome Biology, 2008;9(9):R137.
		else
		    echo "4 - Call the peaks for ${name} with macs3"
		    /home/michael.szalay/anaconda3/envs/macs/bin/macs3 callpeak -t ${inBam} -c ${controlBam} -f BAM -g ${effectiveGenomeSize} -q ${qCutoff} -n ${replicate} --outdir . --tempdir /zssd/scratch/tmp/
		    # Version 3.0.0b1
		    # -c The control or mock data file, e.g. ES_Input.bam
		fi
	    fi
	    
	    if [[ $experimentFlag == "paired-end" ]];
	    then
		if [[ $controlDir == "" ]];
		then
		    echo "4 - Call the peaks for ${name} with macs3"
		    /home/michael.szalay/anaconda3/envs/macs/bin/macs3 callpeak -t ${inBam} -f BAMPE -g ${effectiveGenomeSize} -q ${qCutoff} -n ${replicate} --outdir . --tempdir /zssd/scratch/tmp/
		    # Version 3.0.0b1
		else
		    echo "4 - Call the peaks for ${name} with macs3"
		    /home/michael.szalay/anaconda3/envs/macs/bin/macs3 callpeak -t ${inBam} -c ${controlBam} -f BAMPE -g ${effectiveGenomeSize} -q ${qCutoff} -n ${replicate} --outdir . --tempdir /zssd/scratch/tmp/
		    # Version 3.0.0b1
	    fi
	    fi
	fi
    fi
    
    if [[ ${peakFlag} == "broad" ]];
    then
	if [[ ! -e ${replicate}_peaks.broadPeak ]];	   
	then
	    
	    if [[ $experimentFlag == "single-end" ]];
	    then
		if [[ $controlDir == "" ]];
		then
		    echo "4 - Call the peaks for ${name} with macs3"
		    /home/michael.szalay/anaconda3/envs/macs/bin/macs3 callpeak -t ${inBam} -f BAM -g ${effectiveGenomeSize} -q ${qCutoff} -n ${replicate} --outdir . --broad --broad-cutoff 0.01 --tempdir /zssd/scratch/tmp/
		    # Version 3.0.0b1
		    # --broad If set, MACS will try to call broad peaks using the --broad-cutoff setting. Please tweak '--broad-cutoff' setting to control the peak calling behavior.
		    #         At the meantime, either -q or -p cutoff will be used to define regions with 'stronger enrichment' inside of broad peaks.
		    #         The maximum gap is expanded to 4 * MAXGAP (--max-gap parameter). As a result, MACS will output a 'gappedPeak' and a 'broadPeak' file instead of 'narrowPeak' file.
		    #         Note, a broad peak will be reported even if there is no 'stronger enrichment' inside. DEFAULT: False
                    # --broad-cutoff BROADCUTOFF
                    #         Cutoff for broad region. This option is not available unless --broad is set. If -p is set, this is a pvalue cutoff, otherwise, it's a qvalue cutoff. Please note that in broad
                    #         peakcalling mode, MACS3 uses this setting to control the overall peak calling behavior, then uses -q or -p setting to define regions inside broad region
		    #         as 'stronger' enrichment.
                    #         DEFAULT: 0.1
		else
		    echo "4 - Call the peaks for ${name} with macs3"
		    /home/michael.szalay/anaconda3/envs/macs/bin/macs3 callpeak -t ${inBam} -c ${controlBam} BAM -g ${effectiveGenomeSize} -q ${qCutoff} -n ${replicate} --outdir . --broad --broad-cutoff ${qCutoff} --tempdir /zssd/scratch/tmp/
		    # Version 3.0.0b1
		fi
	    fi
	    
	    if [[ $experimentFlag == "paired-end" ]];
	    then
		if [[ $controlDir == "" ]];
		then
		    echo "4 - Call the peaks for ${name} with macs3"
		    /home/michael.szalay/anaconda3/envs/macs/bin/macs3 callpeak -t ${inBam} -f BAMPE -g ${effectiveGenomeSize} -q ${qCutoff} -n ${replicate} --outdir . --broad --broad-cutoff 0.01 --tempdir /zssd/scratch/tmp/
		    # Version 3.0.0b1
		    # Zhang Y, Liu T, Meyer CA, Eeckhoute J, Johnson DS, Bernstein BE, Nusbaum C, Myers RM, Brown M, Li W, Liu XS. (2008) Model-based Analysis of ChIP-Seq (MACS), Genome Biology, 2008;9(9):R137.
		else
		    echo "4 - Call the peaks for ${name} with macs3"
		    /home/michael.szalay/anaconda3/envs/macs/bin/macs3 callpeak -t ${inBam} -c ${controlBam} -f BAMPE -g ${effectiveGenomeSize} -q ${qCutoff} -n ${replicate} --outdir . --broad --broad-cutoff ${qCutoff} --tempdir /zssd/scratch/tmp/
		    # Version 3.0.0b1
		    # Zhang Y, Liu T, Meyer CA, Eeckhoute J, Johnson DS, Bernstein BE, Nusbaum C, Myers RM, Brown M, Li W, Liu XS. (2008) Model-based Analysis of ChIP-Seq (MACS), Genome Biology, 2008;9(9):R137.		
		fi
	    fi
	fi   
    fi
    cd .. # Exit ${outDir}
    echo ""
    
done # Close cycle over $replicates
