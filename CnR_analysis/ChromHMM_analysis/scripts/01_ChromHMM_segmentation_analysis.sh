assembly=dm6
conditions="NOnub"
chromosomes="chr2L  chr2R  chr3L  chr3R  chr4  chrM  chrX  chrY"
chromSizes="/home/Programs/chrom_sizes_dm6_higlass.txt"

bamDir=bamFiles
bamFlag=0  # 1 is FALSE 0 is TRUE

peakDir=peakFiles
peakFlag=1 # 1 is FALSE 0 is TRUE

for condition in ${conditions} ;
do
    for binning in 200 ; # in bp
    do
	outDir=analysis_${condition}_noH3K9me3_at_${binning}bp
	if [[ -d ${outDir} ]];
	then

	    continue
	fi
	mkdir -p ${outDir}
	    
	#    (
	# 1 - Binarize the .bams
	binarizedBams=binarizedBams_${condition}_at_${binning}bp
	signalDir=signal_${condition}_at_${binning}bp    
	mkdir -p ${binarizedBams}
	mkdir -p ${signalDir}
	
	if [[ ${bamFlag} -eq 0 ]];
	then
	    echo "Generate binarized files from .bam files"       
	    java -mx40G -jar /home/Programs/ChromHMM/ChromHMM.jar BinarizeBam -b ${binning} -t ${signalDir} -p 0.0001 ./scripts/chrom_sizes.txt ./${bamDir}/ cellmarkfiletable.tsv ${binarizedBams} 
	    #-b binsize – The number of base pairs in a bin determining the resolution of the model learning and segmentation. By default this parameter value is set to 200 base pairs.
	    #-p poissonthreshold – This option specifies the tail probability of the poisson distribution that the binarization threshold should correspond to. The default value of this parameter
	    #   is 0.0001.
	    #-t outputsignaldir – If specified signal files will be generated and outputted to the given directory. The files will be named CELL_CHROM_signal.txt. These files could later be
	    #   binarized directly at different thresholds with the BinarizeSignal command. By default no output signal is written.
	fi

	if [[ ${peakFlag} -eq 0 ]];
	then
	    #WD	chrX
	    #H3K27ac	H3K27me3	H3K9me3	Pc
	    #0	0	0	0	    
	    echo "Generate binarized files from peak files"       
	    echo $condition
	    for chrom in ${chromosomes} ;			     
	    do
		binarizedFile=${condition}_${chrom}_binary.txt
		rm -fvr ${binarizedFile}
		chromSize=$(grep -w ${chrom} ${chromSizes} | awk '{print $2}')
		echo ${chromSize}
		
		for file in $(ls -1 ${peakDir}/*${condition}* | grep merge | grep -v PH_ED_PH);
		do			
		    echo $file $chrom

		    target=$(echo ${file##*/} | sed "s/_/ /g" | awk '{print $2"_"$5}')
		    echo $target > ${condition}_${target}_${chrom}_binary.txt
		    
		    awk -v chromSize=${chromSize} -v chrom=${chrom} -v binning=$binning '{if($1==chrom){print $1,int($2/binning)*binning,int($3/binning+1)*binning}}' ${file} | awk -v binning=$binning -v chromSize=${chromSize} 'BEGIN{size=int(chromSize/binning)+1}{s1=int($2/binning); s2=int($3/binning); for(i=s1;i<s2;i++){h[i]=1}}END{for(i=0;i<size;i++){if(h[i]==""){h[i]=0}; print h[i]}}' >> ${condition}_${target}_${chrom}_binary.txt 
		    head ${condition}_${target}_${chrom}_binary.txt
		    
		    if [[ ! -e ${binarizedFile} ]];
		    then
			cp ${condition}_${target}_${chrom}_binary.txt ${binarizedFile}
		    else			
			paste ${binarizedFile} ${condition}_${target}_${chrom}_binary.txt > _tmp_${binarizedFile}
			#paste ${binarizedFile} ${condition}_${chrom}_binary.txt | awk '{for(i=1;i<NF;i++){printf("%s\t",$i)}; printf("%s\n",$NF)}' >> _tmp_${binarizedFile}
			mv _tmp_${binarizedFile} ${binarizedFile}
		    fi
		    head ${binarizedFile}
		done
		cp ${binarizedFile} _tmp_${binarizedFile}
		echo ${condition} ${chrom} | awk '{for(i=1;i<NF;i++){printf("%s\t",$i)}; printf("%s\n",$NF)}' > ${binarizedFile}
		awk '{for(i=1;i<NF;i++){printf("%s\t",$i)}; printf("%s\n",$NF)}' _tmp_${binarizedFile} >> ${binarizedFile}
		head ${binarizedFile}
		mv -v ${binarizedFile} ${binarizedBams}
		rm -fvr _tmp_${binarizedFile} ${condition}_*_${chrom}_binary.txt
	    done
	fi

	# Best guess is 4: active (K27ac); Polycomb (K27me3 and Pc); Repressive (K9me3); Null (None)
	for nstates in 4 ;		       
	do
	    
	    modelDir=model_at_${binning}bp_with_${nstates}_states
	    
	    # 2 - Learn model from the binarized data
	    java -mx40G -jar /home/Programs/ChromHMM/ChromHMM.jar LearnModel -b ${binning} ${binarizedBams} ${modelDir} ${nstates} ${assembly} &> output_at_${binning}bp_with_${nstates}_states.log 
	done
	
	mv -v m*_at_${binning}bp* s*_at_${binning}bp* b*_at_${binning}bp* o*_at_${binning}bp* ${outDir}
	
	#    ) &
    done # Close cycle over binning
done # Close cycle over condition
