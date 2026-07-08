states="PcG Active Het Null"
conditions="SUMOnub NOnub"

analysisTag=withGSH3K9me3

inFile=avgScores_trans1Dinterval_allChrom_all_domains_vs_all_domains.tab

for condition in $conditions ;
do

    outFile=avg_scores_trans1Dinterval_all_domains_vs_all_domains_${condition}_${analysisTag}.tab
    if [[ ! -e ${outFile} ]];    
    then
	echo $states | awk '{for(i=1;i<NF;i++){printf("%s\t",$i)}; printf("%s\n",$NF)}' > ${outFile}
	for state1 in $states ;
	do
	    echo $state1 | awk '{printf("%s",$1)}' >> ${outFile}
	    for state2 in $states ;
	    do
		echo $state1 $state2
		
		grep ${condition} ${inFile} | grep ${state1} | grep ${state2} | awk -v s1=$state1 -v s2=$state2 'BEGIN{l1=length(s1); l2=length(s2)}{if((s1==substr($7,1,l1) && s2==substr($8,1,l2)) || (s1==substr($8,1,l1) && s2==substr($7,1,l2))){s+=$9; cnt++}}END{printf("\t%f",s/cnt)}' >> ${outFile}
	    done
	    echo $state1 | awk '{printf("\n")}' >> ${outFile}
	done
	echo $state1 | awk '{printf("\n")}' >> ${outFile}
    fi

    cp ${outFile} _tmp
    head ${outFile} _tmp
    
    cbmin=-6.5
    cbmax=6.5

    sed -e "s/XXXminXXX/${cbmin}/g" -e "s/XXXmaxXXX/${cbmax}/g" scripts_clean/03_get_chromatinStates_heatmap.R > _tmp.R
    conda run -n hicrep Rscript _tmp.R
    rm -fvr _tmp Rplots.pdf _tmp.R
    
    mv plots.pdf ${outFile%.tab}.pdf

    if [[ $condition == "NOnub" ]];
    then
	inFileNOnub=avg_scores_trans1Dinterval_all_domains_vs_all_domains_NOnub_${analysisTag}.tab
	inFileSUMOnub=avg_scores_trans1Dinterval_all_domains_vs_all_domains_SUMOnub_${analysisTag}.tab
	outFileDiff=avg_scores_trans1Dinterval_all_domains_vs_all_domains_SUMOnubMinusNOnub_${analysisTag}.tab
	
	paste ${inFileSUMOnub} ${inFileNOnub} | awk '{if(NR==1){printf("%s\t%s\t%s\t%s\n",$1,$2,$3,$4)}else{if(NF>0)printf("%-10s\t%-10s\t%-10s\t%-10s\t%-10s\n",$1,$2-$7,$3-$8,$4-$9,$5-$10)};if(NF==0){print $0}}' > ${outFileDiff}

	cp ${outFileDiff} _tmp
	head ${outFileDiff} _tmp

	cbmin=-2.5
	cbmax=2.5
	
	sed -e "s/XXXminXXX/${cbmin}/g" -e "s/XXXmaxXXX/${cbmax}/g" scripts_clean/03_get_chromatinStates_heatmap.R > _tmp.R
	conda run -n hicrep Rscript _tmp.R
	rm -fvr _tmp Rplots.pdf _tmp.R
	
	mv plots.pdf ${outFileDiff%.tab}.pdf
	
    fi
    
done
rm -fvr _tmp
