inFileTADstates=scripts_clean/TopDom_domains_hic_WD_NOnub_merge_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp_with_state.bed 

analysisTag=withGSH3K9me3

inFile=scores_trans1Dinterval_allChrom_all_domains_vs_all_domains_${analysisTag}.tab
if [[ ! -e ${inFile} ]];
then
    cat scores_trans1Dinterval_chr*_all_domains_vs_all_domains_${analysisTag}.tab > ${inFile}
fi

for inFile in $(ls -1 scores_trans1Dinterval_c*_all_domains_vs_all_domains_${analysisTag}.tab);
do     
    outFile=$(echo $inFile | sed -e "s/scores_/avgScores_/g" -e "s/_${analysisTag}//g")
    echo $outFile    
    if [[ -e $outFile ]];
    then
	continue
    fi
    
    echo "inFile ${inFile}"
    echo "outFile ${outFile}"    
    touch ${outFile}
    
    head $inFile
    
    cat <( grep PcG ${inFileTADstates} ) <( grep Active ${inFileTADstates} ) <( grep Het ${inFileTADstates} ) <( grep Null ${inFileTADstates} ) $inFile | awk '{if(NF==4){chrom[$4]=$1; start[$4]=$2; end[$4]=$3; n[$4]++; next}; feature1=chrom[$7]";"start[$7]";"end[$7]";"$7";"$10; feature2=chrom[$8]";"start[$8]";"end[$8]";"$8";"$10; if(n[$7]>n[$8]){feature1=chrom[$8]";"start[$8]";"end[$8]";"$8";"$10; feature2=chrom[$7]";"start[$7]";"end[$7]";"$7";"$10}; sumScore[feature1":"feature2]+=$9; cntScore[feature1":"feature2]++}END{for(i in sumScore) print i,sumScore[i]/cntScore[i]}' | sed -e "s/:/ /g" -e "s/;/ /g" | awk '{print $1,$2,$3,$6,$7,$8,$4,$9,$NF,$5}' | awk '{for(i=1;i<NF;i++){printf("%s\t",$i)}; printf("%s\n",$NF)}' > ${outFile} &

done # Close cycle over $inFile
wait

#cat avgScores_trans1Dinterval_c*_all_domains_vs_all_domains.tab > avgScores_trans1Dinterval_all_domains_vs_all_domains.tab
