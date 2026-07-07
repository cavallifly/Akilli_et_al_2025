inFile1=analysis_withNAH3K9me3/OverlapEnrichment_4EmissionsModel_WD_NOnub_TopDom_domains/TopDom_domains_hic_WD_NOnub_merge_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp_with_state_Original.bed
inFile2=OverlapEnrichment_4EmissionsModel_at_200bp_TopDom_domains_at_5000bp_for_NOnub/TopDom_domains_hic_WD_NOnub_merge_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp_with_state_Original.bed
wc -l $inFile1 $inFile2

inFile=conversion_TADstates_withNAH3K9me3_vs_withGSH3K9me3.txt
paste $inFile1 $inFile2 | awk '{if($2==$8 && $3==$9 && $4==$10){print $0}}' | awk '{print $6""++h1[$6],$NF""++h2[$NF]}' > $inFile
wc -l ${inFile}

inFile=conversion_TADstates_withNAH3K9me3_vs_withGSH3K9me3_for_overlap_analysis.txt
paste $inFile1 $inFile2 | awk '{if($2==$8 && $3==$9 && $4==$10){print $0}}' | awk '{print "TAD"NR,$6""++h1[$6],$NF""++h2[$NF]}' > $inFile
head $inFile
tail $inFile
echo

for category in Active1 PcG1 Null1 Het1 ;
do

    outFile=cluster_${category}_TADs.txt
    
    categoryName=$(echo $category | sed "s/1//g")

    awk '{print $1,$2}' ${inFile} | grep $categoryName | awk '{print $1}' > ${outFile}
    wc -l ${outFile}

done

for category in Active2 PcG2 Null2 Het2;
do

    outFile=cluster_${category}_TADs.txt
    
    categoryName=$(echo $category | sed "s/2//g")

    awk '{print $1,$3}' ${inFile} | grep $categoryName | awk '{print $1}' > ${outFile}
    wc -l ${outFile}

done

Rscript scripts/03_compute_overlaps_between_TADpartitions.R &>> 03_compute_overlaps_between_TADpartitions.log
