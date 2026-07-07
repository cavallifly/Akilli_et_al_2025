assembly=dm6
bamdir=bamFiles
condition=NOnub
binning=200 # nt resolution of the ChromHMM segmentation <- Used in the first submission
#binning=5000 # nt resolution of the ChromHMM segmentation 

analysisTag=withGSH3K9me3 ; nStates=4 # Number of emission used to segment the genome <- Used in the first submission
#analysisTag=withH3K9me3 ; nStates=4 # Number of emission used to segment the genome <- Used in the first submission
#analysisTag=noH3K9me3 ; nStates=3 # Number of emission used to segment the genome

TopDomRes=5000

inDomainFile=$(ls -1 TopDom_domains_hic_WD_NOnub_merge_dm6_NA_window_5_at_${TopDomRes}bp_filtered_refined_at_1000bp.bed)
outDomainFile=$(echo $inDomainFile | sed "s/.bed/_with_state.bed/g")
outDir=OverlapEnrichment_${nStates}EmissionsModel_at_${binning}bp_TopDom_domains_at_${TopDomRes}bp_for_${condition}


# 3 States assignment
# 1 : Active
# 2 : Null
# 3 : PcG

# 4 States assignment
# 1 : Active
# 2 : Null
# 3 : Het
# 4 : PcG

nDomains=$(wc -l ${inDomainFile} | grep -v chrom | awk '{print $1}')
echo $nDomains

if [[ -d ${outDir} ]];
then
    exit
fi

mkdir -p ${outDir}
cd ${outDir}

rm -fvr ${outDomainFile}
for nDomain in $(seq 1 1 ${nDomains});
do

    awk -v nD=${nDomain} '{if(NR==nD){print $0}}' ../${inDomainFile} > _tmp.bed

    chrom=$(awk '{print $1}' _tmp.bed)
    start=$(awk '{print $2}' _tmp.bed)
    end=$(awk '{print $3}' _tmp.bed)
    if [[ ${chrom} == "chrom" ]];
    then
	echo "Wrong domain $chrom $start $end"
	continue
    fi

    outFileName=OverlapEnrichment_${nStates}EmissionsModel_${chrom}_${start}_${end}
    outFile=${outFileName}.txt
    if [[ ! -e ${outFile} ]];
    then
	touch ${outFile}
	#java -mx40G -jar /home/Programs/ChromHMM/ChromHMM.jar OverlapEnrichment ../analysis_${condition}_at_${binning}bp/model_at_${binning}bp_with_${nStates}_states/${condition}_${nStates}_segments.bed _tmp.bed OverlapEnrichment_${nStates}EmissionsModel_${chrom}_${start}_${end}
	java -mx40G -jar /home/Programs/ChromHMM/ChromHMM.jar OverlapEnrichment ../analysis_${condition}_${analysisTag}_at_${binning}bp/model_at_${binning}bp_with_${nStates}_states/${condition}_${nStates}_segments.bed _tmp.bed ${outFileName}
    fi
    rm _tmp.bed

    # 4 states
    echo $nDomain $chrom $start $end $(grep -v "Base\|Genome" OverlapEnrichment_${nStates}EmissionsModel_${chrom}_${start}_${end}.txt | sort -k 3,3n | tail -1 | awk -v state1="Active" -v state2="Null" -v state3="Het" -v state4="PcG" '{if($1==1){state=state1}; if($1==2){state=state2}; if($1==3){state=state3}; if($1==4){state=state4}; print $1,state}') >> ${outDomainFile}

    # 3 states
    #echo $nDomain $chrom $start $end $(grep -v "Base\|Genome" ${outFile} | sort -k 3,3n | tail -1 | awk '{if($1==3){state="PcG"}; if($1==2){state="Null"}; if($1==1){state="Active"}; print $1,state}') >> ${outDomainFile}    

done

cd .. # Exit ${outDir}

