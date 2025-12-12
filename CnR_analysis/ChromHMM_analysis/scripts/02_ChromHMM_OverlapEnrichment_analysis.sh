assembly=dm6
bamdir=bamFiles
condition=NOnub
binning=8000 # nt resolution of the ChromHMM segmentation
nStates=4 # Number of emission used to segment the genome

TopDomRes=${binning}

inDomainFile=$(ls -1 TopDom_domains_microc_contacts_${condition}_merge_all_dm6_NA_window_*_at_${TopDomRes}bp.bed)
#outDomainFile=TopDom_domains_microc_contacts_${condition}_merge_all_dm6_BeS_window_4_at_4000bp_with_state.bed
outDomainFile=$(echo $inDomainFile | sed "s/.bed/_with_state.bed/g")
outDir=OverlapEnrichment_${nStates}EmissionsModel_at_${binning}bp_TopDom_domains_at_${TopDomRes}bp_for_${condition}

# 3 States assignment
# 1 : PcG
# 2 : Null
# 3 : Active

# 4 States assignment
# 1 : PcGUb
# 2 : PcG
# 3 : Null
# 4 : Active

nDomains=$(wc -l ${inDomainFile} | awk '{print $1}')
echo $nDomains

if [[ -d ${outDir} ]];
then
    continue
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
    
    java -mx40G -jar /home/Programs/ChromHMM/ChromHMM.jar OverlapEnrichment ../analysis_PH18_at_${binning}bp/model_at_${binning}bp_with_${nStates}_states/PH18_${nStates}_segments.bed _tmp.bed OverlapEnrichment_${nStates}EmissionsModel_${chrom}_${start}_${end}
    rm _tmp.bed
      
    #echo $nDomain $chrom $start $end $(grep -v "Base\|Genome" OverlapEnrichment_${nStates}EmissionsModel_${chrom}_${start}_${end}.txt | sort -k 3,3n | tail -1 | awk '{if($1==1){state="PcGUb"}; if($1==2){state="PcG"}; if($1==3){state="Null"}; if($1==4){state="Active"}; print $1,state}') >> ${outDomainFile}
    echo $nDomain $chrom $start $end $(grep -v "Base\|Genome" OverlapEnrichment_${nStates}EmissionsModel_${chrom}_${start}_${end}.txt | sort -k 3,3n | tail -1 | awk '{if($1==1){state="PcG"}; if($1==2){state="Null"}; if($1==3){state="Active"}; print $1,state}') >> ${outDomainFile}    
    
done

cd .. # Exit ${outDir}

