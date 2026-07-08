TopDomStatesFile=scripts_clean/TopDom_domains_hic_WD_NOnub_merge_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp_with_state.bed
GeneFile=dmel-all-r6.36_sed_for_rnaseq.gtf
vsdFile=vsd_values_flybase_only_NOnub_SUMOnub.txt

outFile=gene_vsd_values_and_TAD.tab

nTADs=$(wc -l $TopDomStatesFile | awk '{print $1}')
for nTAD in $(seq 1 1 ${nTADs});
do
    TAD=$(  awk -v n=${nTAD} '{if(NR==n){print $4}}' $TopDomStatesFile)
    chrom=$(awk -v n=${nTAD} '{if(NR==n){print $1}}' $TopDomStatesFile | sed "s/chr//g")
    start=$(awk -v n=${nTAD} '{if(NR==n){print $2}}' $TopDomStatesFile)
    end=$(  awk -v n=${nTAD} '{if(NR==n){print $3}}' $TopDomStatesFile)

    nCheck=$(grep -w ${TAD} $outFile | wc -l | awk '{print $1}')
    if [[ $nCheck -lt 0 ]];
    then
	continue
    fi
    
    echo $TAD $chrom $start $end
    awk -v n=${nTAD} '{if(NR==n){print $0}}' $TopDomStatesFile
    
    grep -w gene ${GeneFile} | awk -v c=${chrom} -v s=${start} -v e=${end} -v T=${TAD} '{if($1==c){cnt=0;for(i=$4;i<=$5;i++){if(s<=i && i<=e){cnt++}};if(cnt>0){print $10,"chr"$1,$4,$5,T,cnt}}}' | sed -e "s/\"//g" -e "s/;//g" > _tmp_${TAD}

    for gene in $(awk '{print $1}' _tmp_${TAD});
    do
	grep $gene ${vsdFile} >> _tmp_vsd_${TAD}
    done
    paste _tmp_${TAD} _tmp_vsd_${TAD} | awk '{if($1==$7){for(i=1;i<NF;i++){if(i==7){continue}; printf("%s\t",$i)}; printf("%s\n",$NF)}}' >> ${outFile}
    rm -fr _tmp_${TAD} _tmp_vsd_${TAD}
done

awk '{if(h[$1]!=""){if(overlap[$1]>$6){next}}; h[$1]=$0; overlap[$1]=$6}' ${outFile}
