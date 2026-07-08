#!/bin/sh

#SBATCH --job-name geneLog2ratio
#SBATCH -n 1                    # Number of cores. For now 56 is the number max of core available
#SBATCH -t 4-00:00              # Runtime in D-HH:MM
#SBATCH --mem=16G
#SBATCH -o 05_get_gene_Log2ratio_and_TAD.out # File to which STDOUT will be written
#SBATCH -e 05_get_gene_Log2ratio_and_TAD.out # File to which STDERR will be written 


inTopDomDomains=./scripts_clean/TopDom_domains_hic_WD_NOnub_merge_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp_with_state.bed
inRNAseq=DEseq2_results_WD_SUMOnub_vs_WD_NOnub_flybase_only_NOnub_SUMOnub.tsv
inGenes=dmel-all-r6.36_used_for_rnaseq.gtf

outFileGene=gene_Log2ratio_and_TAD_basedOnGeneLength_preciselyInTAD.tsv
outFileTSS=gene_Log2ratio_and_TAD_basedOnGeneTSS_preciselyInTAD.tsv
if [[ -e ${outFileGene} ]];
then
    sort -k 3,3 -k 4,4n -k 5,5n ${outFileGene} > _tmp_${outFileGene} ; mv _tmp_${outFileGene} ${outFileGene}

    outFileTADs=avgLog2ratio_per_TAD_basedOnGeneLength_preciselyInTAD.tsv
    awk '{if($9==""){next}; chrom[$9]=$10; start[$9]=$11; end[$9]=$12; sum[$9]+=$8; sum2[$9]+=($8*$8); cnt[$9]++}END{for(i in sum){avg=sum[i]/cnt[i]; avg2=sum2[i]/cnt[i]; stddev=sqrt(avg2-avg*avg); print i,chrom[i],start[i],end[i],sum[i],sum2[i],cnt[i],avg,stddev}}' ${outFileGene} | sort -k 2,2 -k 3,3n -k 4,4n > ${outFileTADs}

    outFileTADs=maxLog2ratio_per_TAD_basedOnGeneLength_preciselyInTAD.tsv
    awk '{chrom[$9]=$10; start[$9]=$11; end[$9]=$12; if(max[$9]<$8){max[$9]=$8}; cnt[$9]++}END{for(i in max){print i,chrom[i],start[i],end[i],cnt[i],max[i]}}' ${outFileGene} | sort -k 2,2 -k 3,3n -k 4,4n > ${outFileTADs}    
fi
if [[ -e ${outFileTSS} ]];
then
    sort -k 3,3 -k 4,4n -k 5,5n ${outFileTSS} > _tmp_${outFileTSS} ; mv _tmp_${outFileTSS} ${outFileTSS}

    outFileTADs=avgLog2ratio_per_TAD_basedOnGeneTSS_preciselyInTAD.tsv
    awk '{if($9==""){next}; chrom[$9]=$10; start[$9]=$11; end[$9]=$12; sum[$9]+=$8; sum2[$9]+=($8*$8); cnt[$9]++}END{for(i in sum){avg=sum[i]/cnt[i]; avg2=sum2[i]/cnt[i]; stddev=sqrt(avg2-avg*avg); print i,chrom[i],start[i],end[i],sum[i],sum2[i],cnt[i],avg,stddev}}' ${outFileTSS} | sort -k 2,2 -k 3,3n -k 4,4n > ${outFileTADs}

    outFileTADs=maxLog2ratio_per_TAD_basedOnGeneTSS_preciselyInTAD.tsv
    awk '{chrom[$9]=$10; start[$9]=$11; end[$9]=$12; if(max[$9]<$8){max[$9]=$8}; cnt[$9]++}END{for(i in max){print i,chrom[i],start[i],end[i],cnt[i],max[i]}}' ${outFileTSS} | sort -k 2,2 -k 3,3n -k 4,4n > ${outFileTADs}        
    exit
fi

for gene in $(awk '{print $1}' ${inRNAseq} | grep -vi GENEID);
do
    echo $gene

    
    log2ratio=$(grep -w $gene ${inRNAseq} | awk '{print $3}')

    grep -w gene ${inGenes} | grep -w ${gene} | awk '{tss=$4; if($7=="-"){tss=$5}; print "chr"$1,$4,$5,$7,tss}' > _tmpGene

    TAD=$(cat _tmpGene ${inTopDomDomains} | awk '{if(NF==5){chrom=$1; start=$2; end=$3; tss=$5}else{c=0;if(chrom==$1 && (($2<=start && start<=$3) || ($2<=end && end<=$3))){for(i=start;i<=end;i++){if($2<=i && i<=$3){c++}}; print $0,c}}}' | sort -k 5,5n | tail -1 | awk '{if($1==""){print "NA"}else{print $4,$1,$2,$3}}')

    echo $(grep -w gene ${inGenes} | grep -w ${gene} | awk '{tss=$4; if($7=="-"){tss=$5}; print $10,$12,"chr"$1,$4,$5,$7,tss}' | sed -e "s,\",,g" -e "s,;,,g") ${log2ratio} $TAD | awk '{for(i=1;i<NF;i++){printf("%-10s\t",$i)}; printf("%-10s\n",$NF)}' >> ${outFileGene}
    tail -n 1 ${outFileGene}
    
    TAD=$(cat _tmpGene ${inTopDomDomains} | awk '{if(NF==5){chrom=$1; start=$2; end=$3; tss=$5}else{c=0;if(chrom==$1 && (($2<=start && start<=$3) || ($2<=end && end<=$3))){for(i=tss;i<=tss;i++){if($2<=i && i<=$3){c++}}; print $0,c}}}' | sort -k 5,5n | tail -1 | awk '{if($1!=""){print $4,$1,$2,$3}else{print "NA"}}')
    
    echo $(grep -w gene ${inGenes} | grep -w ${gene} | awk '{tss=$4; if($7=="-"){tss=$5}; print $10,$12,"chr"$1,$4,$5,$7,tss}' | sed -e "s,\",,g" -e "s,;,,g") ${log2ratio} $TAD | awk '{for(i=1;i<NF;i++){printf("%-10s\t",$i)}; printf("%-10s\n",$NF)}' >>  ${outFileTSS}
    tail -n 1 ${outFileTSS}

    rm -fvr _tmpGene
    #exit
done

sort -k 3,3 -k 4,4n -k 5,5n ${outFileGene} > _tmp_${outFileGene} ; mv _tmp_${outFileGene} ${outFileGene}    
sort -k 3,3 -k 4,4n -k 5,5n ${outFileTSS} > _tmp_${outFileTSS} ; mv _tmp_${outFileTSS} ${outFileTSS}
