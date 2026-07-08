#!/bin/sh

#SBATCH --job-name geneLog2ratio
#SBATCH -n 1                    # Number of cores. For now 56 is the number max of core available
#SBATCH -t 4-00:00              # Runtime in D-HH:MM
#SBATCH --mem=16G
#SBATCH -o 05_get_gene_Log2ratio_and_TAD.out # File to which STDOUT will be written
#SBATCH -e 05_get_gene_Log2ratio_and_TAD.out # File to which STDERR will be written 

chrom=$1
analysisTAG=$2
echo $chrom $analysisTAG
maxDist=5000 # bp

inTopDomDomains=./scripts_clean/TopDom_domains_hic_WD_NOnub_merge_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp_with_state.bed
inRNAseq=DEseq2_results_WD_SUMOnub_vs_WD_NOnub_flybase_only_NOnub_SUMOnub.tsv
inGenes=dmel-all-r6.36_used_for_rnaseq.gtf



if [[ ${analysisTAG} == "TSS" ]];
then
    outFileTSS=gene_Log2ratio_and_TAD_basedOnGeneTSS_closestTAD_${chrom}.tsv    
    if [[ -e ${outFileTSS} ]];
    then
	outFileTSSAllChrom=gene_Log2ratio_and_TAD_basedOnGeneTSS_closestTAD.tsv	
	cat gene_Log2ratio_and_TAD_basedOnGeneTSS_closestTAD_chr*.tsv > _tmp
	cat _tmp | sort -k 3,3 -k 4,4n -k 5,5n | awk '{for(i=1;i<NF;i++){printf("%-10s\t",$i)}; printf("%-10s\n",$NF)}' > ${outFileTSSAllChrom}

	echo "Get only pairs closer than ${maxDist}bp"	
	outFileWithinMaxDist=gene_Log2ratio_and_TAD_basedOnGeneTSS_closestTAD_within${maxDist}bp.tsv
	
	echo "Get the pairs in TADs"
	cat ${outFileTSSAllChrom} | awk '{if($9=="in"){printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13)}}' > ${outFileWithinMaxDist}
	echo "Add the pairs outside TADs, but closer than ${maxDist}bp"	
	cat ${outFileTSSAllChrom} | awk '$9 ~ /out/' | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,substr($9,1,length($9)-3),$10,$11,$12,$13}' | awk -v maxDist=${maxDist} '{if($9<=maxDist && NF>9){printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9"out",$10,$11,$12,$13)}else{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9"out","NA","NA","NA","NA")}}' >> ${outFileWithinMaxDist}
	echo "Add the pairs outside TADs for completeness!"
	cat ${outFileTSSAllChrom} | awk '{if(NF<10){printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9"out","NA","NA","NA","NA")}}' >> ${outFileWithinMaxDist}	
	echo "Order the final file"
	cat ${outFileWithinMaxDist} | sort -k 3,3 -k 4,4n -k 5,5n | awk '{for(i=1;i<NF;i++){printf("%-10s\t",$i)}; printf("%-10s\n",$NF)}' | uniq > _tmp ; mv _tmp ${outFileWithinMaxDist}

	inFile=${outFileWithinMaxDist}
	outFile=avgLog2ratio_per_TAD_basedOnGeneTSS_closestTAD_within${maxDist}bp.tsv
	awk '{if($10=="" || $10=="NA"){next}; chrom[$10]=$11; start[$10]=$12; end[$10]=$13; sum[$10]+=$8; sum2[$10]+=($8*$8); cnt[$10]++}END{for(i in sum){avg=sum[i]/cnt[i]; avg2=sum2[i]/cnt[i]; stddev=sqrt(avg2-avg*avg); print i,chrom[i],start[i],end[i],sum[i],sum2[i],cnt[i],avg,stddev}}' ${inFile} | sort -k 2,2 -k 3,3n -k 4,4n > ${outFile}
	
	outFile=maxLog2ratio_per_TAD_basedOnGeneTSS_closestTAD_within${maxDist}bp.tsv
	awk '{if($10=="NA"){next}; chrom[$10]=$11; start[$10]=$12; end[$10]=$13; if(max[$10]<$8){max[$10]=$8}; cnt[$10]++}END{for(i in max){print i,chrom[i],start[i],end[i],cnt[i],max[i]}}' ${inFile} | sort -k 2,2 -k 3,3n -k 4,4n > ${outFile}        
	exit
    fi
    exit
fi

c=$(echo $chrom | sed "s/chr//g")
grep -w gene ${inGenes} | awk -v c=$c '{if($1==c){print $0}}' > _tmpGenes${chrom}${analysisTAG}

nGene=0
for gene in $(cat ${inRNAseq} | grep -vi GENEID | awk '{print $1}'); # | grep "FBgn0266322\|FBgn0263584");
do
    nGene=$(($nGene+1))
    echo $gene $nGene
    
    log2ratio=$(grep -w $gene ${inRNAseq} | awk '{print $3}')

    cat _tmpGenes${chrom}${analysisTAG} | grep -w ${gene} | awk '{tss=$4; if($7=="-"){tss=$5}; print "chr"$1,$4,$5,$7,tss}' > _tmpGene${chrom}${analysisTAG}
    touch _tmpGene${chrom}${analysisTAG}
    if [[ ! -s _tmpGene${chrom}${analysisTAG} ]];
    then
	continue
    fi

    if [[ ${analysisTAG} == "TSS" ]];
    then
	echo "If the TSS of a gene is in TADi, you assign this gene to TADi"
	TAD=$(cat _tmpGene${chrom}${analysisTAG} ${inTopDomDomains} | grep ${chrom} | awk '{if(NF==5){chrom=$1; start=$2; end=$3; tss=$5}else{c=0;if(chrom==$1 && (($2<=start && start<=$3) || ($2<=end && end<=$3))){for(i=tss;i<=tss;i++){if($2<=i && i<=$3){c++}}; print $0,c}}}' | sort -k 5,5n | tail -1 | awk '{if($1!=""){print "in",$4,$1,$2,$3}else{print "NA"}}')
	if [[ $TAD == "" ]];
	then
	    echo "The TSS of gene ${gene} is not in a TAD, we assign it to the TAD closest to the TSS:"	
	    TAD=$(cat _tmpGene${chrom}${analysisTAG} ${inTopDomDomains} | grep ${chrom} | awk '{if(NF==5){chrom=$1; start=$2; end=$3; tss=$5}else{d=0; if(chrom!=$1){next}; ds=sqrt((tss-$2)*(tss-$2)); de=sqrt((tss-$3)*(tss-$3)); d=ds; if(de<ds){d=de}; print $0,d}}' | sort -k 5,5n | head -1 | awk '{if($1==""){print "NA"}else{print $5"out",$4,$1,$2,$3}}')
	fi
	echo $(grep -w gene ${inGenes} | grep -w ${gene} | awk '{tss=$4; if($7=="-"){tss=$5}; print $10,$12,"chr"$1,$4,$5,$7,tss}' | sed -e "s,\",,g" -e "s,;,,g") ${log2ratio} $TAD | awk '{for(i=1;i<NF;i++){printf("%-10s\t",$i)}; printf("%-10s\n",$NF)}' >>  ${outFileTSS}
	tail -n 1 ${outFileTSS}
    fi
    rm -fvr _tmpGene${chrom}${analysisTAG}    
done

if [[ ${analysisTAG} == "TSS" ]];
then
    sort -k 3,3 -k 4,4n -k 5,5n ${outFileTSS} > _tmp_${outFileTSS} ; mv _tmp_${outFileTSS} ${outFileTSS}
fi
rm -fvr _tmpGenes${chrom}${analysisTAG} _tmpGene${chrom}${analysisTAG}

