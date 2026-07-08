if [[ ! -e 00_write_genes_fragments.out ]];
then
    cd scripts

    bash 00_write_genes_fragments.sh &> ../00_write_genes_fragments.out

    cd ..
fi

nGenesPerChunk=500
cellTypes="WD_2801 WD_3601 WD_NOnub WD_SUMOnub"

for chrom in $1 ; #chr2L chr2R chr3L chr3R chrX;
do
    echo $chrom

    if [[ ! -d _tmp_${chrom} ]];
    then
	mkdir _tmp_${chrom}
	cd _tmp_${chrom}	
	grep -w ${chrom} ../scripts/all_gene_fragments_dmel-all-r6.36.gtf > all_gene_fragments_${chrom}_dmel-all-r6.36.gtf
	wc -l all_gene_fragments_${chrom}_dmel-all-r6.36.gtf

	nFragments=$(wc -l ./all_gene_fragments_${chrom}_dmel-all-r6.36.gtf | awk -v nGenesPerChunk=${nGenesPerChunk} '{print int($1/nGenesPerChunk+1)*nGenesPerChunk}')
	echo "Number of fragments on chromosome ${chrom} ${nFragments}"
	nChunks=$(echo ${nFragments} ${nGenesPerChunk} | awk '{print int($1/$2)}')
	echo "Number of chunks to split the work ${nChunks}"

	if [[ ! -e scores_trans1Dinterval_${chrom}_gene_vs_domaina_all_fragments.tab ]];
	then
	    
	    for i in $(seq 1 1 ${nChunks});
	    do
		
		fFragment=$(awk -v i=$i -v nGenesPerChunk=${nGenesPerChunk} 'BEGIN{print (i-1)*nGenesPerChunk+1}')
		lFragment=$(awk -v i=$i -v nGenesPerChunk=${nGenesPerChunk} 'BEGIN{print (i-0)*nGenesPerChunk+1}')	
		echo ${fFragment} ${lFragment}
		awk -v fF=${fFragment} -v lF=${lFragment} '{if(fF<=NR && NR<lF){print $0}}' ./all_gene_fragments_${chrom}_dmel-all-r6.36.gtf > all_gene_fragments_${chrom}_dmel-all-r6.36_1.gtf
		wc -l all_gene_fragments_${chrom}_dmel-all-r6.36_1.gtf
		
		for j in $(seq ${i} 1 ${nChunks});
		do
		    (
			if [[ ! -d _tmp_${i}_${j} ]];
			then
			    echo "Directory _tmp_${i}_${j} exists!"
			    
			    mkdir _tmp_${i}_${j}		    
			    cd _tmp_${i}_${j}
			    
			    fFragment=$(awk -v i=$j -v nGenesPerChunk=${nGenesPerChunk} 'BEGIN{print (i-1)*nGenesPerChunk+1}')
			    lFragment=$(awk -v i=$j -v nGenesPerChunk=${nGenesPerChunk} 'BEGIN{print (i-0)*nGenesPerChunk+1}')	
			    echo ${fFragment} ${lFragment}
			    cp ../all_gene_fragments_${chrom}_dmel-all-r6.36_1.gtf .
			    awk -v fF=${fFragment} -v lF=${lFragment} '{if(fF<=NR && NR<lF){print $0}}' ../all_gene_fragments_${chrom}_dmel-all-r6.36.gtf > all_gene_fragments_${chrom}_dmel-all-r6.36_2.gtf
			    wc -l all_gene_fragments_${chrom}_dmel-all-r6.36_2.gtf

			    Rscript ../../scripts/01_get_scores_between_1D_intervals.R ${chrom} &> 01_get_scores_between_1D_intervals_gene_vs_gene_${chrom}.out
			    #echo "chrom1 start1 end1 chrom2 start2 end2 score interval1 interval2 sample" | awk '{for(i=1;i<NF;i++){printf("%s\t",$i)}; printf("%s\n",$NF)}' > _tmp_scores_trans1Dinterval_${chrom}_gene_vs_domain.tab
			    awk 'BEGIN{n=0}{if(NF==4){n=n+1;chrom[n]=$1;start[n]=$2;end[n]=$3;name[n]=$4"_"$1"_"$2"_"$3}; if(NF==8){for(i=0;i<=n;i++){if(start[i]<=$2 && $3<=end[i]){n1=i; continue}}; for(i=0;i<=n;i++){if(start[i]<=$5 && $6<=end[i]){n2=i; continue}}; print chrom[n1],start[n1],end[n1],chrom[n2],start[n2],end[n2],name[n1],name[n2],$7,$8}}' <( sort -k 1,1 -k 2,2n -k 3,3n all_gene_fragments_${chrom}_dmel-all-r6.36_?.gtf | uniq ) scores_trans1Dinterval_${chrom}_gene_vs_domain.tab | grep -v sample >> _tmp_scores_trans1Dinterval_${chrom}_gene_vs_domain.tab
			    mv _tmp_scores_trans1Dinterval_${chrom}_gene_vs_domain.tab scores_trans1Dinterval_${chrom}_gene_vs_domain.tab
			    if [[ ! -e scores_trans1Dinterval_${chrom}_gene_vs_domain_all.tab ]];
			    then
				mv scores_trans1Dinterval_${chrom}_gene_vs_domain.tab scores_trans1Dinterval_${chrom}_gene_vs_domain_all.tab
			    else		   
				cat scores_trans1Dinterval_${chrom}_gene_vs_domain.tab | grep -v chrom1 >> scores_trans1Dinterval_${chrom}_gene_vs_domain_all.tab
				rm -fvr scores_trans1Dinterval_${chrom}_gene_vs_domain.tab
			    fi
			    cd ..
			fi
		    ) &
		done
		pwd
		wait
	    done
	    rm -fvr all_gene_fragments_${chrom}_dmel-all-r6.36_1.gtf all_gene_fragments_${chrom}_dmel-all-r6.36_2.gtf
	    rm -fvr ./*/all_gene_fragments_${chrom}_dmel-all-r6.36_1.gtf ./*/all_gene_fragments_${chrom}_dmel-all-r6.36_2.gtf
	fi
	cd ..
	rm -fvr all_gene_fragments_${chrom}_dmel-all-r6.36_1.gtf all_gene_fragments_${chrom}_dmel-all-r6.36_2.gtf
	cat _tmp_${chrom}/_tmp_*/scores_trans1Dinterval_${chrom}_gene_vs_domain_all.tab | grep -v chrom1 > scores_trans1Dinterval_${chrom}_gene_vs_domain_all_fragments.tab
    fi
    ls -lrtha scores_trans1Dinterval_${chrom}_gene_vs_domain_all_fragments.tab
    
    # Collect scores per gene and domain
    outFile=scores_trans1Dinterval_${chrom}_all_genes_vs_all_domains.tab
    if [[ ! -e ${outFile} ]];
    then
	head scores_trans1Dinterval_${chrom}_gene_vs_domain_all_fragments.tab
	echo
	cat scores_trans1Dinterval_${chrom}_gene_vs_domain_all_fragments.tab | awk '{count1 = split($7, a1, "_"); count2 = split($8, a2, "_"); for (i=1; i<count1-2; i++){for (j=1; j<count2-2; j++){h[a1[i]"_"a2[j]]++; if(a1[i]!=a2[j]){h[a2[j]"_"a1[i]]++}; if((h[a1[i]"_"a2[j]] == 1 && h[a2[j]"_"a1[i]] == 1)){for(l=1;l<=6;l++){printf("%s\t",$l)}; printf("%s\t%s\t",a1[i],a2[j]); for(l=9;l<NF;l++){printf("%s\t",$l)}; printf("%s\n",$NF)}}}; delete h; delete a1; delete a2}' | awk 'BEGIN{FS=OFS="\t"} {sub(/:.*/,"",$7);sub(/:.*/,"",$8)} 1' > _tmp_${outFile}
	
	awk '{if(NF==4){coords[$4]=$1":"$2":"$3}else{print coords[$7],coords[$8],$7,$8,$9,$10}}' ./scripts/TopDom_domains_hic_WD_NOnub_merge_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp_with_state.bed ./scripts/all_genes_dmel-all-r6.36.gtf _tmp_${outFile} | sed "s/:/ /g" | awk '{for(i=1;i<NF;i++){printf("%s\t",$i)}; printf("%s\n",$NF)}' > ${outFile}

	head ${outFile}
	rm -fvr _tmp_${outFile}
	echo
    fi
    ls -lrtha ${outFile}

done
