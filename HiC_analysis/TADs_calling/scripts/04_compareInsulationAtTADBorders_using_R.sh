for INSnorm in raw 01normalized
do
    for condition in NOnub SUMOnub
    do
	# Top-Dom on hic_WD_${condition}
	inFile1=TopDom_domains_hic_WD_${condition}_merge_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp_with_WD_NOnub_merge_${INSnorm}_insulation.tsv
	inFile2=TopDom_domains_hic_WD_${condition}_merge_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp_with_WD_SUMOnub_merge_${INSnorm}_insulation.tsv
	outFile=TopDom_domains_hic_WD_${condition}_merge_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp_with_WD_NOnub_and_SUMOnub_merge_${INSnorm}_insulation.tsv
	
	head ${inFile1}
	head ${inFile2}
	dos2unix ${inFile1}
	dos2unix ${inFile2}
	
	echo "chrom  start   end     insNOnub1       insNOnub2       insSUMOnub1     insSUMOnub2" | awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7)}' >  ${outFile}
	paste ${inFile1} ${inFile2} | grep -v ins1 | awk '{if($2<$3)printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,int($2/10)*10,int($3/10)*10,$4,$5,$9,$10)}' | awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7)}' >> ${outFile}
	head ${outFile}

	#grep -v "NA\|start" ${outFile} | awk '{if(NR!=1){print "insNOnub",$1"_"int($2/10)*10,$4; print "insNOnub",$1"_"int($3/10)*10,$5; print "insSUMOnub",$1"_"int($2/10)*10,$6; print "insSUMOnub",$1"_"int($3/10)*10,$7;}}' | sort -k 1,1 | uniq | awk '{printf("%s\t%s\t%s\n",$1,$2,$3)}' | sort -k 1,1 > ${outFile%.tsv}_data.tsv
	grep -v "NA\|start" ${outFile} | awk '{if(NR!=1){print "insNOnub",$1"_"int($2/10)*10,$4; print "insNOnub",$1"_"int($3/10)*10,$5; print "insSUMOnub",$1"_"int($2/10)*10,$6; print "insSUMOnub",$1"_"int($3/10)*10,$7;}}' | sort -k 1,1 | uniq | awk '{printf("%s\t%s\n",$1,$3)}' | sort -k 1,1 > ${outFile%.tsv}_data.tsv

	outFileSP=TopDom_domains_hic_WD_${condition}_merge_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp_with_WD_NOnub_and_SUMOnub_merge_${INSnorm}_insulation_scatterData.tsv
	echo "chrom  start   end     insNOnub insSUMOnub" | awk '{printf("%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5)}' >  ${outFileSP}
	paste ${inFile1} ${inFile2} | grep -v ins1 | awk '{if($2<$3)printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,int($2/10)*10,int($3/10)*10,$4,$5,$9,$10)}' | awk '{printf("%s\t%s\t%s\t%s\t%s\n",$1,$2,$2+1000,$4,$6); printf("%s\t%s\t%s\t%s\t%s\n",$1,$3-1000,$3,$5,$7)}' >> ${outFileSP}
	
	#for r in $(awk '{h[$2]++}END{for(i in h){if(h[i]!=2) print i}}' ${outFile%.tsv}_data.tsv);
	#do
	#    echo $r
	#    chrom=$(echo $r | sed "s/_/ /g" | awk '{print $1}')
	#    start=$(echo $r | sed "s/_/ /g" | awk '{print $2}')    
	
	#    grep ${chrom} ${inFile1} ${inFile2} | grep ${start}
	#    grep ${chrom} ${outFile} | grep ${start}
	
	#done
	#exit
	awk '{h[$1]++}END{for(i in h) print i,h[i]}' ${outFile%.tsv}_data.tsv | sort -k 1,1 
    done
done
