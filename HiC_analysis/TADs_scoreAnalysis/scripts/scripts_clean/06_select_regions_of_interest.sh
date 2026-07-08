
criterion=Difference
c=9
part=head
differenceThreshold=10 # Minimum difference
distanceThreshold=1500000 # Distance threshold

outFile=06_select_regions_of_interest_by${criterion}.tab
rm -fr ${outFile} _tmp _tmp1

echo "Get PcG regions that increase the average PcG-PcG score in SUMORNAi and are closer than ${distanceThreshold}bp" > ${outFile}
grep PcG_PcG_Q1 1Ddistances_trans1Dinterval_all_domains_vs_all_domains.tab | sort -k 3,3n | awk -v d=$distanceThreshold '{if($4<$5 && $3<d){printf("%s\t%f\n",$0,$5-$4)}}' | awk -v t=${differenceThreshold} '{if(sqrt($NF*$NF)>t){print $0}}' > _tmp
nLines=$(wc -l _tmp | awk '{print $1}')
for l in $(seq 1 1 $nLines)
do
    tag=$(awk -v l=$l '{if(NR==l){print $6,$3,$4,$5,$7}}' _tmp)
    name1=$(awk -v l=$l '{if(NR==l){print $1}}' _tmp)
    name2=$(awk -v l=$l '{if(NR==l){print $2}}' _tmp)

    region1=$(awk -v r=${name1} '{if($7==r){print $1":"$2"-"$3;exit}}' 1Ddistances_trans1Dinterval_all_domains_vs_all_domains_complete.tab)
    region2=$(awk -v r=${name2} '{if($7==r){print $1":"$2"-"$3;exit}}' 1Ddistances_trans1Dinterval_all_domains_vs_all_domains_complete.tab)
    echo $name1 ${region1} $name2 ${region2} $tag >> _tmp1
done
awk '{for(i=1;i<NF;i++){printf("%s\t",$i)}; printf("%s\n",$NF)}' _tmp1 | sort -k 9,9n | tac >> ${outFile}
rm -fr _tmp1
echo >> ${outFile}



echo "Get PcG regions that decrease the average PcG-PcG score in SUMORNAi and are closer than ${distanceThreshold}bp" >> ${outFile}
grep PcG_PcG_Q4 1Ddistances_trans1Dinterval_all_domains_vs_all_domains.tab | sort -k 3,3n | awk -v d=$distanceThreshold '{if($4>$5 && $3<d){printf("%s\t%f\n",$0,$5-$4)}}' | awk -v t=${differenceThreshold} '{if(sqrt($NF*$NF)>t){print $0}}' > _tmp
nLines=$(wc -l _tmp | awk '{print $1}')
for l in $(seq 1 1 $nLines)
do
    tag=$(awk -v l=$l '{if(NR==l){print $6,$3,$4,$5,$7}}' _tmp)
    name1=$(awk -v l=$l '{if(NR==l){print $1}}' _tmp)
    name2=$(awk -v l=$l '{if(NR==l){print $2}}' _tmp)

    region1=$(awk -v r=${name1} '{if($7==r){print $1":"$2"-"$3;exit}}' 1Ddistances_trans1Dinterval_all_domains_vs_all_domains_complete.tab)
    region2=$(awk -v r=${name2} '{if($7==r){print $1":"$2"-"$3;exit}}' 1Ddistances_trans1Dinterval_all_domains_vs_all_domains_complete.tab)
    echo $name1 ${region1} $name2 ${region2} $tag >> _tmp1
done
awk '{for(i=1;i<NF;i++){printf("%s\t",$i)}; printf("%s\n",$NF)}' _tmp1 | sort -k 9,9n | cat >> ${outFile}
rm -fr _tmp1
echo >> ${outFile}


differenceThreshold=20 # Minimum difference


echo "Get regions that increase the average PcG-Active score in SUMORNAi and are closer than ${distanceThreshold}bp" >> ${outFile}
grep PcG_Active_Q1 1Ddistances_trans1Dinterval_all_domains_vs_all_domains.tab | sort -k 3,3n | awk -v d=$distanceThreshold '{if($4<$5 && $3<d){printf("%s\t%f\n",$0,$5-$4)}}' | awk -v t=${differenceThreshold} '{if(sqrt($NF*$NF)>t){print $0}}' > _tmp
nLines=$(wc -l _tmp | awk '{print $1}')
for l in $(seq 1 1 $nLines)
do
    tag=$(awk -v l=$l '{if(NR==l){print $6,$3,$4,$5,$7}}' _tmp)
    name1=$(awk -v l=$l '{if(NR==l){print $1}}' _tmp)
    name2=$(awk -v l=$l '{if(NR==l){print $2}}' _tmp)

    region1=$(awk -v r=${name1} '{if($7==r){print $1":"$2"-"$3;exit}}' 1Ddistances_trans1Dinterval_all_domains_vs_all_domains_complete.tab)
    region2=$(awk -v r=${name2} '{if($7==r){print $1":"$2"-"$3;exit}}' 1Ddistances_trans1Dinterval_all_domains_vs_all_domains_complete.tab)
    echo $name1 ${region1} $name2 ${region2} $tag >> _tmp1
done
awk '{for(i=1;i<NF;i++){printf("%s\t",$i)}; printf("%s\n",$NF)}' _tmp1 | sort -k 9,9n | tac >> ${outFile}
rm -fr _tmp1
echo >> ${outFile}



echo "Get regions that decrease the average PcG-Active score in SUMORNAi and are closer than ${distanceThreshold}bp" >> ${outFile}
grep PcG_Active_Q4 1Ddistances_trans1Dinterval_all_domains_vs_all_domains.tab | sort -k 3,3n | awk -v d=$distanceThreshold '{if($4>$5 && $3<d){printf("%s\t%f\n",$0,$5-$4)}}' | awk -v t=${differenceThreshold} '{if(sqrt($NF*$NF)>t){print $0}}' > _tmp
nLines=$(wc -l _tmp | awk '{print $1}')
for l in $(seq 1 1 $nLines)
do
    tag=$(awk -v l=$l '{if(NR==l){print $6,$3,$4,$5,$7}}' _tmp)
    name1=$(awk -v l=$l '{if(NR==l){print $1}}' _tmp)
    name2=$(awk -v l=$l '{if(NR==l){print $2}}' _tmp)

    region1=$(awk -v r=${name1} '{if($7==r){print $1":"$2"-"$3;exit}}' 1Ddistances_trans1Dinterval_all_domains_vs_all_domains_complete.tab)
    region2=$(awk -v r=${name2} '{if($7==r){print $1":"$2"-"$3;exit}}' 1Ddistances_trans1Dinterval_all_domains_vs_all_domains_complete.tab)
    echo $name1 ${region1} $name2 ${region2} $tag >> _tmp1
done
awk '{for(i=1;i<NF;i++){printf("%s\t",$i)}; printf("%s\n",$NF)}' _tmp1 | sort -k 9,9n | cat >> ${outFile}
rm -fr _tmp1 _tmp
echo >> ${outFile}

outFileRegions=regions_to_plot_by${criterion}.bed
echo "name chrom start end" | awk '{printf("%s\t%s\t%s\t%s\n",$1,$2,$3,$4)}' > ${outFileRegions}
grep -v Get ${outFile} | sed -e "s/_Q1/_Q1_increase/g" -e "s/_Q4/_Q4_decrease/g" | awk '{if(NF>0){print $0}}' | sed -e "s/:/ /g" -e "s/-/ /g" | awk '{print $1,$5,$9,$2,$3,$8}' | sed "s/_/ /g" | awk '{print $1"_"$2"_"$5"_"$6""++h[$3"_"$4"_"$5"_"$6],$7,$8,$9}' | awk '{printf("%s\t%s\t%s\t%s\n",$1,$2,$3,$4)}' >> ${outFileRegions} 
rsync -avz ${outFileRegions} ../03_TADplots_cHiC_and_ChIPseq/

outFileRegions=interesting_loops_by${criterion}.bedpe
echo "name chrom1 start1 end1 chrom2 start2 end2" | awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7)}' > ${outFileRegions}
cat ${outFileRegions}
grep -v Get ${outFile} | awk '{if(NF>0){print $0}}' | sed -e "s/_Q1/_Q1_increase/g" -e "s/_Q4/_Q4_decrease/g" | sed -e "s/:/ /g" -e "s/-/ /g" -e "s/_/ /g" | awk '{print $1"_"$5"_"$11"_"$12,$2,$3,$4,$6,$7,$8,$9"_"$10"_"$11"_"$12}' | awk '{if(NF>0){print $1""++h[$8],$2,$3,$4,$5,$6,$7}}' | awk '{if(NF>0){printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7)}}' >> ${outFileRegions}
rsync -avz ${outFileRegions} ../03_TADplots_cHiC_and_ChIPseq/

