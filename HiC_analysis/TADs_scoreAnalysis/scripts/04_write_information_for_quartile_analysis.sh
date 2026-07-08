#List of TADs
#chr2L 5001 20000 PcG1
TADs=scripts_clean/TopDom_domains_hic_WD_NOnub_merge_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp_with_state.bed

#interval1 chrom start end sum sum2 cnt avg stddev
#PcG1 chr2L 5001 20000 1.35996 1.84949 1 1.35996 0
log2RatioPerTAD=avgLog2ratio_per_TAD_basedOnGeneTSS_closestTAD_within5000bp.tsv


#interval1	interval2	distance	av_score_NOnub	av_score_SUMOnub	quartile
#PcG1	PcG2	502499.5	-32.3514	-19.1503	PcG_PcG_Q1
distanceAvgScore=1Ddistances_trans1Dinterval_all_domains_vs_all_domains.tab

outFile=avgLog2ratio_1Ddistances_avgScores_trans1Dinterval_all_domains_vs_all_domains_withControlQuartiles.tab
echo "interval1	chrom1 start1 end1 avgLog2Ratio1 interval2	chrom2 start2 end2 avgLog2Ratio2	distance	av_score_NOnub	av_score_SUMOnub	quartile" > ${outFile}
#PcG1	chr2L	5001	20000	1.35996	PcG2	chr2L	489000	541000	1.09032	        502499.5	-32.3514	-19.1503	PcG_PcG_Q1
#PcG1	chr2L	5001	20000	1.35996	PcG3	chr2L	541000	589000	0.208002	552499.5	-40.7376	-18.5372	PcG_PcG_Q1
awk '{if(NF==4){chrom[$4]=$1;start[$4]=$2;end[$4]=$3}; if(NF==9){avgLog2Ratio[$1]=$8}; if(NF==6){if($NF=="quartile"){next}; avgLog2Ratio1=avgLog2Ratio[$1]; if(avgLog2Ratio1==""){avgLog2Ratio1="NA"}; avgLog2Ratio2=avgLog2Ratio[$2]; if(avgLog2Ratio2==""){avgLog2Ratio2="NA"};print $1,chrom[$1],start[$1],end[$1],avgLog2Ratio1,$2,chrom[$2],start[$2],end[$2],avgLog2Ratio2,$3,$4,$5,$6}}' ${TADs} ${log2RatioPerTAD} ${distanceAvgScore} | awk '{for(i=1;i<NF;i++){printf("%s\t",$i)}; printf("%s\n",$NF)}' >> ${outFile}
#awk '{if(NF==4){chrom[$4]=$1;start[$4]=$2;end[$4]; print $0}}' ${TADs} ${log2RatioPerTAD} ${distanceAvgScore} 

### Write files for the violin plots ###
#interval1	chrom1 start1 end1   avgLog2Ratio1 interval2	chrom2 start2 end2    avgLog2Ratio2  distance	av_score_NOnub	av_score_SUMOnub   quartile
#PcG1	        chr2L  5001   20000  1.35996	   PcG2	        chr2L  489000 541000  1.09032	     502499.5	-32.3514	-19.1503	   PcG_PcG_Q1
inFile=avgLog2ratio_1Ddistances_avgScores_trans1Dinterval_all_domains_vs_all_domains_withControlQuartiles.tab
rm -fvr avgScores_trans1Dinterval_all_domains_vs_all_domains_*_data.tab 1Ddistances_trans1Dinterval_all_domains_vs_all_domains_*_data.tab avgLog2ratio_trans1Dinterval_all_domains_vs_all_domains_*_data.tab avgLog2ratio_trans1Dinterval_all_domains_vs_all_domains_*_data_for_*.tab

states="PcG Active Het Null"
conditions="NOnub SUMOnub"
quartiles="All Q1 Q2 Q3 Q4"
for state1 in ${states} ;
do
    for state2 in ${states} ;
    do
	statePair="${state1}_${state2}"
	
	grep ${statePair} ${inFile} > _tmpStatePair
	nLines=$(wc -l _tmpStatePair | awk '{print $1}')
	if [[ ${nLines} -eq 0 ]];
	then
	    continue
	fi
	echo "${statePair} has ${nLines} entries"
	
	for quartile in ${quartiles} ;
	do
	    grep ${quartile} _tmpStatePair > _tmpStatePairQuartile
	    if [[ ${quartile} == "All" ]];
	    then
		cp _tmpStatePair _tmpStatePairQuartile
	    fi
	    nLines=$(wc -l _tmpStatePairQuartile | awk '{print $1}')
	    echo "${statePair} has ${nLines} entries in quartile ${quartile}"    
	    
	    # Average Scores
	    outFileAvgScores=avgScores_trans1Dinterval_all_domains_vs_all_domains_${statePair}_data.tab
	    for condition in $conditions
	    do
		if [[ $condition == "NOnub" ]];		   
		then
		    awk -v c=${condition} -v q=${quartile} '{printf("%s\t%s\n",c""q,$12)}' _tmpStatePairQuartile >> ${outFileAvgScores}
		fi
		if [[ $condition == "SUMOnub" ]];		   
		then
		    awk -v c=${condition} -v q=${quartile} '{printf("%s\t%s\n",c""q,$13)}' _tmpStatePairQuartile >> ${outFileAvgScores}
		fi		
	    done		    	    
	    
	    # 1D distance
	    outFile1Ddistance=1Ddistances_trans1Dinterval_all_domains_vs_all_domains_${statePair}_data.tab
	    cat _tmpStatePairQuartile | awk -v q=${quartile} '{printf("%s\t%s\n",q,$11)}' >> ${outFile1Ddistance}    
	    
	    # Average Log2Ratio
	    outFileAvgLog2ratio=avgLog2ratio_trans1Dinterval_all_domains_vs_all_domains_${statePair}_data.tab
	    cat _tmpStatePairQuartile | awk -v q=${quartile} '{if($5!="NA" && $10!="NA"){printf("%s\t%s\n",q,($5+$10)*0.5)}}' >> ${outFileAvgLog2ratio}	    
	    for state in ${state1} ${state2} ;
	    do
		outFileAvgLog2ratioPerState=avgLog2ratio_trans1Dinterval_all_domains_vs_all_domains_${statePair}_data_for_${state}.tab
		cat _tmpStatePairQuartile | awk '{if($5!="NA"){print $1,$5}; if($10!="NA"){print $6,$10}}' | grep ${state} | sort -k 1,1 -k 2,2n | uniq | awk -v q=${quartile} '{printf("%s\t%s\n",q,$2)}' >> ${outFileAvgLog2ratioPerState}
		continue
	    done	    
	    
	done

    done
done
rm -fvr _tmpStatePair _tmpStatePairQuartile
