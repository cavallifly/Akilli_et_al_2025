inDir=${PWD}/../01_cool_files/balanced_text_matrices/

chrom1=chr3R
start1=30330000
end1=31400000
chrom2=chr3R
start2=30330000
end2=31400000

diagOff=30000

for resolution in 10000 ;
do
    outBins=_bins_${chrom1}_${start1}_${end1}_${chrom2}_${start2}_${end2}
    if [[ ! -e ${outBins} ]];
    then
	for b1 in $(seq ${start1} ${resolution} $((${end1}-${resolution}))) ;
	do
	    for b2 in $(seq ${start2} ${resolution} $((${end2}-${resolution}))) ;
	    do
		echo ${chrom1} ${b1} $((${b1}+${resolution})) ${chrom2} ${b2} $((${b2}+${resolution})) >> ${outBins}
		echo ${chrom2} ${b2} $((${b2}+${resolution})) ${chrom1} ${b1} $((${b1}+${resolution})) >> ${outBins}
	    done
	done
	sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 -k 5,5n -k 6,6n ${outBins} | uniq > _tmp_${outBins} ; mv _tmp_${outBins} ${outBins}
    fi
    #exit

    for refMap in $(ls -1 ${inDir}/*NOnub*_at_${resolution}bp.tab)
    do
	refCondition=$(echo $refMap | sed -e "s,/, ,g" | awk '{print $NF}' | sed -e "s/_/ /g" | awk '{print $4}')
	echo $refMap $refCondition
	tmpFileRefCond=_tmp_${refCondition}_${chrom1}_${start1}_${end1}_${chrom2}_${start2}_${end2}
	if [[ ! -e ${tmpFileRefCond} ]];
	then
	    grep ${chrom1} ${refMap} | grep ${chrom2} | awk -v c1=${chrom1} -v s1=${start1} -v e1=${end1} -v c2=${chrom2} -v s2=${start2} -v e2=${end2} -v diagOff=${diagOff} '{if(($1==c1 && $1==$4) && (NF==8) && (((s1<=$2 && $2<=e1) && (s1<=$3 && $3<=e1) && (s2<=$5 && $5<=e2) && (s2<=$6 && $6<=e2)) || ((s2<=$2 && $2<=e2) && (s2<=$3 && $3<=e2) && (s1<=$5 && $5<=e1) && (s1<=$6 && $6<=e1)))){v=$8; if(sqrt(($2-$5)*($2-$5))<diagOff){v=0} ; print $1,$2,$3,$4,$5,$6,$7,v; print $4,$5,$6,$1,$2,$3,$7,v}}' | sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 -k 5,5n -k 6,6n > ${tmpFileRefCond}
	    awk '{bins=$1"_"$2"_"$3"_"$4"_"$5"_"$6; if(NF==6){h[bins]=0}else{h[bins]=$8}}END{for(i in h){print i,h[i]}}' ${outBins} ${tmpFileRefCond} | sed "s/_/ /g" | sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 -k 5,5n -k 6,6n > _tmp_${tmpFileRefCond} ; mv _tmp_${tmpFileRefCond} ${tmpFileRefCond}	    
	fi
	#exit
	
	
	# Check rows/columns with no inreraction
	#if [[ ! -e _badCols_${refCondition} ]];
	#then
	#    cat ${refMap} | awk '{cnt[$1"_"$2"_"$3]+=$NF; cnt[$4"_"$5"_"$6]+=$NF}END{for(b in cnt){if(cnt[b]==0) print b}}' > _badCols_${refCondition}
	#    exit
	#fi
	for condMap in $(ls -1 ${inDir}/*SUMOnub*_at_${resolution}bp.tab)
	do
	    if [[ $condMap == $refMap ]];
	    then
		continue
	    fi
	    condition=$(echo $condMap | sed -e "s,/, ,g" | awk '{print $NF}' | sed -e "s/_/ /g" | awk '{print $4}')
	    echo $condMap $condition
	    tmpFileCond=_tmp_${condition}_${chrom1}_${start1}_${end1}_${chrom2}_${start2}_${end2}
	    if [[ ! -e ${tmpFileCond} ]];	
	    then
		grep ${chrom1} ${condMap} | grep ${chrom2} | awk -v c1=${chrom1} -v s1=${start1} -v e1=${end1} -v c2=${chrom2} -v s2=${start2} -v e2=${end2} -v diagOff=${diagOff} '{if(($1==c1 && $1==$4) && (NF==8) && (((s1<=$2 && $2<=e1) && (s1<=$3 && $3<=e1) && (s2<=$5 && $5<=e2) && (s2<=$6 && $6<=e2)) || ((s2<=$2 && $2<=e2) && (s2<=$3 && $3<=e2) && (s1<=$5 && $5<=e1) && (s1<=$6 && $6<=e1)))){v=$8; if(sqrt(($2-$5)*($2-$5))<diagOff){v=0} ; print $1,$2,$3,$4,$5,$6,$7,v; print $4,$5,$6,$1,$2,$3,$7,v}}' | sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 -k 5,5n -k 6,6n > ${tmpFileCond}		
		#grep ${chrom1} ${condMap} | grep ${chrom2} | awk -v c1=${chrom1} -v s1=${start1} -v e1=${end1} -v c2=${chrom2} -v s2=${start2} -v e2=${end2} '{if(($1==c1 && $1==$4) && (NF==8) && (((s1<=$2 && $2<=e1) && (s1<=$3 && $3<=e1) && (s2<=$5 && $5<=e2) && (s2<=$6 && $6<=e2)) || ((s2<=$2 && $2<=e2) && (s2<=$3 && $3<=e2) && (s1<=$5 && $5<=e1) && (s1<=$6 && $6<=e1)))){print $0; print $4,$5,$6,$1,$2,$3,$7,$8}}' > ${tmpFileCond}
		awk '{bins=$1"_"$2"_"$3"_"$4"_"$5"_"$6; if(NF==6){h[bins]=0}else{h[bins]=$8}}END{for(i in h){print i,h[i]}}' ${outBins} ${tmpFileCond} | sed "s/_/ /g" | sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 -k 5,5n -k 6,6n > _tmp_${tmpFileCond} ; mv _tmp_${tmpFileCond} ${tmpFileCond}	    		
	    fi
	    #exit
	    
	    #if [[ ! -e _badCols_${condition} ]];
	    #then
		#cat ${condMap} | awk '{cnt[$1"_"$2"_"$3]+=$NF; cnt[$4"_"$5"_"$6]+=$NF}END{for(b in cnt){if(cnt[b]==0)print b}}' > _badCols_${condition}
	    #fi
	    #cat _badCols_${condition} _badCols_${refCondition} | sort | uniq > _badCols_${refCondition}_${condition}
	    outFile=diffZscores_${condition}_vs_${refCondition}_${chrom1}_${start1}_${end1}_${chrom2}_${start2}_${end2}_at_${resolution}bp_diagOff_${diagOff}bp.tsv
	    
            if [[ -e ${outFile} ]];
            then
		sort -k 7,7n <(awk '{if(sqrt(($2-$5)*($2-$5))>300000 && $2>$5){print $0}}' ${outFile}) | head -10
		sort -k 7,7n <(awk '{if(sqrt(($2-$5)*($2-$5))>300000 && $2>$5){print $0}}' ${outFile}) | tail -10
		continue
            fi    
            wc -l ${refMap} ${condMap}
            #(diff - mean(diff))/(standard deviation(diff))	    
            #head $refMap $condMap 
	    
            avg=$(paste ${tmpFileCond} ${tmpFileRefCond}    | awk '{for(i=1;i<=6;i++){if($(i)!=$(i+7)){next}}; if($7==0 && $7==$14){next}; d=$7-$14; s+=d; cnt++}END{print s/cnt}')
            stddev=$(paste ${tmpFileCond} ${tmpFileRefCond} | awk '{for(i=1;i<=6;i++){if($(i)!=$(i+7)){next}}; if($7==0 && $7==$14){next}; d=$7-$14; s+=d; s2+=d*d; cnt++}END{avg=s/cnt; avg2=s2/cnt; stddev=sqrt(avg2-avg*avg); print stddev}')
            cnt=$(paste ${tmpFileCond} ${tmpFileRefCond}    | awk '{for(i=1;i<=6;i++){if($(i)!=$(i+7)){next}}; if($7==0 && $7==$14){next}; cnt++}END{print cnt}')
            tot=$(paste ${tmpFileCond} ${tmpFileRefCond}    | awk '{for(i=1;i<=6;i++){if($(i)!=$(i+7)){next}}; cnt++}END{print cnt}')
            echo "$refMap $condMap mean = $avg std. dev. = $stddev with ${cnt} entries over $tot"

	    awk -v a=$avg -v s=${stddev} '{for(i=1;i<=6;i++){if($(i)!=$(i+7)){next}}; d=$7-$14; v=(d-a)/s; printf("%s\t%s\t%s\t%s\t%s\t%s\t%.10f\n",$1,$2,$3,$4,$5,$6,v)}' <(paste ${tmpFileCond} ${tmpFileRefCond}) > ${outFile}
	    
            head ${outFile}
	    echo
	done # Close cycle over $refMap
    done # Close cycle over $condMap
done #  Close cycle over $resolution
#rm -fvr _*

