#!/bin/sh

#SBATCH --job-name TADoverlaps
#SBATCH -n 1                    # Number of cores. For now 56 is the number max of core available
#SBATCH -t 4-00:00              # Runtime in D-HH:MM
#SBATCH --mem=16G
#SBATCH -o XX_compute_TADs_overlaps.out # File to which STDOUT will be written
#SBATCH -e XX_compute_TADs_overlaps.out # File to which STDERR will be written 

resolution=1000
maxDistGeneTAD=0

### Overlap of Control and SUMORNAi TADs
control=NOnub
replicate1=merge
controlTADs=TopDom_domains_hic_contacts_WD_${control}_${replicate1}_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp_with_state.bed
grep -v chrom ${controlTADs} > _tmp_${controlTADs}
controlTADs=_tmp_${controlTADs}

treatment=SUMOnub
replicate2=merge
treatmentTADs=TopDom_domains_hic_contacts_WD_${treatment}_${replicate2}_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp.bedpe
grep -v chrom ${treatmentTADs} > _tmp_${treatmentTADs}
treatmentTADs=_tmp_${treatmentTADs}

wc -l ${controlTADs} ${treatmentTADs}

outFile=overlaps_of_WD_${control}_${replicate1}_TADs_with_WD_${treatment}_${replicate2}_TADs.txt

if [[ ! -e ${outFile} ]];
then

    nControlTADs=$(wc -l ${controlTADs} | awk '{print $1}')
    for nControlTAD in $(seq 1 1 ${nControlTADs}) ;
    do
	
	chrom=$(awk -v n=${nControlTAD} '{if(NR==n){print $1}}' ${controlTADs})
	start=$(awk -v n=${nControlTAD} '{if(NR==n){print $2}}' ${controlTADs})
	end=$(  awk -v n=${nControlTAD} '{if(NR==n){print $3}}' ${controlTADs})
	state=$(awk -v n=${nControlTAD} '{if(NR==n){print $4}}' ${controlTADs})
	echo ${chrom} ${start} ${end} $state

	rm -fvr _tmp
	awk -v nT=${nControlTAD} -v c=${chrom} -v s=${start} -v e=${end} -v state=${state} -v r=${resolution} '{m++; n=0; t=0; if(c==$1){n=0; for(i=(s+int(0.5*r)); i<=(e-int(0.5*r)); i+=r){if($2<=i && i<=$3){n++}; t++}; if(n>0){print nT,c,s,e,state,"->",m,$1,$2,$3,n,t,n/t*100}}}' ${treatmentTADs} | sort -k 13,13n > _tmp
	touch _tmp
	cat _tmp
	nLines=$(wc -l _tmp | awk '{print $1}')
	if [[ $nLines -gt 0 ]];
	then
	    tail -1 _tmp >> ${outFile}
	else
	    echo "$nControlTAD ${chrom} ${start} ${end} ${state} -> NA NA NA NA NA NA 0" >> ${outFile}
	fi
	
    done
    
    awk '{h[$1]=$0}END{for(i in h){print h[i]}}' ${outFile} | sort -k 1,1n | awk '{printf("TAD%s\t%s\t%s\t%s\t%s\t%s\tTAD%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13)}' > _tmp ; mv _tmp ${outFile}
fi
cat ${outFile} | sort -k 13,13n | awk '{bin=int($13); h[bin]++; t++}END{for(i in h){print i,h[i]/t}}' | awk '{h[$1]=$2;max=$1}END{for(i=0;i<=max;i++){s+=h[i]; print i,s*100}}' > ${outFile%.txt}_cumulative.txt



### Overlap of Control and SUMORNAi TADs
control=SUMOnub
replicate1=merge
controlTADs=TopDom_domains_hic_WD_${control}_${replicate1}_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp.bedpe
grep -v chrom ${controlTADs} > _tmp_${controlTADs}
controlTADs=_tmp_${controlTADs}

treatment=NOnub
replicate2=merge
treatmentTADs=TopDom_domains_hic_WD_${treatment}_${replicate2}_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp_with_state.bed
grep -v chrom ${treatmentTADs} > _tmp_${treatmentTADs}
treatmentTADs=_tmp_${treatmentTADs}

wc -l ${controlTADs} ${treatmentTADs}

outFile=overlaps_of_WD_${control}_${replicate1}_TADs_with_WD_${treatment}_${replicate2}_TADs.txt

if [[ ! -e ${outFile} ]];
then

    nControlTADs=$(wc -l ${controlTADs} | awk '{print $1}')
    for nControlTAD in $(seq 1 1 ${nControlTADs}) ;
    do
	
	chrom=$(awk -v n=${nControlTAD} '{if(NR==n){print $1}}' ${controlTADs})
	start=$(awk -v n=${nControlTAD} '{if(NR==n){print $2}}' ${controlTADs})
	end=$(  awk -v n=${nControlTAD} '{if(NR==n){print $3}}' ${controlTADs})
	state=$(awk -v n=${nControlTAD} '{if(NR==n){print $4}}' ${controlTADs})
	echo ${chrom} ${start} ${end} $state

	rm -fvr _tmp
	awk -v nT=${nControlTAD} -v c=${chrom} -v s=${start} -v e=${end} -v state=${state} -v r=${resolution} '{m++; n=0; t=0; if(c==$1){n=0; for(i=(s+int(0.5*r)); i<=(e-int(0.5*r)); i+=r){if($2<=i && i<=$3){n++}; t++}; if(n>0){print nT,c,s,e,state,"->",m,$1,$2,$3,n,t,n/t*100}}}' ${treatmentTADs} | sort -k 13,13n > _tmp
	touch _tmp
	cat _tmp
	nLines=$(wc -l _tmp | awk '{print $1}')
	if [[ $nLines -gt 0 ]];
	then
	    tail -1 _tmp >> ${outFile}
	else
	    echo "$nControlTAD ${chrom} ${start} ${end} ${state} -> NA NA NA NA NA NA 0" >> ${outFile}
	fi
	
    done
    
    awk '{h[$1]=$0}END{for(i in h){print h[i]}}' ${outFile} | sort -k 1,1n | awk '{printf("TAD%s\t%s\t%s\t%s\t%s\t%s\tTAD%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13)}' > _tmp ; mv _tmp ${outFile}
fi
cat ${outFile} | sort -k 13,13n | awk '{bin=int($13); h[bin]++; t++}END{for(i in h){print i,h[i]/t}}' | awk '{h[$1]=$2;max=$1}END{for(i=0;i<=max;i++){s+=h[i]; print i,s*100}}' > ${outFile%.txt}_cumulative.txt



### Overlap of Control TADs in Rep5 and Rep6
control=NOnub
replicate1=Rep5
controlTADs=TopDom_domains_hic_contacts_WD_${control}_${replicate1}_all_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp.tsv
grep -v chrom ${controlTADs} > _tmp_${controlTADs}
controlTADs=_tmp_${controlTADs}

treatment=NOnub
replicate2=Rep6
treatmentTADs=TopDom_domains_hic_contacts_WD_${treatment}_${replicate2}_all_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp.tsv
grep -v chrom ${treatmentTADs} > _tmp_${treatmentTADs}
treatmentTADs=_tmp_${treatmentTADs}

wc -l ${controlTADs} ${treatmentTADs}

outFile=overlaps_of_WD_${control}_${replicate1}_TADs_with_WD_${treatment}_${replicate2}_TADs.txt
if [[ ! -e ${outFile} ]];
then

    nControlTADs=$(wc -l ${controlTADs} | awk '{print $1}')
    for nControlTAD in $(seq 1 1 ${nControlTADs}) ;
    do
	
	chrom=$(awk -v n=${nControlTAD} '{if(NR==n){print $1}}' ${controlTADs})
	start=$(awk -v n=${nControlTAD} '{if(NR==n){print $2}}' ${controlTADs})
	end=$(  awk -v n=${nControlTAD} '{if(NR==n){print $3}}' ${controlTADs})
	state=$(awk -v n=${nControlTAD} '{if(NR==n){print $4}}' ${controlTADs})
	echo ${chrom} ${start} ${end} $state

	rm -fvr _tmp	
	awk -v nT=${nControlTAD} -v c=${chrom} -v s=${start} -v e=${end} -v state=${state} -v r=${resolution} '{m++; n=0; t=0; if(c==$1){n=0; for(i=(s+int(0.5*r)); i<=(e-int(0.5*r)); i+=r){if($2<=i && i<=$3){n++}; t++}; if(n>0){print nT,c,s,e,state,"->",m,$1,$2,$3,n,t,n/t*100}}}' ${treatmentTADs} | sort -k 13,13n > _tmp
	touch _tmp
	cat _tmp
	nLines=$(wc -l _tmp | awk '{print $1}')
	if [[ $nLines -gt 0 ]];
	then
	    tail -1 _tmp >> ${outFile}
	else
	    echo "$nControlTAD ${chrom} ${start} ${end} ${state} -> NA NA NA NA NA NA 0" >> ${outFile}
	fi
	
    done
    
    awk '{h[$1]=$0}END{for(i in h){print h[i]}}' ${outFile} | sort -k 1,1n | awk '{printf("TAD%s\t%s\t%s\t%s\t%s\t%s\tTAD%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13)}' > _tmp ; mv _tmp ${outFile}
fi
cat ${outFile} | sort -k 13,13n | awk '{bin=int($13); h[bin]++; t++}END{for(i in h){print i,h[i]/t}}' | awk '{h[$1]=$2;max=$1}END{for(i=0;i<=max;i++){s+=h[i]; print i,s*100}}' > ${outFile%.txt}_cumulative.txt



### Overlap of Control TADs in Rep5 and Rep6
control=NOnub
replicate1=Rep6
controlTADs=TopDom_domains_hic_contacts_WD_${control}_${replicate1}_all_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp.tsv
grep -v chrom ${controlTADs} > _tmp_${controlTADs}
controlTADs=_tmp_${controlTADs}

treatment=NOnub
replicate2=Rep5
treatmentTADs=TopDom_domains_hic_contacts_WD_${treatment}_${replicate2}_all_dm6_NA_window_5_at_5000bp_filtered_refined_at_1000bp.tsv
grep -v chrom ${treatmentTADs} > _tmp_${treatmentTADs}
treatmentTADs=_tmp_${treatmentTADs}
wc -l ${controlTADs} ${treatmentTADs}

outFile=overlaps_of_WD_${control}_${replicate1}_TADs_with_WD_${treatment}_${replicate2}_TADs.txt
if [[ ! -e ${outFile} ]];
then

    
    nControlTADs=$(wc -l ${controlTADs} | awk '{print $1}')
    for nControlTAD in $(seq 1 1 ${nControlTADs}) ;
    do
	chrom=$(awk -v n=${nControlTAD} '{if(NR==n){print $1}}' ${controlTADs})
	start=$(awk -v n=${nControlTAD} '{if(NR==n){print $2}}' ${controlTADs})
	end=$(  awk -v n=${nControlTAD} '{if(NR==n){print $3}}' ${controlTADs})
	state=$(awk -v n=${nControlTAD} '{if(NR==n){print $4}}' ${controlTADs})
	echo ${chrom} ${start} ${end} $state

	rm -fvr _tmp	
	awk -v nT=${nControlTAD} -v c=${chrom} -v s=${start} -v e=${end} -v state=${state} -v r=${resolution} '{m++; n=0; t=0; if(c==$1){n=0; for(i=(s+int(0.5*r)); i<=(e-int(0.5*r)); i+=r){if($2<=i && i<=$3){n++}; t++}; if(n>0){print nT,c,s,e,state,"->",m,$1,$2,$3,n,t,n/t*100}}}' ${treatmentTADs} | sort -k 13,13n > _tmp
	touch _tmp
	cat _tmp
	nLines=$(wc -l _tmp | awk '{print $1}')
	if [[ $nLines -gt 0 ]];
	then
	    tail -1 _tmp >> ${outFile}
	else
	    echo "$nControlTAD ${chrom} ${start} ${end} ${state} -> NA NA NA NA NA NA 0" >> ${outFile}
	fi
	
    done
    
    awk '{h[$1]=$0}END{for(i in h){print h[i]}}' ${outFile} | sort -k 1,1n | awk '{printf("TAD%s\t%s\t%s\t%s\t%s\t%s\tTAD%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13)}' > _tmp ; mv _tmp ${outFile}
fi
cat ${outFile} | sort -k 13,13n | awk '{bin=int($13); h[bin]++; t++}END{for(i in h){print i,h[i]/t}}' | awk '{h[$1]=$2;max=$1}END{for(i=0;i<=max;i++){s+=h[i]; print i,s*100}}' > ${outFile%.txt}_cumulative.txt



outFile=overlaps_of_TADs_data.txt
if [[ ! -e ${outFile} ]];
then
    for state in All PcG Active Het Null ;
    do
	control=NOnub
	replicate1=merge
	treatment=SUMOnub
	replicate2=merge    
	grep $state overlaps_of_WD_${control}_${replicate1}_TADs_with_WD_${treatment}_${replicate2}_TADs.txt | awk -v s=${state} '{print s,$NF}' > _tmp
	if [[ $state == "All" ]];
	then
	    awk -v s=${state} '{print s,$NF}' overlaps_of_WD_${control}_${replicate1}_TADs_with_WD_${treatment}_${replicate2}_TADs.txt > _tmp
	    awk -v s=${state} '{print s"Rep1vsRep2",$NF}' overlaps_of_WD_${control}_Rep5_TADs_with_WD_${control}_Rep6_TADs.txt >> _tmp
	fi
	cat _tmp >> ${outFile}
    done
fi
awk '{h[$1]++}END{for(i in h){print i,h[i]}}' ${outFile} | sort -k 1,1n



geneLog2Ratio=gene_Log2ratio_and_TAD_basedOnGeneTSS_closestTAD_within5000bp.tsv
### Number of Down and Up genes vs TAD overlap
for DEG in Up Down NS ;
do
    inFile=${DEG}_DEgenes_WD_SUMOnub_vs_WD_NOnub_padj0.05_FC1.5_flybase_only_NOnub_SUMOnub.tsv
    outFile=overlap_of_${DEG}_DEgenes_WD_SUMOnub_vs_WD_NOnub_padj0.05_FC1.5_flybase_only_NOnub_SUMOnub.txt

    wc -l ${inFile}
    if [[ ! -e ${outFile} ]];
    then
	for gene in $(awk '{print $1}' ${inFile} | grep -v GENEID) ;
	do
	    geneID=$(  grep $gene ${geneLog2Ratio} | awk '{print $2}')
	    TAD=$(     grep $gene ${geneLog2Ratio} | awk '{print $10}')
	    if [[ $TAD != "NA" ]];
	    then
		echo $gene $geneID $TAD $(grep -w $TAD overlaps_of_WD_NOnub_merge_TADs_with_WD_SUMOnub_merge_TADs.txt | awk '{print $NF}') >> ${outFile}
	    else
		echo "$gene $geneID NA NA" >> ${outFile}
	    fi
	done
    fi    
    wc -l ${outFile}
    cat ${outFile} | grep -v NA | sort -k 4,4n | awk '{bin=int($4); h[bin]++; t++}END{for(i in h){print i,h[i]/t}}' | awk '{h[$1]=$2;max=$1}END{for(i=0;i<=max;i++){s+=h[i]; print i,s*100}}' > ${outFile%.txt}_cumulative.txt
done
