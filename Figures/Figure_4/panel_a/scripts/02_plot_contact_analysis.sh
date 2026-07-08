scriptsdir=${PWD}/scripts/

for modMatrix in $(ls -1 d*.tsv 2> /dev/null) ;
do
    start=$(head -1 ${modMatrix} | awk '{print $2}')
    mapResolution=$(head -1 ${modMatrix} | awk '{print $3-$2}')
    size=$(tail -1 ${modMatrix} | awk -v s=${start} -v r=${mapResolution} '{print ($6-s)/r}')
    echo $modMatrix ${start} ${mapResolution} ${size}
    
    modMatrix0=${modMatrix}
    
    awk -v s=$start -v r=$mapResolution '{print ($2-s)/r,($5-s)/r,$NF}' ${modMatrix} | awk '{if($1!=$2) print $0}' > _tmp_mod
    modMatrix=_tmp_mod

    #cbmin=$(awk '{if(NR==1){min=$3}; if($3<min) min=$3}END{print min}' contact_map.tab)
    cbmax=3 #$(awk '{v=sqrt($3*$3); if(v>max) max=v}END{print max}' contact_map.tab)
    cbmaxMid=$(echo ${cbmax} | awk '{printf("%.1f",$1/2.)}')
    echo $cbmin $cbmax
    
    outFile=${modMatrix0%.tsv}_cbmax_${cbmax}.png
    #outFile=matrix_${exp}_rank.png	
    if [[ -e ${outFile} ]];
    then
	ls -lrtha ${outFile}
	rm -fvr _tmp_mod
	continue
    fi
    modFactor=1.0 #$(awk '{d=sqrt(($1-$2)*($1-$2)); if(d==10){sum+=$3;cnt++}}END{print sum/cnt}' ${modMatrix})
    maxRank=$(wc -l $modMatrix | awk '{print $1}')
    
    diagOff=0
    
    # Model
    cp ${modMatrix} contact_map.tab
    
    awk -v size=${size} 'BEGIN{for(i=0;i<size;i++){print i,i,0}}' >> contact_map.tab
    sort -k 1,1n -k 2,2n contact_map.tab > _a.tab ; mv _a.tab contact_map.tab
    
    for maxFactor in 0.1 #1 2 3 4 5 0.001 0.01 0.1 1 ;
    do
	for minFactor in 100 #1 2 3 4 5 0.001 0.01 0.1 1 10 100 ;
	do
	    for scale in 0.75 ; #0.80 0.90 0.95 ;
	    do
		for factor in 1 ;
		do
		    if [[ ! -e ${outFile} ]];
		    then
			wc -l contact_map.tab
			
			sed -e "s/XXXresolutionXXX/${mapResolution}/g" -e "s/XXXstartXXX/${start}/g" -e "s/XXXmaxFactorXXX/${maxFactor}/g" -e "s/XXXminFactorXXX/${minFactor}/g" -e "s/XXXcbminXXX/${cbmin}/g" -e "s/XXXcbmaxMidXXX/${cbmaxMid}/g" -e "s/XXXsizeXXX/${size}/g" -e "s/XXXfactorXXX/${factor}/g" -e "s/XXXcbmaxXXX/${cbmax}/g" -e "s/XXXscaleXXX/${scale}/g" -e "s/XXXconditionXXX/${exp}/g" ${scriptsdir}/02_plot_contact_matrix.gp | gnuplot
			mv contact_map.ps ${outFile%.png}.ps
			
			bash ${scriptsdir}/ps2pdf.sh #&> /dev/null
			
		    fi
		    ls -lrtha ${outFile}
		done
	    done
	done
    done
    #exit
done # Close cycle over $exp
rm -fr contact_map.tab contact_map_diff.tab _tmp_mod _tmp_cHiC _exp_rank _mod_rank
	
cd .. # Exit ${condition}
