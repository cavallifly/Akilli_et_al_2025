#!/bin/bash

#SBATCH --job-name makeQuant
#SBATCH -n 1                    # Number of cores. For now 56 is the number max of core available
#SBATCH -t 4-00:00              # Runtime in D-HH:MM
#SBATCH -o 04_make_violinPlots_for_pilePlot_entries_quantification_betweenConditions.out # File to which STDOUT will be written
#SBATCH -e 04_make_violinPlots_for_pilePlot_entries_quantification_betweenConditions.out # File to which STDERR will be written

conditions="SUMOnub_merge NOnub_merge"

treatment=$(echo $conditions | awk '{print $1}')
treatmentName=$(echo $conditions | awk '{print $1}' | sed "s/_/ /g" | awk '{print $1}')
control=$(echo $conditions | awk '{print $2}')
controlName=$(echo $conditions | awk '{print $2}' | sed "s/_/ /g" | awk '{print $1}')

#window=33
#inFile=$(ls -1 *central_${window}.txt | grep 1000bp)

echo "# Compare merge between ${treatment} and ${control}"
#hic_contacts_WD_SUMOnub_Rep6_all_dm6_NA_at_5000bp_PcGstate_local_rescaled.clpy
samples=$(cat ${inFile} | awk '{print $1}' | sort | uniq | grep "${treatment}\|${control}" | sed -e "s/hic_contacts_WD_/ /g" -e "s/all_dm6_NA_at_//g" -e "s/state_local_rescaled.clpy//g" -e "s/${controlName}_//g" -e "s/${treatmentName}_//g" | sort | uniq)
echo $samples

ls -lrtha enrichmentValues_WD_allConditions_*tab 2> /dev/null
rm -fvr enrichmentValues_WD_allConditions_*.tab

for sample in ${samples};
do
    replicate=$(echo $sample | sed -e "s/_/ /g" | awk '{print $1}')
    resolution=$(echo $sample | sed -e "s/_/ /g" | awk '{print $2}')    
    state=$(echo $sample | sed -e "s/_/ /g" | awk '{print $3}')
    
    for inFile in $(ls -1 **${replicate}*${resolution}*${state}*txt);
    do
	mapRegion=$(echo $inFile | sed -e "s/_/ /g" -e "s/\.txt/ /g" | awk '{print $(NF-1)"_"$NF"pixels_data"}')
	echo $inFile $mapRegion
	
	outFile=enrichmentValues_WD_allConditions_${replicate}_at_${resolution}_${mapRegion}.tab
	echo $outFile

	cat ${inFile} | awk '{if((($2-$3)*($2-$3))>1 && $2>$3){print $0}}' | sed -e "s/hic_contacts_WD_/ /g" -e "s/all_dm6_NA_at_//g" -e "s/state_local_rescaled.clpy//g" -e "s/_${replicate}//g" -e "s/_${resolution}//g" -e "s/${controlName}/Control/g" -e "s/${treatmentName}/SUMORNAi/g" | awk '{print $1,$NF}' >> ${outFile}
    done
    echo $outFile
    awk '{h[$1]++}END{for(i in h){print i,h[i]}}' ${outFile} | sort -k 1,1n    
done

windowINS=25
outFileINS=enrichmentValues_WD_allConditions_merge_at_${resolution}_INSregions_${windowINS}pixels_data.tab
rm -fvr ${outFileINS}
for sample in ${samples};
do
    replicate=$(echo $sample | sed -e "s/_/ /g" | awk '{print $1}')
    resolution=$(echo $sample | sed -e "s/_/ /g" | awk '{print $2}')    
    state=$(echo $sample | sed -e "s/_/ /g" | awk '{print $3}')

    for file1 in $(ls -1 hic_contacts_WD_*_${replicate}_all_dm6_NA_at_${resolution}_*${state}*_bottomRightTAD_66.txt)
    do
	file2=$(echo $file1 | sed "s/_bottomRightTAD_66/_topLeftTAD_66/g")
	head $file1 ${file2}
	
	echo
	wc -l $file1 $file2
	awk -v windowINS=${windowINS} '{if(((0<=$2 && $2<=33) && (0<=$3 && $3<=33)) || ((33<=$2 && $2<=65) && (33<=$3 && $3<=65)) || ((65<=$2 && $2<=99) && (65<=$3 && $3<=99)) || (sqrt(($2-$3)*($2-$3))>=windowINS+1 || sqrt(($2-$3)*($2-$3))<=1 || $2<$3)){next}; print $0}' ${file1} ${file2} | awk '{print $0}' | sed -e "s/hic_contacts_WD_/ /g" -e "s/all_dm6_NA_at_//g" -e "s/state_local_rescaled.clpy//g" -e "s/_${resolution}//g" -e "s/${controlName}/Control/g" -e "s/${treatmentName}/SUMORNAi/g" | awk '{print $1,$NF}' >> ${outFileINS}
    done
    echo $outFileINS
    awk '{h[$1]++}END{for(i in h){print i,h[i]}}' ${outFileINS} | sort -k 1,1n        
done
rm -fvr *topLeftTAD_*pixels_data* *bottomRightTAD_*pixels_data*


Rscript scripts/04_make_violinPlots_for_pilePlot_entries_quantifications.R betweenConditions
exit	
