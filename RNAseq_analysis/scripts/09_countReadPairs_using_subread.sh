
assembly=$1

gtfFile=$(ls -1 /zssd/scratch/DBs/*/UCSC/${assembly}/Annotation/Archives/archive-c*/Genes/genes.gtf)
ls -lrtha ${gtfFile}

checkSingleEnd=$(grep SingleEnd ??_guessStrandedness_using_rseqc.log | wc -l)
if [[ ${checkSingleEnd} -eq 0 ]];
then
    tag=countReadPairs
    tag1=paired-end
else
    tag=countReads
    tag1=single-end
fi

strandedness=$(tail -2 ??_guessStrandedness_using_rseqc.log | awk '{if(NR==1){printf("%s ", $NF)}else{printf("%s\n", $NF)}}' | awk '{if($1>0.80){print 1; next}; if($2>0.80){print 2; next}; print 0}')
#tail -2 ??_guessStrandedness_using_rseqc.log

echo "featureCounts performs strand‐specific read counting."
echo ${strandedness} | awk '{if($1==0){v="unstranded"}; if($1==1){v="stranded"}; if($1==2){v="reversely stranded"}; print "Your experiment is most probably "v", so we will use strandedness value in featureCounts "$1}'

for dir in $(ls -1 | grep Rep | grep all);
do
    if [[ ! -d ${dir} ]];
    then
	echo "${dir} is not a directory for a sample"
	continue
    fi
    
    cd $dir

    bamFile=$(ls -1 *bam)
    echo "Counting reads per gene from ${bamFile} in ${dir} using featureCounts"    
    echo "${strandedness} ${gtfFile} ${bamFile}"

    if [[ -e ${tag}_using_subread_summary.pdf ]]
    then
	cd ..
	continue
    fi
    if [[ -e doing_countReadPairs ]]
    then
	cd ..
	continue
    fi   
    touch doing_countReadPairs
    
    if [[ ! -e ${tag}_using_subread.tab ]]
    then	
	if [[ ${checkSingleEnd} -eq 0 ]];
	then
	    conda run -n subread featureCounts -p --countReadPairs -s ${strandedness} -t exon -g gene_id -a ${gtfFile} -o ${tag}_using_subread.tab ${bamFile}
	    #conda run -n subread featureCounts -p --countReadPairs -s ${strandedness} -t exon -g transcript_id -a ${gtfFile} -o ${tag}_using_subread.tab ${bamFile}	    
	else
	    conda run -n subread featureCounts -s ${strandedness} -t exon -g gene_id -a ${gtfFile} -o ${tag}_using_subread.tab ${bamFile}
	fi
    fi
    
    # Plot the summary statistics as an histogram to check that the assigned reads are more than 90%
    rm ${tag}_using_subread_summary.ps ${tag}_using_subread_summary.pdf 2> /dev/null
    awk '{name[NR]=$1; v[NR]=$2; tot+=$2}END{for(i in v) print name[i],v[i],tot,v[i]/tot*100}' ${tag}_using_subread.tab.summary | grep -v bam | sort -k 1,1d | sed "s/_/-/g" > ${tag}_using_subread.summary

    assignedFraction=$(grep -iw Assigned ${tag}_using_subread.summary | awk '{printf("%.2f", $4)}')
    echo 0  $assignedFraction >  _tmp_line
    echo 15 $assignedFraction >> _tmp_line
    #cat _tmp_line

    goodness=$(echo assignedFraction | awk '{if($1>0.9){print "This is great!"}else{print "WARNING: Possible error in featureCounts -s parameter!!!"}}')

    cp ${tag}_using_subread.summary countReadPairs_using_subread.summary
    sed -e "s/XXXassignedFractionXXX/${assignedFraction}/g" -e "s/XXXgoodnessXXX/${goodness}/g" /home/common_pipelines/rnaseq/09_countReadPairs_using_subread.gp > _tmp.gp
    export GNUPLOT_PS_DIR=/home/michael.szalay/anaconda3/pkgs/gnuplot-5.0.3-4/share/gnuplot/5.0/PostScript/
    conda run -n TADphys gnuplot _tmp.gp ; rm _tmp.gp
    ps2pdf countReadPairs_using_subread_summary.ps
    rm countReadPairs_using_subread_summary.ps
    mv countReadPairs_using_subread_summary.pdf ${tag}_using_subread_summary.pdf
    if [[ ${checkSingleEnd} -ne 0 ]];
    then
	rm _tmp_line countReadPairs_*
    fi
    rm _tmp_line
	
    echo $assignedFraction | awk '{if($1>65){print "The percentage of assigned reads ("$1"%) is higher than 65%. Hence, the count of the read pairs worked!"}else{print "The percentage of Assigned reads ("$1"\%) is lower than 90\%. Hence, the count of the read pairs did not work. Possible problem!!! Please, check the strandedness of your sample!!!"}}'

    rm doing_countReadPairs
    cd .. # Exit $dir

done

# Generate the table of counts per gene
ls -lrtha ./*/*_using_subread.tab
outFile1=${tag}_using_subread.tab
outFile2=${tag}_using_subread_with_geneLength.tab
rm -fvr _tmp ${outFile1} ${outFile2}
touch ${outFile1}
touch ${outFile2}
for dir in $(ls -1 | grep Rep | grep all);
do
    if [[ ! -d ${dir} ]];
    then
	echo "${dir} is not a directory for a sample"
	continue
    fi    
    
    assayName=$(echo ${dir} | sed "s,_, ,g" | awk '{print $1}')
    targetName=$(echo ${dir} | sed "s,_, ,g" | awk '{print $2}')
    cellName=$(echo ${dir} | sed "s,_, ,g" | awk '{print $3}')
    condition=$(echo ${dir} | sed "s,_, ,g" | awk '{print $4}')
    replicateName=$(echo ${dir} | sed "s,_, ,g" | awk '{print $5}')
    laneName=$(echo ${dir} | sed "s,_, ,g" | awk '{print $6}')
    assembly=$(echo ${dir} | sed "s,_, ,g" | awk '{print $7}')
    tag=$(echo ${dir} | sed "s,_, ,g" | awk '{print $8}')
    #echo ${assayName} ${targetName} ${cellName} ${condition} ${replicateName} ${laneName} ${assembly} ${tag}
    
    paste ${outFile1} <(awk '{if(NR==1){next}; printf("%s\t",$1); for(i=7;i<NF;i++){printf("%s\t",$i)}; printf("%s\n",$NF)}' ./${dir}/*_using_subread.tab | sed -e "s/${assayName}_${targetName}_//g" -e "s/_${laneName}_${assembly}_${tag}.bam//g") >> _tmp ; mv _tmp ${outFile1}
    paste ${outFile2} <(awk '{if(NR==1){next}; printf("%s\t",$1); for(i=7;i<=NF;i++){printf("%s\t",$i)}; printf("%s\n",$6)}' ./${dir}/*_using_subread.tab | sed -e "s/${assayName}_${targetName}_//g" -e "s/_${laneName}_${assembly}_${tag}.bam//g") >> _tmp ; mv _tmp ${outFile2}
done    

# Write the tables of counts and counts+length per gene
#awk '{printf("%s",$1); for(i=2;i<=NF;i+=2){printf("\t%s",$i)}; printf("\n")}' ${outFile1} | sed -e "s/-//g"        > _tmp ; mv _tmp ${outFile1}
awk '{printf("%s",$1); for(i=2;i<=NF;i+=2){printf("\t%s",$i)}; printf("\n")}' ${outFile1}        > _tmp ; mv _tmp ${outFile1}
head ${outFile1}
#awk '{printf("%s",$1); for(i=2;i<=NF;i+=3){printf("\t%s",$i)}; printf("\t%s\n",$3)}' ${outFile2} | sed -e "s/-//g" > _tmp ; mv _tmp ${outFile2}
awk '{printf("%s",$1); for(i=2;i<=NF;i+=3){printf("\t%s",$i)}; printf("\t%s\n",$3)}' ${outFile2} > _tmp ; mv _tmp ${outFile2}
head ${outFile2}    

# Generating samplesFile.txt
#id	condition	type
#WD_2801_Rep1	WD_2801	paired-end
awk 'BEGIN{printf("id\tcondition\ttype\n")}' > samplesFiles.txt
head -1 ${outFile1} | awk '{for(i=2;i<=NF;i+=1){printf("%s\n",$i)}}' | awk -v tag=${tag1} '{print $1,$1,tag}' | sed -e "s/_Rep/ /g" -e "s/-//g" |  awk '{print $1"_Rep"$2,$3,$NF}' | sed -e "s/pairedend/paired-end/g" -e "s/singleend/single-end/g" >> samplesFiles.txt
#rm -fvr _tmp_counts    
