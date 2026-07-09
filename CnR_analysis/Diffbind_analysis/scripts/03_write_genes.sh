# genes #
for dir in $(ls -1 | grep genes);
do
    if [[ ! -d ${dir} ]];
    then
	continue
    fi
    cd $dir
    pwd

    inFile=$(ls -1 report_*_WD_SUMOnub_vs_WD_NOnub.tsv 2> /dev/null)
    outFile=${inFile%.tsv}_withGeneNames.tsv
    echo $inFile $outFile

    if [[ $inFile != "" ]];
    then
	head -1 ${inFile} > $outFile
	
	awk '{if(NF==4){name[$1"_"$2]=$4}else{ind=$1"_"$2; printf("%s\t%s\n",name[ind],$0)}}' ../DEseq2_results_WD_SUMOnub_vs_WD_NOnub_genesForDiffBind.tsv ${inFile} >> ${outFile}
    fi
    awk '{if(NF==15){print $0}}' ${outFile}
    echo
    wc -l $inFile $outFile   
    echo
    cd ..   
done

# TSSs #
for dir in $(ls -1 | grep TSSs);
do
    if [[ ! -d ${dir} ]];
    then
	continue
    fi
    cd $dir
    pwd

    inFile=$(ls -1 report_*_WD_SUMOnub_vs_WD_NOnub.tsv 2> /dev/null)
    outFile=${inFile%.tsv}_withGeneNames.tsv
    echo $inFile $outFile

    if [[ $inFile != "" ]];
    then
	head -1 ${inFile} > $outFile
	
	awk '{if(NF==4){name[$1"_"$2]=$4}else{ind=$1"_"$2; printf("%s\t%s\n",name[ind],$0)}}' ../DEseq2_results_WD_SUMOnub_vs_WD_NOnub_TSSsForDiffBind.tsv ${inFile} >> ${outFile}
    fi
    awk '{if(NF==15){print $0}}' ${outFile}
    echo
    wc -l $inFile $outFile   
    echo
    cd ..   
done
