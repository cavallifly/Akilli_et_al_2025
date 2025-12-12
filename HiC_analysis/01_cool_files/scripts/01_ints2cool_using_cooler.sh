
inFiles=$(ls -1 ../01_ints_files/*ints)
#echo $inFiles

assembly=dm6


for inFile in $inFiles ;
do
    (
	#echo $inFile

	outName=$(echo $inFile | sed -e "s/\.ints//g" -e "s/\.noDup//g" -e "s,files/, ,g" -e "s,/,,g" | awk '{print $2}')

	for resolution in 1000 20000 ;
	do
	    
	    outFile=${outName}_${resolution}bp.cool
	    echo $inFile $outFile >> 01_ints2cool_using_cooler_${outName}.out
	    if [[ ! -e ${outFile} ]];
	    then
		touch ${outFile}
		cooler cload pairs --assembly ${assembly} --chrom1 1 --pos1 2 --chrom2 4 --pos2 5 /home/Programs/chrom_sizes_dm6_higlass.txt:${resolution} <( awk '{if($1=="chrom1"){next}; for(i=0;i<$NF;i++){printf("%s\t%d\t%d\t%s\t%d\t%d\t1\n",$1,$2,$3,$4,$5,$6)}}' ${inFile}) ${outFile} &>> 01_ints2cool_using_cooler_${outName}.out
	    fi
	done
    ) &
done
wait

