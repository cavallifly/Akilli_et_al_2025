for coolFile in $(ls -1 *.cool | grep _1000bp)
do
    (
	outFile=${coolFile%_1000bp.cool}.mcool
	if [[ ! -e ${outFile} ]];
	then
	    echo $coolFile $outFile 2>> 02_cool2mcool_${outFile%.mcool}.out	
	    
	    touch ${outFile}
	    conda run -n cooler cooler zoomify -r 1000,2000,5000,10000,15000,20000,40000,50000,100000 ${coolFile} -o ${outFile} --balance 2>> 02_cool2mcool_${outFile%.mcool}.out	
	fi
    ) &
done


