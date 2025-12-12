
assayName=hic
assembly=dm6
author=NA
condition=$1
cellName=WD
blackListFile=/home/Programs/Blacklists/${assembly}-blacklist.v2.bed

for resolution in 5000 ; 
do
    echo "1 - Do the analysis at ${resolution}bp for ${assayName} merge for different values of windows:"
    if [[ ! -e TopDom_domains_${assayName}_${cellName}_${condition}_merge_${assembly}_${author}_window_4_at_${resolution}bp.tsv ]];
    then
	echo "The TopDom analysis at ${resolution}bp is still not done!"
	conda run -n topdom Rscript ./scripts/01_detectDomain_using_TopDom.R ${condition} ${resolution}
    else
	echo "The TopDom analysis at ${resolution}bp is already done!"
    fi
    echo ""

    for file in $(ls -1 TopDom_domains_${assayName}_*_${assembly}_${author}_window_*_at_${resolution}bp.tsv) ;
    do
        echo $file
        echo "Initial number of domains $(cat ${file} | grep chr | awk '{print $2,$3,$4}' | wc -l)"

        echo "2 - Remove TopDom domains in blacklisted regions (e.g., pericentromeric regions)"
        fileFilt=${file%.tsv}_filtered.bed

        if [[ ! -e ${fileFilt} ]];
        then
            echo $fileFilt
            cat ${file} | grep chr | awk '{print $2,$3,$4}' | head
            cat ${file} | grep chr | awk '{print $2,$3,$4}' | tail
            awk '{if($4=="BL"){n++;chrom[n]=$1;start[n]=$2;end[n]=$3;next}else{flag=0;for(i=1;i<=n;i++){if($1==chrom[i]){for(j=start[i];j<=end[i];j+=100){if($2<=j && j<=$3){flag=1}}}};if(flag==0){print $0}}}' <(awk '{print $1,$2,$3,"BL"}' ${blackListFile}) <(cat ${file} | grep chr | awk '{print $2,$3,$4}') > ${fileFilt}
            head ${fileFilt}
        fi
        echo "Number of domains after filtering blacklisted regions $(wc -l ${fileFilt} | awk '{print $1}')"

    done # Close cycle over $file
    echo ""
    continue
    
done
