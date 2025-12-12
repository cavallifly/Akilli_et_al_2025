assayName=hic
assembly=dm6
author=NA
condition=$1
cellName=WD
blackListFile=/home/Programs/Blacklists/${assembly}-blacklist.v2.bed

mkdir -p 01_optimizeW_TopDom
cd 01_optimizeW_TopDom

for resolution in 5000 10000 20000 40000 50000
do
    echo "1 - Do the analysis at ${resolution}bp for ${condition} merge for different values of windows:"
    if [[ ! -e TopDom_domains_${assayName}_${cellName}_${condition}_merge_${assembly}_${author}_window_6_at_${resolution}bp.tsv ]];
    then
	echo "The TopDom analysis at ${resolution}bp is still not done!"
	conda run -n topdom Rscript ../scripts/01_bash_optimizeW_TopDom.R ${condition} ${resolution} #2> 01_detectDomain_using_TopDom_at_${resolution}bp.out
    else
	echo "The TopDom analysis at ${resolution}bp is already done!"
    fi
    echo ""
    
    for file in $(ls -1 TopDom_domains_${assayName}_${cellName}_${condition}_merge_${assembly}_${author}_window_*_at_${resolution}bp.tsv) ;
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
	
	echo "3 - Computing the correlation Pearson correlation coefficients between bins in the same TAD"
	if [[ ! -e  ${fileFilt%.bed}_cisDomainPCC.tsv ]];
	then
	    sed -e "s/XXXresolutionXXX/${resolution}/g" -e "s/XXXregionsFileXXX/${fileFilt}/g" -e "s,XXXmcoolFileXXX,../../01_cool_files/${assayName}_${cellName}_${condition}_merge_${assembly}_${author}.mcool,g" ../scripts/01_cooler_computePCC.py > _tmp.py
	    conda run -n cooler python _tmp.py | awk -v f=${fileFilt%.bed} '{print f,$0}' > ${fileFilt%.bed}_cisDomainPCC.tsv
	    rm _tmp.py
	fi
	echo
    done # Close cycle over $file
    echo ""
    
    echo "4 - Select the best value of window computing the PCC as in the original paper"
    outFile=PCCbins_for_window_optimization_${assayName}_${cellName}_${condition}_merge_${assembly}_${author}_at_${resolution}bp.pdf
    if [[ ! -e ${outFile} ]];
    then
	ls -lrtha TopDom_domains_${assayName}_${cellName}_${condition}_merge_${assembly}_${author}_window_*_at_${resolution}bp_filtered_cisDomainPCC.tsv
	cat TopDom_domains_${assayName}_${cellName}_${condition}_merge_${assembly}_${author}_window_*_at_${resolution}bp_filtered_cisDomainPCC.tsv | sed -e "s/TopDom_domains_${assayName}_${cellName}_${condition}_merge_${assembly}_${author}_window_//g" -e "s/_at_${resolution}bp_filtered//g" | grep chr | awk '{print $1,$(NF-1)}' | sort -k 1,1n > _tmp.tab
	conda run -n DEseq2 Rscript ../scripts/01_plotPCC_for_w_optimization.R
	rm _tmp.tab
	mv PCC.pdf ${outFile}
    fi
done

cd ..
