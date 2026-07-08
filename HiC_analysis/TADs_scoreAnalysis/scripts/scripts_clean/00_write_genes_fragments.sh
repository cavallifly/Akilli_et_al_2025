rm -fvr *~
inDomains=$(ls -1 TopDom_domains_*Original*bed)

inGenes=$(ls -1 ../../../rnaseq/scripts/dmel-all-r6.36.gtf)
echo "File with the dm6 genes ${inGenes}"

outDomains=$(echo $inDomains | sed "s/_Original//g")
echo "File with the Top-Dom domains $inDomains"
awk '{if(NF==6)printf("%s %s %s %s%d\n",$2,$3,$4,$6,++h[$NF])}' ${inDomains} | grep -v "chrom" > ${outDomains}

outGenes=all_genes_dmel-all-r6.36.gtf
rm -fvr ${outGenes}
echo "List of genes"
for chrom in 2L 2R 3L 3R 4 X;
do   
    awk -v c=$chrom '{if($1==c){print $0}}' ${inGenes} | grep -w gene | awk '{print "chr"$1,$4,$5,$12}' | sed -e "s/;//g" -e "s/\"//g" -e "s/:/-/g" >> ${outGenes}
    echo $chrom $(wc -l ${outGenes})
done

outFragments=all_gene_fragments_dmel-all-r6.36.gtf
rm -fvr ${outFragments}
echo "List of genes and domains"
for chrom in chr2L chr2R chr3L chr3R chr4 chrX;
do
    #(
    grep -w ${chrom} ${outDomains} ${outGenes} | awk '{for(i=$2;i<=$3;i++){print i,$NF}}' | awk '{h[$1]++; if(h[$1]==1){name[$1]=$NF}else{name[$1]=name[$1]"_"$NF}}END{for(i in h) print i,name[i]}' | sort -k 1,1n | awk '{if(NR==1){fn=$1;np=$1;gp=$2}else{if($2==gp && $1==np+1){np=$1}else{print fn,np,gp; fn=$1;np=$1;gp=$2}}}END{print fn,np,gp}' | sort -k 1,1n | awk -v c=${chrom} '{if(NF>0 && $1!=$2)print c,$0}' >> ${outFragments}
    echo $chrom $(wc -l ${outFragments})
    #) #&
done

cat ${outFragments} | awk '{h[$NF]++; print $0":"h[$NF]}' > _tmp_${outFragments} ; mv _tmp_${outFragments} ${outFragments}
wc -l ${outFragments}
wait
