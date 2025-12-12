# INPUT: Output file of the subread alignment
# OUTPUT: Statistics on the mapped PE reads in a table with 5 columns:
# Sample                                 	Total_fragments	Mapped         	Perc_mapped    	Uniquely_mapped	Unmapped       	Properly_paired	Not_properly_paired

outFile=MappingStatistics.tsv

if [[ ! -e ${outFile} ]];
then
    # Creating the output table and writing the header
    echo "#Sample Total_fragments Mapped Perc_mapped Uniquely_mapped Unmapped Properly_paired Not_properly_paired" | awk '{printf("%-40s\t",$1); for(i=2;i<NF;i++){printf("%-15s\t",$i)}; printf("%-15s\n",$NF)}' > ${outFile}
fi

grep -i "Total\|Mapped\|Properly\|Analysing" 01_alignReads_using_subread.log | sed -e "s,||,,g" -e "s,Analysing,,g" -e "s,:,,g" -e "s,(,,g" -e "s,),,g" -e "s,%,,g" | awk '{if(NR==1){printf("%-40s",$1); next}; if($1=="Mapped"){printf("%-15s\t%-15s",$2,$3); next} if(NF>1){printf("\t%-15s",$NF)}; if(NF==1){printf("\n%-40s",$1)};}END{printf("\n")}' >> ${outFile}
#Analysing rnaseq_RNA_WD_NOnub_Rep1_L1_dm6_NA
#||             Total fragments : 24,040,899                                   ||
#||                      Mapped : 23,482,999 (97.7%)                           ||
#||             Uniquely mapped : 23,482,999                                   ||
#||                    Unmapped : 557,900                                      ||
#||             Properly paired : 20,330,806                                   ||
#||         Not properly paired : 3,152,193                                    ||

cat ${outFile}
