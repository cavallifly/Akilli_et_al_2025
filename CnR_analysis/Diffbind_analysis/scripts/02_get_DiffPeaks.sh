minFold=1
FDR=0.05
target=$1
outFile=report_${target}_WD_SUMOnub_vs_WD_NOnub_DiffentialPeaks_log2FC_${minFold}_FDR_${FDR}.bed

head -1 report_${target}_WD_SUMOnub_vs_WD_NOnub.tsv > ${outFile}
awk -v FDR=$FDR -v mF=$minFold '{if($11<FDR && sqrt($9*$9)>mF){print $0}}' report_${target}_WD_SUMOnub_vs_WD_NOnub.tsv > ${outFile}
