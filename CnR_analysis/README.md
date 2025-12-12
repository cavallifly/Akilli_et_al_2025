# Dependencies #
We suggest to install [conda](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html) and create an [environment](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html). 

Ensure that you have a working version of R. The scripts in this repository have been tested for version 4.0.5 (2021-03-31).
<!--
Install the following R packages: ggplot2 and DiffBind.
-->
The programs pyGenomeTracks, samtools, sambamba, deeptools, and bowtie2.

# Input data for ChromHMM analysis #
The list of peaks used for the ChromHMM analysis are provided in the peakFiles folder.

<!--
# Input data #
Next, you should download the .fastq files from the GEO entry GSE305456.
From this GEO entry, you can also download directly the .bigWig tracks for all the CnR datasets generated in this work.

Now, you are ready to run the scripts in ./scripts. To do so, you can access the directory CnR_analysis of this repository using
```
cd CnR_analysis
```
and run one after the other the following commands following the messages:
```
scriptsDir=./scripts/

bash    ${scriptsDir}/01_fastqc_analysis.sh                                    &>> 01_fastqc_analysis.out
bash    ${scriptsDir}/02_multiqc_analysis.sh                                   &>> 02_multiqc_analysis.out
bash    ${scriptsDir}/03_align_using_bowtie.sh                                 &>> 03_align_using_bowtie.out
bash    ${scriptsDir}/04_filterReads_using_samtools.sh                         &>> 04_filterReads_using_samtools.out
bash    ${scriptsDir}/05_sortReads_using_samtools.sh                           &>> 05_sortReads_using_samtools.out
bash    ${scriptsDir}/06_dedupReads_using_sambamba.sh                          &>> 06_dedupReads_using_sambamba.out
bash    ${scriptsDir}/07_mergeAndIndexBam_using_samtools.sh                    &>> 07_mergeAndIndexBam_using_samtools.out
bash    ${scriptsDir}/08_generateBigwig_using_deeptools.sh                     &>> 08_generateBigwig_using_deeptools.out
Rscript ${scriptsDir}/09_diffPeak_analysis_on_H3K27me3_tracks_using_diffBind.R &>> 09_diffPeak_analysis_on_H3K27me3_tracks_using_diffBind.out
```


Once the scripts are finished, you will obtain the panels and the data points used to obtain them for all the Figures in the paper.
To obtain the final version of the figures we have visulalized the tracks using [IGV](https://igv.org/) or [pyGenomeTracks](https://pygenometracks.readthedocs.io/en/latest/).
The panels have been assembled using the PowerPoint program.
-->

## Contributions ##
The code in this repository has been developed at the [Cavalli Lab](https://www.igh.cnrs.fr/en/research/departments/genome-dynamics/chromatin-and-cell-biology) with the contributions of Marco Di Stefano, Gonzalo Sabaris, and Giorgio L. Papadopoulos.
