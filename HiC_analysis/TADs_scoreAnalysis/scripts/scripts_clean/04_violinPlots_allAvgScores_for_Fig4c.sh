#!/bin/bash

#SBATCH --job-name violinPlots
#SBATCH -n 1                    # Number of cores. For now 56 is the number max of core available
#SBATCH -t 4-00:00              # Runtime in D-HH:MM
#SBATCH --mem=16G
#SBATCH -o 04_violinPlots_allAvgScores_for_Fig4c.out # File to which STDOUT will be written
#SBATCH -e 04_violinPlots_allAvgScores_for_Fig4c.out # File to which STDERR will be written 

Rscript scripts_clean/04_violinPlots_allAvgScores_for_Fig4c.R


