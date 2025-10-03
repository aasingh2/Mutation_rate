#!/bin/bash
#SBATCH --job-name=mh2d_denovo
#SBATCH -p short


#######################################################################

Rscript ad_table_ompg_functions.R
Rscript hb3_all_gens.R

echo "tables generated"
