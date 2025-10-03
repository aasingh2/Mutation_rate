#!/bin/bash

#SBATCH --job-name=Combine_gvcfs
#SBATCH --ntasks 1
#SBATCH --partition=supermem
#SBATCH --cpus-per-task=4

module load bioinfo/gatk/4.1.4.1

out="/data3/projects/mirage/as-internship-files-26022021/hb3Mdh2/gvcfs"
gvcf_map="/data3/projects/mirage/as-internship-files-26022021/hb3Mdh2/scripts/gvcf_map.txt"
interval="/data3/projects/mirage/as-internship-files-26022021/hb3/reference_hb3/intervals.list"

#############################Merging all the vcf files##################################################
gatk GenomicsDBImport \
        --genomicsdb-workspace-path ${out}/my_database \
        --sample-name-map ${gvcf_map}  --tmp-dir=${out} \
        --batch-size 90 \
        --reader-threads 5 \
        --intervals ${interval}


