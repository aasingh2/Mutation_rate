#!/bin/bash

#SBATCH --job-name=Combine_gvcfs
#SBATCH --ntasks 1
#SBATCH --partition=supermem
#SBATCH --cpus-per-task=4

module load bioinfo/gatk/4.1.4.1

out="/data3/projects/mirage/as-internship-files-26022021/hb3Mdh2/gvcfs"
gvcf_map="/data3/projects/mirage/as-internship-files-26022021/hb3Mdh2/scripts/gvcf_map.txt"
interval="/data3/projects/mirage/as-internship-files-26022021/hb3/reference_hb3/intervals.list"
ref="/data3/projects/mirage/as-internship-files-26022021/hb3/reference_hb3/PlasmoDB-49_PfalciparumHB3_Genome.fasta"
database="data3/projects/mirage/as-internship-files-26022021/hb3Mdh2/gvcfs/my_database"

################################################Genotyping_and_annotation######################################
gatk GenotypeGVCFs\
    -R ${ref} \
    -V gendb:${database} \
    --max-alternate-alleles 2 \
    -O ${out}/calling_GVCF.vcf

#############################################Annotation###########################################################
java -jar /usr/local/snpEff-4.3/snpEff.jar  -c /usr/local/snpEff-4.3/snpEff.config -no-downstream -no-upstream -onlyProtein Plasmodium_falciparum_hb3  \
${out}/calling_GVCF.vcf > ${out}/annotated.vcf
