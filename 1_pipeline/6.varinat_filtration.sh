#!/bin/bash

#SBATCH --job-name=variant_filtration
#SBATCH --ntasks 1
#SBATCH --partition=supermem
#SBATCH --cpus-per-task=4

module load bioinfo/gatk/4.1.4.1
module load bioinfo/bcftools/1.10.2
module load bioinfo/vcftools/0.1.16
module load bioinfo/tabix/0.2.6
module load bioinfo/samtools/1.10
echo "Modules loaded for master script"

out="/data3/projects/mirage/as-internship-files-26022021/hb3Mdh2/gvcfs"
gvcf_map="/data3/projects/mirage/as-internship-files-26022021/hb3Mdh2/scripts/gvcf_map.txt"
interval="/data3/projects/mirage/as-internship-files-26022021/hb3/reference_hb3/intervals.list"
ref="/data3/projects/mirage/as-internship-files-26022021/hb3/reference_hb3/PlasmoDB-49_PfalciparumHB3_Genome.fasta"
#data="data3/projects/mirage/as-internship-files-26022021/hb3Mdh2/gvcfs/annotated.vcf"

###################################################Selecting_variants#################################################
######################################################################################################################
gatk SelectVariants \
    -R ${ref} \
    -V ${out}/annotated.vcf \
    -select-type INDEL\
    --restrict-alleles-to BIALLELIC \
    -O ${out}/all_indels.vcf
echo "SelectVariants for indels PROCESSED"

gatk SelectVariants \
    -R ${ref} \
    -V ${out}/annotated.vcf \
    -select-type SNP \
    --restrict-alleles-to BIALLELIC \
    -O ${out}/all_snps.vcf
echo "SelectVariants for snps PROCESSED"

gatk VariantFiltration \
    -R ${ref} \
    --filter-name Quality -filter "QUAL<100" \
    -V ${out}/all_snps.vcf \
    -O ${out}/quality_snps.vcf
echo "VariantFiltration for snps  PROCESSED"

gatk VariantFiltration \
    -R ${ref} \
    --filter-name Qualtity -filter "QUAL<100" \
    -V ${out}/all_indels.vcf \
    -O ${out}/quality_indels.vcf                         
echo "VariantFiltration for indels PROCESSED"

gatk SelectVariants \
    -R ${ref} \
    -V ${out}/quality_indels.vcf \
    -O ${out}/quality_indels_selected.vcf \
    -select 'vc.isNotFiltered()'
echo "SelectVariants for quality indels PROCESSED"

gatk SelectVariants \
    -R ${ref} \
    -V ${out}/quality_snps.vcf \
    -O ${out}/quality_snps_selected.vcf \
    -select 'vc.isNotFiltered()'
echo "SelectVariants for quality snps PROCESSED"

vcftools --vcf ${out}/quality_snps_selected.vcf --remove-filtered-all --recode --stdout > ${out}/pass_snps.vcf
echo "FIltering PROCESSED for snps"

vcftools --vcf ${out}/quality_indels_selected.vcf --remove-filtered-all --recode --stdout > ${out}/pass_indels.vcf
echo "FIltering PROCESSED for indels"

###################################################VCF Manipulation for AD filtration#####################################
#Extracting chrom ,pos and alt, ref info from pass_snps and pass_indels 
bcftools query -f '%CHROM %POS %REF %ALT\n' ${out}pass_indels.vcf >${out} pass_indels_alt_ref.csv
bcftools query -f '%CHROM %POS %REF %ALT\n' ${out}pass_snps.vcf > ${out}pass_snps_alt_ref.csv

#Extracting Ad info from the vcf
vcftools --vcf ${out}/pass_indels.vcf --extract-FORMAT-info AD --stdout > ${out}/indels_with_ad.vcf
vcftools --vcf ${out}/pass_snps.vcf --extract-FORMAT-info AD --stdout > ${out}/snps_with_ad.vcf

#making a file with the header
cat ${out}/snps_header.csv;head -n1 ${out}/snps_with_ad.vcf| cat >> ${out}/snps_header.csv
cat ${out}/indels_header.csv;head -n1 ${out}/indels_with_ad.vcf| cat >> ${out}/indels_header.csv

#Removing the first line with the header
sed '1d' ${out}/indels_with_ad.vcf > tmpfile; mv tmpfile ${out}/indels_with_ad.vcf_no_head.vcf
sed '1d' ${out}/snps_with_ad.vcf > tmpfile; mv tmpfile ${out}/snps_with_ad.vcf_no_head.vcf

#Sometimes some data points are missin, to correct for that we  add a NA before and after each dot
sed 's/\./NA,NA/g' ${out}/indels_with_ad.vcf_no_head.vcf >> ${out}/indels_with_ad.vcf_no_head_NA.vcf
sed 's/\./NA,NA/g' ${out}/snps_with_ad.vcf_no_head.vcf >> ${out}/snps_with_ad.vcf_no_head_NA.vcf 

#Now we replace each comma with a space
sed 's/,/  /g' ${out}/indels_with_ad.vcf_no_head_NA.vcf >> ${out}/r_ready_indels.vcf
sed 's/,/  /g' ${out}/snps_with_ad.vcf_no_head_NA.vcf >> ${out}/r_ready_snps.vcf
 
mkdir ${out}/r-analysis
mv ${out}/r_ready_indels.vcf ${out}/r-analysis
mv ${out}/r_ready_snps.vcf ${out}/r-analysis
mv ${out}/pass_indels_alt_ref.csv ${out}/r-analysis
mv ${out}/pass_snps_alt_ref.csv ${out}/r-analysis
mv ${out}/snps_header.csv ${out}/r-analysis
mv ${out}/indels_header.csv ${out}/r-analysis

echo "Files are ready for analysis, Good job "
