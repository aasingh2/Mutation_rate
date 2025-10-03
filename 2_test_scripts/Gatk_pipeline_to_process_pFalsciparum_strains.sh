#!/bin/bash
#SBATCH --job-name=gatk_pipeline
#SBATCH -p supermem
#SBATCH --mail-user=aakanksha.singh@etu.umontpellier.fr
#SBATCH --mail-type=ALL
########################################################################################################
##############################################Loading modules############################################
module load bioinfo/gatk/4.1.4.1
module load bioinfo/bcftools/1.10.2
module load bioinfo/vcftools/0.1.16
module load bioinfo/tabix/0.2.6
module load bioinfo/samtools/1.10
echo "Modules loaded for master script"

#############################################Required Global variables####################################
#Path to files TO BE FILLED

#Make sure that the BAMDIR contain just the .bam and its respective indexed file and nothing else
BAMDIR=p_falciparum/hb3/bam_files/bambai
REFERENCE=p_falciparum/hb3/reference_hb3/PlasmoDB-49_PfalciparumHB3_Genome.fasta
REF_FOLDER=p_falciparum/hb3/reference_hb3
SCRIPT_PATH=p_falciparum/hb3/scripts/gatk
#REF_FOLDER contains the path to the folder that has the files from Miles 2016 (Just for 3D7)
#KNOWNSITES=/data3/projects/mirage/p_falciparum/3D7_comp/Known_sites ( Just for 3D7)
GVCF_MAP="$REF_FOLDER"/gvcf_map_hb3.tab #gvcf map file is a tab delimited file containing the file name and the path to the respective .g.vcf.gz file
BATCH_SIZE=83  #batch sixe is the total number of the gvcf files that are to be genotyped
INTERVAL_FILE="$REF_FOLDER"/intervals.list #Path to the interval file containg that list of the chromosomes ( This is nothing but a simple text file containing the list of the chromosomes)
SNPEFF_NAME=Plasmodium_falciparum_hb3

###########################################################################################################
##########################################################################################################
cd "$BAMDIR"
mkdir output_files
OUTPUTDIR="$BAMDIR"/output_files
FILES="$BAMDIR"/*sorted.bam
echo "Path variables defined and the working directory set"
SAMPLES=($(ls -1 $FILES))
source "$SCRIPT_PATH"/producing_a_vcf_from_bam_tested.sh

###########################Starting the loop to call the gatk_function for each file#######################
length_of_array=${#SAMPLES[@]}
start_index=0
number_of_threads=8
sequential_per_process=0
iter=0
while [ $start_index -le $(($length_of_array-$sequential_per_process)) ];do
	sequential_per_process=$(bc <<< $(($length_of_array-$start_index))/$(($number_of_threads-$iter)))
	end_index=$((start_index+$sequential_per_process))
	echo "New batch: ${SAMPLES[$start_index]} to ${SAMPLES[$(($end_index-1))]}"
	for (( i = $start_index; i < $end_index; i++ ));do
		producing_a_vcf_from_bam ${SAMPLES[$i]}
		echo ${SAMPLES[$i]}
	done &
	start_index=$end_index
	((iter++))
done
wait

############################Zipping and Indexing all the GVCFS#########################################
cd "$OUTPUTDIR"
GVCFS="$OUTPUTDIR"/*.bam.g.vcf
for f in $GVCFS
do
 bgzip $f
 tabix -p vcf $f.gz
done
echo "compression and indexation VCF PROCESSED"

#############################Merging all the vcf files##################################################
gatk --java-options "-Xmx4g -Xms4g"  GenomicsDBImport \
        --genomicsdb-workspace-path "$OUTPUTDIR"/my_database \
        --sample-name-map "$GVCF_MAP" --tmp-dir="$OUTPUTDIR" \
        --batch-size $BATCH_SIZE \
        --reader-threads 5 \
        --intervals "$INTERVAL_FILE"
###########################################VARIANT CALLING####################################################
gatk GenotypeGVCFs\
    -R "$REFERENCE" \
    -V gendb://my_database \
    --max-alternate-alleles 2 \
    -O "$OUTPUTDIR"/calling_GVCF.vcf

#############################################Annotation###########################################################
#################################################################################################################

java -jar /usr/local/snpEff-4.3/snpEff.jar  -c /usr/local/snpEff-4.3/snpEff.config -no-downstream -no-upstream -onlyProtein $SNPEFF_NAME \
 "$OUTPUTDIR"/calling_GVCF.vcf > "$OUTPUTDIR"/annotated.vcf

#Just for 3D7
#bcftools annotate \
#   -a "$REF_FOLDER"/regions.tab.gz \
#   -h "$REF_FOLDER"/header.hdr   \
#   -Ov \
#   -o "$OUTPUTDIR"/coreregions.vcf \
#   -c CHROM,FROM,TO,RegionType "$OUTPUTDIR"/annotated.vcf \

###################################################Selecting_variants#################################################
######################################################################################################################
gatk SelectVariants \
    -R "$REFERENCE" \
    -V "$OUTPUTDIR"/annotated.vcf \
    -select-type INDEL\
    --restrict-alleles-to BIALLELIC \
    -O "$OUTPUTDIR"/indels.vcf
echo "SelectVariants for indels PROCESSED"

gatk SelectVariants \
    -R "$REFERENCE" \
    -V "$OUTPUTDIR"/annotated.vcf.vcf \
    -select-type SNP \
    --restrict-alleles-to BIALLELIC \
    -O "$OUTPUTDIR"/snps.vcf
echo "SelectVariants for snps PROCESSED"

gatk VariantFiltration \
    -R "$REFERENCE" \
    --filter-name Quality -filter "QUAL<100" \
    -V "$OUTPUTDIR"/snps.vcf \
    -O "$OUTPUTDIR"/core_snps.vcf
echo "VariantFiltration for snps  PROCESSED"

gatk VariantFiltration \
    -R "$REFERENCE" \
    --filter-name Qualtity -filter "QUAL<100" \
    -V "$OUTPUTDIR"/indels.vcf \
    -O "$OUTPUTDIR"/core_indels.vcf                         
echo "VariantFiltration for indels PROCESSED"

gatk SelectVariants \
    -R "$REFERENCE" \
    -V "$OUTPUTDIR"/core_indels.vcf \
    -O "$OUTPUTDIR"/core_indels_selected.vcf \
    -select 'vc.isNotFiltered()'
echo "SelectVariants for core indels PROCESSED"

gatk SelectVariants \
    -R "$REFERENCE" \
    -V "$OUTPUTDIR"/core_snps.vcf \
    -O "$OUTPUTDIR"/core_snps_selected.vcf \
    -select 'vc.isNotFiltered()'
echo "SelectVariants for core snps PROCESSED"

vcftools --vcf "$OUTPUTDIR"/core_snps_selected.vcf --remove-filtered-all --recode --stdout > "$OUTPUTDIR"/pass_snps.vcf
echo "FIltering PROCESSED for snps"

vcftools --vcf "$OUTPUTDIR"/core_indels_selected.vcf --remove-filtered-all --recode --stdout > "$OUTPUTDIR"/pass_indels.vcf
echo "FIltering PROCESSED for indels"

###################################################VCF Manipulation for AD filtration#####################################
#Extracting chrom ,pos and alt, ref info from pass_snps and pass_indels 
bcftools query -f '%CHROM %POS %REF %ALT\n' pass_indels.vcf > pass_indels_alt_ref.csv
bcftools query -f '%CHROM %POS %REF %ALT\n' pass_snps.vcf > pass_snps_alt_ref.csv

#Extracting Ad info from the vcf
vcftools --vcf "$OUTPUTDIR"/pass_indels.vcf --extract-FORMAT-info AD --stdout > "$OUTPUTDIR"/indels_with_ad.vcf
vcftools --vcf "$OUTPUTDIR"/pass_snps.vcf --extract-FORMAT-info AD --stdout > "$OUTPUTDIR"/snps_with_ad.vcf

#making a file with the header
cat "$OUTPUTDIR"/snps_header.csv;head -n1 "$OUTPUTDIR"/snps_with_ad.vcf| cat >> "$OUTPUTDIR"/snps_header.csv
cat "$OUTPUTDIR"/indels_header.csv;head -n1 "$OUTPUTDIR"/indels_with_ad.vcf| cat >> "$OUTPUTDIR"/indels_header.csv

#Removing the first line with the header
sed '1d' "$OUTPUTDIR"/indels_with_ad.vcf > tmpfile; mv tmpfile "$OUTPUTDIR"/indels_with_ad.vcf_no_head.vcf
sed '1d' "$OUTPUTDIR"/snps_with_ad.vcf > tmpfile; mv tmpfile "$OUTPUTDIR"/snps_with_ad.vcf_no_head.vcf

#Sometimes some data points are missin, to correct for that we  add a NA before and after each dot
sed 's/\./NA,NA/g' "$OUTPUTDIR"/indels_with_ad.vcf_no_head.vcf >> "$OUTPUTDIR"/indels_with_ad.vcf_no_head_NA.vcf
sed 's/\./NA,NA/g' "$OUTPUTDIR"/snps_with_ad.vcf_no_head.vcf >> "$OUTPUTDIR"/snps_with_ad.vcf_no_head_NA.vcf 

#Now we replace each comma with a space
sed 's/,/  /g' "$OUTPUTDIR"/indels_with_ad.vcf_no_head_NA.vcf >> "$OUTPUTDIR"/r_ready_indels.vcf
sed 's/,/  /g' "$OUTPUTDIR"/snps_with_ad.vcf_no_head_NA.vcf >> "$OUTPUTDIR"/r_ready_snps.vcf
 
mkdir "$OUTPUTDIR"/r-analysis
mv r_ready_indels.vcf "$OUTPUTDIR"/r-analysis
mv r_ready_snps.vcf "$OUTPUTDIR"/r-analysis
mv pass_indels_alt_ref.csv "$OUTPUTDIR"/r-analysis
mv pass_snps_alt_ref.csv "$OUTPUTDIR"/r-analysis
mv snps_header.csv "$OUTPUTDIR"/r-analysis
mv indels_header.csv "$OUTPUTDIR"/r-analysis

echo "Files are ready for analysis, Good job "





