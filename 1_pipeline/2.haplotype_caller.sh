#!/bin/bash

#SBATCH --job-name=haplotype_caller
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=supermem
#SBATCH --array=61-90%10


#################Loading_modules###############
module load bioinfo/gatk/4.1.4.1
module load bioinfo/bcftools/1.10.2
module load bioinfo/vcftools/0.1.16
module load bioinfo/tabix/0.2.6
module load bioinfo/samtools/1.10

proc=4
mem=16GB

##################PATH TO FILES ###############
data="/data3/projects/mirage/delly-hb3/sorted-bams"
out="/data3/projects/mirage/as-internship-files-26022021/hb3Mdh2/gvcfs"
ref="/data3/projects/mirage/as-internship-files-26022021/hb3/reference_hb3/PlasmoDB-49_PfalciparumHB3_Genome.fasta"

#####################Slurm array id is the samplename##############
IDX=$SLURM_ARRAY_TASK_ID
name=`sed -n ${IDX}p <samplename.txt`


#Creating Sequence Dictionaries
java -jar /usr/local/picard-tools-2.5.0/picard.jar CreateSequenceDictionary \
     R=${ref} \
     O=${ref}.dict
echo "SequenceDictionary created for sample $INPUTFILE"


#complete the header of each bam file using the AddOrReplaceReadGroups from Picard
java "-Xmx1G" "-Xms1G" "-XX:ConcGCThreads=1" -jar /usr/local/picard-tools-2.5.0/picard.jar AddOrReplaceReadGroups \
        I=${data}/${name}.sorted.bam \
        O=${out}/${name}.fix \
        LB=WGA \
        PL=illumina \
        PU=NA \
        SM=${name}

#Renaming the output file so that it acts as input for the next step
mv ${out}/${name}.fix ${out}/${name}.sorted.bam

#index each fixed bam file using samtools
samtools index ${out}/${name}.sorted.bam
echo "index PROCESSED ${out}/${name}.sorted.bam"

#Using the GATK haplotype Caller
gatk --java-options "-Xmx4G"  HaplotypeCaller \
        -R ${ref} \
        -I ${out}/${name}.sorted.bam \
        --emit-ref-confidence GVCF \
        --pcr-indel-model NONE \
        --sample-ploidy 1 \
        --max-alternate-alleles 2 \
        --output ${out}/${name}.sorted.bam.g.vcf

echo "GVCF produced for ${name}.sorted.bam"
