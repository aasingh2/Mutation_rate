#!/bin/bash

#SBATCH --job-name=hb3MDH2_mapping
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=long
#SBATCH --array=1-90%10


#################Loading_modules###############
module load bioinfo/bwa/0.7.17
module load bioinfo/samtools/1.10

proc=4
mem=16GB

##################PATH TO FILES ###############
data="/data3/projects/mirage/as-internship-files-26022021/hb3Mdh2/rawfastqs"
out="/data3/projects/mirage/as-internship-files-26022021/hb3Mdh2/bams"
ref="/data3/projects/mirage/as-internship-files-26022021/hb3/reference_hb3/PlasmoDB-49_PfalciparumHB3_Genome.fasta"

IDX=$SLURM_ARRAY_TASK_ID
name=`sed -n ${IDX}p <samplename.txt`

###map paired-end read using BWA_MEM (default parameters)
bwa mem -M ${ref} ${data}/${name}_1.fastq.gz ${data}/${name}_2.fastq.gz  -t ${proc} > ${TMPDIR}/${name}.sam

echo "mapping done"

####convert sam to bam and sort it/ keep the sort_bam
samtools view -S -b ${TMPDIR}/${name}.sam  > ${TMPDIR}/${name}.bam
echo " converted sam to bam"
samtools sort ${TMPDIR}/${name}.bam  -o ${out}/${name}.sorted.bam
echo "sorting done"
samtools index ${out}/${name}.sorted.bam
echo "indexing done"

