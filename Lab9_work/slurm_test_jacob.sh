#!/bin/bash
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -o assembly.test.log
#SBATCH --account=jcorey3320      # <-- your account
#SBATCH --partition=silver

# Load required modules
module load biological/samtools_1.23
module load biological/java
module load biological/perl_5.40   # fix for Bowtie2 bug

# Set project directory and SRR ID
export PROJ_DIR=/export/home/bio_class/Lab_9work
cd $PROJ_DIR
export SRR=SRR5324768

# Create output directories if they don't exist
mkdir -p alignment
mkdir -p variants

# Build Bowtie2 genome index (only needed once)
bowtie2-build ncbi_dataset/ncbi_dataset/data/GCA_900604845.1/GCA_900604845.1_TTHNAR1_genomic.fna genome_index

# Align reads and create sorted BAM
bowtie2 -x genome_index \
-1 fastq/${SRR}_pass_1.fastq.gz \
-2 fastq/${SRR}_pass_2.fastq.gz \
--sensitive-local \
--rg-id ${SRR} --rg SM:${SRR} --rg PL:ILLUMINA \
| samtools view -hb - \
| samtools sort -l 5 -o alignment/${SRR}.bam

# Index the BAM file
samtools index alignment/${SRR}.bam

# Generate consensus sequence
samtools consensus -f fasta -o ${SRR}_consensus.fasta alignment/${SRR}.bam
