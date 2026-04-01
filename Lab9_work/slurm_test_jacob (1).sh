#!/bin/bash
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -o assembly.version2test.log
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
#mkdir -p alignment
#mkdir -p variants
'''
java -jar /export/share/software/biological/picard/picard.jar \
	CreateSequenceDictionary \
    REFERENCE=genome/Thermus_thermophilus_TTHNAR1.fa \
    OUTPUT=genome/Thermus_thermophilus_TTHNAR1.dict


# worked
/export/share/software/biological/bowtie2-2.4.2-sra-linux-x86_64/bowtie2-build \
	genome/Thermus_thermophilus_TTHNAR1.fa \
	genome/Thermus_thermophilus_TTHNAR1

# run bowtie2 and samtools separately for testing. bowtie2 failed.
/export/share/software/biological/bowtie2-2.4.2-sra-linux-x86_64/bowtie2 -x \
		genome/Thermus_thermophilus_TTHNAR1 \
        -1 fastq/${SRR}_pass_1.fastq.gz \
        -2 fastq/${SRR}_pass_2.fastq.gz --sensitive-local \
        --rg-id ${SRR} --rg SM:${SRR} --rg PL:ILLUMINA \
        > alignment/${SRR}.sam
'''
# worked
samtools view -hb alignment/${SRR}.sam -o alignment/${SRR}.bam
samtools sort -l 5 alignment/${SRR}.bam -o alignment/${SRR}.bam



# Generate consensus sequence
samtools consensus -f fasta alignment/${SRR}.bam -o ${SRR}_consensus.fasta

