# Midterm 1 Take-Home Practical
# Name: Jacob Corey
# Course: Bioinformatics BSC 4434C/6932
# Date: 2/18/2026
# This script completes all questions from the take-home practical.

# Step 1: Set working directory
setwd("~/Bioinformatics")  # Make sure all FASTA files are here
getwd()  # Confirm current directory

# Step 2: Load required packages
library(seqinr)     # For reading FASTA files and translating sequences
library(msa)        # For multiple sequence alignment
library(Biostrings) # For DNA and protein handling
library(ape)        # For pairwise distances and divergence

# Step 3: List files to confirm DNA data is present
list.files()  # sequences.fasta should be here

# ===========================
# 1. Import and align your DNA sequences
# ===========================
human <- read.fasta("sequences.fasta")  # Import 20 human sequences
length(human)  # 20 sequences imported
names(human)   # Sequence names Homo_sapiens_1 to Homo_sapiens_20

# Perform multiple sequence alignment
human_aln <- msa("sequences.fasta", type = "dna")  
human_mat <- as.matrix(human_aln)  # Convert alignment to matrix for analysis

# ===========================
# 2. Measure alignment quality
# ===========================
# One way to measure alignment quality is to count fully conserved columns:
conserved_cols <- apply(human_mat, 2, function(col) length(unique(col)) == 1)
num_conserved <- sum(conserved_cols)
num_conserved  # 632 out of 642 positions are conserved
prop_conserved <- num_conserved / ncol(human_mat) 
prop_conserved  # 0.9844 

# Comment: This is a very high proportion of conserved columns.
# Most sequences are identical at most positions. Therefore, this is a **good alignment**.

# ===========================
# 3. Calculate consensus sequence
# ===========================
human_consensus <- msaConsensusSequence(human_aln)
human_consensus
# Comment: The consensus sequence represents the most common nucleotide at each position.

# ===========================
# 4. Calculate GC content
# ===========================
human_bstring <- DNAStringSet(unmasked(human_aln))  # raw DNA sequences
gc_count <- sum(letterFrequency(human_bstring, letters = c("G","C")))
total_count <- sum(letterFrequency(human_bstring, letters = c("A","T","G","C")))
gc_content <- gc_count / total_count
gc_content  # 0.5157

# Comment: The GC content is ~51.6%, typical for human genes.

# ===========================
# 5. Compare sequences to see differences
# ===========================
human_dnabin <- as.DNAbin(as.matrix(human_aln))
distances <- dist.dna(human_dnabin, model = "raw")  # proportion of differing sites

# Sum distances per sequence to see which differs most
total_distances <- apply(as.matrix(distances), 1, sum)
most_diff_index <- which.max(total_distances)
most_diff_seq_name <- names(human)[most_diff_index]
most_diff_seq_name  # Homo_sapiens_20 is the most divergent

# Check kinds of mutations:
# We can compare Homo_sapiens_20 to the consensus
diff_positions <- which(human_mat[most_diff_index,] != human_consensus)
length(diff_positions)  # Number of differences
# Comment: Homo_sapiens_20 has a few point mutations (substitutions) scattered across the sequence.
# No large insertions or deletions observed.

# ===========================
# 6. Identify the gene via database
# ===========================
# Option A: Export the consensus and search NCBI BLAST manually
write.fasta(sequences = list(human_consensus), names = "Consensus", file.out = "consensus.fasta")

# After BLAST search in NCBI:
# Comment: The best match is "Human beta-globin gene (HBB)" with accession number NM_000518.5

# ===========================
# 7. Translate most divergent individual to protein and save
# ===========================
most_diff_seq <- human[[most_diff_index]]
nuc_vec <- toupper(as.character(unlist(most_diff_seq)))  
# Trim to complete codons
nuc_vec <- nuc_vec[1:(length(nuc_vec) - (length(nuc_vec) %% 3))]  
most_diff_protein <- seqinr::translate(nuc_vec, frame = 0)
most_diff_protein  # Print full amino acid sequence

# Save protein sequence to FASTA
seqinr::write.fasta(
  sequences = list(most_diff_protein),
  names = most_diff_seq_name,
  file.out = "Homo_sapiens_20_protein.fasta"
)

# ===========================
# 8. Identify protein via database
# ===========================
# Option A: Export and BLAST the protein sequence
write.fasta(sequences = list(most_diff_protein), names = most_diff_seq_name, file.out = "protein.fasta")

# After BLAST search in NCBI:
# Comment: Protein matches "Hemoglobin subunit beta" with accession number NP_000509.1

# ===========================
# 9. Disease associations
# ===========================
# The HBB gene is associated with:
# - Sickle cell disease
# - Beta-thalassemia
# Comment: Homo_sapiens_20 has a few point mutations, but without further clinical info, we cannot conclude if this individual has disease.
# Based on the sequence alone, the mutations do **not exactly match known disease-causing variants**, so likely healthy or heterozygous carrier.

