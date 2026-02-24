getwd()
library(Biostrings)
library(msa)
library(seqinr)
library(ape)
dna <- readDNAStringSet("sequences.fasta")
dna
class(dna)
dna.clean <- replaceAmbiguities(dna,new = "N")
dna.clean
aln <- msa(dna.clean,method = "Muscle")
aln
class(aln) 
?class
human_mat <- as.matrix(aln)
conserved_cols <- apply(human_mat,2,function(col) length(unique(col))==1)
prop_conserved <- sum(conserved_cols) / ncol(human_mat)
prop_conserved
consensus_seq <- msaConsensusSequence(aln)
consensus_seq
gc_count <- sum(letterFrequency(dna.clean,letters = c("G","C")))
?sum
total_count <- sum(letterFrequency(dna.clean,letters = c("A","T","G","C")))
gc_content <- gc_count/total_count
gc_content
gc_count
dna.dnabin <- as.DNAbin(as.matrix(aln))
distances <- dist.dna(dna.dnabin,model = "raw")
?raw
distances
writeXStringSet(dna.clean,"sequences_for_blast.fasta")
total_distances <- apply(as.matrix(distances),1,sum)
most_diff_index <- which.max(total_distances)
most_diff_seq <- dna[[most_diff_index]]
most_diff_seq
most_diff_seq_name
seq_char <- as.character(most_diff_seq)
nuc_vec <- strsplit(seq_char,"")[[1]]
nuc_vec <- toupper(nuc_vec)
nuc_vec <- nuc_vec[1:(length(nuc_vec) - length(nuc_vec) %% 3)]
most_diff_protein <- seqinr::translate(nuc_vec,frame = 0)
most_diff_protein
protein_string <- paste(most_diff_protein,collapse = "")
seqinr::write.fasta(sequences = list(protein_string),names = "mostdivergentprotein",file.out = "most_divergent_protein.fasta")
getwd()
list.files()
