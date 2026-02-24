getwd()
library(Biostrings)
dna_sequences <- readDNAStringSet("sequences.fasta")
dna_sequences
dna.clean <- replaceAmbiguities(dna_sequences, new= "N")
dna.clean
library(msa)
alignment <- msa(dna.clean, method= "Muscle")
alignment
human_mat <- as.matrix(aln)
aln
conserved_cols <- apply(human_mat,2, function(col) length(unique(col)) == 1)
num_conserved <- sum(conserved_cols)
prop_conserved <- num_conserved/ncol(human_mat)
prop_conserved
consensus_seq <- msaConsensusSequence(alignment)
consensus_seq
gc_count <- sum(letterFrequency(DNAStringSet(consensus_seq), letters = c("G","C")))
total_count <- sum(letterFrequency(DNAStringSet(consensus_seq), letters = c("A","T","C","G")))
gc_content <- gc_count/total_count
gc_content
library(ape)
alignmentdnabin <- as.DNAbin(as.matrix(alignment))
distances <- dist.dna(alignmentdnabin, model= "raw")
alignmentdnabin
distances
totaldistances <- apply(as.matrix(distances),1,sum)
totaldistances
mostdiffindex <- which.max(totaldistances)
mostdiffname <- names(dna_sequences)[mostdiffindex]
mostdiffname
mostdiffseq <- dna_sequences[[mostdiffindex]]
mostdiffseq
library(seqinr)
mostdiffseq <- [[mostdiffindex]]
nucvec <- toupper(as.character(unlist(mostdiffseq)))
nucvec
nucvec <- nucvec[1:(length(nucvec) - (length(nucvec) %% 3))]
length(nucvec)
