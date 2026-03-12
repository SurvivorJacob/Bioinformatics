getwd()
# We use this to check our current working directory to know where are files are being read and saved.
.libPaths("C:/Users/jacob/AppData/Local/R/win-library/4.5")
# This informs R to install the proper packages onto my computer.
if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#This allows us to confirm that BiocManager is properly installed which is needed to install Bioconducter packages in the future.
if(!requireNamespace("GenomicAlignments", quietly = TRUE)) BiocManager::install("GenomicAlignments", update = FALSE, ask = FALSE)
# By installing GenomicAlignments we now have the ability to load UnitprotR for our lab.
packages <- c("UniprotR", "protti", "r3dmol", "Biostrings")
# In one row we gather every specific package for this lab so it can all load at once.
for(p in packages){
  if(!requireNamespace(p, quietly = TRUE)) install.packages(p, dependencies = TRUE)
  # This is just a quick check to see if everything so far has been loaded properly.
  library(Biostrings)
  library(UniprotR)
  library(protti)
  library(r3dmol) }
  # These are the commands that we have been installing/loading so now they are properly used in the R system. For example, the 3d model is now alert and ready when we need it. 
dna_seq <- readDNAStringSet("E_coli4.fasta")
# We can properly view/read our Lab 6 DNA data fasta.
dna_char <- as.character(dna_seq)
# We want to check and clean our DNA incase of any errors.
dna_clean <- gsub("[^ATCGatcg]", "N", dna_char)
#We just want to go through any that are not an A,T,C, or G and replace it with an N for transitional error cleanups.
dna_clean <- DNAStringSet(dna_clean)
# We properly put our cleaned DNA back in the active Stringset for reading.
protein_seq <- translate(dna_clean)
# We now can translate our clean DNA into the corresponding Protein sequence.
writeXStringSet(protein_seq, "E_coli4_protein.fasta")
# We want to save our protein sequence into a fasta file for the future of this lab.
protein_check <- readLines("E_coli4_protein.fasta")
# This role is just to check if our protein sequence fasta has been read right.
head(protein_check)
# This shows our first few lines of read for protein sequence output.
accessions <- c("P0A799","P08839")
# After looking through Unitprot we can use our backup sources which we have.
accessions <- as.character(accessions)
# We now can save these into character strings for R.
urls <- paste0("https://www.uniprot.org/uniprot/", accessions, ".fasta")
# We now have active urls from Unitprot for viewing in Github or R itself.
print(urls)
# We want to check the urls are correct.
seq_list <- list()
# We want to create an empty hub for our proteins to be stored this allows us to do so.
for(i in 1:length(urls)) {
  file_name <- paste0(accessions[i], ".fasta")
  download.file(urls[i], destfile = file_name, mode = "wb")
  seq_list[[accessions[i]]] <- readAAStringSet(file_name)
  cat("Downloaded:", file_name, "\n")
}
# Every fasta file from Unitprot will be read right here in the R environment.
go_info <- GetProteinGOInfo(accessions)
# This retrieves Gene Ontology (GO) annotations for each protein accession from the UniProt database.
head(go_info)
# This displays the first few rows of GO term data to confirm it was retrieved successfully.
PlotGoInfo(go_info)
# This generates a visualization showing the distribution of GO terms for the proteins.
PlotGOAll(
  GOObj = go_info,
  Top = 10,
  directorypath = getwd(),
  width = 8,
  height = 5
)
# This is very important it creates a visual for every GO term of the proteins.
path_info <- GetPathology_Biotech(accessions)
# This retrieves information about any known pathology or disease associations linked to the proteins.
path_info
# This displays the pathology information retrieved from UniProt.
disease_info <- Get.diseases(path_info)
# This extracts disease information from the pathology dataset.
disease_info
# This displays any disease associations identified.
uniprot_info <- fetch_uniprot(accessions)
# This retrieves a large dataset of metadata for each protein accession from UniProt.
View(uniprot_info)
# This opens the dataframe so the information can be visually seen.
uniprot_info$xref_pdb
# This shows the Protein Data Bank (PDB) IDs associated with these proteins if structures exist.
pdb_info <- fetch_pdb(c("1ZMR","2HWG"))
# This downloads structural information for the selected PDB protein structures.
pdb_info
# This displays the PDB structural data.
af_info <- fetch_alphafold_prediction(accessions)
# This retrieves AlphaFold predicted 3D structure information for the proteins.
af_info
# This displays the AlphaFold structure data for the proteins.