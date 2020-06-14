#Install BioManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.11")

# Load the BiocInstaller package
library(BiocInstaller)

# Using Bioconductor version 3.6
BiocManager::install("BSgenome")

# Load the BSgenome package
library(BSgenome)

# Check the version of the BSgenome package
packageVersion("BSgenome")

# Check the BSgenome class results and find its parent classes (Extends) and the classes that inherit from it (Subclasses).
showClass("BSgenome")

# Learn more about BSgenome
??BSgenome


# What genomes are currently installed:
installed.genomes()

# What genomes are available:
available.genomes()
summary(available.genomes())

# Split the package names in parts:
av_gen <- available.genomes(splitNameParts=TRUE)
table(av_gen$organism) # several genomes available including Ecoli, Scerevisiae, Mmusculus and Hsapiens
table(av_gen$provider) # most of the genomes available come from UCSC provider

# Assign data to the Ecoli_genome object
Ecoli_genome <- BSgenome.Ecoli.NCBI.20080805::BSgenome.Ecoli.NCBI.20080805



# Investigate the Ecoli_genome using show()
# Ecoli_genome contains the genome of E.coli with 13 sequences
Ecoli_genome
show(Ecoli_genome)

# Investigate some other accesors
organism(Ecoli_genome)  #E.coli organism
provider(Ecoli_genome)  # NCBI provider
seqinfo(Ecoli_genome)  # 13 sequences (13 circular) ==> the table contains the following columns: seqnames(ID sequence), seqlengths(base pair length), iscircular(true) and genome(date)



#dowload a genome who is available but not already installed: download the genome of drosophila (Dmelanogaster)
# Load and install BiocInstaller
install.packages("BiocInstaller",repos="http://bioconductor.org/packages/3.4/bioc")
library(BiocInstaller)

#choose a genome available : e.g Dmelanogaster
biocLite("BSgenome.Dmelanogaster.UCSC.dm2")

# Load the package and display the index of sequences for this genome:
library(BSgenome.Dmelanogaster.UCSC.dm2)
Dmelanogaster # same as BSgenome.Scerevisiae.UCSC.sacCer1, Scerevisiae has 13 sequences





#Analysis using  the object previously created Ecoli_genome - 13 circular sequences
#Summary
summary(Ecoli_genome)
#name of the 13 sequences
seqnames(Ecoli_genome)
#lengths of the 13 sequences (base pair length)
seqlengths(Ecoli_genome)

# Print chromosome NC_008253 ==> DNA string
Ecoli_genome$NC_008253

# Count characters of the chromosome NC_008253 
nchar(Ecoli_genome$NC_008253)

# Get the first 30 bases of each chromosome
getSeq(Ecoli_genome, start = 1, end = 30)





#------------------------------------
#STUDY CASE :BAT CORONAVIRUS FROM ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/bat_coronavirus_1a_uid29247/NC_010437.ffn
#------------------------------------------
#Load packages
#Biostrings package ==> efficient manipulation of biological strings
library(Biostrings)


dir.create("data", showWarnings = FALSE)
system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/bat_coronavirus_1a_uid29247/NC_010437.ffn")
download.file("ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/bat_coronavirus_1a_uid29247/NC_010437.ffn", "data/NC_010437.ffn")

#Import FASTA file with readDNAStringSet
myseq <- readDNAStringSet("data/NC_010437.ffn")
#viruslize all the 7 sequences coronavirus from the link from NCBI
myseq
#number of coronavirus sequences from file
length(myseq)
#visualize only the 3 first coronavirus sequences
myseq[1:3]
corona1seq <- myseq[1:1]
corona1seq 
# "tricky part" ==> Subset sequences with regular expression on sequence name field
#sub <- myseq[grep( "",names(myseq))]
#length(sub) # 7 coronavirus sequences from file

#Export subsetted sequences to FASTA file
#writeXStringSet(sub, file="./data/NC_010437sub.ffn")




# Unlist the set, select the first 30 letters, and assign to dna_seq
dna_seq <- subseq(unlist(corona1seq), end = 30)
dna_seq

# Transcribe dna_seq into an RNAString object and print it
rna_seq <- RNAString(dna_seq) 
rna_seq

# Translate rna_seq into an AAString object and print it
aa_seq <- translate(rna_seq)
aa_seq























#------------------------------------------------------------------------------
#working with coronavirus genome
#coronavirus FATA seq from NCBI https://www.ncbi.nlm.nih.gov/nuccore/1798174254
#------------------------------------------------------------------------------

#Load packages
#Biostrings package ==> efficient manipulation of biological strings
library(Biostrings)

readDNAStringSet("coronavirus.fa")



#BiocManager::install("seqinr")
#install.packages("seqinr")

#BiocManager::install("GenomicTools.fileHandler")
#install.packages("GenomicTools.fileHandler")

# Define here the location on HDD for the example file
#fpath <- system.file("coronavirus_fa", package="GenomicTools.fileHandler")
# Import the example fasta file
#fastaFile <- importFA(file=fpath)

















