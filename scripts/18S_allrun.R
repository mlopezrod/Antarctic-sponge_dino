
##Load libraries 

library("dada2")
library("phyloseq")
library("Biostrings")
library("ggplot2")
library("dplyr")
library("tidyr")
library("readxl")
library("readr")
library("stringr")
library("tibble")
library("DECIPHER")
library("ape")
library("phangorn")

#Set the current working directory
setwd("/home/mlopez/amplicons/allrun")


#Assign storage directories
raw_dir <- getwd()
results_dir <- "./results/" 

#Create storage directories
dir.create("results")



# Read the RDS files and convert them to data frames
out.1<-as.data.frame(readRDS("./out.1.rds"))
out.2<-as.data.frame(readRDS("./out.2.rds"))
out.3<-as.data.frame(readRDS("./out.3.rds"))
# Combine the data frames into a single data frame
out<-bind_rows(out.1,out.2,out.3)

# Read the RDS files
sample.names.1<-readRDS("./sample.names.1.rds")
sample.names.2<-readRDS("./sample.names.2.rds")
sample.names.3<-readRDS("./samples.names.3.rds")

# Combine the vectors into a single vector
sample.names<- c(sample.names.1,sample.names.2,sample.names.3)

# Read the RDS files
dada.R1.1<-readRDS("./dada.R1.1.rds")
dada.R2.1<-readRDS("./dada.R2.1.rds")
dada.R1.2<-readRDS("./dada.R1.2.rds")
dada.R2.2<-readRDS("./dada.R2.2.rds")
dada.R1.3<-readRDS("./dada.R1.3.rds")
dada.R2.3<-readRDS("./dada.R2.3.rds")
# Combine the list into a single list
dada.R1<- c(dada.R1.1,dada.R1.2,dada.R1.3)
dada.R2<- c(dada.R2.1,dada.R2.2,dada.R2.3)

# Read the RDS files
merge.1<-readRDS("./mergers.1.rds")
merge.2<-readRDS("./mergers.2.rds")
merge.3<-readRDS("./mergers.3.rds")
#Combine the merged data from into a single list
merge<-c(merge.1,merge.2,merge.3)

# Read the RDS files
st1 <- readRDS("./seqtab.1.rds")
st2 <- readRDS("./seqtab.2.rds")
st3 <- readRDS("./seqtab.3.rds")
# Merge the sequence tables into a single sequence table and print its dimensions
st.all <- mergeSequenceTables(st1,st2,st3)
dim(st.all)

# Sequence lenghts
# Calculate the distribution of sequence lengths
table(nchar(getSequences(st.all)))
#Save a histogram of sequence lengths
png("histograma.png", units = "cm", width = 15, height = 10, res = 300)
hist(nchar(getSequences(st.all)))
dev.off()

# Collapse identical sequences and print dimensions of the table sequence
st.all<-collapseNoMismatch(st.all,verbose=T)
dim(st.all)

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(st.all, method = "consensus", 
                                    multithread = 20, verbose = TRUE)

t_seqtab.nochim <- as.data.frame(t(seqtab.nochim)) %>% rownames_to_column(var = "sequence") %>% mutate(sequence = str_replace_all(sequence, "(-|\\.)", ""))


#Extracts sequences and saves them in a FASTA file.
seqs <- getSequences(seqtab.nochim)
unfiltered.seqs <- DNAStringSet(seqs)
names(unfiltered.seqs) <- paste0("ASV", sprintf("%03d", seq_along(seqs)))
writeXStringSet(unfiltered.seqs, "./UnfilteredASVsequence.fa", append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")

#Track reads
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dada.R1, getN), sapply(merge, getN), rowSums(st.all), 
               rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
write.table(track,file = "./results/summarize-stats.tsv", sep="\t",  na = "NA", 
                        row.names = T,col.names = T, quote = F)

#Assign taxonomy
set.seed(100)
# Define a vector of taxonomic levels used in the PR2 database
PR2_tax_levels <- c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", 
                    "Genus", "Species")
ref_fasta <- "/datos2/ecogenomicalab/BasesDeDatos/databases_PR2/pr2_version_4.14.0_SSU_dada2.fasta"
taxa <- assignTaxonomy(seqtab.nochim, refFasta = ref_fasta, taxLevels = PR2_tax_levels, minBoot = 80, 
                             outputBootstraps = F, verbose = T, multithread = 20)
taxa_ta<- as.data.frame(taxa) %>% rownames_to_column(var = "sequence") %>% mutate(sequence = str_replace_all(sequence, "(-|\\.)", ""))


#Remove terrestrial plants, macroalgae,metazoa, mitochondria, chloroplasts, bacteria and ASVs not assigned the Kingdom, Supergroup, Division, and Class levels.
is.clmina <- taxa_ta[,"Class"] %in% c("Embryophyceae","Florideophyceae", "Phaeophyceae", "NA") | 
  taxa_ta[,"Division"] %in% c("Metazoa","NA") | 
  taxa_ta[,"Kingdom"] %in% c("Bacteria", "Eukaryota:mito", "Eukaryota:plas", "NA") |
  taxa_ta[,"Supergroup"] %in% "NA" 

seqtab.nochim.ok <- seqtab.nochim[,!is.clmina]
t_seqtab.nochim.ok <- t(seqtab.nochim.ok)
taxa.ok<- taxa_ta[!is.clmina,]



# Save files

write.table(t_seqtab.nochim, file = "./results/UnfilteredASVcount.tsv", sep="\t", 
            na = "NA", row.names =F, col.names = T, quote = F)

write.table(taxa_ta,file = "./results/UnfilteredASVtaxonomy.tsv", sep="\t", na = "NA", 
            row.names =F, col.names = T, quote = F)

write.table(t_seqtab.nochim.ok, "./results/seqtable_nochim.ok.tsv", sep="\t", na = "NA", 
            row.names =F, col.names = T,quote = F)
write.table(taxa.ok, "./results/taxonomy.ok.tsv", sep="\t", na = "NA", 
            row.names =F, col.names = T, quote = F)

# Phylogenetic tree
seqs <- getSequences(seqtab.nochim.ok)
names(seqs) <- seqs 
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)

phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) 
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

detach("package:phangorn", unload=TRUE)
tree_18S<-fitGTR$tree
saveRDS(tree_18S, "tree_18S.rds")