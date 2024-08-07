# Bash script to use cutadapt
# ```{bash}
# for i in *_R1.fastq.gz; do
#     SAMPLE=$(echo ${i} | sed "s/_R1\.fastq\.gz//")
#     echo "Processing: ${SAMPLE}_R1.fastq.gz and ${SAMPLE}_R2.fastq.gz"
#     cutadapt -a ^CYGCGGTAATTCCAGCTC...CRAAGAYGATYAGATACCRT \ 
#     -A AYGGTATCTRATCRTCTTYG...GAGCTGGAATTACCGCRG \
#     -o run3/${SAMPLE}_R1_trimmed.fastq.gz -p run3/${SAMPLE}_R2_trimmed.fastq.gz \
#     ${SAMPLE}_R1.fastq.gz ${SAMPLE}_R2.fastq.gz \
#     >> run3/cutadapt_primer_trimming_stats.txt 2>&1
# done
# 
# ```
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

##Set the current working directory
setwd("/home/mlopez/amplicons/run3")

##Assign storage directories
raw_dir <- getwd()
quality_dir <- "./quality/"  
filtered_dir <- "./filtered/"  
results_dir <- "./results/" 

##Create storage directories
dir.create("quality")
dir.create("filtered")
dir.create("results")

#List fastq files
fns <- sort(list.files(raw_dir, pattern = ".fastq", full.names = FALSE))
fns_R1 <- fns[str_detect(basename(fns), "R1")]
fns_R2 <- fns[str_detect(basename(fns), "R2")]

##Extract and save sample names
sample.names <- str_split(basename(fns_R1), pattern = "_", simplify = TRUE)
sample.names <- sample.names[, 1]


##Calculate the number of reads
df <- data.frame()
for (i in 1:length(fns_R1)) {
  geom <- fastq.geometry(fns_R1[i])
  df_one_row <- data.frame(n_seq = geom[1], file_name = basename(sample.names[i]))
  df <- bind_rows(df, df_one_row)
}

#Plot a histogram of the number of reads
pdf("nreads.pdf")
ggplot(df, aes(x = n_seq)) + 
  geom_histogram(alpha = 0.5, position = "identity", binwidth = 100)
dev.off()

##Quality plots of the reads per sample

for (i in 1:length(fns)) {
  qplots <- plotQualityProfile(fns[i])
  qplots_png<- paste0(quality_dir, basename(fns[i]), "-qual.png")
  ggsave(plot = qplots, filename = qplots_png, device = "png", width = 15, 
         height = 15, scale = 1, units = "cm")
}

##Assign names to the FASTQ files that will contain the filtered reads
filt_R1 <- str_c(filtered_dir, sample.names, "_R1_filt.fastq")
filt_R2 <- str_c(filtered_dir, sample.names, "_R2_filt.fastq")

##Filter and trimming
out <- filterAndTrim(fns_R1, filt_R1, fns_R2, filt_R2, truncLen = c(262, 195), 
                     maxN = 0, maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE, 
                     compress = FALSE, multithread = 20)

#Save table
write.table(out, file = './results/reads_filtered.tsv', sep='\t', 
            row.names = TRUE, col.names = TRUE, quote = FALSE)

##Denoising
err_R1 <- learnErrors(filt_R1, multithread = 20)
err_R2 <- learnErrors(filt_R2, multithread = 20)
#Plot estimated error rates
pdf("eF.18S.3.pdf")
plotErrors(err_R1, nominalQ = TRUE)
dev.off()
pdf("eR.18S.3.pdf")
plotErrors(err_R2, nominalQ = TRUE)
dev.off()

## Inference of ASVs
#Deduplication
derep_R1 <- derepFastq(filt_R1, verbose = FALSE)
derep_R2 <- derepFastq(filt_R2, verbose = FALSE)
#Rename the files with the deduplicated reads
names(derep_R1) <- sample.names
names(derep_R2) <- sample.names
#Infer ASVs
dada_R1 <- dada(derep_R1, err = err_R1, multithread = 20, pool = FALSE)
dada_R2 <- dada(derep_R2, err = err_R2, multithread = 20, pool = FALSE)

## Merge paired sequences

mergers <- mergePairs(dada_R1, derep_R1, dada_R2, derep_R2, verbose = TRUE)

## Create a sequence table and print its dimensions
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

## Sequence lenghts
# Calculate the distribution of sequence lengths
table(nchar(getSequences(seqtab)))
#Save a histogram of sequence lengths
pdf("hisrun3.pdf")
hist(nchar(getSequences(seqtab)))
dev.off()

##Save objects to RDS files
saveRDS(sample.names, "sample.names.3.rds")
saveRDS(out,"out.3.rds")
saveRDS(dada_R1, "dada.R1.3.rds")
saveRDS(dada_R2,"dada.R2.3.rds")
saveRDS(mergers,"mergers.3.rds")
saveRDS(seqtab,"seqtab.3.rds")

