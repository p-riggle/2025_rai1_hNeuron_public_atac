# this code reads mouse cortical bulk rnaseq data from ENCODE, wildtype
# and combines and prepares data for downstream PCA plot

# load all libraries
library(tidyr)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggrepel)
library(tibble)
library(readr)

# read in the ENCODE data (rep1)
data0 <- read.delim("data/huang/p10_ctx.tsv", sep = "\t", header = TRUE)
data <- read.delim("data/huang/p14_ctx.tsv", sep = "\t", header = TRUE)
data1 <- read.delim("data/huang/p25_ctx.tsv", sep = "\t", header = TRUE)
data2 <- read.delim("data/huang/p36_ctx.tsv", sep = "\t", header = TRUE)
data3 <- read.delim("data/huang/p60_ctx.tsv", sep = "\t", header = TRUE)

# read in rep 2
rep0 <- read.delim("data/huang/p10_ctx_rep2.tsv", sep = "\t", header = TRUE)
rep <- read.delim("data/huang/p14_ctx_rep2.tsv", sep = "\t", header = TRUE)
rep1 <- read.delim("data/huang/p25_ctx_rep2.tsv", sep = "\t", header = TRUE)
rep2 <- read.delim("data/huang/p36_ctx_rep2.tsv", sep = "\t", header = TRUE)
rep3 <- read.delim("data/huang/p60_ctx_rep2.tsv", sep = "\t", header = TRUE)


# the gene_id is by NCBI identifiers
colnames(data)

# only want to keep the gene_id, expected_count, and FPKM columns
data_subset0 <- data0[, c("gene_id", "expected_count", "FPKM")]
data_subset0 <- rename(data_subset0, FPKM_p10 = FPKM, count_p10 = expected_count)

data_subset <- data[, c("gene_id", "expected_count", "FPKM")]
data_subset <- rename(data_subset, FPKM_p14 = FPKM, count_p14 = expected_count)
head(data_subset)

data_subset1 <- data1[, c("gene_id", "expected_count", "FPKM")]
data_subset1 <- rename(data_subset1, FPKM_p25 = FPKM, count_p25 = expected_count)

data_subset2 <- data2[, c("gene_id", "expected_count", "FPKM")]
data_subset2 <- rename(data_subset2, FPKM_p36 = FPKM, count_p36 = expected_count)

data_subset3 <- data3[, c("gene_id", "expected_count", "FPKM")]
data_subset3 <- rename(data_subset3, FPKM_p60 = FPKM, count_p60 = expected_count)


rep0_subset0 <- rep0[, c("gene_id", "expected_count", "FPKM")]
rep0_subset0 <- rename(rep0_subset0, FPKM_p10_rep2 = FPKM, count_p10_rep2 = expected_count)

rep_subset <- rep[, c("gene_id", "expected_count", "FPKM")]
rep_subset <- rename(rep_subset, FPKM_p14_rep2 = FPKM, count_p14_rep2 = expected_count)

rep1_subset1 <- rep1[, c("gene_id", "expected_count", "FPKM")]
rep1_subset1 <- rename(rep1_subset1, FPKM_p25_rep2 = FPKM, count_p25_rep2 = expected_count)

rep2_subset2 <- rep2[, c("gene_id", "expected_count", "FPKM")]
rep2_subset2 <- rename(rep2_subset2, FPKM_p36_rep2 = FPKM, count_p36_rep2 = expected_count)

rep3_subset3 <- rep3[, c("gene_id", "expected_count", "FPKM")]
rep3_subset3 <- rename(rep3_subset3, FPKM_p60_rep2 = FPKM, count_p60_rep2 = expected_count)


#combine the data_subset into one file by gene_id
combined <- data_subset0 %>%
  full_join(data_subset, by = "gene_id") %>%
  full_join(data_subset1, by = "gene_id") %>%
  full_join(data_subset2, by = "gene_id") %>%
  full_join(data_subset3, by = "gene_id") %>%
  full_join(rep0_subset0, by = "gene_id") %>%
  full_join(rep_subset, by = "gene_id") %>%
  full_join(rep1_subset1, by = "gene_id") %>%
  full_join(rep2_subset2, by = "gene_id") %>%
  full_join(rep3_subset3, by = "gene_id")

head(combined)

# now, prepare the data for PCA later

# filter to include rows with ENSMUSG ids only
library(dplyr)
library(stringr)

combined <- combined %>%
  filter(str_starts(gene_id, "ENSMUSG"))

head(combined)

# remove the version number from ENSEMBL ids to match with the Huang data for downstream subsetting
combined <- combined %>%
  mutate(ensembl_id = sub("\\.[0-9]+$", "", gene_id))
head(combined)
# make gene_id as rownames
rownames(combined) <- combined$ensembl_id
combined$gene_id <- NULL
combined$ensembl_id <- NULL
head(combined)

######################################################################################
# save a version with FPKM at this stage for use later 
#take rowmeans of replicates

combined_FPKM <- combined %>%
  select(starts_with("FPKM"), starts_with("ensembl"))

combined_FPKM$FPKM_p10 <- rowMeans(combined_FPKM[, c(1,6)])
combined_FPKM$FPKM_p14 <- rowMeans(combined_FPKM[, c(2,7)])
combined_FPKM$FPKM_p25 <- rowMeans(combined_FPKM[, c(3,8)])
combined_FPKM$FPKM_p36 <- rowMeans(combined_FPKM[, c(4,9)])
combined_FPKM$FPKM_p60 <- rowMeans(combined_FPKM[, c(5,10)])

combined_FPKM[, c(6:10)] <- NULL

# move ensembl_id to first column
combined_FPKM <- combined_FPKM %>% 
  select(ensembl_id, everything())

# now save
write_tsv(combined_FPKM, "data/huang/combined_FPKM_encode_2reps.tsv")
###################################################################################

# next, delete the FPKM column 
combined$FPKM_p10 <- NULL
combined$FPKM_p14 <- NULL
combined$FPKM_p25 <- NULL
combined$FPKM_p36 <- NULL
combined$FPKM_p60 <- NULL

combined$FPKM_p10_rep2 <- NULL
combined$FPKM_p14_rep2 <- NULL
combined$FPKM_p25_rep2 <- NULL
combined$FPKM_p36_rep2 <- NULL
combined$FPKM_p60_rep2 <- NULL

# round numbers to make integers
combined <- round(combined) # now the data is all ready!

# save "combined" for use in PCA analysis later. 
write.table(combined, "data/huang/combined_ENCODE_data.tsv", row.names = TRUE, quote = FALSE, header = TRUE, sep = "\t")

# end of script