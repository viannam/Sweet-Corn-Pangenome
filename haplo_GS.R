rm(list=ls())

setwd("/red/mresende/viannam/PHG")

# Load required packages
library(BGLR)
library(AGHmatrix)
library(dplyr)
library(ggplot2)
library(data.table)
library(superheat)

## Phenotyping information
## Pheno file should have geno ID column (same order as in the Gmatrix) + trait columns

## Read phenotypic data
BLUEs_trait <- read.csv("/blue/mresende/share/viannam/GS/SC_Blues_Multy-Year_GXE.csv", sep = ",")

# Display the first few rows of the phenotypic data
head(BLUEs_trait)

# Enssure that the genotypes are in the same order on both files
##1341488 haplotypes
#reads_names <- fread("./haplotype_PAV_matrix.txt", sep = "\t", header = TRUE)
reads_names <- fread("./haplotype_PAV_matrix_filtered.txt", sep = "\t", header = TRUE)
reads_names[1:10,1:10]
##294153 haplotypes
#df = fread("./haplotype_PAV_matrix.txt", sep = "\t", header = TRUE, select = c("Sample", "000001eea716c7f13e3fbd69f162afa1", "00000b40f29fbde0d7ea44e660d27615"))
df = fread("./haplotype_PAV_matrix_filtered.txt", sep = "\t", header = TRUE, select = c("Sample", "00000ad150c404d8d7e50a8e0799b6f9", "0000398a9414d91d8d72f998025fa42a"))
t_names <- as.data.frame(df$Sample)
t_names <- data.frame(Sample = df$Sample)

t_names <- as.data.frame(reads_names$Sample)
t_names <- data.frame(Sample = reads_names$Sample)
t_names$Sample <- gsub("_G1", "", t_names$Sample)

# Display the first few rows of haplo_names
head(t_names)

# Reorder rows of BLUEs_trait based on genotype IDs
BLUEs <- BLUEs_trait[match(t_names$Sample, BLUEs_trait$VCF_ID), ]

# Display the first few rows of BLUEs
head(BLUEs)

#haplotypes <- as.matrix(reads_names[,-1])
#haplotypes_transpose <- t(haplotypes)

#haplotypes <- scale(haplotypes, scale = F, center = TRUE)
#write.table(haplotypes_transpose, "haplotypes_matrix_filtered_input.txt", row.names = TRUE, col.names = FALSE, quote = FALSE, sep = "\t")
#haplo <- read.table("haplotypes_gemma_input.txt")
# haplo <- read.table("haplotype_PAV_matrix_filtered.txt")
# haplo [1:10,1:10]

# # Check if there are any non-numeric or NA values
#any(!is.numeric(haplo[,-1]))
#any(!is.numeric(haplotypes[,-1]))
# This should return FALSE
# Remove rows with non-numeric entries
# Ensure SNP IDs are in the first column, genotypes in subsequent columns

# haplo <- haplo[apply(haplo[,-1], 1, function(x) all(is.numeric(x))), ]
# haplo[,-1] <- apply(haplo[,-1], 2, function(x) as.numeric(as.character(x)))
# haplo[is.na(haplo)] <- 0
# haplo <- cbind(SNP_ID = haplo[, 1], haplo[, -1])

# Save the corrected file
#write.table(haplo, "haplotypes_gemma_input_corrected.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
#haplo <- read.table("haplotypes_gemma_input_corrected.txt")

# #calculate the kinship matrix
#G_mat <-  ((haplotypes %*% t(haplotypes))/(ncol(haplotypes)))
#save(G_mat, file = "./G_mat_haplo.RData")
# load( "./G_mat_haplo.RData")
write.table(G_mat, file = "haplo_kinship_matrix_filtered.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
# superheat::superheat(G_mat)
# dim(G_mat)
# range(G_mat)
#G_mat <- fread("./haplo_kinship_matrix.txt")
G_mat <- fread("./haplo_kinship_matrix_filtered.txt")

# Extract phenotypic data for traits
y <- (BLUEs[, 3:22])

# Define a function to run GBLUP with cross-validation
run_GBLUP_CV <- function(data, folds, kinship, nIter, burnIn, nReps) {
  predM1_all <- data.frame()
  traitNames = colnames(data)

  for (trait in 1:ncol(data)) {
    predM1 <- data.frame()
    Y <- data[, trait]

    for (Rep in 1:nReps) {
      nFolds= sample(1:folds, size = length(Y), replace = TRUE)

      for (i in 1:max(nFolds)) {
        tst <- which( nFolds == i)
        yNA <- Y
        yNA[tst] <- NA

        # Fit GBLUP model
        fm <- BGLR(y = as.matrix(yNA),
                   ETA = list(list(K = as.matrix(kinship), model = 'RKHS')),
                   nIter = nIter,
                   burnIn = burnIn,
                   verbose = FALSE)

        # Predicted values
        PC <- cor(Y[tst], fm$yHat[tst], use = "pairwise.complete.obs")
        predM1 <- rbind(predM1, data.frame(Trait = traitNames[trait],
                                           k_Fold = i,
                                           Rep = Rep,
                                           Acc = PC))
      }
    }
    predM1_all <- rbind(predM1, predM1_all)  # Store results for each trait in a list element
  }
  return(predM1_all)
}


# Run GBLUP with cross-validation
set.seed(1234)
predM1_all <- run_GBLUP_CV(data = y,  # Extract phenotypic data
                           folds = 5,
                           kinship = G_mat,
                           nIter = 10000,
                           burnIn = 1000,
                           nReps = 10)

# Save model output for each trait
save(predM1_all, file = "./Accuracy_GBLUP_haplo_Multi-year.RData")
