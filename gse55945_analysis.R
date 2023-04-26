# my_paths <- c("C:/Program Files/R/R-4.2.3/library", .libPaths())
# .libPaths(my_paths)
# 
# new_paths <- .libPaths()[.libPaths() != "C:/Program Files/R/R-4.2.3/library"]
# .libPaths(new_paths)

# install.packages("BiocManager")
# BiocManager::install("affy", force=TRUE)
# BiocManager::install("affyio", force=TRUE)
# BiocManager::install("GEOquery", force=TRUE)

# Install packages
if (!require("pacman"))
  install.packages("pacman")

pacman::p_load(dplyr, affy, affyio)

# Define file structure
data_path <- "C:/Users/user/Documents/Addon Projects/GSE55945/data/"
src_file <- "GSE55945_RAW.tar"
src_file_path <- paste(data_path,src_file, sep="")

# Extract the contents of the tar file
# Only for the first time
# untar(src_file_path, exdir = data_path)

# 2 Corrupted files removed
# 1. GSM1348937_110807_HGU133_PLUS_2.0_MS_36A6
# 2. GSM1348948_011508_HGU133_PLUS_2.0_MS_36D2

# Load the CEL files
list.celfiles(data_path)
affyData <- ReadAffy(celfile.path = data_path)
# normalize data
eset.mas5 = mas5(affyData)
# Expression Values
exprSet_nologs = exprs(eset.mas5)
# Get the column names
new_colnames <- gsub("^(.*?)_.*$", "\\1", colnames(exprSet_nologs))
# Assign the new column names to the matrix
colnames(exprSet_nologs) <- new_colnames
# List the column (chip) names
colnames(exprSet_nologs)

# Create heatmap
# heatmap(exprSet_nologs,main="Normalized ME matrix for brain, liver, N=19")

# log transform for further normalization
exprSet <- log(exprSet_nologs,2)

# Run the Affy A/P call algorithm on the CEL files we processed above
mas5callsAP <- mas5calls(affyData)

# Get the actual A/P calls
# Generates matrix for P = Present and A = Absent
mas5callsAP_calls <- exprs(mas5callsAP)

# Create 3 lists for 3 comparison groups
grA_ergneg <- c('GSM1348933','GSM1348934','GSM1348935','GSM1348936','GSM1348938','GSM1348939')
grB_ergpos <- c('GSM1348940','GSM1348941','GSM1348942','GSM1348943','GSM1348944','GSM1348945')
grC_normal <- c('GSM1348946','GSM1348947','GSM1348949','GSM1348950','GSM1348951','GSM1348952','GSM1348953')
grD_AandB <- c(grA_ergneg, grB_ergpos)

# Calculate Expression Ratios of genes between different groups
# Step1:Calculate mean
grA_mean <- apply(exprSet[,grA_ergneg],1,mean)
grB_mean <- apply(exprSet[,grB_ergpos],1,mean)
grC_mean <- apply(exprSet[,grC_normal],1,mean)
grD_mean <- apply(exprSet[,grD_AandB],1,mean)

# Compare between log transformed sets
grC_vs_grA <- grC_mean - grA_mean
grC_vs_grB <- grC_mean - grB_mean
grC_vs_grD <- grC_mean - grD_mean
