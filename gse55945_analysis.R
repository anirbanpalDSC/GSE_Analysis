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
exprSet.nologs = exprs(eset.mas5)
# List the column (chip) names
colnames(exprSet.nologs)

# Create heatmap
heatmap(exprSet.nologs,main="Normalized ME matrix for brain, liver, N=19")
