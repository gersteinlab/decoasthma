# quantile-normalize.R
#
# This script reads in the exceRpt outputs exceRpt_gencode_ReadCounts.txt and 
# exceRpt_gencode_ReadsPerMillion.txt and clinical data that define the sample
# set and then outputs an RData file with the expression matrices reduced
# to only the sample set in the clinical data, excluding outliers, 
# with quantile normalization, protein subsetting, and rpkm normalization.
#
# Dan Spakowicz
#

## Load required packages
if (!require("install.load")) {
  install.packages("install.load")
  library(install.load)
}

# CRAN packages
list.of.packages <- c("tidyverse", "preprocessCore")
install.load::install_load(list.of.packages)


# Load exxpression data
date <- "2017-01-21"
counts <- read.table(paste("../exceRpt", date,
                            "exceRpt_gencode_ReadCounts.txt", sep = "/"),
                      stringsAsFactors = FALSE)
rpm <- read.table(paste("../exceRpt", date,
                           "exceRpt_gencode_ReadsPerMillion.txt", sep = "/"),
                     stringsAsFactors = FALSE)

expr.mats <- list(counts = counts,
                  rpm = rpm)

# Load clinical data
# Saved from "exceRpt_summary.Rmd" where I subsetted the clinical data by 
# removing COPD and ensuring one visit per person
load("../exceRpt/clinical_subset.Rdata")
  
expr.mats.sub <- lapply(expr.mats, function(x) {
  x[, colnames(x) %in% make.names(df$filename)]
})

# Function for quantile normalization
QuantileNormalize <- function(edata) {
  # Convert to matrix
  mat <- edata %>% as.matrix
  # Quantile normalize
  nedata <- preprocessCore::normalize.quantiles(mat)
  # Replace same and gene names
  row.names(nedata) <- row.names(edata)
  colnames(nedata) <- colnames(edata)
  
  return(nedata)
}

# Quantile normalize the two matrices
qnorm.mats <- lapply(expr.mats.sub, QuantileNormalize)

# Find outliers by PCA of quantile normalized count data
pca.check <- prcomp(t(qnorm.mats[["counts"]]))
plot(pca.check$x[,1], pca.check$x[,2])
outliers <- rownames(pca.check$x)[pca.check$x[,1] > 1e6]

# Function to remove any outliers
RemoveOutliers <- function(exprs){
  exprs <- exprs[, -which(colnames(exprs) %in% outliers)]
  return(exprs)
}

# Remove the outliers
expr.mats.sub.o <- lapply(expr.mats.sub, RemoveOutliers)
qnorm.mats.o <- lapply(qnorm.mats, RemoveOutliers)

# Exclude the outlier from the clinical data
clinical <- df[!make.names(df$filename) %in% outliers,]

# Create an object of all the expression matrices
all.mats <- list(counts = expr.mats.sub.o[["counts"]],
                rpm = expr.mats.sub.o[["rpm"]],
                qnorm.counts = qnorm.mats.o[["counts"]],
                qnorm.rpm = qnorm.mats.o[["rpm"]])

# Reduce to only protein coding genes
onlyProteinCoding <- function(x) {
  x[grep(":protein_coding$", row.names(x)),] 
}

all.mats.protein <- lapply(all.mats, onlyProteinCoding)

# Function to convert protein rpm to rpkm
# Read in key of protein lengths
x <- read.table("../exceRpt/gencode_v24_protein_lengths.txt", header = TRUE)

# Find the row matches in the key
matches <- match(rownames(all.mats.protein[["rpm"]]), x$name)

# rpkm calculation
rpkm.mats <- list(protein.rpkm = sweep(all.mats.protein[["rpm"]], 1, 
                                       x$length[matches], "/"),
                  qnorm.protein.rpkm = sweep(all.mats.protein[["qnorm.rpm"]], 1,
                                             x$length[matches], "/"))


##
## Save 
##

# Save the matrices and the clinical data in a single RData file
save(all.mats, all.mats.protein, rpkm.mats, clinical, 
     file = "../exceRpt/counts.rpm.protein.rpkm.clinical.Rdata")

# md5checksums to ensure accurate upload
system({
  "md5 ../exceRpt/counts.rpm.protein.rpkm.clinical.Rdata > ../exceRpt/counts.rpm.protein.rpkm.clinical.Rdata.md5"
})

# Save each as separate csv files and compress them
mats <- ls(pattern = "all\\.mats.*|rpkm.mats")
for (i in mats) {
  for (j in names(get(i))) {
    write.csv(get(i)[[j]], file = paste("../exceRpt/", i, "_", j, ".csv", 
                                        sep = ""))
  }
}
write.csv(clinical, "../exceRpt/clinical.csv")
# Compress csv's into a single tarball
system({
  "tar -cvzf ../exceRpt/expression_clinical_csvs.tar.gz ../exceRpt/*.csv"
})
# Delete the uncompressed csvs
system({
  "rm ../exceRpt/*.csv"
})

# md5checksums to ensure accurate upload
system({
  "md5 ../exceRpt/expression_clinical_csvs.tar.gz > ../exceRpt/expression_clinical_csvs.tar.gz.md5"
})
