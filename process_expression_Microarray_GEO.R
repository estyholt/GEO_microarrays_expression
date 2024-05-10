#!/usr/bin/env Rscript

library(GEOquery)
library(affy)
library(dplyr)
library(tidyverse)
library(biomaRt)
library (data.table)

## check which platform was used for Microarray (Affy....)

## change my_id to be the desired dataset:
my_id <- "GSE16879" #[[Here]]
gse <- getGEO(my_id)

###############################################################
# Download  expression data

## get supplementary files/;
getGEOSuppFiles(my_id)
untar(paste(my_id, "/", my_id, "_RAW.tar",sep =""), exdir = 'data/')

## reading in .cel files
raw_data <- ReadAffy(celfile.path = "data/")

## performing RMA normalization
normalized_data <- affy::rma(raw_data)

## get expression estimates
normalized_expr <- as.data.frame(exprs(normalized_data))
##change names of sample columns which end with .CEL.gz:
names(normalized_expr) = gsub(pattern = "*\\.CEL\\.gz", replacement = "", x = names(normalized_expr))

## map probe IDs to gene symbols
gse <- getGEO(my_id, GSEMatrix = TRUE)

###############################################################
# Process patient data

## fetch sample data to get disease type mapping
sample_info <- pData(gse[[1]])
sample_info_modified = sample_info %>%
  dplyr::select(c("disease:ch1"))  %>%
  rownames_to_column (var="sample") %>%
  dplyr::rename(disease = "disease:ch1")  %>%
  dplyr::mutate(sample_disease = paste(disease, "_", sample))
View(sample_info_modified)

## going to change column names of samples in exp data to disease_sample by looking at 2 dfs:
## ex:
## GSM364627 -> Control_ GSM364627
## GSM364628 -> Control_GSM364628
## GSM364641 -> UC_GSM364641
for (i in 1:nrow(sample_info_modified)){
  sample <- sample_info_modified[i,1] 
  disease_sample <- sample_info_modified[i,3]
  disease_sample <- gsub(" ", "", disease_sample)
  normalized_expr <- normalized_expr %>% data.table::setnames(old = sample, new = disease_sample)
}

## fetch feature data to get ID - gene symbol mapping
feature_data <- fData(gse[[1]]) ## Usually there will only be 1 platform and the dataset we want to analyse will be the 1st obj in the list, hence [[]]


###############################################################
## Retrieving gene information from _______ Biomart ______

ensembl.con <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
attr <- listAttributes(ensembl.con)
biomart_info <- getBM(attributes = c('affy_hg_u133_plus_2', 'external_gene_name', 'description' ), 
                  filters = 'affy_hg_u133_plus_2', 
                  values = feature_data$ID,
                  mart = ensembl.con)
unique_biomart_info <- biomart_info[!duplicated(biomart_info$affy_hg_u133_plus_2),] # Biomart sometimes gives more than 1 line for a probe

#Joining expression matrix with gene info from biomart:
joined <- left_join(rownames_to_column(normalized_expr) , unique_biomart_info, by=c("rowname" = "affy_hg_u133_plus_2"))


###############################################################
# Write data to files

final_normalised_expr <- joined %>% dplyr::rename("affy_hg_u133_plus_2" = rowname)
write.table(final_normalised_expr, "exp_matrix.csv", quote=FALSE, sep = ',', col.names=NA) 
write.table(sample_info, "patient_info.csv", quote=FALSE, sep = '\t', col.names=NA)
