# TCGA-BRCA data download and preprocessing
# This script downloads BRCA transcriptome data from GDC and saves prepared datasets.

if (!dir.exists("data")) {
  dir.create("data")
}
if (!dir.exists("results")) {
  dir.create("results")
}

options(timeout = 600)

library(TCGAbiolinks)
library(dplyr)

message("Querying TCGA-BRCA transcriptome profiling data...")
query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

message("Downloading TCGA data. This can take a while depending on the connection.")
GDCdownload(query, directory = "data", method = "api", files.per.chunk = 50)

message("Preparing the downloaded data into a SummarizedExperiment object...")
brca_se <- GDCprepare(query, directory = "data")

message("Saving raw SummarizedExperiment object to disk...")
saveRDS(brca_se, file = "data/BRCA_SE_raw.rds")

message("Extracting clinical and count data...")
clinical <- colData(brca_se) %>% as.data.frame()
counts <- assay(brca_se)

saveRDS(clinical, file = "data/BRCA_clinical_raw.rds")
saveRDS(counts, file = "data/BRCA_counts_raw.rds")

message("Filtering primary tumor and normal samples...")
clinical <- clinical %>%
  mutate(sample_type = as.character(sample_type)) %>%
  filter(sample_type %in% c("Primary Tumor", "Solid Tissue Normal"))

sample_ids <- rownames(clinical)
counts <- counts[, sample_ids]

saveRDS(clinical, file = "data/BRCA_clinical.rds")
saveRDS(counts, file = "data/BRCA_counts.rds")

message("Downloaded and preprocessed TCGA-BRCA count data and clinical metadata.")

