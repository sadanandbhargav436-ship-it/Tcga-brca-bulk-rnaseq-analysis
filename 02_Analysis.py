# Differential expression, visualization, and survival analysis for TCGA-BRCA.
library(TCGAbiolinks)
library(dplyr)
library(ggplot2)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(survival)
library(survminer)

if (!dir.exists("results")) {
  dir.create("results")
}

message("Loading data...")
counts <- readRDS("data/BRCA_counts.rds")
clinical <- readRDS("data/BRCA_clinical.rds")

message("Preparing sample groups...")
clinical <- clinical %>%
  mutate(
    sample_type = factor(sample_type, levels = c("Solid Tissue Normal", "Primary Tumor")),
    group = ifelse(sample_type == "Primary Tumor", "Tumor", "Normal")
  )

# Keep only samples present in both clinical and count matrices
common_samples <- intersect(colnames(counts), rownames(clinical))
counts <- counts[, common_samples]
clinical <- clinical[common_samples, , drop = FALSE]

message("Filtering low-count genes...")
gene_filter <- rowSums(counts >= 10) >= 5
counts <- counts[gene_filter, ]

message("Running differential expression analysis...")
de_results <- TCGAanalyze_DEA(
  mat1 = counts[, clinical$group == "Tumor"],
  mat2 = counts[, clinical$group == "Normal"],
  Cond1type = "Tumor",
  Cond2type = "Normal",
  method = "glmLRT"
)

message("Saving DEA results...")
write.csv(de_results, file = "results/BRCA_DEA_results.csv", row.names = FALSE)

if (!"genes" %in% colnames(de_results)) {
  if ("gene" %in% colnames(de_results)) {
    de_results$genes <- de_results$gene
  } else {
    de_results$genes <- rownames(de_results)
  }
}

if (nrow(de_results) > 0) {
  top20 <- de_results %>%
    arrange(FDR) %>%
    slice_head(n = 20)
  write.csv(top20, file = "results/BRCA_DEA_top20.csv", row.names = FALSE)
}

message("Creating volcano plot...")
if ("logFC" %in% colnames(de_results) && "FDR" %in% colnames(de_results)) {
  png("results/BRCA_volcano.png", width = 1400, height = 1000, res = 150)
  
  # Note: Explicit print() required for ggplot/EnhancedVolcano inside graphic device
  p_volcano <- EnhancedVolcano(
    de_results,
    lab = de_results$genes,
    x = "logFC",
    y = "FDR",
    title = "TCGA-BRCA Tumor vs Normal",
    pCutoff = 0.05,
    FCcutoff = 1.0,
    pointSize = 2.0,
    labSize = 3.0
  )
  print(p_volcano)
  dev.off()
}

message("Selecting top genes for heatmap...")
if (nrow(de_results) >= 50) {
  key_genes <- de_results %>%
    arrange(FDR) %>%
    slice_head(n = 50) %>%
    pull(genes)
  
  key_genes <- unique(key_genes)
  
  # Log2 transform (CPM / normalized) and scale across rows (Z-score)
  log_counts <- log2(counts[key_genes, , drop = FALSE] + 1)
  scaled_matrix <- t(scale(t(log_counts)))
  
  clinical_annotation <- HeatmapAnnotation(
    Group = clinical$group,
    col = list(Group = c("Tumor" = "#E41A1C", "Normal" = "#377EB8"))
  )
  
  png("results/BRCA_heatmap_top50.png", width = 1600, height = 1200, res = 150)
  
  Heatmap(
    scaled_matrix,
    name = "Z-Score",
    top_annotation = clinical_annotation,
    show_row_names = TRUE,
    show_column_names = FALSE,
    cluster_rows = TRUE,
    cluster_columns = TRUE
  )
  
  dev.off()
}

message("Running pathway enrichment on significant DE genes...")
sig_genes <- de_results %>%
  filter(FDR <= 0.05, abs(logFC) >= 1) %>%
  pull(genes)

if (length(sig_genes) >= 10) {
  tryCatch({
    enrich <- TCGAanalyze_EAcomplete(
      geneList = sig_genes,
      organism = "hsa",
      databases = c("KEGG", "GO_Biological_Process", "GO_Molecular_Function")
    )
    write.csv(enrich, file = "results/BRCA_enrichment_results.csv", row.names = FALSE)
  }, error = function(e) {
    message("Enrichment analysis failed: ", conditionMessage(e))
  })
}

message("Preparing survival analysis...")
req_surv_cols <- c("vital_status", "days_to_death", "days_to_last_follow_up")

if (all(req_surv_cols %in% colnames(clinical))) {
  surv_data <- clinical %>%
    mutate(
      survival_time = dplyr::coalesce(
        suppressWarnings(as.numeric(days_to_death)),
        suppressWarnings(as.numeric(days_to_last_follow_up))
      ),
      status = ifelse(toupper(vital_status) == "DEAD", 1, 0)
    ) %>%
    filter(!is.na(survival_time) & survival_time > 0 & group == "Tumor")
  
  # Stratify by top DE gene log2 expression rather than unweighted mean
  top_gene <- de_results %>% arrange(FDR) %>% pull(genes) %>% .[1]
  
  if (nrow(surv_data) >= 20 && !is.na(top_gene)) {
    gene_expr <- log2(counts[top_gene, rownames(surv_data)] + 1)
    surv_data$expression_score <- gene_expr
    median_score <- median(surv_data$expression_score, na.rm = TRUE)
    surv_data$risk_group <- ifelse(surv_data$expression_score >= median_score, "High", "Low")
    
    fit <- survival::survfit(
      survival::Surv(survival_time, status) ~ risk_group,
      data = surv_data
    )
    
    png("results/BRCA_survival_plot.png", width = 1400, height = 1000, res = 150)
    
    p_surv <- ggsurvplot(
      fit,
      data = surv_data,
      risk.table = TRUE,
      pval = TRUE,
      title = paste0("Survival Stratified by ", top_gene, " Expression")
    )
    print(p_surv)
    dev.off()
    
    summary_text <- capture.output(summary(fit))
    writeLines(summary_text, "results/BRCA_survival_summary.txt")
  }
} else {
  message("Required clinical survival columns not present.")
}

message("Analysis complete. Review results in the results/ directory.")
