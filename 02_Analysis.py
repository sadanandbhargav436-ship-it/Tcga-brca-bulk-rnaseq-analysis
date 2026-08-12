# Differential expression, visualization, and survival analysis for TCGA-BRCA.
# Robust script with placeholders and fallbacks so it runs for testing.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(survival)
  library(survminer)
})

# Try to load optional packages (don't stop if missing)
has_tcgabiolinks <- requireNamespace("TCGAbiolinks", quietly = TRUE)
if (has_tcgabiolinks) library(TCGAbiolinks)
has_limma <- requireNamespace("limma", quietly = TRUE)
if (has_limma) library(limma)
has_edgeR <- requireNamespace("edgeR", quietly = TRUE)
if (has_edgeR) library(edgeR)
has_EnhancedVolcano <- requireNamespace("EnhancedVolcano", quietly = TRUE)
if (has_EnhancedVolcano) library(EnhancedVolcano)
has_ComplexHeatmap <- requireNamespace("ComplexHeatmap", quietly = TRUE)
if (has_ComplexHeatmap) library(ComplexHeatmap)
has_pheatmap <- requireNamespace("pheatmap", quietly = TRUE)
if (has_pheatmap) library(pheatmap)

# Ensure results/data directories exist
if (!dir.exists("results")) dir.create("results", recursive = TRUE)
if (!dir.exists("data")) dir.create("data", recursive = TRUE)

# Create small placeholder data if missing (for testing)
counts_path <- "data/BRCA_counts.rds"
clinical_path <- "data/BRCA_clinical.rds"

if (!file.exists(counts_path) || !file.exists(clinical_path)) {
  message("Placeholder data files not found. Creating small synthetic data for testing at data/*.rds")
  set.seed(42)
  n_genes <- 200
  n_samples <- 12
  gene_ids <- paste0("Gene", seq_len(n_genes))
  sample_ids <- paste0("Sample", seq_len(n_samples))
  # count matrix: integers, moderate range
  counts <- matrix(rpois(n_genes * n_samples, lambda = 50), nrow = n_genes, ncol = n_samples,
                   dimnames = list(gene_ids, sample_ids))
  # clinical: rows must match sample names
  sample_type <- rep(c("Primary Tumor", "Solid Tissue Normal"), length.out = n_samples)
  vital_status <- sample(c("Alive", "Dead"), n_samples, replace = TRUE)
  days_to_death <- ifelse(vital_status == "Dead", sample(100:2000, n_samples, replace = TRUE), NA_integer_)
  days_to_last_follow_up <- ifelse(is.na(days_to_death), sample(100:2000, n_samples, replace = TRUE), NA_integer_)
  clinical <- data.frame(sample_type = sample_type,
                         vital_status = vital_status,
                         days_to_death = days_to_death,
                         days_to_last_follow_up = days_to_last_follow_up,
                         stringsAsFactors = FALSE,
                         row.names = sample_ids)
  saveRDS(counts, counts_path)
  saveRDS(clinical, clinical_path)
  message("Synthetic data created. Replace them with real data for real results.")
}

# Read inputs
message("Loading data...")
counts <- readRDS(counts_path)
clinical <- readRDS(clinical_path)

# Basic checks
if (is.null(colnames(counts))) stop("Counts matrix must have column names (sample IDs).")
if (is.null(rownames(counts))) rownames(counts) <- paste0("Gene", seq_len(nrow(counts)))
if (is.null(rownames(clinical))) stop("Clinical data must have row.names matching sample IDs.")

# Prepare sample groups
message("Preparing sample groups...")
clinical <- clinical %>%
  mutate(sample_type = as.character(sample_type)) %>%
  mutate(sample_type = ifelse(is.na(sample_type), "Primary Tumor", sample_type)) %>%
  mutate(sample_type = factor(sample_type, levels = c("Solid Tissue Normal", "Primary Tumor"))) %>%
  mutate(group = ifelse(sample_type == "Primary Tumor", "Tumor", "Normal"))

# Keep only samples present in both clinical and counts
common_samples <- intersect(colnames(counts), rownames(clinical))
if (length(common_samples) == 0) stop("No common samples between counts and clinical (check row/col names).")
counts <- counts[, common_samples, drop = FALSE]
clinical <- clinical[common_samples, , drop = FALSE]

# Filter low-count genes
message("Filtering low-count genes...")
gene_filter <- rowSums(counts >= 10) >= 2  # relaxed for small placeholders
counts <- counts[gene_filter, , drop = FALSE]

# Normalization: prefer TCGAbiolinks function if available, else simple library-size normalization (log2 CPM)
message("Normalizing count data for DEA...")
norm_counts <- NULL
if (has_tcgabiolinks && exists("TCGAanalyze_Normalization")) {
  # If geneInfoHT not available, fall back to simple normalization
  if (exists("geneInfoHT")) {
    tryCatch({
      norm_counts <- TCGAanalyze_Normalization(tabDF = counts, geneInfo = geneInfoHT, method = "geneLength")
    }, error = function(e) {
      message("TCGAanalyze_Normalization failed: ", conditionMessage(e), " - falling back to log2 CPM.")
    })
  } else {
    message("geneInfoHT object not found; falling back to log2 CPM normalization.")
  }
}
if (is.null(norm_counts)) {
  # Simple CPM normalization followed by log2
  libsizes <- colSums(counts)
  cpm <- t(t(counts) / libsizes) * 1e6
  norm_counts <- log2(cpm + 1)
}

# Differential expression analysis
message("Running differential expression analysis...")
de_results <- NULL

# Try TCGAbiolinks DEA first
if (has_tcgabiolinks && exists("TCGAanalyze_DEA")) {
  tryCatch({
    # TCGAanalyze_DEA expects raw or normalized matrices depending on method; use counts if available
    de_results <- TCGAanalyze_DEA(mat1 = counts[, clinical$group == "Tumor", drop = FALSE],
                                  mat2 = counts[, clinical$group == "Normal", drop = FALSE],
                                  Cond1type = "Tumor",
                                  Cond2type = "Normal",
                                  method = "glmLRT")
    # Ensure consistent column names
    if ("adj.P.Val" %in% colnames(de_results) && !"FDR" %in% colnames(de_results)) {
      de_results$FDR <- de_results$adj.P.Val
    }
    if (!"logFC" %in% colnames(de_results) && "log2FoldChange" %in% colnames(de_results)) {
      de_results$logFC <- de_results$log2FoldChange
    }
  }, error = function(e) {
    message("TCGAanalyze_DEA failed: ", conditionMessage(e), " - will try limma or t-test fallback.")
    de_results <<- NULL
  })
}

# If TCGAbiolinks DEA not available or failed, try limma (recommended)
if (is.null(de_results) && has_limma) {
  tryCatch({
    message("Using limma for DEA (on log2 CPM normalized counts).")
    # ensure group factor has Normal as reference
    clinical$group <- factor(clinical$group, levels = c("Normal", "Tumor"))
    design <- model.matrix(~0 + group, data = clinical)
    colnames(design) <- sub("^group", "", colnames(design))
    # fit
    fit <- lmFit(norm_counts, design)
    if (ncol(design) == 2 && all(c("Normal", "Tumor") %in% colnames(design))) {
      contr <- makeContrasts(contrast = "Tumor - Normal", levels = design)
      fit2 <- contrasts.fit(fit, contr)
    } else {
      # fallback: use first vs second column
      contr <- diag(ncol(design))
      fit2 <- fit
    }
    fit2 <- eBayes(fit2)
    tt <- topTable(fit2, number = Inf, sort.by = "P")
    tt$genes <- rownames(tt)
    # standardize column names
    if (!"FDR" %in% colnames(tt) && "adj.P.Val" %in% colnames(tt)) tt$FDR <- tt$adj.P.Val
    de_results <- tt
  }, error = function(e) {
    message("limma DEA failed: ", conditionMessage(e), " - will use t-test fallback.")
    de_results <<- NULL
  })
}

# t-test fallback
if (is.null(de_results)) {
  message("Falling back to per-gene t-tests (quick-and-dirty).")
  group_samples <- split(colnames(norm_counts), clinical$group)
  if (!("Tumor" %in% names(group_samples) && "Normal" %in% names(group_samples))) {
    stop("Need both 'Tumor' and 'Normal' groups for DEA.")
  }
  tumor_idx <- group_samples$Tumor
  normal_idx <- group_samples$Normal
  gene_list <- rownames(norm_counts)
  pvals <- vapply(gene_list, function(g) {
    x <- as.numeric(norm_counts[g, tumor_idx])
    y <- as.numeric(norm_counts[g, normal_idx])
    if (length(unique(c(x,y))) <= 1) return(1)
    t_res <- try(t.test(x, y), silent = TRUE)
    if (inherits(t_res, "try-error")) return(1)
    t_res$p.value
  }, numeric(1))
  logFC <- vapply(gene_list, function(g) {
    mean(as.numeric(norm_counts[g, tumor_idx]), na.rm = TRUE) - mean(as.numeric(norm_counts[g, normal_idx]), na.rm = TRUE)
  }, numeric(1))
  df <- data.frame(logFC = logFC, P.Value = pvals, FDR = p.adjust(pvals, method = "BH"), genes = gene_list, stringsAsFactors = FALSE, row.names = NULL)
  de_results <- df %>% arrange(FDR)
}

# Ensure genes column exists and results are a data frame
if (!is.data.frame(de_results)) de_results <- as.data.frame(de_results)
if (!"genes" %in% colnames(de_results)) de_results$genes <- rownames(de_results)
if (!"FDR" %in% colnames(de_results) && "adj.P.Val" %in% colnames(de_results)) de_results$FDR <- de_results$adj.P.Val
if (!"logFC" %in% colnames(de_results) && "log2FoldChange" %in% colnames(de_results)) de_results$logFC <- de_results$log2FoldChange

# Save DEA results
message("Saving DEA results...")
write.csv(de_results, file = "results/BRCA_DEA_results.csv", row.names = FALSE)
if (nrow(de_results) > 0) {
  top20 <- de_results %>% arrange(FDR) %>% slice_head(n = 20)
  write.csv(top20, file = "results/BRCA_DEA_top20.csv", row.names = FALSE)
}

# Volcano plot: use EnhancedVolcano if available, else simple ggplot
message("Creating volcano plot...")
if ("logFC" %in% colnames(de_results) && "FDR" %in% colnames(de_results)) {
  volcano_png <- "results/BRCA_volcano.png"
  png(volcano_png, width = 1400, height = 1000, res = 150)
  if (has_EnhancedVolcano) {
    EnhancedVolcano::EnhancedVolcano(de_results,
                                    lab = de_results$genes,
                                    x = "logFC",
                                    y = "FDR",
                                    title = "TCGA-BRCA: Tumor vs Normal",
                                    pCutoff = 0.05,
                                    FCcutoff = 1.0,
                                    pointSize = 2.0,
                                    labSize = 3.0)
  } else {
    # basic ggplot volcano
    dfv <- de_results %>% mutate(logP = -log10(FDR + 1e-300))
    ggplot(dfv, aes(x = logFC, y = logP)) +
      geom_point(alpha = 0.6) +
      theme_minimal() +
      xlab("logFC") + ylab("-log10(FDR)") +
      ggtitle("Volcano: Tumor vs Normal")
  }
  dev.off()
  message("Volcano saved to ", volcano_png)
} else {
  message("Skipping volcano: required columns not present.")
}

# Heatmap of top genes
message("Selecting top genes for heatmap...")
if (nrow(de_results) >= 10) {
  top_n <- min(50, nrow(de_results))
  key_genes <- de_results %>% arrange(FDR) %>% slice_head(n = top_n) %>% pull(genes) %>% unique()
  present_genes <- intersect(key_genes, rownames(norm_counts))
  if (length(present_genes) >= 2) {
    log_counts <- norm_counts[present_genes, , drop = FALSE]
    # remove zero variance
    var_filter <- apply(log_counts, 1, var, na.rm = TRUE) > 0
    log_counts <- log_counts[var_filter, , drop = FALSE]
    scaled_matrix <- t(scale(t(log_counts)))
    heat_png <- "results/BRCA_heatmap_top50.png"
    png(heat_png, width = 1600, height = 1200, res = 150)
    if (has_ComplexHeatmap) {
      clinical_annotation <- ComplexHeatmap::HeatmapAnnotation(Group = clinical$group,
                                                               col = list(Group = c("Tumor" = "#E41A1C", "Normal" = "#377EB8")))
      ComplexHeatmap::Heatmap(scaled_matrix,
                              name = "Z-Score",
                              top_annotation = clinical_annotation,
                              show_row_names = TRUE,
                              show_column_names = FALSE,
                              cluster_rows = TRUE,
                              cluster_columns = TRUE)
    } else if (has_pheatmap) {
      annotation_col <- data.frame(Group = clinical$group)
      rownames(annotation_col) <- rownames(clinical)
      pheatmap::pheatmap(scaled_matrix, annotation_col = annotation_col, show_rownames = TRUE, show_colnames = FALSE)
    } else {
      # Basic heatmap as fallback
      heatmap(scaled_matrix, labRow = rownames(scaled_matrix), labCol = colnames(scaled_matrix), main = "Top genes heatmap")
    }
    dev.off()
    message("Heatmap saved to ", heat_png)
  } else {
    message("Not enough top genes present in normalized matrix for heatmap.")
  }
} else {
  message("Not enough DE results to create heatmap.")
}

# Pathway enrichment: try TCGAanalyze_EAcomplete if available, else skip
message("Running pathway enrichment on significant DE genes...")
sig_genes <- de_results %>% filter(FDR <= 0.05 & abs(logFC) >= 1) %>% pull(genes)
if (length(sig_genes) >= 5) {
  if (has_tcgabiolinks && exists("TCGAanalyze_EAcomplete")) {
    tryCatch({
      enrich <- TCGAanalyze_EAcomplete(geneList = sig_genes, organism = "hsa")
      write.csv(enrich, file = "results/BRCA_enrichment_results.csv", row.names = FALSE)
      message("Enrichment results saved.")
    }, error = function(e) {
      message("Enrichment function failed: ", conditionMessage(e))
    })
  } else {
    message("TCGAanalyze_EAcomplete not available; skipping enrichment.")
  }
} else {
  message("Not enough significant genes for enrichment analysis. Found: ", length(sig_genes))
}

# Survival analysis
message("Preparing survival analysis...")
req_surv_cols <- c("vital_status", "days_to_death", "days_to_last_follow_up")
if (any(req_surv_cols %in% colnames(clinical))) {
  surv_data <- clinical %>%
    mutate(
      survival_time = dplyr::coalesce(
        suppressWarnings(as.numeric(days_to_death)),
        suppressWarnings(as.numeric(days_to_last_follow_up))
      ),
      status = ifelse(toupper(vital_status) %in% c("DEAD", "DECEASED"), 1, 0)
    ) %>%
    filter(!is.na(survival_time))
  # Focus on Tumor samples
  surv_data <- surv_data %>% filter(group == "Tumor")
  
  if (nrow(surv_data) >= 6 && length(sig_genes) >= 1) {
    # select a representative gene (top by FDR)
    top_gene <- de_results %>% arrange(FDR) %>% pull(genes) %>% .[1]
    if (!is.na(top_gene) && top_gene %in% rownames(norm_counts)) {
      surv_data$expression_score <- as.numeric(norm_counts[top_gene, rownames(surv_data)])
      median_score <- median(surv_data$expression_score, na.rm = TRUE)
      surv_data$risk_group <- ifelse(surv_data$expression_score >= median_score, "High", "Low")
      
      fit <- survival::survfit(survival::Surv(survival_time, status) ~ risk_group, data = surv_data)
      png("results/BRCA_survival_plot.png", width = 1400, height = 1000, res = 150)
      try({
        ggsurvplot(fit, data = surv_data, risk.table = TRUE, pval = TRUE,
                   title = paste0("BRCA Survival by ", top_gene, " expression"))
      }, silent = TRUE)
      dev.off()
      writeLines(capture.output(summary(fit)), "results/BRCA_survival_summary.txt")
      message("Survival plot and summary saved.")
    } else {
      message("Top gene not found in normalized counts; skipping survival plotting.")
    }
  } else {
    message("Not enough survival data to fit a model. Need at least 6 tumor samples with survival info and >=1 significant gene.")
  }
} else {
  message("Required clinical survival columns not present. Skipping survival analysis.")
}

# Save session info
writeLines(capture.output(sessionInfo()), "results/session_info.txt")
message("Analysis complete. Review outputs in the results/ directory.")
