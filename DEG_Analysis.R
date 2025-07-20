###############################################################################
# geo_deg_qc.R
#
# One-stop R pipeline to:
#   1. Fetch raw-count matrices from one or more GEO Series (GSE IDs)
#   2. Merge them (union of genes, fill missing with zero)
#   3. Run DESeq2 with a simple batch-aware design
#   4. Perform **core QC only**:
#        • Library-size & detected-gene bar plot
#        • Variance-stabilising transform (VST) PCA
#        • DESeq2 dispersion plot
#        • MA plot
#
# All QC figures are saved as PNG files in the directory supplied to `qc_dir`.
#
# Author: Abtin Tondar — 2025-07-19
###############################################################################

suppressPackageStartupMessages({
  library(GEOquery)
  library(DESeq2)
  library(data.table)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

# --------------------------------------------------------------------------- #
# UTILITIES
# --------------------------------------------------------------------------- #
save_png <- function(plot_obj, path, w = 7, h = 5) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  png(path, width = w, height = h, units = "in", res = 300)
  print(plot_obj)
  dev.off()
}

# --------------------------------------------------------------------------- #
# 1. DOWNLOAD A COUNT MATRIX FROM GEO
# --------------------------------------------------------------------------- #
fetch_geo_counts <- function(gse_id,
                             destdir  = "geo_downloads",
                             pattern  = "(count|matrix).*\\.(txt|csv|tsv)$",
                             gene_col = 1) {

  message(sprintf(">>> Fetching supplementary files for %s", gse_id))
  GEOquery::getGEOSuppFiles(gse_id,
                            makeDirectory = TRUE,
                            baseDir       = destdir,
                            fetch_files   = TRUE,
                            filter_regex  = pattern)

  files <- list.files(file.path(destdir, gse_id),
                      pattern = pattern,
                      full.names = TRUE, ignore.case = TRUE)

  if (!length(files))
    stop(sprintf("No count file found for %s", gse_id))

  f <- files[1]
  sep <- ifelse(grepl("\\.csv$", f, ignore.case = TRUE), ",", "\t")
  mat <- fread(f, sep = sep, data.table = FALSE)

  names(mat)[gene_col] <- "gene"
  rownames(mat) <- mat$gene
  mat$gene <- NULL

  message(sprintf("    Imported %d genes × %d samples", nrow(mat), ncol(mat)))
  as.matrix(mat)
}

# --------------------------------------------------------------------------- #
# 2. MERGE MULTIPLE COUNT MATRICES (UNION OF GENES)
# --------------------------------------------------------------------------- #
merge_counts <- function(mats) {
  all_genes <- Reduce(union, lapply(mats, rownames))
  merged <- do.call(cbind, lapply(mats, function(m) {
    missing <- setdiff(all_genes, rownames(m))
    if (length(missing))
      m <- rbind(m,
                 matrix(0L, nrow = length(missing), ncol = ncol(m),
                        dimnames = list(missing, colnames(m))))
    m[all_genes, , drop = FALSE]
  }))
  storage.mode(merged) <- "integer"
  merged
}

# --------------------------------------------------------------------------- #
# 3. CORE QC PLOTS
# --------------------------------------------------------------------------- #
qc_library_barplot <- function(counts, limit = 1e6) {
  df <- data.frame(sample = colnames(counts),
                   library_size = colSums(counts),
                   detected_genes = colSums(counts > 0))
  ggplot(df, aes(x = reorder(sample, library_size),
                 y = library_size / 1e6,
                 fill = library_size < limit)) +
    geom_col(colour = "black") +
    coord_flip() +
    scale_fill_manual(values = c(`FALSE` = "grey70", `TRUE` = "red"),
                      guide = "none") +
    labs(y = "Library size (millions)", x = NULL) +
    theme_minimal(base_size = 12)
}

qc_pca_plot <- function(dds, intgroup = "condition") {
  vsd <- vst(dds, blind = FALSE)
  pca <- prcomp(t(assay(vsd)))
  pc_df <- data.frame(pca$x[, 1:2],
                      colData(dds)[, intgroup, drop = FALSE] |> as.data.frame())
  ggplot(pc_df, aes(PC1, PC2, colour = .data[[intgroup]])) +
    geom_point(size = 3) +
    theme_minimal(base_size = 12)
}

# --------------------------------------------------------------------------- #
# 4. MASTER FUNCTION
# --------------------------------------------------------------------------- #
run_geo_deg <- function(gse_ids,
                        sample_sheet,
                        design      = ~ dataset + condition,
                        contrast    = c("condition", "case", "control"),
                        min_counts  = 10,
                        min_lib     = 1e6,
                        qc_dir      = "QC") {

  # -------- download & merge counts
  mats <- setNames(lapply(gse_ids, fetch_geo_counts), gse_ids)
  counts <- merge_counts(mats)
  counts <- counts[rowSums(counts) >= min_counts, ]
  message(sprintf(">>> %d genes retained after low-count filter", nrow(counts)))

  # -------- sample metadata check
  if (!all(colnames(counts) %in% sample_sheet$sample_id))
    stop("Column names of count matrix do not match sample_sheet$sample_id")
  sample_sheet <- sample_sheet |>
                  slice(match(colnames(counts), sample_id))
  rownames(sample_sheet) <- sample_sheet$sample_id

  # -------- library-size QC
  save_png(qc_library_barplot(counts, limit = min_lib),
           file.path(qc_dir, "01_library_sizes.png"),
           h = max(4, ncol(counts) / 5))

  # -------- DESeq2
  dds <- DESeqDataSetFromMatrix(counts, sample_sheet, design)
  dds <- DESeq(dds)
  res <- lfcShrink(dds, contrast = contrast, type = "ashr")

  # -------- PCA + diagnostics
  save_png(qc_pca_plot(dds, intgroup = "condition"),
           file.path(qc_dir, "02_PCA.png"))
  png(file.path(qc_dir, "03_dispersion.png"), 7, 5, "in", res = 300)
  plotDispEsts(dds); dev.off()
  png(file.path(qc_dir, "04_MA.png"), 7, 5, "in", res = 300)
  plotMA(res); dev.off()

  # -------- tidy DEG output
  deg <- as.data.frame(res) |>
         rownames_to_column("gene") |>
         arrange(padj)

  list(dds = dds,
       deg = deg)
}

###############################################################################

