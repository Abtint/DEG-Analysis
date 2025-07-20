
# DEG-Analysis

A minimal R pipeline that **downloads raw-count matrices from any number of GEO Series**, merges them, runs **DESeq2** with batch correction, and produces the **four most common RNA-seq QC plots**—all from a single function call.

**But, hold a second. What is DEG Analysis and why do we run it?**
Differentially expressed genes (DEGs) are genes whose activity levels rise or fall significantly when you compare one biological condition to another. For example, tumor tissue versus healthy tissue or drug-treated cells versus untreated controls. Pinpointing these expression shifts is powerful because it highlights the molecular circuitry most altered in a disease or response to a therapy. DEGs often identify key drivers of pathology, reveal early biomarkers that can enhance diagnosis or prognosis, and highlight pathways that are promising targets for new drugs. In short, DEG analysis converts raw transcriptome data into an actionable list of genes that advances our understanding of disease mechanisms and accelerates the search for more precise and effective treatments.

---

## Clone the repository

Choose any method below, then switch into the project folder:

```bash
# HTTPS
git clone https://github.com/Abtint/DEG-Analysis.git

# SSH  (requires a configured SSH key)
git clone git@github.com:Abtint/DEG-Analysis.git

# GitHub CLI
gh repo clone Abtint/DEG-Analysis

cd DEG-Analysis
````

---

## Script in this repo

| File               | Purpose                                                                           |
| ------------------ | --------------------------------------------------------------------------------- |
| **`DE_Analysis.R`** | All-in-one pipeline: download counts → merge → DESeq2 → QC plots → tidy DEG table |

---

## Installation

```r
install.packages(c(
  "GEOquery",   # fetch GEO data
  "DESeq2",     # differential expression
  "data.table",
  "dplyr",
  "tibble",
  "ggplot2"
))
```

---

## Quick start

```r
# 1. Load the script
source("geo_deg_qc.R")

# 2. List the GEO Series you want
gse_ids <- c("GSE111111", "GSE222222")      # replace with real IDs

# 3. Build a sample metadata data.frame
meta <- data.frame(
  sample_id = c("GSM0001","GSM0002","GSM0003","GSM0004"),
  condition = c("case","case","control","control"),  # factor to test
  dataset   = c("GSE111111","GSE111111","GSE222222","GSE222222")
)

# 4. Run the pipeline
res <- run_geo_deg(
         gse_ids      = gse_ids,
         sample_sheet = meta,
         contrast     = c("condition","case","control")  # adjust labels as needed
       )

# 5. Explore results
head(res$deg)                       # tidy DEG table
write.csv(res$deg, "DEGs.csv", row.names = FALSE)

# QC figures are saved automatically under QC/
```

---

## What the pipeline does

| Step            | Action                                                         | Output                    |
| --------------- | -------------------------------------------------------------- | ------------------------- |
| **Download**    | Fetch supplemental count files for each GSE ID with `GEOquery` | `geo_downloads/…`         |
| **Merge**       | Union of gene lists, zero-padded where missing                 | consolidated count matrix |
| **Gene filter** | Drop genes with < 10 total reads (default)                     | cleaner matrix            |
| **DESeq2**      | Design `~ dataset + condition`; shrinkage LFCs                 | `res$deg`                 |
| **QC #1**       | Library-size bar plot (red < 1 M reads)                        | `QC/01_library_sizes.png` |
| **QC #2**       | VST-based PCA for outlier scan                                 | `QC/02_PCA.png`           |
| **QC #3**       | Dispersion fit plot                                            | `QC/03_dispersion.png`    |
| **QC #4**       | MA plot                                                        | `QC/04_MA.png`            |

---

## Input requirements

* **GSE IDs** – any GEO Series that provides raw-count matrices.
* **sample\_sheet** – data frame with three columns:

  * `sample_id` exactly matches column names in count matrix
  * `condition` factor you want to contrast (e.g., *case* vs *control*)
  * `dataset`  the GSE each sample belongs to (enables batch term)

---

## Output

* **`QC/`** Four PNG QC plots.
* **`res$deg`** Tidy DESeq2 results (gene, log2FC, padj, …).
* **`res$dds`** Full DESeq2 object for downstream analyses.

---

## Common problems

| Problem                                  | Symptom                                                                   | Resolution                                                                                                                                                                                                                                            |
| ---------------------------------------- | ------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Count file not detected                  | `fetch_geo_counts` throws “No count file found”                           | Check the supplementary files on the GEO page. If the counts are inside a compressed archive extract them first or change the `pattern` argument to match the real file name, for example `"counts_matrix.csv"`.                                      |
| Mixed file formats                       | Script stops with read errors or wrong column count                       | Inspect each downloaded file. Set the `pattern` argument to capture the correct file or adjust the separator inside `fetch_geo_counts` to comma or tab as needed.                                                                                     |
| Sample identifiers do not match          | Function stops with “column names do not match sample\_sheet\$sample\_id” | Open the count matrix and your metadata and be sure the column names are identical, including upper and lower case. Rename columns or rows so they match exactly.                                                                                     |
| Libraries too small                      | Library size plot shows red bars and dispersion plot looks noisy          | Remove or resequence samples with fewer than one million reads, or raise `min_lib` if your project demands deeper coverage. Rerun the pipeline after dropping weak libraries from the metadata.                                                       |
| Strong batch effect                      | PCA clusters by dataset not by condition                                  | Increase biological replicates if possible, or include additional covariates in the design, for example `~ dataset + sex + condition`. Consider running a batch correction method such as limma removeBatchEffect on VST counts before visualization. |
| Inconsistent gene identifiers            | Many blank rows in results or zero counts after merge                     | Make sure every dataset reports gene symbols or every dataset reports Ensembl identifiers, not a mix. Convert with biomaRt or another gene mapping tool before running the merge step.                                                                |
| Memory limit reached                     | R session crashes or freezes when merging very large matrices             | Run the script on a system with more memory, or pre filter each matrix to protein coding genes only, or work chromosome by chromosome and combine results later.                                                                                      |
| Internet or firewall blocks GEO download | `getGEOSuppFiles` fails with curl or timeout message                      | Download the files manually through a browser, place them inside `geo_downloads/GSEXXXXX`, and skip the automatic download by commenting out the call to `getGEOSuppFiles`.                                                                           |
| Blank or empty QC figures                | PNG files created but plots appear white                                  | Verify that the count matrix contains more than one sample per condition and that the low count filter did not remove all genes. Lower `min_counts` if necessary.                                                                                     |

---

## License

MIT.

---

## Citation

Tondar A. 2025. https://github.com/Abtint/DEG-Analysis

```
```

