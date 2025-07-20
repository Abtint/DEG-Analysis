
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

## License

MIT.

---

## Citation

Tondar A. 2025. https://github.com/Abtint/DEG-Analysis

```
```

