# 🧬 RAD-Seq Analysis Pipeline using STACKS and R

This repository contains a comprehensive bioinformatics pipeline for processing and analyzing **Restriction site-associated DNA sequencing (RAD-Seq)** data.  
The workflow takes raw FASTQ files, processes them through the **STACKS `denovo_map.pl`** pipeline across various parameter sets, and performs downstream population genetics analyses and visualizations in **R**, primarily using the `dartR` ecosystem.

---

## 📖 Overview

The pipeline automates the analysis of genetic data for **population structure**, **diversity**, and **differentiation**.  
It is designed to test different STACKS parameters and identify the optimal assembly configuration for a given dataset.

### Workflow structure
The pipeline is divided into three main parts:

1. **Data Pre-processing and STACKS Execution**  
   Runs the core STACKS pipeline (`process_radtags` and `denovo_map.pl`) across user-defined parameters.
2. **Downstream Genetic Analysis in R**  
   Performs population genetic analyses using R (`dartR`, `LEA`, `poppr`, etc.).
3. **Run Metrics Extraction**  
   Extracts and aggregates key metrics from STACKS logs and VCF files for parameter comparison.

---

## ⚙️ Pipeline Workflow

### **Part 1 – STACKS Data Processing (Main STACKS Pipe DArT.R – Bash section)**

Automates the core data processing steps.

**Steps:**
1. **Concatenate FASTQ files** (optional) – merges technical replicates using a mapping CSV.
2. **Define datasets & parameters** – configure datasets, directories, popmaps, and parameter grids (`-M`, `-m`, `-n`).
3. **Clean input files** – removes carriage returns from popmap/catalog files.
4. **Run `process_radtags`** – demultiplex and clean FASTQ data.
5. **Execute `denovo_map.pl`** – runs the full STACKS pipeline iteratively for each parameter combination, producing `populations.snps.vcf`.

---

### **Part 2 – Downstream Analysis (Main STACKS Pipe DArT.R – R section)**

Analytical core of the pipeline, executed after STACKS completion.

**Steps:**

- **Load and filter data:** recursively load all `populations.snps.vcf` files as `genlight` objects using `vcfR`.
- **Data quality control:**
  - Assess read depth and call rates.  
  - Manage technical replicates and calculate mismatch rates.  
  - Apply filters for call rate, monomorphic loci, and HWE.
- **Population structure (SNMF):**
  - Run `LEA::snmf` for K = 1–30.
  - Determine optimal K (lowest mean cross-entropy).
  - Generate admixture plots and cross-entropy curves.
- **Hierarchical analysis (AMOVA):**
  - Quantify variance between SNMF clusters and populations.  
  - For UK datasets, perform additional AMOVA by region.
- **Genetic differentiation (Fst):**
  - Compute pairwise Fst values between populations.  
  - Export genetic distance matrices for tools like Core Hunter.
- **Diversity metrics:**
  - Calculate observed heterozygosity (Ho).  
  - Count private alleles per population.
- **Final summary table:**  
  Compile key metrics (best K, loci counts, heterozygosity, replicate distance) into  
  `all_runs_summary.csv`.

---

### **Part 3 – STACKS Run Metrics Extraction (`extract_STACKS_run_metrics.md`)**

Collects and aggregates run statistics from STACKS output.

**Steps:**
1. Iterate through parameter sets.
2. Parse `gstacks.log.distribs` for coverage and phasing stats.
3. Parse `populations.sumstats_summary.tsv` for SNP summaries.
4. Use `vcftools` for depth, quality, and missingness.
5. Use `bcftools gtcheck` to compute sample distances.
6. Combine all metrics into CSV summaries (`combined_gstacks.csv`, `combined_loci_metrics.csv`, etc.).

---

## 🧰 Prerequisites

### **Software**
- [STACKS (v2.x)](https://catchenlab.life.illinois.edu/stacks/)
- [VCFtools](https://vcftools.github.io/)
- [BCFtools](http://samtools.github.io/bcftools/)
- [R ≥ 4.0](https://www.r-project.org/)

### **R Packages**
```r
required_packages <- c(
  "dplyr", "readr", "tidyr", "vcfR", "dartRverse", "poppr",
  "pegas", "hierfstat", "SNPRelate", "LEA", "ggplot2", "patchwork",
  "ragg", "scales", "RColorBrewer", "ggraph", "mapmixture",
  "rnaturalearth", "rnaturalearthdata", "sf", "devtools"
)
```

Install them in R with:
```r
install.packages(required_packages)
```

---

## 📁 Setup & Configuration

### **Directory Structure**
```
/path/to/project/
├── fastq_files/                 # Raw .fastq.gz files
├── stacks_results/              # STACKS output
├── scripts/                     # Pipeline scripts
└── metadata/
    ├── popmap_uk.txt
    ├── popmap_all.txt
    └── master_id_map.csv
```

### **Input Files**
- **FASTQ files:** gzipped raw reads per individual.  
- **Population maps:** tab-separated `sample_id → population` files.  
- **Catalog files:** optional list for STACKS catalog construction.  
- **Coordinate file:** `Master_popmap.csv` with `Site, Longitude, Latitude`.  
- **Region map:** `UK_population_regions.csv` for AMOVA region grouping.

---

## 🧾 Script Configuration

### **Bash Script**
Edit variables at the top of the main script:
- `datasets`, `fastq_dirs`, `popmaps`, `catalogs`
- `truncate_sizes`, `parameter_sets`
- `BASE_DIR`

### **R Script**
In the **ANALYSIS CONTROL PANEL** section:
- `DATASET_FILTER`, `TRUNCATION_FILTER`, `PARAMS_FILTER`
- Ensure paths like `BASE_DIR` match your directory layout.

---

## 🚀 How to Run

### **Part 1 – STACKS**
```bash
chmod +x ./main_script.sh
./main_script.sh
```

> ⚠️ This step can take several hours depending on dataset size and parameter combinations.

### **Part 2 – R Analysis**
```bash
Rscript ./downstream_analysis.R
```

Or run interactively in **RStudio**.

### **Part 3 – Metrics Extraction**
```bash
chmod +x ./extract_metrics.sh
./extract_metrics.sh
```

---

## 📊 Key Output Files

Each parameter combination creates a directory:
```
ustacks_out_[dataset]_t[size]_[label]/
```

**Contains:**
- `populations.snps.vcf` – filtered SNPs.  
- `*_SNMF_Composite_Mapmixture.png` – admixture map & cross-entropy plot.  
- `*_best_K_results.csv` – best K and locus count.  
- `*_UK_regional_AMOVA_summary.csv` – AMOVA results for UK runs.  
- `*_replicate_QC_analysis.csv` – mismatch rates between replicates.

**Global outputs:**
- `all_runs_summary.csv` – summary across all parameter sets.  
- `run_summary_output/` – combined metrics (`combined_gstacks.csv`, `combined_loci_metrics.csv`, etc.).

---

## 🧭 Citation
If you use this workflow, please cite:
- Catchen *et al.* (2013) **STACKS: Building and Genotyping Loci from Short-Read Sequences**.  
- Jombart *et al.* (2010) **Adegenet: A R Package for the Multivariate Analysis of Genetic Markers**.  
- Gruber *et al.* (2018) **dartR: An R package to facilitate analysis of SNP data generated from reduced representation genome sequencing**.

---

## 🧑‍💻 Author
**James Baker**
*Ecology and Evolution Group, UK Centre for Ecology & Hydrology (UKCEH)*
**Susheel Bhanu Busi**  
*Molecular Ecology Group, UK Centre for Ecology & Hydrology (UKCEH)*  
