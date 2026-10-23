---
description: Guide for preparing the three required input data.frames (data, genos, genos.pos) for the Qploidy2 R package from VCF, Axiom, or Illumina array files.
---

# Qploidy2 Input Preparation

This skill guides you through preparing the three required input objects for `Qploidy2`:
- **`data`**: allele intensities or read counts per marker and sample
- **`genos`**: genotype dosages for reference samples
- **`genos.pos`**: genomic positions for each marker

## 1. Prerequisites — Check Before Starting

**Qploidy2 WORKS with:**
- Axiom or Illumina genotyping arrays
- Targeted sequencing (DArTag, GTseq, AgriSeq)
- Whole-genome sequencing (with sufficient coverage and filtering)
- Same library preparation protocol for all samples
- At least ~60 samples with known ploidy OR a known predominant ploidy
- Heterozygous samples

**Qploidy2 does NOT work with:**
- RADseq / GBS data (sparse, inconsistent marker coverage)
- Mixed sequencing batches in a single run — analyze each batch separately, then merge with `merge_qploidy_datas()`
- Exclusively inbred/homozygous lines
- First-generation colchicine-treated individuals (use their progeny instead)

## 2. Installation

```r
devtools::install_github("Breeding-Insight/Qploidy2")
library(Qploidy2)
```

## 3. Create the `data` Object

### From VCF (requires AD format field)

```r
data <- qploidy_read_vcf(vcf_file = "your_file.vcf")
head(data)
# Columns: MarkerName, SampleName, X (ref allele), Y (alt allele), R (X+Y), ratio (Y/R)
```

### From Axiom Array

```r
# summary file from apt-probeset-genotype (APT platform)
data <- read_axiom("AxiomGT1_summary.txt")
```

### From Illumina Array

```r
# File must contain [Data] section with columns: SNP Name, Sample ID, Theta, X, Y, etc.
data <- read_illumina_array(
  "data/illumina/Set1.txt",
  "data/illumina/Set2.txt"  # merge multiple files
)
```

All three functions return the same structure: `MarkerName`, `SampleName`, `X`, `Y`, `R`, `ratio`.

## 4. Create the `genos` Object — Choose Path (1) or (2)

### Path (1): Subset with Known Ploidy (preferred)

Filter your samples to those with known ploidy for use as the standardization reference:

```r
known_tetraploids <- c("Sample1", "Sample2", ...)  # your known-ploidy sample names
genos <- genos_full[which(genos_full$SampleName %in% known_tetraploids), ]
```

### Path (2): Use All Samples (most common ploidy unknown)

Run Qploidy2 on all samples first, then after identifying euploid individuals, re-standardize:

```r
# Use entire dataset as reference for first run
data_reference <- data
```

## 5. Call Dosages

### From VCF with GT field already present

```r
genos     <- qploidy_read_vcf("your_file.vcf", geno = TRUE)
genos.pos <- qploidy_read_vcf("your_file.vcf", geno.pos = TRUE)

# If MarkerName column has NAs (missing VCF ID field):
genos.pos$MarkerName <- paste0(genos.pos$Chromosome, "_", genos.pos$Position)

# Filter to known-ploidy samples (Path 1):
genos <- genos[which(genos$SampleName %in% known_tetraploids), ]
```

### From VCF without GT — Use updog (sequencing data)

```r
library(tidyr)
library(updog)

ref  <- pivot_wider(data[, 1:3], names_from = SampleName, values_from = X)
size <- pivot_wider(data[, c(1,2,5)], names_from = SampleName, values_from = R)
# convert to matrix...

multidog_obj <- multidog(
  refmat = ref, sizemat = size,
  model  = "norm",  # adjust to your population
  ploidy = 4,       # most common/known ploidy
  nc     = 2
)

genos <- data.frame(
  MarkerName = multidog_obj$inddf$snp,
  SampleName = multidog_obj$inddf$ind,
  geno       = multidog_obj$inddf$geno,
  prob       = multidog_obj$inddf$maxpostprob
)
```

### From Axiom — Use fitPoly

```r
library(fitPoly)
saveMarkerModels(
  ploidy     = 4,
  data       = data_reference,
  filePrefix = "fitpoly_output/prefix",
  ncores     = 2
)

library(vroom)
fitpoly_scores <- vroom("fitpoly_output/prefix_scores.dat")

genos <- data.frame(
  MarkerName = fitpoly_scores$MarkerName,
  SampleName = fitpoly_scores$SampleName,
  geno       = fitpoly_scores$maxgeno,
  prob       = fitpoly_scores$maxP
)
```

### Genomic positions for array data (not in summary files)

Provide a separate file with chromosome and position per marker:

```r
genos.pos_file <- read.table("marker_positions.txt", header = TRUE)
genos.pos <- data.frame(
  MarkerName = genos.pos_file$probes,
  Chromosome = genos.pos_file$chr,
  Position   = genos.pos_file$pos
)
```

## 6. Final Check

By the end you should have three objects ready:

| Object      | Required columns |
|-------------|------------------|
| `data`      | MarkerName, SampleName, X, Y, R, ratio |
| `genos`     | MarkerName, SampleName, geno, prob |
| `genos.pos` | MarkerName, Chromosome, Position |

Proceed to the Qploidy2 analysis workflow (`qploidy2-analysis` skill).

## Key Tips

- **Batch effect**: Never mix samples from different library preparation batches in one run.
- **Path (2) users**: After the first run, go back and re-run dosage calling using only the identified euploids, then re-standardize.
- **Compute tip**: Dosage calling (updog / fitPoly) is slow — run on HPC with multiple cores.
- **Recommended dosage callers**: `updog` or `polyRAD` for sequencing; `fitPoly` for arrays.
