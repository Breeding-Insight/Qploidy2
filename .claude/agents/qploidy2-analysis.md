---
description: Step-by-step workflow for running Qploidy2: batch-effect check, standardization, BAF model selection, HMM copy number estimation, filtering, visualization, and export.
---

# Qploidy2 Analysis Workflow

This skill covers the full analysis pipeline after input preparation (see `qploidy2-input-prep` for creating `data`, `genos`, and `genos.pos`).

## 1. Raw Data Quality Check

Inspect raw allele intensities and heterozygosity before standardization:

```r
library(Qploidy2)

samples <- unique(data$SampleName)
plot_raw(data = data, geno.pos = geno.pos, sample = samples[1])
```

The `plot_raw` output shows (top to bottom):
- Heterozygosity per genomic window
- Allele ratio (B/(A+B)) vs. genomic position
- Per-chromosome ratio histograms
- Overall ratio histogram
- Total read depth / intensities

**If histograms already show clear peaks** → skip to BAF Model Selection.
**If peaks are unclear** → proceed with standardization.

**Check heterozygosity:** Samples with < 5% heterozygosity are treated as haploid (CN = 1) and cannot be distinguished from inbred lines.

## 2. Check for Batch Effects

Batch effects arise from different DNA extraction dates, plate runs, or library preparation protocols. **Always check before standardizing.**

```r
# Build or load a passport file with a batch-identifying column
passport_file <- data.frame(
  SampleName = unique(data$SampleName),
  batch       = c(...)  # e.g., plate number, extraction year
)

p <- pca_plot(
  input           = data,
  geno.pos        = geno.pos,
  passport_file   = passport_file,
  group_column    = "batch",
  sampleID_column = "SampleName",
  col2use         = "R"
)
p$PC1xPC2
```

Also inspect a single marker across all samples:

```r
plot_xy_with_ploidy_guides(data, ploidy = 4, marker = 1, samples = passport_file, palette = "auto")
```

**If samples cluster by batch** → run `standardize()` separately for each batch, then merge with `merge_qploidy_datas()`.

## 3. Standardization

Inspect total depth distribution to set depth filters:

```r
hist(data$R)
summary(data$R)
```

Run standardization:

```r
qploidy_std <- standardize(
  data                   = data,
  genos                  = genos,         # reference samples (known ploidy or all)
  geno.pos               = geno.pos,
  ploidy.standardization = 4,             # known or predominant ploidy
  threshold.n.clusters   = 4,            # min dosage clusters per marker (max = ploidy+1)
  threshold.missing.geno = 0.9,           # max missing fraction per marker
  min.depth              = 30,            # for sequencing data
  max.depth              = 600,
  n.cores                = 2,
  out_filename           = "standardization.tsv.gz",
  verbose                = TRUE
)
qploidy_std  # view filtering summary
```

**Key parameters:**

| Argument | Guidance |
|----------|----------|
| `ploidy.standardization` | Set to known or most frequent ploidy |
| `threshold.n.clusters` | Start with `ploidy+1`; relax (e.g. 4) if too many markers drop |
| `min.depth` / `max.depth` | Based on `hist(data$R)`; center around the mean |
| `threshold.missing.geno` | 0.9 = allow up to 90% missing per marker |

Visually inspect standardization quality:

```r
plot_geno_by_marker(qploidy_std, alpha = 0.1, marker = c(1, 10, 20, 30))
```

Reload a saved standardization in a new session:

```r
qploidy_std <- read_qploidy_standardization("standardization.tsv.gz")
```

### Visualize standardization effect per sample

```r
plot_qploidy_standardization(x = qploidy_std, sample = "SampleName",
                             type = c("Ratio_hist_overall", "BAF_hist_overall"))
plot_qploidy_standardization(x = qploidy_std, sample = "SampleName",
                             type = c("zscore", "R"))
```

### Optional: Chromosome-level Z-score screen

```r
chromosome_level_test_plot_qploidy(input = qploidy_std,
                                   selected_samples = samples[1:12],
                                   col2test = "z")
```

## 4. BAF Model Selection

Select the best statistical distribution model for BAF for one reference sample:

```r
selected_model <- select_best_baf_model(
  qploidy_standardization = qploidy_std,
  sample                  = "SampleName",
  cn_grid                 = c(2, 3, 4, 5, 6),  # candidate copy numbers
  bw                      = c(0.02, 0.04, 0.06, 0.08, 0.1),  # bandwidth grid
  plot                    = TRUE
)
selected_model$plot
head(selected_model$grid_results)
```

Use a sample with known ploidy to validate model selection.

## 5. HMM Copy Number Estimation — Single Sample

```r
res <- hmm_estimate_CN(
  qploidy_standarize_result = qploidy_std,
  sample_id                 = "SampleName",
  cn_grid                   = c(2, 3, 4, 5, 6)
)
res  # tabular summary
```

**Visualize karyotype:**

```r
plot_karyotype(res, sample_name = "SampleName", nrow = 2)
plot_karyotype(res, sample_name = "SampleName", nrow = 2, color_by = "post_max_CN")
plot_karyotype(res, sample_name = "SampleName", nrow = 2, color_by = "w_baf")
```

**Detailed diagnostic plot:**

```r
p <- plot_cn_track(hmm_CN = res, data = data, geno.pos = geno.pos,
                   sample_id = "SampleName", show_window_lines = FALSE, z_by_mean = TRUE)
p$arranged  # all panels
p$baf; p$z; p$cn; p$baf_dist  # individual panels
```

**Understanding `CN_reliability`:** ranges 0–1; combines HMM posterior × BAF weight. Calls relying only on Z-score (no heterozygous markers in window) are capped at 0.25 × post_max. Low values warrant manual inspection.

## 6. HMM Copy Number Estimation — All Samples

```r
multi_esti <- hmm_estimate_CN_multi(
  qploidy_standarize_result = qploidy_std,
  sample_ids                = "all",
  n_cores                   = 2,
  cn_grid                   = c(2, 3, 4, 5, 6)
)

compare_cn_track(hmm_CN = multi_esti, samples_to_plot = samples[1:30], facet_ncol = 8)
```

## 7. Filter by Quality Metrics

```r
multi_esti_filtered <- filter_hmm_CN(
  multi_esti,
  min_CN_reliability = 0.25,  # joint BAF × HMM confidence
  min_post_max       = 0.6    # raw HMM posterior threshold
)
```

- `min_CN_reliability = 0`: no filter; `1`: maximum stringency
- `min_post_max`: useful when keeping z-score-only calls but filtering poor HMM fits

## 8. Population Summary

```r
# Classify samples as euploid / aneuploid / segmental
types_df <- count_types(multi_esti_filtered)
str(types_df)  # named list of sample IDs per category

# Karyotype matrix overview
p_summ <- compare_cn_track_summary(
  multi_esti_filtered,
  samples_to_plot = c(types_df$aneuploid, types_df$segmental)
)
p_summ
```

## 9. Optional: Re-standardization

After identifying euploid samples via HMM, re-standardize using them as the reference for improved BAF resolution:

```r
re_std <- re_standardize(
  data                   = data,
  geno.pos               = geno.pos,
  hmm_CN_multi           = multi_esti_filtered,
  use_estimated_dosages  = TRUE,
  ploidy.standardization = 4,
  threshold.n.clusters   = 4,
  threshold.missing.geno = 0.9,
  min.depth              = 30, max.depth = 600,
  n.cores                = 2,
  out_filename           = "re_standardization.tsv.gz"
)

# Re-run HMM with new standardization
re_multi <- hmm_estimate_CN_multi(qploidy_standarize_result = re_std, sample_ids = "all",
                                  n_cores = 2, cn_grid = c(2,3,4,5,6))
re_multi_filtered <- filter_hmm_CN(re_multi, min_CN_reliability = 0.25, min_post_max = 0.6)
```

## 10. Refine a Single Sample

```r
res_tweak <- hmm_estimate_CN(
  qploidy_standarize_result = re_std,
  sample_id                 = "SampleName",
  min_snps_per_window       = 40,  # increase for more smoothing
  cn_grid                   = c(2, 3, 4, 5, 6)
)

# Update the multi-sample object
final_res <- update_hmm_CN_multi(re_multi_filtered, res_tweak, rm_sample = FALSE)
# Set rm_sample = TRUE to remove inconclusive samples
```

## 11. Dosage Calling

```r
dosages <- call_hmm_dosages(hmm_CN = final_res, n.cores = 2)

plot_karyotype(dosages, sample_name = "SampleName", nrow = 2, color_by = "dosage")
plot_karyotype(dosages, sample_name = "SampleName", nrow = 2, color_by = "post_max_dosage")
```

## 12. Export Results

### Save / reload HMM results

```r
write_hmm_CN(final_res, prefix = "my_project")
# Creates: my_project_by_marker.csv.gz, my_project_by_window.csv.gz, my_project_params.rds

final_res <- read_hmm_CN(
  by_marker_file = "my_project_by_marker.csv.gz",
  by_window_file = "my_project_by_window.csv.gz",
  params_file    = "my_project_params.rds"
)
```

### Summary tables

```r
summ_chr    <- summarize_cn_mode(final_res, level = "chromosome")
summ_sample <- summarize_cn_mode(final_res, level = "sample")

write.csv(summ_chr,    "CN_by_chromosome.csv")
write.csv(summ_sample, "CN_by_sample.csv")
write.csv(final_res$by_window, "CN_by_window.csv")
write.csv(final_res$by_marker, "CN_by_marker.csv")
```

### Karyotype notation

```r
k_not <- karyotype_notation(dosages)
head(k_not)
# Columns: SampleName, Ploidy, ChrID, ChrCN, Segment (ISCN-inspired notation)
```

### VCF export

```r
export_VCF(dosages, file = "output_CN_dosages.vcf")
# GT field: e.g. 0/0/1/1 = dosage 2 in a tetraploid
# PMD field: probability of the most likely dosage call
```

### MAPpoly export

```r
for_mappoly <- export_mappoly(final_res)
# Pass to MAPpoly::filter_aneuploidy() to remove aneuploid chromosomes before mapping
```

## 13. Interactive Curation with GenoBrew

Upload `write_hmm_CN()` output files to [GenoBrew](https://github.com/Breeding-Insight/GenoBrew) Shiny app → CNV profiles tab for interactive inspection and per-sample HMM re-runs.

## Quick Reference: Key Function Sequence

```
plot_raw()                          # raw data QC
pca_plot()                          # batch effect check
standardize()                       # standardize intensities
plot_geno_by_marker()               # inspect standardization
select_best_baf_model()             # choose BAF model (one sample)
hmm_estimate_CN_multi()             # run HMM for all samples
filter_hmm_CN()                     # quality filter
count_types()                       # summarize euploid/aneuploid/segmental
re_standardize()  [optional]        # improve with euploid reference
call_hmm_dosages()                  # call per-locus dosages
export_VCF() / write_hmm_CN()       # export results
```
