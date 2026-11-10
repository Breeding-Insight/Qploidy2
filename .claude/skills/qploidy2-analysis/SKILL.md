---
name: qploidy2-analysis
description: "Step-by-step workflow for running Qploidy2: batch-effect check, standardization, BAF model selection, HMM copy number estimation, filtering, visualization, and export."
---

# Qploidy2 Analysis Workflow

This skill covers the full analysis pipeline after input preparation (see `qploidy2-input-prep` for creating `data`, `genos`, and `genos.pos`).

See also the full narrative tutorials: `docs/Qploidy_alfalfa_tutorial.Rmd` (grid-based BAF model selection + HMM) and `docs/Qploidy2_nQuack_tutorial.Rmd` (nQuack EM-based alternative, sequencing data only).

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

### Special case: Sex chromosomes

Sex chromosomes are handled in the HMM and downstream functions through `sex_system` and `sex_chr` (see **Section 5b**). Supported systems:

| `sex_system` | Females | Males | Sex inferred from |
|---|---|---|---|
| `"XY"` | XX | XY | Y present (CN >= 1) → male |
| `"ZW"` (birds, Lepidoptera) | ZW | ZZ | W present → female |
| `"X0"` | XX | X0 | X CN = ploidy → female; ploidy/2 → male |
| `"haplodiploid"` (bees, ants, wasps) | 2x | 1x | Sample ploidy (no `sex_chr` needed) |

Before running, check `plot_raw()` on the sex-chromosome markers to confirm the X (or Z) dosage difference between sexes is visible.

**Alternative for sequencing data with clear raw peaks:** sex-chromosome markers can also be analysed without standardization, using `use_values = c("ratio", "R")` in `hmm_estimate_CN_multi()`. Only do this when `plot_raw` shows well-defined peaks; it is **not appropriate for array data**, where standardization is always required.

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
| `filter_R` | `TRUE` = compute z-scores only for markers that passed the BAF filters (missing data + number of clusters) |

**Reference samples:** only samples present in `genos` (after filtering) are used as reference — for BAF cluster centers and for the per-marker mean/SD of `R` behind the z-score. Every sample is then transformed against that reference, so z-scores are relative to the known-ploidy samples.

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

## 4b. Alternative: nQuack EM-Based Route (sequencing data only)

For sequencing-based data (WGS, DArTag, GTseq, AgriSeq) with allele read counts (`X`, `Y`, `R`), `hmm_estimate_CN_nQuack()` replaces the grid-based model selection + HMM (steps 4–5) with nQuack's EM mixture-model fitting (Normal/Beta/Beta-Binomial), selecting the ploidy by BIC for **one** `(nQuack_dists, nQuack_type, nQuack_uniform)` combination per call. BIC is not comparable across distributions, types, or uniform settings — choose the combination on samples with known ploidy first, following the nQuack tutorials. **Not usable on intensity-based array data** (Axiom/Illumina have no read counts). Requires `install.packages("nQuack")` (>= 1.0.5, on CRAN).

```r
set.seed(42)  # allele sampling inside nQuack is stochastic — set a seed for reproducibility
res_nq <- hmm_estimate_CN_nQuack(
  qploidy_standarize_result = qploidy_std,
  sample_id                 = "SampleName",
  cn_grid                   = c(2, 3, 4, 5, 6),  # subset of 1:6; nQuack models 2:6 (see CN 1 notes below)
  max_cn_dev                = NULL,      # e.g. 2 to restrict a wide grid around exp_ploidy (see "Choosing cn_grid")
  nQuack_dists              = "normal",  # "normal" (default) | "beta" | "beta-binomial" — exactly one
  nQuack_type               = "fixed",   # "fixed" (default, recommended) | "fixed_2" | "fixed_3"
  nQuack_uniform            = 1,         # 1 (default, include a uniform mixture component) | 0
  verbose                   = TRUE
)

plot_karyotype(res_nq, sample_name = "SampleName", nrow = 2, color_by = "post_max_CN")
p <- plot_cn_track(hmm_CN = res_nq, data = data, geno.pos = geno.pos, sample_id = "SampleName",
                   show_window_lines = FALSE, z_by_mean = TRUE)
p$arranged
```

For all samples, `hmm_estimate_CN_multi()` auto-dispatches to the nQuack route whenever `nQuack_dists`, `nQuack_type`, and/or `nQuack_uniform` are passed via `...`:

```r
multi_esti_nq <- hmm_estimate_CN_multi(
  qploidy_standarize_result = qploidy_std,
  sample_ids                = "all",
  n_cores                   = 3,             # parallelism across samples
  cn_grid                   = c(2, 3, 4, 5, 6),
  nQuack_dists              = "normal",
  nQuack_type               = "fixed",
  nQuack_uniform            = 1
)
```

**CN 1 in the nQuack route** (`cn_grid` may include it, e.g. `1:6`; CN 0 is not allowed — it only makes sense for array background signal, use `hmm_estimate_CN()`):
- **CN 1** has no nQuack model: windows with heterozygous fraction ≥ `min_het_frac` rule it out; otherwise BAF is neutral and z decides. If the whole sample has heterozygous fraction ≤ `min_het_frac` (and `exp_ploidy` isn't given), it is treated as haploid (`exp_ploidy = 1`) and the nQuack global fit is skipped — useful for haplodiploid males.

**Grid vs. nQuack:**

| Feature | Grid (`hmm_estimate_CN()`) | nQuack (`hmm_estimate_CN_nQuack()`) |
|---|---|---|
| Platform | Arrays and sequencing | Sequencing only (read counts X, Y, R) |
| Speed | Faster | Slower |
| CN range | 0–10 (CN 1 has a BAF template; CN 0 for arrays) | 1–6; nQuack models 2–6, CN 1 via `min_het_frac` rule; no CN 0 |
| Reproducibility | Deterministic | Stochastic allele sampling — set `set.seed()` |

To inspect or export the nQuack input matrix directly (e.g. to call `nQuack::bestquack()` yourself), use `convert2nQuack(qploidy_standardization_object = qploidy_std)`.

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

### Interpreting `plot_cn_track()` and `plot_karyotype()`

Use this section to explain the two diagnostic figures to a user (e.g. for a figure legend). Describe what the code actually draws; do not improvise.

**`plot_cn_track()` layout.** The `$arranged` figure has a top row with three parts, then three tracks stacked below, faceted with one column per chromosome:

1. **BAF distribution (top left).** Density of the BAF (B-allele frequency, the proportion of reads or signal carrying the alternative allele) of all markers in the sample, pooled across chromosomes. Peak positions show the sample's overall dosage classes, e.g. 0/0.5/1 for a diploid and 0/0.25/0.5/0.75/1 for a tetraploid.
2. **Legends (top centre).** Colour scales for the BAF weight (`w_baf`, used in the BAF and z tracks) and for the CN reliability (used in the CN track). The CN legend reads `P(CN call)` instead when the object has no `CN_reliability` column (older objects).
3. **Parameters (top right).** Sample name; `Mean total depth/Z CN` (the `exp_ploidy`, i.e. the CN taken as the sample baseline); the CN grid tested; the minimum SNPs per window; the emission distribution.
4. **BAF track.** One dot per marker (not per window) at its genomic position, y = BAF.
   - Heterozygous markers cluster at dosage/CN, e.g. 0.5 in a diploid region. A trisomic region has bands near 0.33 and 0.67 and none at 0.5. A monosomic, hemizygous or LOH region has no heterozygous band.
   - The grey shape behind the dots is that chromosome's BAF density, drawn sideways.
   - Dots are coloured by the window's `w_baf`.
5. **Z-score track.** One dot per marker. The z-score is the marker's depth/intensity `R` compared with the same marker in the reference samples: `z = (R - mean(R_ref)) / sd(R_ref)` (`get_zscore()`). It is per marker against the reference panel, not normalised across autosomes.
   - Black horizontal segments are the window mean z (`z_by_mean = TRUE`; a loess line otherwise).
   - Dots are coloured by `w_baf`. Outliers, and with `depth_zero_as_x = TRUE` markers with `R == 0`, are drawn as `x`.
   - In a euploid sample, the window means stay flat at one level along the genome, with only random scatter. Do not call this a "uniform distribution".
   - A whole-chromosome gain or loss shifts all windows of that chromosome up or down together. A segmental event is a step change at the breakpoint.
6. **HMM CN track.** y = `CN_call` per window (Viterbi), drawn as horizontal segments. Single-marker windows are drawn as points. Colour = `CN_reliability`, 0–1, calculated as `post_max × baf_info × zbaf_agree × w_baf_norm`:
   - `post_max`: HMM posterior probability of the called state.
   - `baf_info`: 1 − normalised Shannon entropy of a 10-bin BAF histogram for the window. Sharp heterozygous peaks → 1; dispersed BAF → 0.
   - `zbaf_agree`: 1 − total variation distance between the per-state likelihoods from z and from BAF (softmax-normalised).
   - `w_baf_norm`: BAF weight normalised to 0–1.
   - Windows with no BAF support (`w_baf == 0`, depth-only calls) use the floor `z_no_baf_scale` (default 0.25), so they can reach at most 0.25 × `post_max`. They always look low-reliability. Low values warrant manual inspection or `filter_hmm_CN()`.

**`w_baf` (BAF weight).** `sqrt(min(1, n_het / q))`, where `n_het` is the number of heterozygous dosage calls in the window and `q` is a quantile (`het_quantile`) of `n_het` across windows. It tempers the BAF emission in the HMM (`ll = llz + w_baf * llb`), so windows with few heterozygous markers are called mostly from depth (z). When `exp_ploidy` ≈ 1, `w_baf = 1`.

**Sex chromosomes in `plot_cn_track()`** (when `sex_system`/`sex_chr` are set). The X/Z-specific and PAR groups appear as separate facets.
- The reference panel defines CN = baseline. When the reference contains only the homogametic sex (e.g. female-only standardization in an XY species), X-specific markers sit at baseline in females.
- In XY males, the X-specific region has about half the reference depth (negative z) and no heterozygous markers, so it is called CN 1 in a diploid. In polyploids, it is called the number of X copies (XXYY → 2, XXXY → 3).
- Y/W-specific markers have no depth in a homogametic reference, so their z is undefined and they get no CN calls. They do not appear in the CN track; Y presence is reported by `plot_karyotype()` / `karyotype_notation()` instead.

**`plot_karyotype()` layout.** One facet per chromosome. Within a facet, one vertical column is drawn per chromosome copy (x-axis = copy number; y-axis = position in Mb, reversed so the chromosome start is at the top). Each chromosome is not drawn as a single bar.
- **Segments:** consecutive regions with the same `CN_call` form a segment that shows as many columns as its CN.
- **Markers:** horizontal ticks on every copy of their segment. Ticks are black when `color_by = NULL`; otherwise they are coloured by `dosage`, `post_max_CN`, `post_max_dosage`, `w_baf` or `CN_reliability`.
- **Grey areas:** light grey (`grey80`) = gap with no markers; darker grey (`grey60`) = windows filtered out (`CN_call` set to `NA`, e.g. by `filter_hmm_CN()`).
- **Dashed outlines (ghost copies, `add_ghost = TRUE`):** copies that are expected but absent. Expected = max(sample ploidy mode, highest CN on that chromosome), or the expected sex-chromosome complement when `sex_system` is set (no ghost second X in an XY male). Dashed outlines therefore mark losses, not all aneuploid events:
  - Whole-chromosome loss → one dashed column along the whole chromosome.
  - Segmental loss → a dashed column only over the lost region.
  - Segmental gain → an extra solid column over the gained region, with a dashed outline where that copy is absent elsewhere on the chromosome.
  - Whole-chromosome gain → an extra solid column with no dashed outline.
- **Labels (`notation = TRUE`):** strips show the chromosome CN (`CN: X`), plus deviating segments if `show_segments = TRUE`. The subtitle shows `Ploidy: X`. With `sex_system`, it also shows the sex, the sex-chromosome karyotype and Y/W marker coverage, e.g. `Sex: male (XXYY), Y markers covered 5/9`, and `[check sex chromosomes]` when `SexChr` is `?` or inconsistent.

**Suggested figure pairing.** Show `plot_cn_track()` (upper) for the evidence behind the calls, and `plot_karyotype()` (lower) for the resulting karyotype. Read them together: a CN change in the karyotype should be supported by a z-shift and a matching BAF band pattern in the track, with reasonable `CN_reliability`. Depth-only changes (`w_baf` ≈ 0) are weaker evidence.

**Window size (`min_snps_per_window`):** sets the minimum segment length for z-score changepoint detection. Windows are never removed — a chromosome with fewer markers than the minimum becomes a single window. Set it per chromosome with a named vector, e.g. `min_snps_per_window = c(Chr01 = 300, Chr02 = 300, ChrX = 20)`; unnamed chromosomes use the automatic default (printed when verbose).

**Mixed-ploidy populations (`max_cn_dev`):** see **Choosing `cn_grid`** below.

### Choosing `cn_grid` (range of tested CN)

**Allowed values:**

| Route | Allowed `cn_grid` | Notes |
|---|---|---|
| Grid (`hmm_estimate_CN()`) | `0:10` | CN 1 scored with a BAF template (haploid/LOH). CN 0 (nullisomy) detected from z-scores only, through a z-only pre-pass, and returned as `NA` — meaningful for **array data** (probes keep background signal; tune `cn0_background_frac`). |
| nQuack (`hmm_estimate_CN_nQuack()`) | `1:6` | nQuack models 2–6; CN 1 via the `min_het_frac` rule. **CN 0 is not allowed** (sequencing data has no CN 0 background signal). |

**Larger range → more false positives.** The initial z mean of each CN state is a ramp spread across the sample's observed z range, so every extra state in `cn_grid` shrinks the z distance between neighbouring states. Small depth fluctuations or noise can then be called as a different CN. Also, each extra state is one more alternative the HMM can choose. Recommendations:
- Use the narrowest `cn_grid` that covers the copy numbers you actually expect.
- Include `0` only when nullisomy is plausible and the data are from arrays; include `1` only when haploids or LOH are expected.
- In **mixed-ploidy populations** that need a wide grid (e.g. `1:6`), set `max_cn_dev`. `exp_ploidy` is still estimated on the full grid, then each sample's grid is restricted to states within `max_cn_dev` of it, e.g. `max_cn_dev = 2` tests `1:4` for a 2x sample and `3:6` for a 5x sample. On the grid route CN 0 is kept whenever it is in `cn_grid`. The restricted grid is printed (verbose) and stored in `params$cn_grid`.
- States whose z means are not monotonic with CN, or that collapse onto a neighbour, are dropped automatically during EM; check `params$cn_grid` for the final grid.

## 5b. Sex Chromosomes

Pass `sex_system` and `sex_chr` to the HMM functions. They don't change the CN calling; they are stored in `params`, and `karyotype_notation()`, `plot_karyotype()`, and `count_types()` read them automatically.

```r
multi_esti <- hmm_estimate_CN_multi(
  qploidy_standarize_result = qploidy_std,
  sample_ids                = "all",
  n_cores                   = 2,
  cn_grid                   = 0:4,  # grid route on arrays: include 0 to detect Y absence (CN 0) in females
  min_snps_per_window       = c(ChrX = 20, ChrY = 10, ChrX_PAR = 20),  # small contigs
  sex_system                = "XY",
  sex_chr                   = list(X = "ChrX", Y = "ChrY", PAR = "ChrX_PAR")
)
```

`sex_chr` rules:
- Names: `X`, `Y`, `PAR` for `"XY"`; `Z`, `W`, `PAR` for `"ZW"`; `X` for `"X0"`; none for `"haplodiploid"`. Each element may hold several chromosome IDs.
- `X`/`Z` and `Y`/`W` hold markers specific to that chromosome; `PAR` holds pseudoautosomal markers.
- For `"XY"`/`"ZW"`, give `Y` (or `W`), or both `X` (or `Z`) and `PAR`. A missing group is derived from the PAR (Y = PAR − X; X = PAR − Y). When all three are given, X + Y = PAR is checked (`SexChrConsistent`).
- Detecting Y/W absence from Y/W-specific markers requires `0` in `cn_grid`, which is only available on the grid route (array data). On sequencing data, Y-specific markers have no reads in females and drop out, so also give `PAR` to derive Y = PAR − X.

Downstream effects:
- `karyotype_notation()` adds `Sex`, `SexChr` (e.g. `"XY"`, `"XXXY"`, `"ZW"`, `"X0"`), `SexChrConsistent`, `ChrType`, `ExpectedCN`. Sex chromosomes are excluded from `Ploidy` and compared to the expected complement (X + Y = ploidy), so one X in an XY male is not a deviation but an XXY diploid is.
- `plot_karyotype()` draws ghost copies from the expected complement, labels sex-chromosome strips, and adds sex to the subtitle.
- `count_types()` does not count males as aneuploid and adds a `$sex` table plus sex counts in the printout.
- `compare_cn_track_summary()` uses the same sex-aware classification and adds a `Sex` column (e.g. `M (XY)`, `F (XX)`, `M (XXXY)`, `?` if unknown) before the block labels. Set `sex_system = NULL` there to disable it.
- `print()` of a single-sample result shows the sex system; of a multi-sample result, the sex counts.
- Passing `sex_system`/`sex_chr` explicitly to these functions overrides the stored values.

## 6. HMM Copy Number Estimation — All Samples

```r
multi_esti <- hmm_estimate_CN_multi(
  qploidy_standarize_result = qploidy_std,
  sample_ids                = "all",
  n_cores                   = 2,
  cn_grid                   = c(2, 3, 4, 5, 6)
)

# Mixed-ploidy population: wide grid, restricted per sample around its exp_ploidy
multi_esti_mixed <- hmm_estimate_CN_multi(
  qploidy_standarize_result = qploidy_std,
  sample_ids                = "all",
  n_cores                   = 2,
  cn_grid                   = 1:6,
  max_cn_dev                = 2
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

`types_df$by_cn` groups sample IDs by sample-level CN and category (e.g. `types_df$by_cn$\`4x\`$euploid`), useful for picking representative samples to inspect with `plot_karyotype()`. Pass `add_ghost = FALSE` to drop a fully-missing chromosome placeholder from the plot. When the HMM was run with `sex_system`, `types_df$sex` lists `Sample`, `Sex`, `SexChr`.

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
  min_snps_per_window       = 40,  # increase for more smoothing; or a named vector per chromosome
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
# With sex_system set (stored in the HMM params or passed here): Sex, SexChr,
# SexChrConsistent, ChrType, ExpectedCN
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
select_best_baf_model()             # choose BAF model (one sample); or hmm_estimate_CN_nQuack() for seq-data EM route
hmm_estimate_CN_multi()             # run HMM for all samples (narrow cn_grid, or max_cn_dev; sex_system/sex_chr if needed)
plot_karyotype()                    # visualize per-chromosome CN calls for biological interpretation
plot_cn_track()                     # detailed BAF/z-score/CN diagnostic track for a single sample
compare_cn_track()                  # multi-sample CNV profile grid
filter_hmm_CN()                     # quality filter
count_types()                       # summarize euploid/aneuploid/segmental
compare_cn_track_summary()          # population-level karyotype matrix (samples x chromosomes)
re_standardize()  [optional]        # improve with euploid reference
call_hmm_dosages()                  # call per-locus dosages
export_VCF() / write_hmm_CN()       # export results
```

Among the visualization functions, these four are particularly useful and complementary:
- `plot_karyotype()` gives a genome-wide, per-chromosome view of called copy number for one sample (colorable by `post_max_CN`, `w_baf`, `CN_reliability`, or `dosage`) — the primary way to biologically interpret aneuploidy/segmental results.
- `plot_cn_track()` overlays BAF, z-score, CN calls, and BAF distribution per window for one sample — the primary diagnostic to judge *why* the HMM made a call and whether it should be trusted or refined (see `CN_reliability`).
- `compare_cn_track()` tiles CNV profiles for many samples at once — useful for spotting outliers across the population.
- `compare_cn_track_summary()` condenses the whole population into a samples × chromosomes karyotype matrix, grouped by ploidy category — useful for a high-level population overview (with a `Sex` column when a sex system is set).