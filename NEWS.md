# Qploidy2 1.22.0

* Added `nQuack_EM`, `nQuack_mixtures`, and `nQuack_cores` parameters to `hmm_estimate_CN()`. When `nQuack_EM = TRUE`, the three nQuack mixture-model distributions (`quackNormal`, `quackBeta`, `quackBetaBinom`) are fitted to whole-sample heterozygous-site coverage data and the BIC winner is used to set `exp_ploidy` — the expected ploidy that initialises the HMM state distribution and z-scale — instead of relying on the template BAF model. Per-window BAF likelihoods continue to use the existing template/histogram method. `nQuack_mixtures` controls which ploidal-level names are passed to `quackit()` (default: derived from `cn_grid` over 2:6); `nQuack_cores` passes the core count to the nQuack fitting functions. Requires `nQuack` installed via `remotes::install_github("mgaynor1/nQuack")`.

* Added `from_standardized_baf` parameter (default `TRUE`) to `convert2nQuack()`. When `TRUE`, allele counts are derived by inverting the standardization-corrected `baf` column back into reference/alternative counts (`alt = round(baf * R)`, `ref = R - alt`), the same approach already used internally by `.nquack_estimate_ploidy()`, instead of using the raw `X`/`Y` allele counts.

* Added the `Qploidy2_nQuack_tutorial.Rmd` vignette, covering the nQuack EM-based BAF model selection route (`nQuack_EM`/`nQuack_exp_ploidy`) and exporting standardized BAF data to nQuack's input format via `convert2nQuack()`.

# Qploidy2 1.21.0

* Added two-signal consistency gating to `hmm_estimate_CN()`: for each chromosome, the chromosome-level z-score consensus (sum of per-window z log-likelihoods) is compared to the chromosome-level BAF consensus. Gating operates at the chromosome level — not per window — because individual window z-scores are too noisy to reliably indicate disagreement, especially with wide CN grids. When the chromosome mean BAF weight is meaningful (>= 0.1) and the two signals disagree on the best CN state, the BAF contribution is suppressed for all windows on that chromosome and emissions are set to z-score only. When they agree, a reinforcement bonus (`log(3) × mean_w_baf`) is added to the agreed state. Chromosomes with sparse BAF (mean weight < 0.1) are unaffected.

* Added post-EM collapsed-state pruning to `hmm_estimate_CN()`: after EM convergence, adjacent CN states whose means differ by less than 10% of sigma are removed before decoding. This prevents arbitrary HMM state-flipping when the data spread is small relative to the CN grid spacing (e.g. a diploid sample tested with `cn_grid=1:6`), where EM can collapse two CN means to virtually identical values. The `exp_ploidy` state is never dropped regardless of gap. Only truly degenerate pairs (gap < 0.10σ) are removed; genuine low-z_range separations are preserved.

* Added post-decoding z-weakness override to `hmm_estimate_CN()`: after HMM decoding, chromosomes whose called CN differs from the sample-mode CN are overridden back to the sample mode when the chromosome mean z is within `hom_z_sigma_inflate × sigma` of the expected-ploidy z mean. This fires regardless of BAF weight — when the depth signal is too weak to support a non-modal CN call, a marginal BAF signal is also insufficient to sustain it. Prevents small z fluctuations (with or without low-weight BAF) from driving spurious aneuploid calls.

* Redesigned `CN_reliability` in `hmm_estimate_CN()` as a product of four independent evidence components: (1) `post_max` — HMM posterior certainty for the called state; (2) `baf_info` — BAF informativeness computed as 1 minus the normalised Shannon entropy of a 10-bin BAF histogram on [0,1] per window. Measures how structured/peaked the BAF distribution is within a window — a dispersed or uniform BAF (high background noise, impossible to assign to copy-number peaks) yields a low score even when many heterozygous SNPs are present, while BAF clearly clustered at expected heterozygous positions (e.g. 0, 1/3, 1/2, 2/3, 1) yields a high score. BAF values are clipped to [0,1] before binning so `ratio` input columns are handled safely. No-BAF windows use `z_no_baf_scale` as the floor; (3) `zbaf_agree` — per-window z vs BAF agreement computed as 1 minus the total variation distance between the softmax of the per-window z log-likelihood vector and the softmax of the per-window BAF log-likelihood vector. 1 = z and BAF posterior distributions are identical; 0 = completely opposite. Windows with no BAF data default to 1 (no penalty). Degrades gracefully to 1 when gating is skipped (`z_only`, `baf_only`, or `K==1`); (4) `w_baf_norm` — existing BAF-weight scaling (unchanged). `CN_reliability` is now low for calls where BAF is dispersed/uninformative (even when `post_max` is high) and for calls where z and BAF favour different CN states.

# Qploidy2 1.20.0

* Add Claude AI Qploidy2 skills 
* Display number of individuals for print standardization object
* Fixed `hmm_estimate_CN()`: windows with no valid BAF values were assigned a BAF log-likelihood of `0` (i.e. likelihood of 1.0) via `is.nan()` replacement, which could mask non-finite values. These windows are now explicitly excluded from BAF scoring (`w_baf` forced to `0`) instead of relying on incidental cancellation.
* Fixed a bug in `hmm_estimate_CN()` where the `for (chr in unique(win_df$Chr))` loop variable used in the `baf_only` per-chromosome refit and in per-chromosome decoding silently overwrote the `chr` function argument (used for chromosome subsetting), causing the `rerun_overall_ploidy` recursive call to receive the wrong chromosome(s). Loop variables renamed to `chr_i`.
* Fixed `pi0 <- rep((1 - initial_prob) / (K - 1), K)` in `hmm_estimate_CN()` producing `NaN` when only one CN state (`K == 1`) is tested; `pi0` is now set to `1` in that case.
* Fixed the non-monotonic z-mean check in `hmm_estimate_CN()` (`while (any(mu != sort(mu)) ...)`), which used exact floating-point comparison and could spuriously trigger on near-equal `mu` values; it now uses a `1e-6` tolerance.
* Fixed `em_hmm_cn()`: dividing the BAF log-likelihood by `n_baf` could produce `Inf`/`-Inf` for windows with zero valid BAF markers, which the previous `is.nan()` guard did not catch. Now guarded with `ifelse(n_baf > 0, llb / n_baf, 0)` and `!is.finite()`.
* Changed the combined z/BAF emission in `hmm_estimate_CN()` and `em_hmm_cn()` from the heuristic weighted average `(1 - w_baf) * llz + w_baf * llb` to a tempered log-product `llz + w_baf * llb`, so the z-emission always contributes fully and BAF evidence is scaled by `w_baf` rather than both terms being down-weighted together. This changes some CN posterior probabilities relative to previous versions.
* Documented in `hmm_estimate_CN()` that a single `sigma` is shared across all CN states (not state-specific), and that the combined emission is a tempered log-product rather than a statistically derived mixture.
* Improved performance of `hmm_estimate_CN()`: window-level summary statistics are now computed with `dplyr::group_by()`/`summarize()` instead of `by()` + `do.call(rbind, ...)`, and per-window BAF values are extracted via a single `split()` keyed lookup instead of rescanning the full marker table for every window (O(N) instead of O(W x N)). Also simplified the `n_baf` computation.


# Qploidy2 1.19.0

* Added `min_samples_by_cluster` to `standardize()`, `re_standardize()`, and `get_centers()`. Default `1` keeps all clusters. When set higher, clusters with fewer non-NA samples are dropped before center estimation, preventing weak clusters from anchoring BAF standardization.
* Added `split_mk_type()` for segmental allopolyploid species analysis to separate sub-genome-specific markers from markers not specific to a sub-genome using reference diploid samples.
* Improved segmental aneuploidy detection in `hmm_estimate_CN()`: the per-chromosome forward pass in `fb_smooth()` is now initialized from a `pi0` derived from the summed BAF log-likelihoods across windows on the chromosome, rather than a uniform prior. This prevents a high z-score in the first window from anchoring the whole chromosome in the wrong CN state.
* Fixed the raw-data input handling in `hmm_estimate_CN()`: when `qploidy_standarize_result = NULL`, the BAF/z signal columns are resolved from `use_values` using the correct `ratio`/`R` and `BAF`/`z` mappings, and the function now rejects invalid or missing column names with clearer checks.
* Improved chromosome filtering and NA handling in `hmm_estimate_CN()`: explicit chromosome subsets now ignore missing/NA chromosome labels, invalid chromosome requests stop cleanly, and rows with `NA` chromosome values are dropped before analysis so uninformative markers cannot leak into downstream segmentation or HMM fitting.
* Fixed asymmetric copy-number gain/loss detection in `define_z_limits()`. The initial directional steps now calibrate the mean for CN states below and above `exp_ploidy` separately, avoiding misplacement when gains and losses are unequal.
* Added `cn0_background_frac` to `hmm_estimate_CN()` (default `0`). When `0` is included in `cn_grid`, this parameter shifts the CN=0 emission mean upward by the expected background fraction so array probes with residual signal at CN=0 are modeled correctly. `CN = 0` is also now allowed in `cn_grid`.
* Fixed a bug in `re_standardize()`: markers in windows called as CN=0 were not being masked when `CN_call` was `NA`, because `NA != ploidy.standardization` evaluates to `NA` instead of `TRUE`. Both code paths now explicitly exclude `is.na(CN_call)` before cluster standardization.
* Update `transition_jump` default value to 0.9 (before was 0.995)

# Qploidy2 1.18.3

* Reverted the unconditional mu update in `em_hmm_cn` back to the conditional guard, as the unconditional change introduced downstream bugs.
* Changed default `min_het_frac` from `0.05` to `0.01` and default `het_range` from `c(0.2, 0.8)` to `c(0.05, 0.95)` in `compute_baf_likelihoods`, `select_best_baf_model`, `hmm_estimate_CN`, `call_hmm_dosages`, `select_best_raw_model`, and `re_standardize`. The wider interval and lower threshold are more appropriate for high-ploidy species where inner heterozygous peaks sit close to the boundaries, and produce fewer spurious CN=1 exclusions.
* Allelic ratio plots now use the symbol `θ` (rendered via `expression(theta)`) instead of "Ratio" as the y-axis / x-axis label in `plot_raw` and `plot_standardization`, with facet strip labels parsed via `label_parsed`.
* Bugfix in `re_standardize`: when `ploidy.standardization` is provided but differs from the mode CN call in `hmm_CN_multi`, a warning is issued and `ploidy.standardization` is overridden to the mode CN to avoid unexpected standardization results.
* Fixed R CMD check NOTEs: added `CN_mode` to `globalVariables` in `export_results.R` and `CN_call_filled` to `globalVariables` in `karyotype_notation.R`.


# Qploidy2 1.18.2

* Revert mu stability guard introduced in 1.18.0: `em_hmm_cn` now always updates mu unconditionally (standard EM update with `pmax` denominator guard). The while-loop in `hmm_estimate_CN` already handles non-monotonic mu.
* Add `hom_z_sigma_inflate` parameter to `hmm_estimate_CN` (default `1.5`). Applies two corrections for all-homozygous windows (`w_baf == 0`): (1) sigma is inflated by this factor during decoding to reduce z-emission discriminability; (2) post-decoding, chromosomes where every window has `w_baf == 0` are overridden to the sample-mode CN if their mean z-score is within `hom_z_sigma_inflate * sigma` of the baseline — suppressing false positive CNV calls driven by reference-homozygous depth bias while preserving genuine CN changes with larger z-deviations.
* CN=1 samples are now labelled "1x or inbred" instead of "1x" in `compare_cn_track_summary` right-side block annotations, `count_types` printed summary, and `print.count_types` per-type rows.
* Bugfix on `re_standardize` - when `use_estimated_dosages = TRUE` and majority of CN estimated doesn't match `ploidy.standardization` defined, replace the value by the majority estimated and warn the user with message

# Qploidy2 1.18.1

* Bugfix on `fb_smooth` when there is only one window at the chromosome
* Adapt `plot_karyotype` to objects coming from `filter_hmm_cn`. Now, a gray bar is plotted on filtered regions
* Add color option by `CN_reliability` on `plot_karyotype`
* Add `CN_reliability` to results by window 
* Add generalized Rmd document to reproduce alfalfa tutorial with any `qploidy_standardized` initial file
* Update alfalfa tutorial adding new functions

# Qploidy2 1.18.0

## HMM stability improvements

### Issue 1 – First window of the first chromosome had higher CN variation

**Changes that remained in the EM (`em_hmm_cn.R`)**:

* `update_pi0 = FALSE` (default)
* Mu stability guard

**Changes in the decoder (`hmm_main.R`)**:

* Per-chromosome forward-backward decoding(`fb_smooth`)
* Prior transition matrix for decoding (`A_dec`)

### Issue 2 – Higher variation in all-homozygous chromosomes

* New output metric (`CN_reliability`): `hmm_estimate_CN` now returns a `CN_reliability` column in `by_window`
* hmm results plots updated with this new metric

### Issue 3 – CN estimation results differed based on chromosome ordering / naming

**Changes in `hmm_main.R`**:

* Sort windows by chromosome z-deviation before the EM
* Recompute `ll_em` in original window order after the EM

## New functions

* `compare_cn_track_summary` — karyotype-style matrix (samples × chromosomes),
  ordered by sample-level CN and ploidy category (euploid / aneuploid /
  segmental), with colour-coded tiles, separating lines, and right-side labels.
  Accepts `filter_type` and `filter_cn` arguments.
* `filter_hmm_CN` — masks CN calls (sets to `NA`) in windows that fail one or
  more quality thresholds: `min_CN_reliability`, `min_post_max`, `min_w_baf`,
  `min_n_snps`, `min_n_het`, `min_window_size`, `max_CN_call`.
* `count_types` — classifies every sample as euploid, aneuploid, or segmental
  aneuploidy; returns a `count_types` object (a named list of sample-ID vectors
  plus a `$by_cn` breakdown by sample-level CN) and auto-prints a summary table
  via `print.count_types`.  

## New tutorial template
 * QuickStart_template.Rmd was created to allow reprodutibility of the alfalfa tutorial
 with any initial standardized file input.


# Qploidy2 1.17.0

* Change `plot_cn_track` labels from `Est overall ploidy` to `Mean total depth/Z CN` and y-axis label `Copy Number` to `HMM-based CN`
* Add argument `threshold.missing.samples` to `standardize` and `re_standardize`. It filters samples with fraction of markers containig missing information higher than defined by the argument. Default set to `1` (filter not applied).
* Function `export_mappoly` to produce input file for `filter_aneuploidy` function in `MAPpoly`
* Alfalfa tutorial updates

# Qploidy2 1.16.2

* Add different VCF required fields for when `qploidy_read_vcf` is `geno = TRUE` or `FALSE`. When it is `FALSE`, `GT` is not required

# Qploidy2 1.16.1

* Fix README link to new tutorial
* Adapt Alfalfa tutorial 
* Add cn_grid limits 1-10

# Qploidy2 1.16.0

* Modify `call_hmm_dosages` and `call_BAF_dosages` to output the BAF weights and likelihood of each CN and dosage called
* Adapt `export_VCF` to export likelihoods information
* Add functions `plot_karyotype` and `karyotype_notation`
* Update alfalfa tutorial

# Qploidy2 1.15.0

* Create arguments `rerun_overall_ploidy` and `recycled_obj_rerun_overall_ploidy` in `hmm_estimate_CN`. When `rerun_overall_ploidy` is set to TRUE it re-run the model selection and HMM removing markers estimated having CN difference than the mode. This improves overall ploidy estimation and corrects the variation of total depth. `recycled_obj_rerun_overall_ploidy` is of exclusve use of internal process
* Fix and update `hmm_estimate_CN` documentation
* Add more tests to `hmm_estimage_CN` function
* Reduce RAM consumption of parallel dosage call
* Create convert2nQuack function based on Michelle Gaynor tutorial (added as author)
* Modify defaults after testing in different scenarios (selecting values that work best in majority): 
  - remove beta distribution default testing
  - add higher variances to be tested in the grid
  - `rerun_overall_ploidy` is set to TRUE
* `re_standardize` new `use_estimated_dosages` parameter — added with default FALSE (uses original dosages for non-circular re-standardization); replaces the previous genos argument approach
* hmm_estimate_CN bug fixes:
    - z_only = TRUE caused a crash (n_baf was never defined)
    - pi0 didn't sum to 1 — best CN state was hardcoded to 0.85 instead of using initial_prob
    - Loop counter idx was overwritten by inner for (idx in keep_lower/higher) loops, breaking the non-monotonic correction guard
* Add `plot_xy_with_ploidy_guides` to tutorial and add arguments to be able to use standardization object as input
* Make available Prepare Inputs tutorial
* Update Alfalfa tutorial
* **warning**: this version change functions default values. Therefore, results may differ from previous versions  

# Qploidy2 1.14.0

* Parallelize dosage call in `call_hmm_dosages`
* Change error by warning in read_hmm_CN when one or more samples in multi sample hmm object failed to be estimated and lack parameters information
* Fix conflict of out_filename while merging qploidy_standardization objects in `merge_qploidy_data`
* Function `rename_samples` to rename sample names in qploidy_standardization and hmm_CN objects

# Qploidy2 1.13.0

* Improve dosage call. Now mixed model parameters are defined by sample and copy number. Parameters carried by the hmm object is passed to `call_hmm_dosage`
* The change above remove the selected_model requirement from `re-standardize` function once the models parameters are being used from the hmm object o, if absent (rare), estimated by sample again via grid and select_baf_model
* Add tests for `re_standardize`
* Updates vignette

# Qploidy2 1.12.0

* Report on print.standardization object the number of markers containing z-score metric
* add function merge_qploidy_datas

# Qploidy2 1.11.0

* Add count breakpoints functions
* Fix notes in CMD checks

# Qploidy2 1.10.0

* Changing default value of `add_uniform_grid` in `hmm_estimate_CN` to FALSE for better performance
* Add option to filter min and max total depth (R) on `standardize` function

# Qploidy2 1.9.2

* Specific error messages for standardization when required inputs are not provided
* threshold.n.clusters made as optional. Default is ploidy.standardization + 1

# Qploidy2 1.9.1

* Add `GenoBrew` links 

# Qploidy2 1.9.0

* Fork from Cristianetaniguti/Qploidy@development
* Re-branding with Qploidy2
* Adding updated vignettes 
* Allow `hmm_estimate_CN_multi` to run with raw data (not standardized)

# Qploidy 1.8.7

* Adding input checks for read_hmm_CN function
* bugfix update_hmm_CN_multi

# Qploidy 1.8.6

* Order columns of hmm_CN$by_window object when results from different samples are merged
* Add update_hmm_CN to drop a specific sample results from a multi sample hmm_CN and/or replace by the results in a single sample hmm_CN
* R CMD checks okay

# Qploidy 1.8.5

* Remove defunct code for Shiny interface - migrating the Breeding-Insight/QploidyApp
* remove NonASCII characters
* Better address results when number of windows is 1
* Removing old vignette

# Qploidy 1.8.4

* Bugfix in `pca_plot()`: markers where `sd(na.rm = TRUE)` returns `NA` (e.g. only one non-NA sample) were silently retained as all-NA columns due to R's `NA > 0 = NA` subsetting behaviour; these are now correctly excluded before PCA.
* Bugfix in `pca_plot()`: `NaN`/`Inf` values introduced when imputing columns that are entirely `NA` are now replaced and removed, preventing the `prcomp()` error "cannot rescale a constant/zero column to unit variance".


# Qploidy 1.8.3

* `hmm_estimate_CN()` gains a `use_values` argument supporting all combinations of BAF/ratio and z/R as input signals.
* `plot_cn_track()` now accepts raw `data` + `geno.pos` as an alternative to a `qploidy_standardization` object, enabling visualization of non-standardized ratio/R data directly.
* Y-axis labels in `plot_cn_track()` now adapt to the input type: "BAF" or "Ratio" for the BAF panel, and "z" or "R" for the z-score panel.
* BAF density polygon in `plot_cn_track()` is now normalized per chromosome, ensuring consistent visibility across chromosomes with varying marker densities.
* `plot_cn_track()` with `summarized = TRUE` now renders the BAF density polygon in black for improved contrast, and correctly scales the polygon per chromosome.
* The distribution summary panel in `plot_cn_track()` now uses the correct x-axis limits and label based on the input data type (BAF or Ratio).
* `compare_cn_track()` gains an `add_het` argument (with `hmm_dosage_calls`) to display a per-sample heterozygosity sidebar alongside the CN track.
* `compare_cn_track()` gains an `interactive` argument to produce a plotly figure with per-segment tooltips.
* `compare_cn_track()` gains a `gray_CN` argument to manually specify the baseline copy-number value for color scaling.
* `compare_cn_track()` gains a `facet_nrow` argument for controlling facet layout.
* New function `plot_heterozygosity()` for visualizing per-sample heterozygosity as a heatmap grid, with optional plotly interactivity.
* `pca_plot()` gains `samples` and `palette` arguments for subsetting samples and customizing colors.

# Qploidy 1.8.2

* BUGfix on read vcf functions when there are duplicated markers and lack of marker ID
* Adjustments in the plot_standardization

# Qploidy 1.8.1

* Added cross-parameter validation in `standardize()`: `threshold.n.clusters` is now checked to not exceed `ploidy.standardization + 1`, with an informative error message.
* `select_best_baf_model()` now supports dual-input mode: BAF values can be supplied directly via `baf_vec` (with an optional `chr_vec`) or extracted from a `qploidy_standardization` object and a sample name. Both paths are mutually exclusive and validated.
* Added diagnostic warning in `hmm_estimate_CN()` when all heterozygous markers in one or more windows are discarded by the `dosage_threshold` filter, reporting the affected windows by chromosome and window ID.
* Added `parallel_type` argument to `hmm_estimate_CN_multi()` (default `"auto"`): automatically selects `"FORK"` on Unix/macOS (faster, no symbol re-export needed) and `"PSOCK"` on Windows. Can be set explicitly to `"FORK"` or `"PSOCK"`.
* Warnings emitted inside parallel workers are now captured and re-issued on the main R session after `parLapply`, preventing diagnostic messages from being silently lost in PSOCK clusters.
* `select_best_baf_model()` now emits a warning when the top three BIC-ranked models disagree on the best CN estimate, alerting users to ambiguous model selection.
* `hmm_estimate_CN()` now excludes CN = 1 as a candidate state for windows that contain BAF values within `het_range` (default `c(0.2, 0.8)`), preventing spurious CN = 1 calls driven by a low heterozygous-to-homozygous ratio.

# Qploidy 1.8.0

* Enhanced console messages for standardization and HMM steps, improving clarity and user feedback.
* Improved CN grid selection: unlike CN values are now removed after the EM loop by checking z-score means. The CN grid is reordered from lower to higher ploidy, and EM is rerun if disruptive CN values are detected.
* The min_snp_per_window argument now defaults to 10% of the smallest chromosome size, with a minimum of 5 SNPs, for more adaptive windowing.
* Added an argument to define z-score intervals by subtracting z_range (instead of adding), resulting in smaller intervals.
* Refactored EM loop into a dedicated function for better modularity and maintainability.
* Updated dosage call plot colors for improved visual distinction and interpretability.
* Exception created in HMM when expected ploidy = 1. In this case, heterozygous are not expected so the BAF weight is set to 0.5 regardless if heterozygous are present or not. Parameter can still be controled by baf_weight argument (default 1).

# Qploidy 1.7.1

* New function depth_pca_plot to check for batch effects using a pca based on total depth
* Allow input data for function chrom_ttest, window_ttest and depth_pca_plot to have only MarkerName, SampleName, and R columns if geno.pos data.frame is provided with Chromosome and Position information.
* Bugfix on DESCRIPTION authors 
* Improve documentation for new functions

# Qploidy 1.7.0

* New functions chrom_ttest and window_ttest
* Adding to plot_geno_by_marker:
  - transparency argument alpha
  - identify decreasing ratios by dosage for expected dot plot

# Qploidy 1.6.9

* Add filter_R argument on standardize. If FALSE (default) filters defined are not applied for z score calculation
* Add depth_zero_as_x to plot_cn_track. Total counts (R) of missing genotypes are counted as 0 and converted to z score. depth_zero_as_x = TRUE mark on the graphic which z values refer to the 0 total counts

# Qploidy 1.6.8

* Increase speed to export VCF
* Add re_standardize function to run standardization based on HMM results
* Reduce text on plot_cn_track

# Qploidy 1.6.7

* Make initial probability 0.95 for the expected ploidy
* Bugfix plot_cn_track

# Qploidy 1.6.6

* Add correction factor for Z and BAF likelihood to have same weight if BAF_weight = 0.5 (new default)
* Filter low probability heterozygous for counts on determining the BAF_weight
* If correct_scale = TRUE (default) BAF likelihood is corrected by the number of markers (with BAF values) used for the distribution
* plot_cn_track now contains the sample-level BAF histogram and parameters descriptions written on the figure

# Qploidy 1.6.5

* Number of heterozygous to define the BAF weight now is counted using dosage call based only on BAF
* Dosage call functions added based only on BAF and based on BAF according to HMM CN call result
* Export called dosages in a VCF format
* Add co-pilot-instructions to the repository

# Qploidy 1.6.2

* Avoid markers without chromosome information in the standardized dataset
* Add plot to compare ratios and standardazed ratios (BAF) in a geno vs value format
* Add new type of plot for plot_qploidy_standardized for raw total depth (R)
* Add checks for standardize input files. It now requires specific column names and ordering.


# Qploidy 1.6.1

* write and read function for hmm_CN object
* print function summarizing estimated paramters in hmm_CN object
* hmm now return two data.frames, one with information by window and another by marker
* hmm returns list of parameters defined by user and estimated by function for reproducibility

# Qploidy 1.6.0 

* Major refactor and improvements to BAF likelihood/model selection workflow and plotting 
* Added support for testing different emission distributions for BAF (nQuack idea), including negative binomial, with automated model selection via BIC
* Refactoring of BAF template generation and likelihood computation into dedicated functions
* Improved plotting functions:
    - plot_cn_track has now a fallback to hmm_CN$updated_data if qploidy_standarize_result is NULL
    - function compare_cn_track added
    - Improved color palette for CN plots
    - Chromosome sorting and x-axis labeling improved for clarity
    - Legends and color mapping clarified and improved
    - Subsetting chromosomes on plot - not only on hmm
* HMM segmentation and CN calling:
    - Implementation of changepoint-based z-score segmentation as an alternative to fixed-window approaches
    - Single HMM across all chromosomes (not per-chromosome)
    - Robust outlier handling for z-scores (rm_outlier supports z argument and outlier column)
    - Exception handling for single-window case (assigns CN by BAF likelihood only)
    - exp_ploidy argument added and integrated throughout (including Shiny UI)
    - Improved verbose output and summary messages
    - By default, z_range is now automatically calculated based on the inter-quantile range of z-scores and the number of CN states tested, for more robust and adaptive HMM initialization
* User experience and error handling:
    - Improved error messages and fallback logic in plotting and HMM functions
    - Documentation expanded and clarified for all major functions
* Shiny UI:
    - Sample-level BAF distributions plot added 
    - exp_ploidy is now an advanced option for z_only=TRUE
    - Argument passing and advanced options printing improved

# Qploidy 1.5.2

* Z-score calculation is now only applied to markers not filtered by genotype probabilities and missing data
* Markers filtered out during standardization are not considered to define the window boundaries for HMM CNV estimation 
* Last window of each chromosome are never smaller than the defined window size
* CN HMM plots follow same scale for easier comparison between metrics by window

# Qploidy 1.5.1

* Fix HMM CNV estimation when some of the windows don't have heterozygous genotypes (use only z-score for them)
* Improve Shiny interface for HMM CNV estimation

# Qploidy 1.5.0

* Add HMM to estimate CN (beta version)
* Shiny interface improvements
* Add CSV or TSV as possible input files

# Qploidy 1.1.0

* Beta version of Qploidy + BIGapp interface
* Fix Bioconductor dependency call
* Fix %>% call

# Qploidy 1.0.0

* First release
* Functions improvements including comprehensive documentation
* Vignette with functional simulated example
* Testthat tests included
* CRAN submission

# Qploidy 0.0.0.9000

* This is a beta version
* Initial addition of functions
* First version of the vignette



