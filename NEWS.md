# Qploidy2 1.17.0

## HMM stability improvements

### Issue 1 – First window of the first chromosome had higher CN variation

The root cause was a **circular feedback loop in the EM M-step**: at each iteration
the initial-state distribution was updated as `pi0 ← gamma[1,]`.  When the very
first window was slightly biased toward a wrong CN state (e.g., because the initial
z-score happened to be closer to the wrong mu), the M-step would progressively
increase `pi0` for that wrong state.  By convergence, `mu` for the wrong state had
collapsed to `z[1]` (driven by the `pmax(sum(w), 1e-12)` denominator guard when
gamma was near-zero), making the emission maximally favorable for the wrong call
at window 1.  Any subsequent decoder — including bidirectional Viterbi — would
then also pick the wrong state because `ll_em[window1, wrong_state]` was at its
peak.

**Changes that remained in the EM (`em_hmm_cn.R`)**:

* `update_pi0 = FALSE` (default): the initial-state distribution `pi0` is **no
  longer updated** in the M-step.  `pi0` is fixed at the value set from
  `initial_prob` and `exp_ploidy` before the EM starts, breaking the feedback
  loop that caused window 1 to lock into the wrong state.
* **Mu stability guard**: `mu[k]` is only updated when the total posterior weight
  for state `k` exceeds `0.5`.  States that are genuinely absent from the data
  (total weight ≈ 0 due to the `pmax` denominator guard) keep their physically
  motivated initial estimate from `define_z_limits`, preventing mu from
  collapsing to `mean(z)` and making emissions degenerate.

**Changes in the decoder (`hmm_main.R`)**:

* **Per-chromosome forward-backward decoding** (`fb_smooth`): instead of running
  a single HMM chain across all chromosomes, `fb_smooth` is now called
  independently for each chromosome with a fresh uniform initial distribution.
  This makes the result independent of chromosome ordering — no chromosome is
  disadvantaged by its position in the chain — and gives window 1 of every
  chromosome proper backward context from the rest of its own chromosome.
* **Prior transition matrix for decoding** (`A_dec`): the decoder uses the
  _prior_ transition matrix (strong diagonal = `transition_jump`, rare
  off-diagonal ≈ `1e-3/(K-1)`) rather than the EM-converged `A`.  After EM
  convergence on a sample dominated by one CN state, states never visited receive
  uniform rows in `A` (due to the same `pmax` guard), giving them a "free jump"
  advantage of only ~1.4 log-units.  The prior matrix imposes a ~8 log-unit
  penalty for any off-diagonal transition, making it virtually impossible for a
  degenerate state to win at window 1 through the backward pass alone.

### Issue 2 – Higher variation in all-homozygous chromosomes

When an entire chromosome has no heterozygous markers (`w_baf = 0` for all
windows), the HMM emission reduces to z-score only.  Before the fixes in
Issue 1, two compounding problems caused wrong CN calls on these chromosomes:

1. The mu-collapse bug (Issue 1) could pull `mu[wrong_state]` toward the
   chromosome's mean z-score, making the emission flat and removing any
   discriminating power.
2. Degenerate rows in the EM-learned `A` (uniform for never-visited states)
   gave those states a ~1.4 log-unit "free jump" advantage in the backward
   pass, easily overcoming the weak emission signal.

The fixes from Issue 1 — the mu stability guard, `update_pi0 = FALSE`, the
prior `A_dec` (8 log-unit penalty per off-diagonal transition), and
per-chromosome decoding — together eliminate both compounding problems and are
the primary solution for all-homozygous chromosomes as well.

**New output metric (`CN_reliability`)**:

* `hmm_estimate_CN` now returns a `CN_reliability` column in `by_window`:

  ```
  CN_reliability = post_max × max(z_no_baf_scale, min(1, w_baf / baf_weight))
  ```

  This metric does **not** affect the HMM emission or the CN call itself.  It
  is a post-hoc quality indicator that is low whenever the HMM posterior is
  low _or_ BAF evidence is absent.  The `z_no_baf_scale` parameter (default
  `0.25`) sets the floor: a window with no BAF support can reach at most
  `0.25 × post_max` in reliability, signalling to users that the call is
  z-score only and should be treated with more caution.  Old `hmm_CN` objects
  without this column automatically fall back to `post_max` in both plot
  functions.

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

## Other changes

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



