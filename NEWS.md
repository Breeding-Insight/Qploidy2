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



