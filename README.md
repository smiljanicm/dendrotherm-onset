# dendrotherm-onset

An R analysis pipeline for estimating temperature thresholds associated with the onset of dendrometer-derived stem growth in 14 temperate tree species. It processes dendrometer and temperature records, fits binomial generalized additive models (GAMs), estimates threshold distributions with a hierarchical bootstrap on SLURM, and produces sensitivity analyses and manuscript tables and figures.

This README documents the current code on `main`. The manuscript title, authors, data provenance, and final methodological wording should be checked against the paper before publication.

## What the analysis estimates

The response `occ` indicates whether a day's derived maximum growth rate (`GRO_max_rate`) is greater than zero. Models relate this binary response to temperature during the preceding 24-hour growth window. The three temperature summaries used in the main model family are the mean (`Temperature`), minimum (`TempOverallMin`), and maximum (`TempMax`). Additional models include day of year (DOY) or growing degree days (GDD) calculated with base temperatures 0–6 °C.

For each fitted probability curve, the code reports the lowest temperature on its prediction grid at which the probability reaches 5%, 10%, 25%, or 50%. `onset_p*` uses an absolute probability target (for example, `onset_p50` means probability 0.5); `onset_rel*` uses the same fraction of the maximum fitted probability on that grid. A threshold can be missing if a target is not attained. These are operational thresholds derived from fitted curves, rather than direct observations of cambial activity.

## Repository layout

| Path | Role |
| --- | --- |
| `R/data/` | Retrieve, clean, filter, and create full and restricted analysis datasets. |
| `data/` | Intermediate and analysis CSV files, including committed datasets. |
| `hpc/R/` | GAM specifications, hierarchical bootstrap, threshold extraction, and consolidation. |
| `hpc/scripts/` | SLURM submission and compilation scripts. |
| `R/posthpc/` | Aggregate bootstrap thresholds and model diagnostics into tables. |
| `R/analysis/` | Freeze–thaw summary, descriptive statistics, and leave-one-year-out Brier evaluation. |
| `R/summary/` | Meta-analysis, manuscript numbers, and figure generation. |
| `output/tables/` | Aggregated CSV results; `german/` contains semicolon-delimited variants. |
| `paper/` | Generated manuscript figures and tables. |

All relative paths in the scripts assume that commands are run from the repository root.

## Requirements

- R with the packages used by the selected stage. The scripts use `tidyverse` (including `dplyr`, `purrr`, and `readr`), `lubridate`, `hms`, `mgcv`, `qgam`, `duckdb`, `DBI`, `metafor`, and `zip`. Raw-data retrieval additionally uses `dendrogreif` and `treenetproc`.
- A SLURM cluster with `sbatch` for the supplied bootstrap launcher. Its scripts assume a `micromamba` environment named `r-env`, a `snowball` partition, and access to the required R packages. Adjust those local settings for your cluster.
- Credentials and access to the source dendrometer database **only if regenerating** `data/initial_data.csv` from raw records. `R/data/00_get_data.R` currently calls `dendrogreif::con_init("../../../creds_local")`; set this path for your environment and keep credentials outside version control.

There is no package lockfile or single installation script in the repository. Package versions and the database snapshot used for the manuscript should be recorded separately for an exact rerun.

## Reproduce the pipeline

The repository includes the derived CSV inputs under `data/`, so readers with no database access can start from the appropriate committed dataset. To regenerate from the source database, first configure the connection in `R/data/00_get_data.R` and run that stage. It writes `data/initial_data.csv` and currently requests data through 2026-12-31; reproducing a historical manuscript dataset requires the corresponding source snapshot and date cutoff.

### 1. Prepare analysis data

```bash
Rscript R/data/01_clean_data.R
Rscript R/data/02_growth_season_filter.R
Rscript R/data/03_restricted_data.R
```

The first script filters tree-years with more than 60 observed days, a growth-day proportion above 0.1, and a first observation before DOY 185; it also excludes 2013. The second calculates GDD from DOY 31 and writes three spring windows ending at DOY 172: `gs_data.csv` (DOY 60–172), `wide_data.csv` (DOY 32–172), and `narrow_data.csv` (DOY 91–172). The third makes `_restricted.csv` versions for years after 2022.

If using the committed `data/gs_data.csv`, `data/wide_data.csv`, and `data/narrow_data.csv`, these preparation commands are unnecessary for running the models.

### 2. Run the bootstrap on SLURM

From the repository root, create the working directories and submit the arrays:

```bash
mkdir -p hpc/logs hpc/temp hpc/output
bash hpc/scripts/submit_all.sh
```

The launcher submits one array per species for each of the `full` and `restricted` datasets. Each array has 100 tasks; each task processes 20 replicate indices for each of the three spring windows, yielding 2,000 bootstrap replicates per species, window, and period when all tasks finish successfully. The bootstrap resamples observed years and then tree-years within each selected year. The dependent compile job gathers segment and diagnostic CSV files into `hpc/output/seg_{full,restricted}.csv` and `hpc/output/diag_{full,restricted}.csv`; processed chunks move to `hpc/temp/compiled_archive/`.

The model set includes temperature-only GAMs with `k = 10` or `20` and `k = 10` temperature models with DOY or one of seven GDD covariates. Fitting uses a binomial logit GAM with REML. The scripts record AIC, BIC, `k.check`, and concurvity diagnostics. `hpc/output/` is ignored by Git, so a checkout may contain no compiled bootstrap results.

### 3. Summarize bootstrap outputs

Run the postprocessing scripts in numerical order for each period:

```bash
for period in full restricted; do
  Rscript R/posthpc/01-threshold_table.R "$period"
  Rscript R/posthpc/02-delta_table.R "$period"
  Rscript R/posthpc/03-extract_aic.R "$period"
  Rscript R/posthpc/04-extract_concurvity.R "$period"
  Rscript R/posthpc/05-kcheck.R "$period"
done
```

The threshold table groups bootstrap estimates by species, model, and spring window. Its `*_med`, `*_low`, and `*_high` columns use the median, 5th percentile, and 95th percentile of available estimates, respectively; these bounds form a 90% central percentile interval. A separate `n` column counts rows in the group, which should not be assumed to equal the number of nonmissing estimates for every threshold.

### 4. Optional analyses and paper outputs

```bash
Rscript R/analysis/01_freeze_thaw.R
Rscript R/analysis/02_summary_stats.R
Rscript R/analysis/03-loyo_brier.R
Rscript R/summary/03_ms_numbers.R
Rscript R/summary/04_loyo_meta_analysis.R
```

The leave-one-year-out (LOYO) script evaluates species with more than six observed years on the wide window and writes full and summary Brier-score tables. The summary scripts consume those and the postprocessed tables. `R/summary/02_anatomy_meta_analysis.R` defines a callable `run_anatomy_meta_analysis()` function; it needs `R/summary/01_meta.R` and explicit output filenames for its four result types. The remaining scripts in `R/summary/` create manuscript figures and related summaries; inspect their input and output paths before rerunning individual figures.

## Main outputs and interpretation

| Output | Contents |
| --- | --- |
| `output/tables/thresholds_full.csv`, `thresholds_restricted.csv` | Bootstrap threshold summaries for each species, temperature metric, model, and window. |
| `output/tables/delta_table_*.csv` | Differences from temperature-only models under DOY and GDD adjustments. |
| `output/tables/aic_comparison_all_sweep_*.csv` | AIC comparisons across the GDD base-temperature sweep. |
| `output/tables/concurvity_all_sweep_*.csv`, `kcheck_summary_*.csv` | GAM diagnostics. |
| `output/tables/loyo_brier_full.csv`, `loyo_brier_summary.csv` | Held-out-year predictive Brier scores. |
| `paper/figures/`, `paper/tables/` | Selected manuscript-ready figures and tables. |

The `full` and `restricted` suffixes refer to the available period and the post-2022 subset. Window labels in output tables use `60-172`, `32-172`, and `91-172`. Some generated tables are committed; rerunning scripts may overwrite them. Preserve a copy of any publication snapshot before regeneration.

## Citation and license

The repository is licensed under GPL-3.0; see `LICENSE`. Software citation metadata is in `CITATION.cff`. Its manuscript reference and some author metadata are placeholders pending final publication. Please use the final paper citation once available.

## Reproducibility notes for maintainers

- Record the source-data snapshot, actual end date, R/package versions, and cluster environment associated with a published release.
- Check the manuscript's stated interval coverage against the code's 5th–95th percentile aggregation, and clarify how boundary and missing thresholds were treated.
- Confirm the exact manuscript title, species/site description, author list, DOI, and final table/figure numbering before releasing the README and completing `CITATION.cff`.
