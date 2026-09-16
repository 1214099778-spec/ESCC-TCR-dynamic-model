# ESCC-TCR-dynamic-model

Analysis code accompanying the paper:

**Dynamic Peripheral TCR Remodeling Reflects Intratumoral Cytotoxic Immunity and Predicts Response to NICT in ESCC**

This repository contains six R scripts for longitudinal peripheral TCR feature calculation, discovery-only model development, independent validation, cross-center external validation, TCR + MTV integration, and ROC/survival analyses. It contains code and documentation only. Study data, fitted model objects, final fitted coefficients, and generated outputs are excluded.

## Analysis design

- A discovery-defined GLIPH2 cluster dictionary maps exact CDR3b/TRBV/TRBJ keys to clusters. Dynamic features are log2((Post + 1)/(Pre + 1)). Dictionary generation is an upstream prerequisite, outside these scripts.
- Feature filtering, Mann-Whitney ranking, Top10 selection, and L1-penalized logistic model development use the discovery cohort only.
- Fixed-Top10 LOOCV fixes the discovery-selected features before leave-one-out fitting; coefficients and inner-CV regularization are refitted within each training fold.
- Independent validation applies the locked discovery model without refitting.
- Cross-center external validation projects repertoires onto the discovery dictionary and applies the locked model without refitting. The external response endpoint is **TRS0/1 versus TRS2/3** (1 versus 0).
- TCR + MTV integration uses post-treatment metabolic tumor volume (`MTV_post`). Candidate PET comparisons and final logistic models use discovery data; independent validation uses those fitted models. Integration-layer LOOCV keeps the previously calculated TCR score fixed and is not end-to-end nested LOOCV.
- ROC curves use binormal smoothing and bootstrap AUC intervals. Paired DeLong tests use empirical ROC objects. Survival cutpoints are estimated within the analyzed pooled or external data; external survival analysis is exploratory.

## Data availability

Processed bulk TCR repertoire data have been deposited in OMIX under accession **OMIX020653**.

Single-cell RNA/TCR data: **PRJCA028740**.

These accessions do not imply that access has been publicly released. Access conditions must be checked with the data repositories. Single-cell analysis is not implemented by these six scripts. No accession-associated files, raw reads, or participant-level supplementary files are included here.

## Files and execution

Run scripts from the repository root in numerical order. `R/00_utils.R` is sourced by the analysis scripts.

| Script | Purpose |
| --- | --- |
| `R/00_utils.R` | Shared preprocessing, fitting, prediction, and LOOCV helpers |
| `R/01_tcr_dynamic_feature_calculation.R` | Discovery and independent-validation dynamic features |
| `R/02_tcr_score_calculation.R` | Discovery feature selection, TCR model, and fixed-Top10 LOOCV |
| `R/03_independent_and_external_validation.R` | Locked-model independent and external predictions |
| `R/04_tcr_mtv_model.R` | Discovery PET comparisons and TCR + MTV integration |
| `R/05_roc_survival_analyses.R` | ROC, DeLong, thresholds, and exploratory survival analyses |

Required R packages: `data.table`, `dplyr`, `tidyr`, `glmnet`, `pROC`, `survival`, `survminer`, and `ggplot2`. Package versions have not been locked. Syntax checks use R 4.6.0; no data-dependent analysis has been executed for this repository preparation.

See [input specifications](data/README.md), [output handling](results/README.md), and [analysis workflow and review notes](docs/analysis_workflow.md) before execution. The repository is not a standalone executable reproduction without authorized study inputs.

## Review status

This is a private review copy. Only double-bracket syntax corrections were made in scripts 04 and 05. Analytical logic was preserved. Known data-join and dictionary-mapping questions are documented in the workflow and require author review before claiming reproducibility. The code uses the MIT License (2026); the copyright holder and unconfirmed citation details are explicitly marked TODO in `LICENSE` and `CITATION.cff` and require review before public release.
