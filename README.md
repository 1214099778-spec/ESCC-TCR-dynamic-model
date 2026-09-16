# ESCC-TCR-dynamic-model

Analysis code accompanying the paper:

**Dynamic Peripheral TCR Remodeling Reflects Intratumoral Cytotoxic Immunity and Predicts Response to NICT in ESCC**

This repository contains six R scripts for longitudinal peripheral TCR feature calculation, discovery-cohort model development, independent validation, cross-center external validation, TCR + MTV integration, and ROC/survival analyses.

## Analysis design

- A discovery-defined GLIPH2 cluster dictionary maps exact CDR3b/TRBV/TRBJ keys to TCR clusters. Dynamic features are calculated as log2((Post + 1)/(Pre + 1).
- Feature filtering, Mann–Whitney ranking, Top10 selection, and L1-penalized logistic model development are performed using the discovery cohort only.
- Fixed-Top10 LOOCV evaluates the discovery-selected feature set by leave-one-out refitting within the discovery cohort and does not represent fully nested feature-selection cross-validation.
- Independent validation applies the locked discovery-derived model without feature reselection, coefficient refitting, or recalibration.
- Cross-center external validation applies the locked discovery-derived model to repertoires mapped to the discovery-defined TCR cluster dictionary. The external pathological-response endpoint is TRS0/1 versus TRS2/3.
- TCR + MTV integration uses post-NICT metabolic tumor volume (`MTV_post`). Candidate imaging variables are evaluated in the discovery cohort, and the fitted TCR + MTV model is subsequently assessed in the independent validation cohort.
- ROC curves are smoothed using the binormal method, with bootstrap confidence intervals for AUC estimation. Paired model comparisons are performed using DeLong tests on empirical ROC objects.
- Survival analyses use Kaplan–Meier estimates and data-driven cutpoints. External-cohort survival analyses are exploratory.

## Data availability

Processed bulk TCR repertoire data have been deposited in OMIX under accession **OMIX020653**.

Single-cell RNA/TCR sequencing data are available under BioProject **PRJCA028740**.

Access to these datasets is subject to the policies and access conditions of the corresponding repositories. No raw sequencing data or participant-level study data are included in this repository.

## Files and execution

Run scripts from the repository root in numerical order. `R/00_utils.R` contains shared functions used by the analysis scripts.

| Script | Purpose |
| --- | --- |
| `R/00_utils.R` | Shared preprocessing, fitting, prediction, and LOOCV functions |
| `R/01_tcr_dynamic_feature_calculation.R` | Longitudinal TCR dynamic-feature calculation |
| `R/02_tcr_score_calculation.R` | Discovery-cohort feature selection, TCR model development, and fixed-Top10 LOOCV |
| `R/03_independent_and_external_validation.R` | Locked-model independent and external validation |
| `R/04_tcr_mtv_model.R` | PET/CT candidate-variable assessment and TCR + MTV integration |
| `R/05_roc_survival_analyses.R` | ROC, DeLong, threshold, and survival analyses |

Required R packages include `data.table`, `dplyr`, `tidyr`, `glmnet`, `pROC`, `survival`, `survminer`, and `ggplot2`. The code was syntax-checked under R 4.6.0.

See [input specifications](data/README.md), [output handling](results/README.md), and [analysis workflow](docs/analysis_workflow.md) for additional details. Execution requires access to the study inputs described in the Data Availability section.
