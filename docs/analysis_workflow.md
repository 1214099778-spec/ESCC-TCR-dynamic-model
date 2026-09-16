# Analysis workflow

1. Prepare approved local input tables and the upstream discovery-defined GLIPH2 dictionary described in `data/README.md`.
2. Run script 01 to calculate paired cluster abundance changes for discovery and independent-validation cohorts.
3. Run script 02 to apply discovery-only prevalence/variance filtering, Mann-Whitney ranking and Top10 selection, then fit the L1 logistic model. Fixed-Top10 LOOCV refits on each training fold without repeating feature selection. Locally save the locked model and generated coefficients.
4. Run script 03 to apply the locked model to independent-validation features and to cross-center external repertoires. External response is TRS0/1 versus TRS2/3. No predictive model is refitted on validation or external outcomes.
5. After resolving the clinical-join review item below, run script 04 for discovery PET AIC comparisons and final TCR + MTV integration. The final PET variable is MTV_post. The integration-layer LOOCV holds the precomputed TCR score fixed; it is not a nested validation of the entire upstream pipeline. Independent predictions use discovery-fitted models.
6. Run script 05 for binormal-smoothed ROC/AUC analysis and empirical paired DeLong comparisons. Youden thresholds are estimated on discovery fitted predictions. Pooled and external survival cutpoints are estimated within their respective analyzed data; external survival results are exploratory.

## Syntax verification and scope

R parse checks validate syntax only and do not execute analyses or establish numerical reproducibility. The private copy repairs two split double-bracket expressions in script 04 and seven in script 05 to valid `[[...]]` expressions. No analytical parameters, model coefficients, selection logic, endpoint definitions, or original source files were changed.

## Author review required

- **Clinical join in script 04:** script 01 requires `Target` in the clinical tables, and script 02/03 outputs also contain `Target`. Joining these full tables by `Patient_ID` can produce `Target.x` and `Target.y`, leaving no plain `Target` for subsequent operations. Confirm which outcome column is authoritative and approve a logic-level correction before execution. This has deliberately not been changed by syntax repair.
- **External dictionary ambiguity:** script 01 removes ambiguous clonotype keys across the full dictionary. Script 03 restricts to model features before removing ambiguous keys. A key shared with an unselected cluster may therefore be treated differently. Confirm intended projection rules before altering the analysis.
- **Input quality and folds:** verify complete binary outcomes, valid nonnegative counts, unique clinical linkage rows, consistent cluster names, both paired timepoints, and adequate outcome-class counts for inner cross-validation. These assumptions are not fully enforced in the supplied scripts.
- **Survival interpretation:** predictive-model refitting is absent in external validation, but survival cutpoints are estimated on the analyzed external cohort. Do not describe these survival groups as externally validated locked cutpoints.
- **Metadata:** MIT License (2026) was selected for this code. Confirm the TODO copyright holder, manuscript authors and citation details before public release. No manuscript DOI or publication status is asserted.

## Data boundary

This repository contains source code and explanatory text only. Actual input data, derived tables, spreadsheets, raw sequence files, fitted models, fitted coefficients and generated figures are excluded. References to local input/output filenames and generic schema fields in the code are not embedded study records.
