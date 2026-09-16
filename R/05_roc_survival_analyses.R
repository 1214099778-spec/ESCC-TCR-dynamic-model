# ============================================================
# 05. ROC and survival analyses
# ============================================================

source("R/00_utils.R")

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(ggplot2)
})

ensure_dir("results/private")
ensure_dir("figures")

BOOT_N <- 2000


# ------------------------------------------------------------
# ROC helper
# ------------------------------------------------------------

calculate_smoothed_roc <- function(
    outcome,
    predictor,
    boot_n = 2000
) {
  
  empirical_roc <- roc(
    outcome,
    predictor,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )
  
  
  smoothed_roc <- smooth(
    empirical_roc,
    method = "binormal"
  )
  
  
  set.seed(GLOBAL_SEED)
  
  
  ci <- ci.auc(
    smoothed_roc,
    boot.n = boot_n,
    boot.stratified = TRUE
  )
  
  
  list(
    empirical =
      empirical_roc,
    
    smooth =
      smoothed_roc,
    
    AUC =
      as.numeric(
        auc(
          smoothed_roc
        )
      ),
    
    CI_low =
      as.numeric(
        ci[1]
      ),
    
    CI_high =
      as.numeric(
        ci[3]
      )
  )
}


# ============================================================
# 1. TCR ROC
# ============================================================

discovery_tcr <- fread(
  "derived_data/discovery_tcr_scores.csv",
  data.table = FALSE
)

validation_tcr <- fread(
  "derived_data/validation_tcr_scores.csv",
  data.table = FALSE
)


roc_discovery_tcr <-
  calculate_smoothed_roc(
    discovery_tcr$Target,
    discovery_tcr$TCR_LOOCV_Probability,
    BOOT_N
  )


roc_validation_tcr <-
  calculate_smoothed_roc(
    validation_tcr$Target,
    validation_tcr$TCR_Probability,
    BOOT_N
  )


# ============================================================
# 2. External response ROC
# ============================================================

if (
  file.exists(
    "derived_data/external_tcr_scores.csv"
  )
) {
  
  external_tcr <- fread(
    "derived_data/external_tcr_scores.csv",
    data.table = FALSE
  )
  
  
  roc_external_tcr <-
    calculate_smoothed_roc(
      external_tcr$External_Response,
      external_tcr$TCR_Probability,
      BOOT_N
    )
}


# ============================================================
# 3. TCR + MTV ROC
# ============================================================

discovery_combo <- fread(
  "derived_data/discovery_tcr_mtv_scores.csv",
  data.table = FALSE
)

validation_combo <- fread(
  "derived_data/validation_tcr_mtv_scores.csv",
  data.table = FALSE
)


# Use the LOOCV linear predictor for discovery combined-model
# discrimination to avoid numerical probability saturation.

roc_discovery_combo <-
  calculate_smoothed_roc(
    discovery_combo$Target,
    discovery_combo$Combo_LOOCV_Score,
    BOOT_N
  )


roc_validation_combo <-
  calculate_smoothed_roc(
    validation_combo$Target,
    validation_combo$Combo_Probability,
    BOOT_N
  )


# ============================================================
# 4. DeLong tests
#
# DeLong tests use empirical ROC objects.
# ============================================================

roc_tcr_validation <- roc(
  validation_combo$Target,
  validation_combo$TCR_Model_Probability,
  quiet = TRUE
)

roc_mtv_validation <- roc(
  validation_combo$Target,
  validation_combo$MTV_Model_Probability,
  quiet = TRUE
)

roc_combo_validation <- roc(
  validation_combo$Target,
  validation_combo$Combo_Probability,
  quiet = TRUE
)


delong_tcr_combo <- roc.test(
  roc_tcr_validation,
  roc_combo_validation,
  paired = TRUE,
  method = "delong"
)


delong_mtv_combo <- roc.test(
  roc_mtv_validation,
  roc_combo_validation,
  paired = TRUE,
  method = "delong"
)


# ============================================================
# 5. Discovery Youden thresholds
# ============================================================

get_youden <- function(
    outcome,
    probability
) {
  
  roc_obj <- roc(
    outcome,
    probability,
    quiet = TRUE
  )
  
  as.numeric(
    coords(
      roc_obj,
      x = "best",
      best.method = "youden",
      ret = "threshold",
      transpose = FALSE
    )[[1]]
  )
}


threshold_tcr <- get_youden(
  discovery_combo$Target,
  discovery_combo$TCR_Model_Probability
)

threshold_mtv <- get_youden(
  discovery_combo$Target,
  discovery_combo$MTV_Model_Probability
)

threshold_combo <- get_youden(
  discovery_combo$Target,
  discovery_combo$Combo_Probability
)


# ============================================================
# 6. Survival cutpoint
# ============================================================

derive_cutpoint <- function(
    data,
    time_col,
    event_col,
    score_col
) {
  
  survival_data <- data %>%
    filter(
      !is.na(
        .data[[time_col]]
      ),
      !is.na(
        .data[[event_col]]
      ),
      !is.na(
        .data[[score_col]]
      )
    )
  
  
  cut_result <- surv_cutpoint(
    survival_data,
    time = time_col,
    event = event_col,
    variables = score_col
  )
  
  
  as.numeric(
    cut_result$cutpoint[
      score_col,
      "cutpoint"
    ]
  )
}


# ============================================================
# 7. Kaplan–Meier function
#
# Higher TCR score represents the favorable / low-risk group.
# ============================================================

run_km <- function(
    data,
    time_col,
    event_col,
    score_col,
    title
) {
  
  cutoff <- derive_cutpoint(
    data,
    time_col,
    event_col,
    score_col
  )
  
  
  km_data <- data %>%
    filter(
      !is.na(
        .data[[time_col]]
      ),
      !is.na(
        .data[[event_col]]
      ),
      !is.na(
        .data[[score_col]]
      )
    )
  
  
  km_data$Risk_Group <- ifelse(
    km_data[[score_col]] >= cutoff,
    "Low risk",
    "High risk"
  )
  
  
  km_data$Risk_Group <- factor(
    km_data$Risk_Group,
    levels = c(
      "High risk",
      "Low risk"
    )
  )
  
  
  fit <- survfit(
    as.formula(
      paste0(
        "Surv(",
        time_col,
        ",",
        event_col,
        ") ~ Risk_Group"
      )
    ),
    data = km_data
  )
  
  
  ggsurvplot(
    fit,
    data = km_data,
    risk.table = TRUE,
    pval = TRUE,
    legend.title = "",
    legend.labs = c(
      "High risk",
      "Low risk"
    ),
    title = title
  )
}


# ============================================================
# 8. Pooled discovery + validation survival
# ============================================================

pooled <- bind_rows(
  discovery_combo %>%
    mutate(
      Cohort =
        "Discovery"
    ),
  
  validation_combo %>%
    mutate(
      Cohort =
        "Validation"
    )
)


if (
  all(
    c(
      "PFS_time",
      "PFS_event"
    ) %in%
    colnames(pooled)
  )
) {
  
  run_km(
    pooled,
    "PFS_time",
    "PFS_event",
    "TCR_Score",
    "TCR model PFS"
  )
  
  
  run_km(
    pooled,
    "PFS_time",
    "PFS_event",
    "Combo_Score",
    "Combined model PFS"
  )
}


if (
  all(
    c(
      "CSS_time",
      "CSS_event"
    ) %in%
    colnames(pooled)
  )
) {
  
  run_km(
    pooled,
    "CSS_time",
    "CSS_event",
    "TCR_Score",
    "TCR model CSS"
  )
  
  
  run_km(
    pooled,
    "CSS_time",
    "CSS_event",
    "Combo_Score",
    "Combined model CSS"
  )
}


# ============================================================
# 9. External exploratory survival analysis
# ============================================================

if (
  exists(
    "external_tcr"
  )
) {
  
  if (
    all(
      c(
        "PFS_time",
        "PFS_event"
      ) %in%
      colnames(
        external_tcr
      )
    )
  ) {
    
    run_km(
      external_tcr,
      "PFS_time",
      "PFS_event",
      "TCR_Score",
      "External cohort TCR model PFS"
    )
  }
  
  
  if (
    all(
      c(
        "OS_time",
        "OS_event"
      ) %in%
      colnames(
        external_tcr
      )
    )
  ) {
    
    run_km(
      external_tcr,
      "OS_time",
      "OS_event",
      "TCR_Score",
      "External cohort TCR model OS"
    )
  }
}