# ============================================================
# 02. TCR dynamic score development
# ============================================================

source("R/00_utils.R")

ensure_dir("derived_data")
ensure_dir("results/private")

DISCOVERY_FEATURE_FILE <-
  "derived_data/discovery_dynamic_features.csv"

MIN_NONZERO_PATIENTS <- 3
TOP_N <- 10
INNER_CV_FOLDS <- 5


# ------------------------------------------------------------
# 1. Load discovery cohort
# ------------------------------------------------------------

discovery <- fread(
  DISCOVERY_FEATURE_FILE,
  data.table = FALSE,
  check.names = FALSE
)

assert_columns(
  discovery,
  c(
    "Patient_ID",
    "Target"
  ),
  "discovery feature matrix"
)

discovery$Target <-
  clean_target(
    discovery$Target
  )


feature_cols <- setdiff(
  colnames(discovery),
  c(
    "Patient_ID",
    "Target"
  )
)


# ------------------------------------------------------------
# 2. Remove patients with all feature values missing
# ------------------------------------------------------------

all_na_patient <- apply(
  discovery[
    ,
    feature_cols,
    drop = FALSE
  ],
  1,
  function(x) all(is.na(x))
)

discovery <-
  discovery[
    !all_na_patient,
    ,
    drop = FALSE
  ]


# ------------------------------------------------------------
# 3. Discovery-only unsupervised pre-filtering
# ------------------------------------------------------------

x_all <- discovery[
  ,
  feature_cols,
  drop = FALSE
]

nonzero_n <- colSums(
  x_all != 0,
  na.rm = TRUE
)

variance <- apply(
  x_all,
  2,
  var,
  na.rm = TRUE
)

eligible_features <- names(
  which(
    nonzero_n >= MIN_NONZERO_PATIENTS &
      variance > 0
  )
)


# ------------------------------------------------------------
# 4. Mann–Whitney ranking in discovery cohort
# ------------------------------------------------------------

mw_p <- sapply(
  eligible_features,
  function(feature) {
    
    group_pcr <- discovery[
      discovery$Target == 1,
      feature,
      drop = TRUE
    ]
    
    group_nonpcr <- discovery[
      discovery$Target == 0,
      feature,
      drop = TRUE
    ]
    
    suppressWarnings(
      wilcox.test(
        group_pcr,
        group_nonpcr,
        exact = FALSE
      )$p.value
    )
  }
)


mw_table <- data.frame(
  Feature = eligible_features,
  P_value = mw_p
) %>%
  arrange(
    P_value,
    Feature
  ) %>%
  mutate(
    FDR =
      p.adjust(
        P_value,
        method = "BH"
      )
  )


top10_features <-
  head(
    mw_table$Feature,
    TOP_N
  )


# ------------------------------------------------------------
# 5. L1-penalized logistic regression
# ------------------------------------------------------------

tcr_model <- fit_l1_logistic(
  df = discovery,
  features = top10_features,
  nfolds = INNER_CV_FOLDS,
  lambda_choice = "lambda.min"
)


# ------------------------------------------------------------
# 6. Fit full-discovery TCR dynamic score
# ------------------------------------------------------------

discovery$TCR_Score <-
  predict_l1(
    tcr_model,
    discovery,
    type = "link"
  )

discovery$TCR_Probability <-
  predict_l1(
    tcr_model,
    discovery,
    type = "response"
  )


# ------------------------------------------------------------
# 7. Identify retained non-zero model features
#
# Coefficients are NOT hard-coded or printed.
# ------------------------------------------------------------

coef_matrix <- as.matrix(
  coef(
    tcr_model$cvfit,
    s = tcr_model$lambda
  )
)

coef_table <- data.frame(
  term =
    rownames(
      coef_matrix
    ),
  
  coefficient =
    as.numeric(
      coef_matrix
    )
)

nonzero_features <- coef_table %>%
  filter(
    term != "(Intercept)",
    coefficient != 0
  ) %>%
  pull(term)


# ------------------------------------------------------------
# 8. Fixed-Top10 LOOCV
#
# Top10 features are fixed before LOOCV.
# LASSO coefficients are refitted in each leave-one-out fold.
#
# This evaluates the selected feature set and is not a fully
# nested re-selection of the entire discovery pipeline.
# ------------------------------------------------------------

loocv <- fixed_top10_loocv(
  df = discovery,
  features = top10_features,
  nfolds_inner = INNER_CV_FOLDS
)


discovery <- discovery %>%
  left_join(
    loocv %>%
      select(
        Patient_ID,
        
        TCR_LOOCV_Score =
          LOOCV_Score,
        
        TCR_LOOCV_Probability =
          LOOCV_Probability
      ),
    by = "Patient_ID"
  )


# ------------------------------------------------------------
# 9. Save locked model object locally
# ------------------------------------------------------------

model_bundle <- list(
  
  model =
    tcr_model,
  
  top10_features =
    top10_features,
  
  nonzero_features =
    nonzero_features,
  
  min_nonzero_patients =
    MIN_NONZERO_PATIENTS,
  
  seed =
    GLOBAL_SEED
)


saveRDS(
  model_bundle,
  "results/private/tcr_dynamic_score_model.rds"
)


# Coefficients are generated locally but should not be
# committed to the public GitHub repository.

fwrite(
  coef_table,
  "results/private/tcr_dynamic_score_coefficients.csv"
)

fwrite(
  mw_table,
  "results/private/tcr_feature_ranking.csv"
)


# ------------------------------------------------------------
# 10. Save discovery scores
# ------------------------------------------------------------

fwrite(
  discovery %>%
    select(
      Patient_ID,
      Target,
      TCR_Score,
      TCR_Probability,
      TCR_LOOCV_Score,
      TCR_LOOCV_Probability
    ),
  "derived_data/discovery_tcr_scores.csv"
)