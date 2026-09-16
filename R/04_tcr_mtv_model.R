# ============================================================
# 04. TCR + MTV multimodal model
# ============================================================

source("R/00_utils.R")

ensure_dir("derived_data")
ensure_dir("results/private")


DISCOVERY_TCR_FILE <-
  "derived_data/discovery_tcr_scores.csv"

VALIDATION_TCR_FILE <-
  "derived_data/validation_tcr_scores.csv"

DISCOVERY_CLINICAL_FILE <-
  "data/discovery_clinical.csv"

VALIDATION_CLINICAL_FILE <-
  "data/validation_clinical.csv"


PET_CANDIDATES <- c(
  "SUVmax_post",
  "SUVmean_post",
  "SUVTBR_post",
  "TLG_post",
  "MTV_post"
)

FINAL_PET_VAR <-
  "MTV_post"


# ------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------

discovery_tcr <- fread(
  DISCOVERY_TCR_FILE,
  data.table = FALSE
)

validation_tcr <- fread(
  VALIDATION_TCR_FILE,
  data.table = FALSE
)

discovery_clinical <- fread(
  DISCOVERY_CLINICAL_FILE,
  data.table = FALSE
)

validation_clinical <- fread(
  VALIDATION_CLINICAL_FILE,
  data.table = FALSE
)


# ------------------------------------------------------------
# 2. Candidate PET/CT AIC comparison
# ------------------------------------------------------------

candidate_data <- discovery_tcr %>%
  select(
    Patient_ID,
    Target,
    TCR_Score
  ) %>%
  left_join(
    discovery_clinical %>%
      select(
        Patient_ID,
        all_of(
          PET_CANDIDATES
        )
      ),
    by = "Patient_ID"
  )


aic_table <- bind_rows(
  
  lapply(
    PET_CANDIDATES,
    function(pet_variable) {
      
      analysis_data <- candidate_data %>%
        select(
          Target,
          TCR_Score,
          all_of(
            pet_variable
          )
        ) %>%
        drop_na()
      
      model_formula <- as.formula(
        paste(
          "Target ~ TCR_Score +",
          pet_variable
        )
      )
      
      fit <- suppressWarnings(
        glm(
          model_formula,
          data = analysis_data,
          family = binomial()
        )
      )
      
      data.frame(
        PET_variable =
          pet_variable,
        
        n =
          nrow(
            analysis_data
          ),
        
        AIC =
          AIC(fit)
      )
    }
  )
) %>%
  arrange(AIC)


fwrite(
  aic_table,
  "results/private/pet_candidate_aic.csv"
)


# ------------------------------------------------------------
# 3. Build final TCR + MTV datasets
# ------------------------------------------------------------

discovery_joint <- discovery_tcr %>%
  left_join(
    discovery_clinical,
    by = "Patient_ID"
  ) %>%
  mutate(
    MTV =
      .data[[FINAL_PET_VAR]]
  ) %>%
  drop_na(
    Target,
    TCR_Score,
    MTV
  )


validation_joint <- validation_tcr %>%
  left_join(
    validation_clinical,
    by = "Patient_ID"
  ) %>%
  mutate(
    MTV =
      .data[[FINAL_PET_VAR]]
  ) %>%
  drop_na(
    Target,
    TCR_Score,
    MTV
  )


# ------------------------------------------------------------
# 4. Fit discovery models
# ------------------------------------------------------------

fit_tcr <- suppressWarnings(
  glm(
    Target ~ TCR_Score,
    data = discovery_joint,
    family = binomial()
  )
)


fit_mtv <- suppressWarnings(
  glm(
    Target ~ MTV,
    data = discovery_joint,
    family = binomial()
  )
)


fit_combo <- suppressWarnings(
  glm(
    Target ~
      TCR_Score +
      MTV,
    data = discovery_joint,
    family = binomial()
  )
)


# Do not expose coefficients in public source code.

saveRDS(
  list(
    TCR = fit_tcr,
    MTV = fit_mtv,
    TCR_MTV = fit_combo
  ),
  "results/private/tcr_mtv_models.rds"
)


# ------------------------------------------------------------
# 5. Discovery integration-layer LOOCV
# ------------------------------------------------------------

loocv_combo <- loocv_glm(
  discovery_joint,
  Target ~
    TCR_Score +
    MTV
)


discovery_joint <- discovery_joint %>%
  left_join(
    loocv_combo %>%
      select(
        Patient_ID,
        
        Combo_LOOCV_Score =
          LOOCV_Score,
        
        Combo_LOOCV_Probability =
          LOOCV_Probability
      ),
    by = "Patient_ID"
  )


# ------------------------------------------------------------
# 6. Discovery fitted probabilities
# ------------------------------------------------------------

discovery_joint$TCR_Model_Probability <-
  predict(
    fit_tcr,
    newdata = discovery_joint,
    type = "response"
  )

discovery_joint$MTV_Model_Probability <-
  predict(
    fit_mtv,
    newdata = discovery_joint,
    type = "response"
  )

discovery_joint$Combo_Score <-
  predict(
    fit_combo,
    newdata = discovery_joint,
    type = "link"
  )

discovery_joint$Combo_Probability <-
  predict(
    fit_combo,
    newdata = discovery_joint,
    type = "response"
  )


# ------------------------------------------------------------
# 7. Independent validation
#
# Models are not refitted.
# ------------------------------------------------------------

validation_joint$TCR_Model_Probability <-
  predict(
    fit_tcr,
    newdata = validation_joint,
    type = "response"
  )

validation_joint$MTV_Model_Probability <-
  predict(
    fit_mtv,
    newdata = validation_joint,
    type = "response"
  )

validation_joint$Combo_Score <-
  predict(
    fit_combo,
    newdata = validation_joint,
    type = "link"
  )

validation_joint$Combo_Probability <-
  predict(
    fit_combo,
    newdata = validation_joint,
    type = "response"
  )


fwrite(
  discovery_joint,
  "derived_data/discovery_tcr_mtv_scores.csv"
)

fwrite(
  validation_joint,
  "derived_data/validation_tcr_mtv_scores.csv"
)