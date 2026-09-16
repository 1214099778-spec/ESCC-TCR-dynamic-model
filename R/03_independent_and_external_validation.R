# ============================================================
# 03. Independent and external validation
# ============================================================

source("R/00_utils.R")

ensure_dir("derived_data")

MODEL_FILE <-
  "results/private/tcr_dynamic_score_model.rds"

VALIDATION_FEATURE_FILE <-
  "derived_data/validation_dynamic_features.csv"

EXTERNAL_TCR_FILE <-
  "data/external_clonotypes.tsv"

EXTERNAL_CLINICAL_FILE <-
  "data/external_clinical.csv"

CLUSTER_DICTIONARY_FILE <-
  "data/discovery_cluster_dictionary.tsv"

PSEUDOCOUNT <- 1


# ============================================================
# Part A
# Independent validation cohort
# ============================================================

model_bundle <-
  readRDS(
    MODEL_FILE
  )

tcr_model <-
  model_bundle$model

model_features <-
  model_bundle$top10_features


validation <- fread(
  VALIDATION_FEATURE_FILE,
  data.table = FALSE,
  check.names = FALSE
)


missing_features <- setdiff(
  model_features,
  colnames(validation)
)

if (length(missing_features) > 0) {
  
  stop(
    "Validation cohort is missing locked model features."
  )
}


# No refitting.

validation$TCR_Score <-
  predict_l1(
    tcr_model,
    validation,
    type = "link"
  )

validation$TCR_Probability <-
  predict_l1(
    tcr_model,
    validation,
    type = "response"
  )


fwrite(
  validation %>%
    select(
      Patient_ID,
      Target,
      TCR_Score,
      TCR_Probability
    ),
  "derived_data/validation_tcr_scores.csv"
)


# ============================================================
# Part B
# Cross-center external cohort
# ============================================================

external_tcr <- fread(
  EXTERNAL_TCR_FILE,
  data.table = FALSE
)

external_clinical <- fread(
  EXTERNAL_CLINICAL_FILE,
  data.table = FALSE
)

cluster_dict <- fread(
  CLUSTER_DICTIONARY_FILE,
  data.table = FALSE
)


external_tcr <- external_tcr %>%
  mutate(
    Patient_ID =
      as.character(
        Patient_ID
      ),
    
    Timepoint =
      standardize_timepoint(
        Timepoint
      ),
    
    CDR3b =
      trimws(
        as.character(
          CDR3b
        )
      ),
    
    TRBV =
      clean_vj(TRBV),
    
    TRBJ =
      clean_vj(TRBJ),
    
    count =
      as.numeric(count)
  )


cluster_dict <- cluster_dict %>%
  mutate(
    CDR3b =
      trimws(
        as.character(
          CDR3b
        )
      ),
    
    TRBV =
      clean_vj(TRBV),
    
    TRBJ =
      clean_vj(TRBJ)
  ) %>%
  
  filter(
    Cluster_Tag %in%
      model_features
  ) %>%
  
  distinct(
    Cluster_Tag,
    CDR3b,
    TRBV,
    TRBJ
  )


# Remove ambiguous matching keys.

cluster_dict <- cluster_dict %>%
  group_by(
    CDR3b,
    TRBV,
    TRBJ
  ) %>%
  filter(
    n_distinct(
      Cluster_Tag
    ) == 1
  ) %>%
  ungroup()


# ------------------------------------------------------------
# Identify paired external patients
# ------------------------------------------------------------

paired_ids <- external_tcr %>%
  distinct(
    Patient_ID,
    Timepoint
  ) %>%
  count(
    Patient_ID,
    Timepoint
  ) %>%
  pivot_wider(
    names_from = Timepoint,
    values_from = n,
    values_fill = 0
  ) %>%
  filter(
    Pre > 0,
    Post > 0
  ) %>%
  pull(Patient_ID)


# ------------------------------------------------------------
# Exact match:
# CDR3b + TRBV + TRBJ
# ------------------------------------------------------------

matched <- external_tcr %>%
  
  filter(
    Patient_ID %in%
      paired_ids
  ) %>%
  
  inner_join(
    cluster_dict,
    by = c(
      "CDR3b",
      "TRBV",
      "TRBJ"
    )
  ) %>%
  
  distinct(
    Patient_ID,
    Timepoint,
    Cluster_Tag,
    CDR3b,
    TRBV,
    TRBJ,
    .keep_all = TRUE
  )


# ------------------------------------------------------------
# Cluster-level abundance
# ------------------------------------------------------------

cluster_abundance <- matched %>%
  group_by(
    Patient_ID,
    Timepoint,
    Cluster_Tag
  ) %>%
  summarise(
    abundance =
      sum(
        count,
        na.rm = TRUE
      ),
    .groups = "drop"
  )


skeleton <- expand_grid(
  Patient_ID = paired_ids,
  Cluster_Tag = model_features
)


external_long <- skeleton %>%
  left_join(
    cluster_abundance %>%
      pivot_wider(
        names_from = Timepoint,
        values_from = abundance,
        values_fill = 0
      ),
    by = c(
      "Patient_ID",
      "Cluster_Tag"
    )
  ) %>%
  mutate(
    Pre =
      ifelse(
        is.na(Pre),
        0,
        Pre
      ),
    
    Post =
      ifelse(
        is.na(Post),
        0,
        Post
      ),
    
    log2FC =
      log2(
        (Post + PSEUDOCOUNT) /
          (Pre + PSEUDOCOUNT)
      )
  )


external_features <- external_long %>%
  select(
    Patient_ID,
    Cluster_Tag,
    log2FC
  ) %>%
  pivot_wider(
    names_from = Cluster_Tag,
    values_from = log2FC,
    values_fill = 0
  )


# ------------------------------------------------------------
# Apply locked discovery model
# ------------------------------------------------------------

external_features$TCR_Score <-
  predict_l1(
    tcr_model,
    external_features,
    type = "link"
  )

external_features$TCR_Probability <-
  predict_l1(
    tcr_model,
    external_features,
    type = "response"
  )


external_results <- external_features %>%
  left_join(
    external_clinical,
    by = "Patient_ID"
  )


# ------------------------------------------------------------
# Pathological responder endpoint:
# TRS0/1 vs TRS2/3
# ------------------------------------------------------------

if ("TRS" %in%
    colnames(
      external_results
    )) {
  
  trs <- toupper(
    gsub(
      "\\s+",
      "",
      external_results$TRS
    )
  )
  
  external_results$External_Response <-
    ifelse(
      trs %in%
        c(
          "TRS0",
          "TRS1"
        ),
      1,
      ifelse(
        trs %in%
          c(
            "TRS2",
            "TRS3"
          ),
        0,
        NA
      )
    )
}


fwrite(
  external_results,
  "derived_data/external_tcr_scores.csv"
)