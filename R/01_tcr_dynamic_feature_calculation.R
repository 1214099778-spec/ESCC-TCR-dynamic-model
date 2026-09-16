# ============================================================
# 01. TCR dynamic feature calculation
# ============================================================

source("R/00_utils.R")

ensure_dir("derived_data")

PSEUDOCOUNT <- 1

DISCOVERY_TCR_FILE <-
  "data/discovery_clonotypes.tsv"

VALIDATION_TCR_FILE <-
  "data/validation_clonotypes.tsv"

CLUSTER_DICTIONARY_FILE <-
  "data/discovery_cluster_dictionary.tsv"

DISCOVERY_CLINICAL_FILE <-
  "data/discovery_clinical.csv"

VALIDATION_CLINICAL_FILE <-
  "data/validation_clinical.csv"


# ------------------------------------------------------------
# 1. Read discovery-defined GLIPH2 cluster dictionary
# ------------------------------------------------------------

cluster_dict <- fread(
  CLUSTER_DICTIONARY_FILE,
  data.table = FALSE
)

assert_columns(
  cluster_dict,
  c(
    "Cluster_Tag",
    "CDR3b",
    "TRBV",
    "TRBJ"
  ),
  "cluster dictionary"
)

cluster_dict <- cluster_dict %>%
  mutate(
    Cluster_Tag = as.character(Cluster_Tag),
    CDR3b = trimws(as.character(CDR3b)),
    TRBV = clean_vj(TRBV),
    TRBJ = clean_vj(TRBJ)
  ) %>%
  distinct(
    Cluster_Tag,
    CDR3b,
    TRBV,
    TRBJ
  )


# ------------------------------------------------------------
# 2. Remove ambiguous clonotype keys
#
# A clonotype key mapping to multiple clusters is excluded
# before external projection.
# ------------------------------------------------------------

cluster_dict <- cluster_dict %>%
  group_by(
    CDR3b,
    TRBV,
    TRBJ
  ) %>%
  filter(
    n_distinct(Cluster_Tag) == 1
  ) %>%
  ungroup()


# ------------------------------------------------------------
# 3. Function to calculate patient-level dynamic features
# ------------------------------------------------------------

calculate_dynamic_matrix <- function(
    tcr_file,
    clinical_file,
    cluster_dict,
    pseudocount = 1
) {
  
  tcr <- fread(
    tcr_file,
    data.table = FALSE
  )
  
  assert_columns(
    tcr,
    c(
      "Patient_ID",
      "Timepoint",
      "CDR3b",
      "TRBV",
      "TRBJ",
      "count"
    ),
    "TCR clonotype table"
  )
  
  tcr <- tcr %>%
    mutate(
      Patient_ID = as.character(Patient_ID),
      
      Timepoint =
        standardize_timepoint(Timepoint),
      
      CDR3b =
        trimws(as.character(CDR3b)),
      
      TRBV =
        clean_vj(TRBV),
      
      TRBJ =
        clean_vj(TRBJ),
      
      count =
        as.numeric(count)
    ) %>%
    filter(
      Timepoint %in% c(
        "Pre",
        "Post"
      ),
      !is.na(CDR3b),
      !is.na(TRBV),
      !is.na(TRBJ),
      !is.na(count)
    )
  
  
  # ----------------------------------------------------------
  # Keep paired patients
  # ----------------------------------------------------------
  
  paired_ids <- tcr %>%
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
  
  
  tcr <- tcr %>%
    filter(
      Patient_ID %in% paired_ids
    )
  
  
  # ----------------------------------------------------------
  # Exact clonotype-to-cluster mapping
  # ----------------------------------------------------------
  
  matched <- tcr %>%
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
  
  
  # ----------------------------------------------------------
  # Aggregate abundance within each GLIPH2 cluster
  # ----------------------------------------------------------
  
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
  
  
  all_clusters <-
    sort(
      unique(
        cluster_dict$Cluster_Tag
      )
    )
  
  
  skeleton <- expand_grid(
    Patient_ID = paired_ids,
    Cluster_Tag = all_clusters
  )
  
  
  feature_long <- skeleton %>%
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
          (Post + pseudocount) /
            (Pre + pseudocount)
        )
    )
  
  
  # ----------------------------------------------------------
  # Patient × cluster dynamic feature matrix
  # ----------------------------------------------------------
  
  feature_matrix <- feature_long %>%
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
  
  
  # ----------------------------------------------------------
  # Add clinical outcome
  # ----------------------------------------------------------
  
  clinical <- fread(
    clinical_file,
    data.table = FALSE
  )
  
  assert_columns(
    clinical,
    c(
      "Patient_ID",
      "Target"
    ),
    "clinical table"
  )
  
  clinical$Patient_ID <-
    as.character(
      clinical$Patient_ID
    )
  
  clinical$Target <-
    clean_target(
      clinical$Target
    )
  
  
  feature_matrix <- feature_matrix %>%
    left_join(
      clinical %>%
        select(
          Patient_ID,
          Target
        ),
      by = "Patient_ID"
    ) %>%
    relocate(
      Patient_ID,
      Target
    )
  
  
  list(
    feature_matrix = feature_matrix,
    feature_long = feature_long
  )
}


# ------------------------------------------------------------
# 4. Discovery cohort
# ------------------------------------------------------------

discovery_features <-
  calculate_dynamic_matrix(
    DISCOVERY_TCR_FILE,
    DISCOVERY_CLINICAL_FILE,
    cluster_dict,
    PSEUDOCOUNT
  )


# ------------------------------------------------------------
# 5. Independent validation cohort
# ------------------------------------------------------------

validation_features <-
  calculate_dynamic_matrix(
    VALIDATION_TCR_FILE,
    VALIDATION_CLINICAL_FILE,
    cluster_dict,
    PSEUDOCOUNT
  )


fwrite(
  discovery_features$feature_matrix,
  "derived_data/discovery_dynamic_features.csv"
)

fwrite(
  validation_features$feature_matrix,
  "derived_data/validation_dynamic_features.csv"
)