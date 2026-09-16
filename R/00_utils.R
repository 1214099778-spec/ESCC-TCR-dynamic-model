# ============================================================
# Shared utilities
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(glmnet)
  library(pROC)
})

GLOBAL_SEED <- 20260717L

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

assert_columns <- function(df, required, object_name = "data") {
  missing_cols <- setdiff(required, colnames(df))
  
  if (length(missing_cols) > 0) {
    stop(
      object_name,
      " is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
}

clean_target <- function(x) {
  
  if (is.numeric(x) || is.integer(x)) {
    return(as.numeric(x))
  }
  
  x <- tolower(trimws(as.character(x)))
  
  out <- rep(NA_real_, length(x))
  
  out[x %in% c(
    "1", "pcr", "yes", "responder", "positive"
  )] <- 1
  
  out[x %in% c(
    "0", "non-pcr", "nonpcr",
    "no", "non-responder",
    "nonresponder", "negative"
  )] <- 0
  
  out
}

standardize_timepoint <- function(x) {
  
  x0 <- tolower(trimws(as.character(x)))
  
  case_when(
    x0 %in% c(
      "pre", "baseline", "before",
      "pretreatment", "pre-treatment"
    ) ~ "Pre",
    
    x0 %in% c(
      "post", "after",
      "post-treatment", "posttreatment"
    ) ~ "Post",
    
    TRUE ~ as.character(x)
  )
}

clean_vj <- function(x) {
  
  x <- trimws(as.character(x))
  
  # Remove allele suffix:
  # TRBV12-3*01 -> TRBV12-3
  x <- sub("\\*.*$", "", x)
  
  x[x %in% c("", "NA", "NaN", "NULL", "<NA>")] <- NA_character_
  
  x
}

balanced_weights <- function(y) {
  
  y <- as.numeric(y)
  
  n  <- length(y)
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  
  if (n1 == 0 || n0 == 0) {
    stop("Both outcome classes are required.")
  }
  
  ifelse(
    y == 1,
    n / (2 * n1),
    n / (2 * n0)
  )
}

make_stratified_foldid <- function(
    y,
    nfolds = 5,
    seed = GLOBAL_SEED
) {
  
  set.seed(seed)
  
  foldid <- rep(NA_integer_, length(y))
  
  for (cl in sort(unique(y))) {
    
    idx <- which(y == cl)
    
    foldid[idx] <- sample(
      rep(seq_len(nfolds), length.out = length(idx))
    )
  }
  
  foldid
}

fit_l1_logistic <- function(
    df,
    features,
    nfolds = 5,
    lambda_choice = "lambda.min",
    seed = GLOBAL_SEED
) {
  
  x <- as.matrix(
    df[, features, drop = FALSE]
  )
  
  storage.mode(x) <- "double"
  
  y <- as.numeric(df$Target)
  
  nfolds_use <- min(
    nfolds,
    min(table(y))
  )
  
  foldid <- make_stratified_foldid(
    y,
    nfolds = nfolds_use,
    seed = seed
  )
  
  cvfit <- cv.glmnet(
    x = x,
    y = y,
    family = "binomial",
    alpha = 1,
    weights = balanced_weights(y),
    foldid = foldid,
    type.measure = "deviance",
    standardize = TRUE
  )
  
  list(
    cvfit = cvfit,
    lambda = cvfit[[lambda_choice]],
    features = features,
    lambda_choice = lambda_choice
  )
}

predict_l1 <- function(
    model,
    new_df,
    type = c("link", "response")
) {
  
  type <- match.arg(type)
  
  x <- as.matrix(
    new_df[, model$features, drop = FALSE]
  )
  
  storage.mode(x) <- "double"
  
  as.numeric(
    predict(
      model$cvfit,
      newx = x,
      s = model$lambda,
      type = type
    )
  )
}

fixed_top10_loocv <- function(
    df,
    features,
    nfolds_inner = 5
) {
  
  n <- nrow(df)
  
  result <- data.frame(
    Patient_ID = df$Patient_ID,
    Target = df$Target,
    LOOCV_Score = NA_real_,
    LOOCV_Probability = NA_real_
  )
  
  for (i in seq_len(n)) {
    
    train_fold <- df[-i, , drop = FALSE]
    test_fold  <- df[i, , drop = FALSE]
    
    model_i <- fit_l1_logistic(
      train_fold,
      features = features,
      nfolds = nfolds_inner,
      seed = GLOBAL_SEED + i
    )
    
    result$LOOCV_Score[i] <- predict_l1(
      model_i,
      test_fold,
      type = "link"
    )
    
    result$LOOCV_Probability[i] <- predict_l1(
      model_i,
      test_fold,
      type = "response"
    )
  }
  
  result
}

loocv_glm <- function(
    df,
    formula
) {
  
  n <- nrow(df)
  
  output <- data.frame(
    Patient_ID = df$Patient_ID,
    Target = df$Target,
    LOOCV_Score = NA_real_,
    LOOCV_Probability = NA_real_
  )
  
  for (i in seq_len(n)) {
    
    train_fold <- df[-i, , drop = FALSE]
    test_fold  <- df[i, , drop = FALSE]
    
    fit_i <- suppressWarnings(
      glm(
        formula,
        data = train_fold,
        family = binomial()
      )
    )
    
    output$LOOCV_Score[i] <- as.numeric(
      predict(
        fit_i,
        newdata = test_fold,
        type = "link"
      )
    )
    
    output$LOOCV_Probability[i] <- as.numeric(
      predict(
        fit_i,
        newdata = test_fold,
        type = "response"
      )
    )
  }
  
  output
}