phecap_set_rng_seed <- function(seed)
{
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    oldseed <- .GlobalEnv$.Random.seed
  } else {
    oldseed <- NULL
  }
  set.seed(seed)
  return(oldseed)
}

phecap_restore_rng_seed <- function(oldseed)
{
  if (!is.null(oldseed)) {
    .GlobalEnv$.Random.seed <- oldseed
  } else {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  return(invisible(NULL))
}


phecap_get_roc_auc_with_splits <- function(
  x, y, penalty_weight = NULL, 
  method = "lasso_bic", 
  train_percent = 0.7, num_splits = 200L, 
  start_seed = 1L, verbose = 0L)
{
  oldseed <- phecap_set_rng_seed(start_seed)
  on.exit(phecap_restore_rng_seed(oldseed))
  
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  n_total <- length(y)
  n_train <- as.integer(n_total * train_percent)
  n_test <- n_total - n_train
  
  if (is.character(method)) {
    if (method == "plain") {
      fit_function <- fit_plain
      predict_function <- predict_linear
    } else if (method == "lasso_cv") {
      fit_function <- fit_lasso_cv
      predict_function <- predict_linear
    } else if (method == "lasso_bic") {
      fit_function <- fit_lasso_bic
      predict_function <- predict_linear
    } else if (method == "svm") {
      fit_function <- fit_svm
      predict_function <- predict_svm
    } else if (method == "rf") {
      fit_function <- fit_rf
      predict_function <- predict_rf
    } else {
      stop("Unknown value for method")
    }
  } else if (is.list(method) && all(
    c("fit", "predict") %in% names(method))) {
    fit_function <- method$fit
    predict_function <- method$predict
  } else {
    stop("Unrecognize specification for method")
  }
  
  train_beta <- fit_function(x, y, penalty_weight)
  train_prob <- predict_function(train_beta, x)
  train_roc <- get_roc(y, train_prob)
  train_auc <- get_auc(y, train_prob)
  
  split_roc <- vector("list", num_splits)
  split_auc <- vector("list", num_splits)
  for(sp in seq_len(num_splits)) {
    if (verbose > 0L && (sp %% verbose == 0L || 
                         sp == 1L || sp == num_splits)) {
      cat(sprintf("Split %d/%d\n", sp, num_splits))
    }
    set.seed(start_seed + sp)
    i_train <- sort.int(sample.int(n_total, n_train, FALSE))
    i_test <- setdiff(seq_len(n_total), i_train)
    split_beta <- fit_function(
      x[i_train, , drop = FALSE], y[i_train], penalty_weight)
    split_prob <- predict_function(
      split_beta, x[i_test, , drop = FALSE])
    split_roc[[sp]] <- get_roc(y[i_test], split_prob)
    split_auc[[sp]] <- get_auc(y[i_test], split_prob)
  }
  split_roc <- Reduce(`+`, split_roc) / num_splits
  split_auc <- Reduce(`+`, split_auc) / num_splits
  
  return(list(
    coefficients = train_beta, method = method,
    fit_function = fit_function,
    predict_function = predict_function,
    train_roc = train_roc, train_auc = train_auc,
    split_roc = split_roc, split_auc = split_auc))
}


phecap_read_or_set_frame <- function(source)
{
  data <- NULL
  if (is.character(source)) {
    if (endsWith(source, ".csv")) {
      data <- read.csv(
        source, header = TRUE, sep = ",",
        stringsAsFactors = FALSE)
    } else {
      data <- read.table(
        source, header = TRUE, sep = "\t",
        stringsAsFactors = FALSE)
    }
  } else {
    data <- as.data.frame(source)
  }
  if (nrow(data) == 0L || ncol(data) == 0L) {
    stop("An empty dataset encountered")
  }
  
  return(data)
}


#' Define or Read Datasets for Phenotyping
#' 
#' Specify the data to be used for phenotyping.
#' 
#' @usage PhecapData(
#' data, hu_feature, label, validation,
#' patient_id = NULL, seed = 12300L, feature_transformation = log1p)
#' 
#' @param data A data.frame consisting of all the variables needed for phenotyping,
#' or a character scalar of the path to the data,
#' or a list consisting of either character scalar or data.frame.
#' If a list is given, patient_id cannot be NULL.
#' All the datasets in the list will be joined into a single dataset
#' according to the columns specified by patient_id.
#' @param hu_feature A character scalar or vector specifying the names of
#' one of more healthcare utilization (HU) variables.
#' There variables are always included in the phenotyping model.
#' @param label A character scalar of the column name that gives the phenotype status
#' (1 or TRUE: present, 0 or FALSE: absent).
#' If label is not ready yet, just put a column filled with NA in data.
#' In such cases only the feature extraction step can be done.
#' @param validation A character scalar, a real number strictly between 0 and 1,
#' or an integer not less than 2.
#' If a character scalar is used, it is treated as the column name
#' in the data that specifies whether this observation
#' belongs to the validation samples
#' (1 or TRUE: validation, 0 or FALSE: training).
#' If a real number strictly between 0 and 1 is used, it is treated as
#' the proportion of the validation samples. The actual validation samples
#' will be drawn from all labeled samples.
#' If an integer not less than 2 is used, it is treated as
#' the size of the validation samples. The actual validation samples
#' will be drawn from all labeled samples.
#' @param patient_id A character vector for the column names, if any, that uniquely identifies
#' each patient. Such variables must appear in the data.
#' patient_id can be NULL if such fields are not contained in the data.
#' @param seed If validation samples need to be drawn from all labeled samples,
#' seed specifies the random seed for sampling.
#' @param feature_transformation A function that will be applied to all the features.
#' Since count data are typically right-skewed,
#' by default \code{log1p} will be used.
#' feature_transformation can be NULL, in which case
#' no transformation will be done on any of the feature.
#' 
#' @return An object of class \code{PhecapData}.
#' 
#' @export
PhecapData <- function(
  data, hu_feature, label, validation, 
  patient_id = NULL,  seed = 12300L, 
  feature_transformation = log1p)
{
  if (is.list(data) && !is.data.frame(data)) {
    if (length(data) == 0L) {
      stop("'data' is an empty list")
    }
    if (is.null(patient_id)) {
      stop("'patient_id' cannot be NULL when 'data' is a list")
    }
    patient_id <- as.character(patient_id)
    data <- lapply(data, phecap_read_or_set_frame)
    for (i in seq_along(data)) {
      columns <- names(data[[i]])
      if (!all(patient_id %in% columns)) {
        stop(sprintf(
          "Some 'patient_id' are not found in the %d-th dataset: %s", 
          i, paste0(setdiff(patient_id, columns), collapse = ", ")))
      }
    }
    dt <- data[[1L]]
    for (i in seq_along(data)[-1L]) {
      dt <- merge(dt, data[[i]], by = patient_id, all = TRUE)
    }
    data <- dt
    patient_id_part <- data[patient_id]
    data[patient_id] <- NULL
  } else {
    data <- phecap_read_or_set_frame(data)
    patient_id <- as.character(patient_id)
    patient_id_part <- data[intersect(names(data), patient_id)]
    data <- data[setdiff(names(data), patient_id)]
  }
  patient_id <- patient_id_part
  
  label_part <- data[[label]]
  eps <- sqrt(.Machine$double.eps)
  label_part[abs(label_part) < eps] <- 0
  label_part[abs(label_part - 1) < eps] <- 1
  data[is.na(data)] <- 0
  data[[label]] <- label_part
  
  bad <- !sapply(data, is.numeric)
  if (any(bad)) {
    warning(sprintf(
      "Ignoring non-numeric columns: %s",
      paste0(names(bad)[bad], collapse = ", ")))
    data <- data[names(bad)[!bad]]
  }
  columns <- names(data)
  
  if (!is.character(hu_feature)) {
    stop("'hu_feature' should be of type character ")
  }
  if (!all(hu_feature %in% columns)) {
    stop(sprintf(
      "Some 'hu_feature' are not found in the data: %s",
      paste0(setdiff(hu_feature, columns), collapse = ", ")))
  }
  
  if (!is.character(label)) {
    stop("'label' should be of type character")
  }
  if (length(label) != 1L) {
    stop("'label' should specify a single column")
  }
  if (!all(label %in% columns)) {
    stop(sprintf("%s is not found in the data", label))
  }
  
  index_labeled <- which(!is.na(data[[label]]))
  if (is.character(validation)) {
    index_masked <- which(as.logical(dt[[validation]]))
    training_set <- sort.int(setdiff(index_labeled, index_masked))
    validation_set <- sort.int(intersect(index_labeled, index_masked))
    data[[validation]] <- NULL
  } else {
    oldseed <- phecap_set_rng_seed(seed)
    on.exit(phecap_restore_rng_seed(oldseed))
    num_labeled <- length(index_labeled)
    if (length(validation) != 1L || !all(is.numeric(validation))) {
      stop("Unrecognized value for 'validation'")
    }
    if (validation > 0.0 && validation < 1.0) {
      size <- as.integer(num_labeled * validation + 0.5)
    } else {
      size <- as.integer(validation + 0.5)
    }
    if (size > num_labeled) {
      stop("The size of validation samples exceeds that of labeled samples")
    }
    size <- min(max(size, 0L), num_labeled)
    index_masked <- sample(index_labeled, size, replace = FALSE)
    training_set <- sort.int(setdiff(index_labeled, index_masked))
    validation_set <- sort.int(intersect(index_labeled, index_masked))
  }
  
  if (is.null(feature_transformation)) {
    feature_transformation <- function(x) x
  } else if (!is.function(feature_transformation)) {
    stop("'feature_transformation' should be a function or NULL")
  }
  
  data <- list(
    frame = data,
    hu_feature = hu_feature,
    label = label,
    patient_id = patient_id,
    training_set = training_set,
    validation_set = validation_set,
    feature_transformation = feature_transformation)
  class(data) <- "PhecapData"
  
  return(data)
}


#' Define a Surrogate Variable used inSurrogate-Assisted Feature Extraction (SAFE)
#' 
#' Define a surrogate varible from existing features, and specify
#' associated lower and upper cutoffs. 
#' 
#' This function only stores the definition. No calculation is done.
#' 
#' @usage PhecapSurrogate(variable_names, lower_cutoff = 1L, upper_cutoff = 10L)
#' 
#' @param variable_names a character scalar or vector consisting of variable names.
#' If a vector is given, the value of the surrogate is defined as the sum of
#' the values of each variable.
#' @param lower_cutoff a numeric scalar. If the surrogate value of a patient is less than 
#' or equal to this cutoff, then this patient is treated as a control in SAFE.
#' @param upper_cutoff a numeric scalar. If the surrogate value of a patient is greater than 
#' or equal to this cutoff, then this patient is treated as a case in SAFE.
#' 
#' @return An object of class \code{PhecapSurrogate}.
#' 
#' @export
PhecapSurrogate <- function(
  variable_names, 
  lower_cutoff = 1L, upper_cutoff = 10L)
{
  result <- list(
    variable_names = variable_names,
    lower_cutoff = lower_cutoff,
    upper_cutoff = upper_cutoff)
  class(result) <- "PhecapSurrogate"
  
  return(result)
}


phecap_check_surrogates <- function(
  surrogates, variable_list)
{
  for (surrogate in surrogates) {
    if (!all(surrogate$variable_names %in% variable_list)) {
      stop(sprintf("Variable(s) %s not found", paste0(
        setdiff(surrogate$variable_names, variable_list), 
        collapse = ", ")))
    }
  }
  return(invisible())
}


#' Run Surrogate-Assisted Feature Extraction (SAFE)
#' 
#' Run surrogate-assisted feature extraction (SAFE) using unlabeled data and subsampling. 
#' 
#' In this unlabeled setting, the extremes of each surrogate are used to define cases and 
#' controls. The variables selected are those selected in at least half of the subsamples.
#' 
#' @usage phecap_run_feature_extraction(
#' data, surrogates,
#' subsample_size = 1000L, num_subsamples = 200L,
#' start_seed = 45600L, verbose = 0L)
#' 
#' @param data An object of class PhecapData, obtained by calling PhecapData(...)
#' @param surrogates A list of objects of class PhecapSurrogate, obtained by something like
#' list(PhecapSurrogate(...), PhecapSurrogate(...))
#' @param subsample_size An integer scalar giving the size of each subsample
#' @param num_subsamples The number of subsamples drawn for each surrogate
#' @param start_seed in the i-th subsample, the seed is set to start_seed + i.
#' @param verbose print progress every \code{verbose} subsample if \code{verbose} is positive,
#' or remain quiet if \code{verbose} is zero
#' 
#' @return A character vector consisting of the names of the variables selected,
#' with attribute frequency, which consists of the proportion of being selected
#' among all subsamples for each variable.
#' 
#' @export
phecap_run_feature_extraction <- function(
  data, surrogates, 
  subsample_size = 1000L, num_subsamples = 200L, 
  start_seed = 45600L, verbose = 0L)
{
  oldseed <- phecap_set_rng_seed(start_seed)
  on.exit(phecap_restore_rng_seed(oldseed))
  
  variable_list <- setdiff(names(data$frame), data$label)
  variable_matrix <- as.matrix(data$frame[, variable_list, drop = FALSE])
  phecap_check_surrogates(surrogates, variable_list)
  
  extremes <- lapply(
    surrogates, function(surrogate) {
      x <- rowSums(variable_matrix[, surrogate$variable_names, 
                                   drop = FALSE])
      cases <- which(x >= surrogate$upper_cutoff)
      controls  <- which(x <= surrogate$lower_cutoff)
      if (length(cases) <= subsample_size %/% 2L) {
        stop(sprintf("'%s' has too few cases; %s or %s",
                     paste0(surrogate$variable_names, collapse = "&"),
                     "decrease upper_cutoff", "decrease subsample_size"))
      }
      if (length(controls) <= subsample_size %/% 2L) {
        stop(sprintf("'%s' has too few controls; %s or %s", 
                     paste0(surrogate$variable_names, collapse = "&"),
                     "increase lower_cutoff", "decrease subsample_size"))
      }
      list(cases = cases, controls = controls)
    })
  if (verbose > 0L) {
    message <- sapply(
      extremes, function(z) 
        c(NumCases = length(z$cases),
          NumControls = length(z$controls)))
    colnames(message) <- sapply(surrogates, function(surrogate)
      paste0(surrogate$variable_names, collapse = "&"))
    print(message)
  }
  
  variable_matrix <- data$feature_transformation(variable_matrix)
  selection <- lapply(seq_along(surrogates), function(k) {
    surrogate <- surrogates[[k]]
    exclusion <- match(surrogate$variable_names, variable_list)
    cases <- extremes[[k]]$cases
    controls <- extremes[[k]]$controls
    half_size <- subsample_size %/% 2L
    if (verbose > 0L) {
      cat("Using surrogate", 
          paste0(surrogate$variable_names, collapse = "&"), "\n")
    }
    
    # subsampling
    nonzero <- t(sapply(seq_len(num_subsamples), function(ss) {
      if (verbose > 0L && (ss %% verbose == 0L || 
                           ss == 1L || ss == num_subsamples)) {
        cat(sprintf("Subsample %d/%d\n", ss, num_subsamples))
      }
      set.seed(start_seed + ss)
      ipos <- sort.int(sample(cases, subsample_size %/% 2L, FALSE))
      ineg <- sort.int(sample(controls, subsample_size %/% 2L, FALSE))
      
      y <- c(rep(1.0, length(ipos)), rep(0.0, length(ineg)))
      x <- variable_matrix[c(ipos, ineg), -exclusion, drop = FALSE]
      
      alpha <- fit_lasso_bic(x, y)
      alpha <- alpha[-1L]  # drop intercept
      as.integer(alpha != 0.0)
    }))
    
    nonzero_final <- matrix(1, nrow(nonzero), ncol(variable_matrix))
    nonzero_final[, -exclusion] <- nonzero
    nonzero_final
  })
  
  selection <- do.call("rbind", selection)
  frequency <- colMeans(selection)
  names(frequency) <- variable_list
  selected <- variable_list[frequency >= 0.5]
  
  result <- list(selected = selected, frequency = frequency)
  class(result) <- "PhecapFeatureExtraction"
  return(result)
}


phecap_generate_feature_matrix <- function(
  data, surrogates, feature_selected)
{
  if (is.list(feature_selected) &&
      "selected" %in% names(feature_selected)) {
    feature_selected <- feature_selected$selected
  } else {
    feature_selected <- as.character(feature_selected)
  }
  
  frame <- data$frame
  surrogate_matrix <- sapply(surrogates, function(surrogate) {
    rowSums(frame[, surrogate$variable_names, drop = FALSE])
  })
  surrogate_matrix <- data$feature_transformation(surrogate_matrix)
  colnames(surrogate_matrix) <- sapply(surrogates, function(surrogate) {
    paste0(surrogate$variable_names, collapse = "&")
  })
  
  hu_matrix <- as.matrix(frame[, data$hu_feature, drop = FALSE])
  hu_matrix <- data$feature_transformation(hu_matrix)
  
  other_matrix <- as.matrix(
    frame[, setdiff(feature_selected, c(
      colnames(surrogate_matrix), colnames(hu_matrix))),
      drop = FALSE])
  other_matrix <- data$feature_transformation(other_matrix)
  other_matrix <- qr.resid(qr(cbind(
    1.0, surrogate_matrix, hu_matrix)), other_matrix)
  
  result <- cbind(
    surrogate_matrix,
    hu_matrix,
    other_matrix)
  attr(result, "free") <- ncol(surrogate_matrix) + ncol(hu_matrix)
  return(result)
}


#' Train Phenotyping Model using the Training Labels
#' 
#' Train the phenotyping model on the training dataset,
#' and evaluate its performance via random splits of the training dataset.
#' 
#' @usage phecap_train_phenotyping_model(
#' data, surrogates, feature_selected,
#' method = "lasso_bic", train_percent = 0.7, num_splits = 200L,
#' start_seed = 78900L, verbose = 0L)
#' 
#' @param data an object of class \code{PhecapData}, 
#' obtained by calling \code{PhecapData(...)}.
#' 
#' @param surrogates a list of objects of class \code{PhecapSurrogate}, obtained by 
#' something like \code{list(PhecapSurrogate(...), PhecapSurrogate(...))}.
#' 
#' @param feature_selected a character vector of the features that should be included 
#' in the model, probably returned by \code{phecap_run_feature_extraction}.
#' 
#' @param method Either a character scalar or a list of two components.
#' If a character scalar is used, possible values are
#' \code{'plain'} (logistic regression without penalty),
#' \code{'lasso_cv'} (logistic regression with lasso penalty and CV tuning),
#' \code{'lasso_bic'} (logistic regression with lasso penalty and BIC tuning),
#' \code{'svm'} (support vector machine with CV tuning, package \code{e1071} needed), and
#' \code{'rf'} (random forest with default parameters, package \code{randomForest} needed).
#' If a list is used, it should contain two named components:
#' \code{fit} --- a function for model fitting, and
#' \code{predict} ---- a function for prediction.
#' 
#' @param train_percent The percentage (between 0 and 1) of labels that are used for 
#' model training during random splits
#' 
#' @param num_splits The number of random splits.
#' 
#' @param start_seed in the i-th split, the seed is set to start_seed + i.
#' 
#' @param verbose print progress every verbose splits if verbose is positive,
#' or remain quiet if verbose is zero
#' 
#' @return An object of class \code{PhecapModel}, with components
#' \item{coefficients}{the fitted object}
#' \item{method}{the method used for model training}
#' \item{feature_selected}{the feature selected by SAFE}
#' \item{train_roc}{ROC on training dataset}
#' \item{train_auc}{AUC on training dataset}
#' \item{split_roc}{average ROC on random splits of training dataset}
#' \item{split_auc}{average AUC on random splits of training dataset}
#' \item{fit_function}{the function used for fitting}
#' \item{predict_function}{the function used for prediction}
#' 
#' @export
phecap_train_phenotyping_model <- function(
  data, surrogates, feature_selected, 
  method = "lasso_bic",
  train_percent = 0.7, num_splits = 200L,
  start_seed = 78900L, verbose = 0L)
{
  if (length(data$training_set) < 2L) {
    stop("Too few training samples")
  }
  if (is.list(feature_selected) &&
      "selected" %in% names(feature_selected)) {
    feature_selected <- feature_selected$selected
  } else {
    feature_selected <- as.character(feature_selected)
  }
  
  feature <- phecap_generate_feature_matrix(
    data, surrogates, feature_selected)
  label <- data$frame[, data$label]
  ii <- data$training_set
  x <- feature[ii, , drop = FALSE]
  y <- label[ii]
  penalty_weight <- c(
    rep.int(0.0, attr(feature, "free")),
    rep.int(1.0, ncol(feature) - attr(feature, "free")))
  
  result <- phecap_get_roc_auc_with_splits(
    x, y, penalty_weight, method = method,
    train_percent = train_percent, num_splits = num_splits,
    start_seed = start_seed, verbose = verbose)
  if (is.numeric(result$coefficients)) {
    names(result$coefficients) <- c("(Intercept)", colnames(x))
  }
  result$surrogates <- surrogates
  result$feature_selected <- feature_selected
  
  class(result) <- "PhecapModel"
  return(result)
}


#' Validate the Phenotyping Model using the Validation Labels
#' 
#' Apply the trained model to all patients in the validation dataset,
#' and measure the prediction accuracy via ROC and AUC.
#' 
#' @usage phecap_validate_phenotyping_model(data, model)
#' 
#' @param data an object of class \code{PhecapData}, 
#' obtained by calling \code{PhecapData(...)}
#' @param model an object of class \code{PhecapModel}, obtained by calling
#' \code{phecap_train_phenotyping_model}.
#' 
#' @return An object of class \code{PhecapValidation}, with components
#' \item{method}{the method used for model training}
#' \item{train_roc}{ROC on training dataset}
#' \item{train_auc}{AUC on training dataset}
#' \item{split_roc}{average ROC on random splits of training dataset}
#' \item{split_auc}{average AUC on random splits of training dataset}
#' \item{valid_roc}{ROC on validation dataset}
#' \item{valid_auc}{AUC on validation dataset}
#' 
#' @export
phecap_validate_phenotyping_model <- function(
  data, model)
{
  if (length(data$validation_set) < 2L) {
    stop("Too few validation samples")
  }
  
  surrogates <- model$surrogates
  feature_selected <- model$feature_selected
  
  feature <- phecap_generate_feature_matrix(
    data, surrogates, feature_selected)
  label <- data$frame[, data$label]
  ii <- data$validation_set
  x <- feature[ii, , drop = FALSE]
  y <- label[ii]
  
  prediction <- model$predict_function(model$coefficients, x)
  valid_roc <- get_roc(y, prediction)
  valid_auc <- get_auc(y, prediction)
  
  result <- list(
    coefficients = model$coefficients,
    method = model$method,
    train_roc = model$train_roc,
    train_auc = model$train_auc,
    split_roc = model$split_roc,
    split_auc = model$split_auc,
    valid_roc = valid_roc,
    valid_auc = valid_auc)
  class(result) <- "PhecapValidation"
  
  return(result)
}


#' Predict Phenotype
#' 
#' Compute predicted probability of having the phenotype
#' for each patient in the dataset.
#' 
#' @usage phecap_predict_phenotype(data, model)
#' 
#' @param data an object of class \code{PhecapData}, obtained by calling \code{PhecapData(...)}.
#' @param model an object of class \code{PhecapModel}, probably returned from
#' \code{phecap_train_phenotyping_model}.
#' 
#' @return A \code{data.frame} with two columns:
#' \item{patient_index}{patient identifier},
#' \item{prediction}{predicted phenotype}.
#' 
#' @export
phecap_predict_phenotype <- function(
  data, model)
{
  surrogates <- model$surrogates
  feature_selected <- model$feature_selected
  
  feature <- phecap_generate_feature_matrix(
    data, surrogates, feature_selected)
  prediction <- model$predict_function(
    model$coefficients, feature)
  
  result <- data.frame(
    data$patient_id,
    prediction = prediction)
  return(result)
}


#' @export
print.PhecapData <- function(x, ...)
{
  cat("PheCAP Data\n")
  cat(sprintf(
    "Feature: %d observations of %d variables\n",
    nrow(x$frame), ncol(x$frame) - 1L))
  cat(sprintf(
    "Label: %d yes, %d no, %d missing\n",
    sum(x$frame[[x$label]] == 1, na.rm = TRUE),
    sum(x$frame[[x$label]] == 0, na.rm = TRUE),
    sum(is.na(x$frame[[x$label]]))))
  cat(sprintf(
    "Size of training samples: %d\n", 
    length(x$training_set)))
  cat(sprintf(
    "Size of validation samples: %d\n", 
    length(x$validation_set)))
}


#' @export
print.PhecapFeatureExtraction <- function(x, ...)
{
  cat("Feature(s) selected by",
      "surrogate-assisted feature extraction (SAFE)\n")
  print(x$selected, ...)
}


#' @export
print.PhecapModel <- function(x, ...)
{
  cat("Phenotyping model:\n")
  print(x$coefficients, ...)
  cat("AUC on training data:", 
      format(x$train_auc, digits = 3L), "\n")
  cat("Average AUC on random splits:", 
      format(x$split_auc, digits = 3L), "\n")
}


#' @export
print.PhecapValidation <- function(x, ...)
{
  cat("AUC on validation data:", 
      format(x$valid_auc, digits = 3L), "\n")
  cat("AUC on training data:", 
      format(x$train_auc, digits = 3L), "\n")
  cat("Average AUC on random splits:", 
      format(x$split_auc, digits = 3L), "\n")
}


#' Plot ROC and Related Curves for Phenotyping Models
#' 
#' Plot ROC-like curves to illustrate phenotyping accuracy.
#' 
#' @usage phecap_plot_roc_curves(
#' x, axis_x = "1 - spec", axis_y = "sen",
#' what = c("training", "random-splits", "validation"),
#' ggplot = TRUE, ...)
#' 
#' @param x either a single object of class PhecapModel or PhecapValidation 
#' (returned from \code{phecap_train_phenotyping_model} or
#' \code{phecap_validate_phenotyping_model}), or a named list of such objects
#' 
#' @param axis_x 
#' an expression that leads to the \code{x} coordinate.
#' Recognized quantities include:
#' \code{cut} (probability cutoff),
#' \code{pct} (percent of predicted cases),
#' \code{acc} (accuracy),
#' \code{tpr} (true positive rate),
#' \code{fpr} (false positive rate),
#' \code{tnr} (true negative rate),
#' \code{ppv} (positive predictive value),
#' \code{fdr} (false discovery rate),
#' \code{npv} (negative predictive value),
#' \code{sen} (sensitivity),
#' \code{spec} (specificity),
#' \code{prec} (precision),
#' \code{rec} (recall),
#' \code{f1} (F1 score).
#' @param axis_y an expression that leads to the \code{y} coordinate.
#' Recognized quantities are the same as those in \code{axis_x}.
#' @param what TThe curves to be included in the figure.
#' @param ggplot if TRUE and ggplot2 is installed, ggplot will be used for the figure.
#' Otherwise, the base R graphics functions will be used.
#' @param \dots arguments to be ignored.
#' 
#' @export
phecap_plot_roc_curves <- function(
  x, axis_x = "1 - spec", axis_y = "sen", 
  what = c("training", "random-splits", "validation"),
  ggplot = TRUE,
  ...)
{
  object <- x
  if (is(x, "PhecapModel") || is(x, "PhecapValidation")) {
    object <- list(x)
    names(object) <- deparse(substitute(x))
  } else if (!is.list(x)) {
    stop("Not a PhecapModel / PhecapValidation object or a list of them")
  } else if (is.null(names(x))) {
    stop("List should be named")
  }
  
  df <- vector("list", length(object) * 3L)
  ii <- 1L
  for (kk in seq_along(object)) {
    oo <- object[[kk]]
    for (ww in what) {
      if (ww == "training" && "train_roc" %in% names(oo)) {
        df[[ii]] <- data.frame(
          cut = oo$train_roc$cut,
          value_x = eval(parse(text = axis_x), oo$train_roc),
          value_y = eval(parse(text = axis_y), oo$train_roc),
          kk = names(object)[kk], ww = ww)
      } else if (ww == "random-splits" && "split_roc" %in% names(oo)) {
        df[[ii]] <- data.frame(
          cut = oo$split_roc$cut,
          value_x = eval(parse(text = axis_x), oo$split_roc),
          value_y = eval(parse(text = axis_y), oo$split_roc),
          kk = names(object)[kk], ww = ww)
      } else if (ww == "validation" && "valid_roc" %in% names(oo)) {
        df[[ii]] <- data.frame(
          cut = oo$valid_roc$cut,
          value_x = eval(parse(text = axis_x), oo$valid_roc),
          value_y = eval(parse(text = axis_y), oo$valid_roc),
          kk = names(object)[kk], ww = ww)
      }
      ii <- ii + 1L
    }
  }
  
  df <- do.call("rbind", df)
  df <- aggregate(cbind(cut, value_y) ~ kk + ww + value_x, df, 
                  max)
  df <- df[order(df$kk, df$ww, df$cut, df$value_x, df$value_y), ]
  if (length(unique(df$kk)) > 1L) {
    if (length(unique(df$ww)) > 1L) {
      df$type <- paste(df$kk, df$ww, sep = ":")
    } else {
      df$type <- df$kk
    }
  } else if (length(unique(df$ww)) > 1L) {
    df$type <- df$ww
  }
  
  if (ggplot && requireNamespace("ggplot2", quietly = TRUE)) {
    pp <- ggplot2::ggplot(df)
    if ("type" %in% names(df)) {
      pp <- pp + ggplot2::geom_path(ggplot2::aes(
        x = df$value_x, y = df$value_y, color = type))
    } else {
      pp <- pp + ggplot2::geom_path(ggplot2::aes(
        x = df$value_x, y = df$value_y))
    }
    pp <- pp + ggplot2::xlab(axis_x) + ggplot2::ylab(axis_y)
    print(pp)
  } else {
    plot(NULL, NULL, 
         xlim = range(df$value_x), ylim = range(df$value_y),
         xlab = axis_x, ylab = axis_y)
    if ("type" %in% names(df)) {
      col <- 1L
      for (type in unique(df$type)) {
        lines(df[df$type == type, "value_x"],
              df[df$type == type, "value_y"], col = col)
        col <- col + 1L
      }
    } else {
      lines(df[, "value_x"], df[, "value_y"])
    }
    legend("bottomright", 
           legend = unique(df$type),
           lty = 1,
           col = seq_along(unique(df$type)))
  }
}