## model fitting ----

fit_plain <- function(
  x, y, subject_weight, penalty_weight = NULL, ...)
{
  family <- binomial()
  if (length(unique(y)) > 2L) {
    family <- quasibinomial()
  }
  model <- glm(y ~ x, family = family, weights = subject_weight)
  as.double(coef(model))
}


fit_ridge_cv <- function(
  x, y, subject_weight, penalty_weight = NULL, ...)
{
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  if (is.null(penalty_weight)) {
    penalty_weight <- rep(1.0, ncol(x))
  }
  if (sum(penalty_weight) < 1e-4 || ncol(x) == 1L) {
    return(fit_plain(x, y, subject_weight))
  }
  if (!is.matrix(y) && length(unique(y)) > 2L) {
    y <- cbind(1.0 - y, y)
  }
  penalty_weight <- pmin(pmax(penalty_weight, 1e-4), 1e4)
  model <- cv.glmnet(
    x, y, weights = subject_weight, alpha = 0.05,
    family = "binomial", thresh = 5e-7,
    penalty.factor = penalty_weight,
    nlambda = 100L, lambda.min.ratio = 1e-4)
  as.double(coef(model, s = "lambda.min"))
}


fit_lasso_cv <- function(
  x, y, subject_weight, penalty_weight = NULL, ...)
{
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  if (is.null(penalty_weight)) {
    penalty_weight <- rep(1.0, ncol(x))
  }
  if (sum(penalty_weight) < 1e-4 || ncol(x) == 1L) {
    return(fit_plain(x, y, subject_weight))
  }
  if (!is.matrix(y) && length(unique(y)) > 2L) {
    y <- cbind(1.0 - y, y)
  }
  penalty_weight <- pmin(pmax(penalty_weight, 1e-4), 1e4)
  model <- cv.glmnet(
    x, y, weights = subject_weight,
    family = "binomial", thresh = 5e-7,
    penalty.factor = penalty_weight,
    nlambda = 100L, lambda.min.ratio = 1e-3)
  as.double(coef(model, s = "lambda.min"))
}


fit_lasso_bic <- function(
  x, y, subject_weight, penalty_weight = NULL, ...)
{
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  if (is.null(penalty_weight)) {
    penalty_weight <- rep(1.0, ncol(x))
  }
  if (sum(penalty_weight) < 1e-4 || ncol(x) == 1L) {
    return(fit_plain(x, y, subject_weight))
  }
  if (!is.matrix(y) && length(unique(y)) > 2L) {
    y <- cbind(1.0 - y, y)
  }
  penalty_weight <- pmin(pmax(penalty_weight, 1e-4), 1e4)
  model <- glmnet(
    x, y, weights = subject_weight,
    family = "binomial", thresh = 5e-7,
    penalty.factor = penalty_weight,
    nlambda = 100L, lambda.min.ratio = 1e-3)
  n <- length(y)
  dev <- deviance(model)
  nonzero <- sapply(predict(model, type = "nonzero"), length) + 1L
  complexity <- nonzero * log(n)
  lambda <- model$lambda[which.min(dev + complexity)]
  as.double(predict(model, s = lambda, type = "coefficients"))
}


fit_alasso_cv <- function(
  x, y, subject_weight, penalty_weight = NULL, gamma = 1, ...)
{
  init <- fit_ridge_cv(x, y, subject_weight, penalty_weight, ...)
  alasso_weight <- 1.0 / abs(init[-1L]) ** gamma
  alasso_weight <- pmin(pmax(alasso_weight, 1e-4), 1e4)
  if (is.null(penalty_weight)) {
    penalty_weight <- alasso_weight
  } else {
    penalty_weight <- penalty_weight * alasso_weight
  }
  fit_lasso_cv(x, y, subject_weight, penalty_weight, ...)
}


fit_alasso_bic <- function(
  x, y, subject_weight, penalty_weight = NULL, gamma = 1, ...)
{
  init <- fit_ridge_cv(x, y, subject_weight, penalty_weight, ...)
  alasso_weight <- 1.0 / abs(init[-1L]) ** gamma
  alasso_weight <- pmin(pmax(alasso_weight, 1e-4), 1e4)
  if (is.null(penalty_weight)) {
    penalty_weight <- alasso_weight
  } else {
    penalty_weight <- penalty_weight * alasso_weight
  }
  fit_lasso_bic(x, y, subject_weight, penalty_weight, ...)
}


predict_linear <- function(beta, x, ...)
{
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  plogis(beta[1L] + drop(x %*% beta[-1L]))
}


fit_svm <- function(
  x, y, subject_weight, ...)
{
  if (!(
    missing(subject_weight) || is.null(subject_weight) ||
    sd(subject_weight) < 1e-8)) {
    warning("'subject_weight' not supported in SVM")
  }
  if (requireNamespace("e1071", quietly = TRUE)) {
    y1 <- factor(y, c(0, 1))
    tuning <- e1071::tune.svm(
      x, y1, gamma = c(0.2, 1, 5) / ncol(x), cost = 4.0 ** (-5L : 5L),
      kernel = "radial", type = "C-classification",
      probability = TRUE)
    return(tuning$best.model)
  } else {
    stop("Package e1071 not found")
  }
}


predict_svm <- function(beta, x, ...)
{
  if (requireNamespace("e1071", quietly = TRUE)) {
    return(attr(predict(
      beta, x, probability = TRUE), "probabilities")[, "1"])
  } else {
    stop("Package e1071 not found")
  }
}


fit_rf <- function(
  x, y, subject_weight, ...)
{
  if (requireNamespace("randomForestSRC", quietly = TRUE)) {
    y <- factor(y, c(0, 1))
    return(randomForestSRC::rfsrc(
      y ~ ., data = data.frame(y = y, x = x),
      case.wt = subject_weight))
  } else {
    stop("Package randomForestSRC not found")
  }
}


predict_rf <- function(beta, x, ...)
{
  if (requireNamespace("randomForestSRC", quietly = TRUE)) {
    return(as.numeric(
      predict(beta, data.frame(x = x))$predicted[, "1"]))
  } else {
    stop("Package randomForestSRC not found")
  }
}


fit_xgb <- function(
  x, y, subject_weight, ...)
{
  if (requireNamespace("xgboost", quietly = TRUE)) {
    return(xgboost::xgboost(
      data = x, label = y,
      weight = subject_weight,
      nrounds = 20, verbose = 0,
      objective = "binary:logistic"))
  } else {
    stop("Package xgboost not found")
  }
}


predict_xgb <- function(beta, x, ...)
{
  if (requireNamespace("xgboost", quietly = TRUE)) {
    return(as.numeric(predict(beta, x)))
  } else {
    stop("Package xgboost not found")
  }
}


## multiple models

fit_multiple <- function(method, ...)
{
  all <- list(
    plain = fit_plain,
    ridge_cv = fit_ridge_cv,
    lasso_cv = fit_lasso_cv,
    lasso_bic = fit_lasso_bic,
    alasso_cv = fit_alasso_cv,
    alasso_bic = fit_alasso_bic,
    svm = fit_svm,
    rf = fit_rf,
    xgb = fit_xgb)
  relevant <- all[method]
  beta <- lapply(relevant, function(f) f(...))
  beta
}


predict_multiple <- function(beta_list, ...)
{
  all <- list(
    plain = predict_linear,
    ridge_cv = predict_linear,
    lasso_cv = predict_linear,
    lasso_bic = predict_linear,
    alasso_cv = predict_linear,
    alasso_bic = predict_linear,
    svm = predict_svm,
    rf = predict_rf,
    xgb = predict_xgb)
  predictions <- lapply(names(beta_list), function(method)
    all[[method]](beta_list[[method]], ...))
  predictions <- do.call("cbind", predictions)
  rowMeans(predictions)
}


## roc and auc ----

get_raw_roc_auc <- function(y_true, y_score, subject_weight)
{
  ties <- aggregate(
    data.frame(
      w0 = (1 - y_true) * subject_weight,
      w1 = y_true * subject_weight),
    list(y_score = y_score), sum)
  ties <- ties[order(ties$y_score, decreasing = TRUE), ]

  fact1 <- sum(ties$w1)
  fact0 <- sum(ties$w0)
  total <- fact1 + fact0

  true_positive <- cumsum(ties$w1)
  false_positive <- cumsum(ties$w0)
  true_negative <- fact0 - false_positive

  prediction1 <- true_positive + false_positive
  prediction0 <- total - prediction1

  thr <- ties$y_score
  pct <- prediction1 / total
  acc <- (true_positive + true_negative) / total

  tpr <- true_positive / fact1
  fpr <- false_positive / fact0
  tnr <- true_negative / fact0

  ppv <- true_positive / prediction1
  fdr <- false_positive / prediction1
  npv <- true_negative / prediction0
  npv[length(npv)] <- 1

  sen <- tpr
  rec <- tpr
  spec <- tnr
  prec <- ppv
  f1 <- prec * rec / (prec + rec)

  term1 <- diff(c(0, 1 - spec, 1))
  term2 <- c(0, sen) + c(sen, 1)
  auc <- sum(0.5 * term1 * term2)

  roc <- data.frame(
    thr = thr, pct = pct, acc = acc,
    tpr = tpr, fpr = fpr, tnr = tnr,
    ppv = ppv, fdr = fdr, npv = npv,
    sen = sen, rec = rec, spec = spec, prec = prec,
    f1 = f1)

  list(auc = auc, roc = roc)
}


get_roc <- function(
  y_true, y_score, subject_weight,
  cut = seq(0.001, 0.999, 0.001))
{
  df <- get_raw_roc_auc(y_true, y_score, subject_weight)$roc

  f <- function(x, y) {
    approxfun(x, y, method = "constant", rule = 2L)(cut)
  }
  df <- data.frame(
    cut = cut,
    pct = f(df$thr, df$pct),
    acc = f(df$thr, df$acc),
    tpr = f(df$thr, df$tpr),
    fpr = f(df$thr, df$fpr),
    tnr = f(df$thr, df$tnr),
    ppv = f(df$thr, df$ppv),
    fdr = f(df$thr, df$fdr),
    npv = f(df$thr, df$npv),
    sen = f(df$thr, df$sen),
    spec = f(df$thr, df$spec),
    prec = f(df$thr, df$prec),
    rec = f(df$thr, df$rec),
    f1 = f(df$thr, df$f1))

  return(df)
}


get_auc <- function(
  y_true, y_score, subject_weight)
{
  get_raw_roc_auc(y_true, y_score, subject_weight)$auc
}
