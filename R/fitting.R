## model fitting ----

fit_plain <- function(x, y, ...)
{
  family <- binomial()
  if (length(unique(y)) > 2L) {
    family <- quasibinomial()
  }
  as.double(coef(glm(y ~ x, family = family)))
}


fit_lasso_cv <- function(
  x, y, penalty_weight = NULL, ...)
{
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  if (is.null(penalty_weight)) {
    penalty_weight <- rep(1.0, ncol(x))
  }
  if (sum(penalty_weight) < 1e-4 || ncol(x) == 1L) {
    return(fit_plain(x, y))
  }
  if (!is.matrix(y) && length(unique(y)) > 2L) {
    y <- cbind(1.0 - y, y)
  }
  penalty_weight <- pmin(pmax(penalty_weight, 1e-4), 1e4)
  model <- cv.glmnet(
    x, y, family = "binomial", thresh = 5e-7,
    penalty.factor = penalty_weight, 
    nlambda = 50L, lambda.min.ratio = 5e-3)
  as.double(coef(model, s = "lambda.min"))
}


fit_lasso_bic <- function(
  x, y, penalty_weight = NULL, ...)
{
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  if (is.null(penalty_weight)) {
    penalty_weight <- rep(1.0, ncol(x))
  }
  if (sum(penalty_weight) < 1e-4 || ncol(x) == 1L) {
    return(fit_plain(x, y))
  }
  if (!is.matrix(y) && length(unique(y)) > 2L) {
    y <- cbind(1.0 - y, y)
  }
  penalty_weight <- pmin(pmax(penalty_weight, 1e-4), 1e4)
  model <- glmnet(
    x, y, family = "binomial", thresh = 5e-7,
    penalty.factor = penalty_weight, 
    nlambda = 50L, lambda.min.ratio = 5e-3)
  n <- length(y)
  dev <- deviance(model)
  nonzero <- sapply(predict(model, type = "nonzero"), length) + 1L
  complexity <- nonzero * log(n)
  lambda <- model$lambda[which.min(dev + complexity)]
  as.double(predict(model, s = lambda, type = "coefficients"))
}


predict_linear <- function(beta, x, ...)
{
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  plogis(beta[1L] + drop(x %*% beta[-1L]))
}


fit_svm <- function(x, y, ...)
{
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


fit_rf <- function(x, y, ...)
{
  if (requireNamespace("randomForest", quietly = TRUE)) {
    return(randomForest::randomForest(x, factor(y, c(0, 1))))
  } else {
    stop("Package randomForest not found")
  }
}


predict_rf <- function(beta, x, ...)
{
  if (requireNamespace("randomForest", quietly = TRUE)) {
    return(as.numeric(predict(beta, x, type = "prob")[, "1"]))
  } else {
    stop("Package randomForest not found")
  }
}


## roc and auc ----

get_roc <- function(
  y_true, y_score, 
  cut = seq(0.001, 0.999, 0.001),
  tol = sqrt(.Machine$double.eps))
{
  i <- order(y_score, decreasing = TRUE)
  y_score <- y_score[i]
  y_true <- y_true[i]
  
  n <- length(y_true)
  t1 <- sum(y_true)
  t0 <- n - t1
  p1 <- seq_len(n)
  p0 <- n - p1
  
  true_positive <- cumsum(y_true)
  false_positive <- p1 - true_positive
  true_negative <- t0 - false_positive
  
  thr <- y_score
  pct <- p1 / n
  acc <- (true_positive + true_negative) / n
  
  tpr <- true_positive / t1
  fpr <- false_positive / t0
  tnr <- true_negative / t0
  
  ppv <- true_positive / p1
  fdr <- false_positive / p1
  npv <- true_negative / p0
  
  sen <- tpr
  rec <- tpr
  spec <- tnr
  prec <- ppv
  f1 <- prec * rec / (prec + rec)
  
  i <- which(diff(thr) < -tol)
  df <- data.frame(
    thr = thr[i], 
    pct = pct[i], 
    acc = acc[i],
    tpr = tpr[i], 
    fpr = fpr[i], 
    tnr = tnr[i],
    ppv = ppv[i], 
    fdr = fdr[i], 
    npv = npv[i],
    sen = sen[i], 
    spec = spec[i], 
    prec = prec[i], 
    rec = rec[i],
    f1 = f1[i])
  df <- df[order(df$thr), ]
  
  f <- function(x, y) {
    approxfun(x, y, method = "constant", rule = 2L)(cut)
  }
  df <- data.frame(
    cut = cut,
    pct = f(df$thr, df$pct),
    acc = f(df$thr, df$acc),
    tpr = f(df$thr, df$tor),
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


get_auc <- function(y_true, y_score) 
{
  t1 <- sum(y_true)
  t0 <- length(y_true) - t1
  r <- rank(y_score, ties.method = "average")
  auc <- (sum(r[y_true > 0.5]) - t1 * (t1 + 1) / 2) / t1 / t0
  
  return(auc)
}
