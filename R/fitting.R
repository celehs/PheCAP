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
ROC.Est.FUN=function(Di,yyi,yy0,fpr0=NULL,wgti=NULL,yes.smooth=F)
{
  out.yy <- out.pp <- out.AUC <- out.TPR <- out.FPR <- out.PPV <- out.NPV <- out.F1 <- NULL
  if(is.null(wgti)){wgti=rep(1,length(Di))}; yyi = as.matrix(yyi); pp=ncol(as.matrix(yyi));  
  mu0 = sum(wgti*(1-Di))/sum(wgti); mu1 = 1-mu0  
  for(k in 1:pp)
  {
    yy = yy0; 
    if(!is.null(fpr0)){
      tpr.all = S.FUN(yyi[,k],Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth); 
      fpr.all = S.FUN(yyi[,k],Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth);
      TPR = approx(c(0,fpr.all,1),c(0,tpr.all,1),fpr0,method="linear",rule=2)$y; 
      TPR = c(S.FUN(yy0,Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth), TPR); 
      yy = c(yy,Sinv.FUN(fpr0,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth))           
      FPR = S.FUN(yy,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth)
    }else{
      TPR = S.FUN(yy,Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth); 
      FPR = S.FUN(yy,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth)
    }
    out.yy = cbind(out.yy, yy)
    out.pp = cbind(out.pp, S.FUN(yy,Yi=yyi[,k],wgti,yes.smooth=yes.smooth))
    out.TPR = cbind(out.TPR,  TPR);  out.FPR  <- cbind(out.FPR,  FPR)
    PPV <- 1/(1+FPR*mu0/(TPR*mu1)); NPV <- 1/(1+(1-TPR)*mu1/((1-FPR)*mu0))
    out.PPV <- cbind(out.PPV, PPV); out.NPV <- cbind(out.NPV, NPV)
    AUC = sum(S.FUN(yyi[,k],Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth)*(1-Di)*wgti)/sum((1-Di)*wgti)
    F1=2*TPR*PPV/(TPR+PPV)
    F1[is.na(F1)==1]=0
    out.AUC <- c(out.AUC, AUC)
    out.F1=c(out.F1, F1)
  }
  out = c(out.AUC,out.yy,out.pp,out.FPR,out.TPR,out.PPV,out.NPV,out.F1)
  out
}
get_roc<- function(y_true, y_score, subject_weight)
{
  junk=ROC.Est.FUN(y_true,y_score,yy0=0.5,fpr0=seq(0,1,0.01),wgti=subject_weight,yes.smooth=F)[-1]
  df=matrix(junk, ncol=7)[-1,]
  colnames(df)=c("cutoff","pos.rate","FPR","TPR","PPV","NPV","F1")
  return(df)
}
get_auc<- function(y_true, y_score, subject_weight)
{
  ROC.Est.FUN(y_true,y_score,yy0=0.5,fpr0=seq(0,1,0.01),wgti=subject_weight,yes.smooth=F)[1]
}

S.FUN=function(yy,Yi,Di,yes.smooth=F)
{
  if(yes.smooth){
    Y1i = Yi[Di==1]; n1 = sum(Di); bw = bw.nrd(Y1i)/n1^0.6
    c(t(rep(1/n1,n1))%*%pnorm((Y1i-VTM(yy,n1))/bw))
  }else{
    return((sum.I(yy,"<",Yi,Vi=Di)+sum.I(yy,"<=",Yi,Vi=Di))/sum(Di)/2)
  }
  ##sum.I(yy,"<=",Yi,Vi=Di)/sum(Di)
}

sum.I <- function(yy,FUN,Yi,Vi=NULL)
  ## sum_i I(yy FUN Yi)Vi
  # Vi weight
{
  if (FUN=="<"|FUN==">=") { yy <- -yy; Yi <- -Yi}
  # for each distinct ordered failure time t[j], number of Xi < t[j]
  pos <- rank(c(yy,Yi),ties.method='f')[1:length(yy)]-rank(yy,ties.method='f')    
  if (substring(FUN,2,2)=="=") pos <- length(Yi)-pos # number of Xi>= t[j]
  if (!is.null(Vi)) {
    ## if FUN contains '=', tmpind is the order of decending
    if(substring(FUN,2,2)=="=") tmpind <- order(-Yi) else  tmpind <- order(Yi)
    ##Vi <- cumsum2(as.matrix(Vi)[tmpind,])
    Vi <- apply(as.matrix(Vi)[tmpind,,drop=F],2,cumsum)
    return(rbind(0,Vi)[pos+1,])
  } else return(pos)
}
Sinv.FUN <- function(uu,Yi,Di,yes.smooth=F)
{
  yy0<-unique(sort(Yi,decreasing=T)); ss0 <- S.FUN(yy0,Yi,Di,yes.smooth=yes.smooth) 
  return(approx(ss0[!duplicated(ss0)],yy0[!duplicated(ss0)],uu,method="linear",f=0,rule=2)$y)
}

VTM<-function(vc, dm){
  matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
}

