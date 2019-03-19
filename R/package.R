#' High-Throughput Phenotyping with Electronic Health Records 
#' using a Common Automated Pipeline
#' 
#' Implement surrogate-assisted feature extraction (SAFE) and common 
#' machine learning approaches to train and validate phenotyping models.
#' Background and details about the methods can be found at 
#' <doi:10.1093/jamia/ocw135> and <doi:10.1136/bmj.h1885>.
#' 
#' PheCAP provides a straightforward interface for conducting
#' phenotyping on eletronic health records. One can specify the
#' data via \code{\link{PhecapData}}, define surrogate using
#' \code{\link{PhecapSurrogate}}. Next, one may run
#' surrogate-assisted feature extraction (SAFE) by calling
#' \code{\link{phecap_run_feature_extraction}}, and then
#' train and validate phenotyping models via
#' \code{\link{phecap_train_phenotyping_model}} and
#' \code{\link{phecap_validate_phenotyping_model}}.
#' The predictive performance can be visualized using
#' \code{\link{phecap_plot_roc_curves}}.
#' Predicted phenotype is provided by
#' \code{\link{phecap_predict_phenotype}}.
#' 
#' @name PheCAP-package
#' 
#' @aliases PheCAP-package
#' @aliases PheCAP
#' 
#' @docType package
#' 
#' @keywords package
NULL


#' @importFrom glmnet glmnet
#' @importFrom glmnet cv.glmnet
#' @importFrom graphics legend
#' @importFrom graphics lines
#' @importFrom graphics plot
#' @importFrom methods is
#' @importFrom stats aggregate
#' @importFrom stats approxfun
#' @importFrom stats binomial
#' @importFrom stats coef
#' @importFrom stats deviance
#' @importFrom stats glm
#' @importFrom stats plogis
#' @importFrom stats predict
#' @importFrom stats quasibinomial
#' @importFrom utils read.csv
#' @importFrom utils read.table
NULL

