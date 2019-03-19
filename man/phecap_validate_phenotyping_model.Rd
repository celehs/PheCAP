\name{phecap_validate_phenotyping_model}
\alias{phecap_validate_phenotyping_model}
\title{
Validate the Phenotyping Model using the Validation Labels
}
\description{
Apply the trained model to all patients in the validation dataset,
and measure the prediction accuracy via ROC and AUC.
}
\usage{
phecap_validate_phenotyping_model(data, model)
}
\arguments{
  \item{data}{
an object of class \code{PhecapData}, obtained by calling \code{PhecapData(...)
}}
  \item{model}{
an object of class \code{PhecapModel}, obtained by calling
\code{phecap_train_phenotyping_model}.
}
}
\value{
An object of class \code{PhecapValidation}, with components
  \item{method}{the method used for model training}
  \item{train_roc}{ROC on training dataset}
  \item{train_auc}{AUC on training dataset}
  \item{split_roc}{average ROC on random splits of training dataset}
  \item{split_auc}{average AUC on random splits of training dataset}
  \item{valid_roc}{ROC on validation dataset}
  \item{valid_auc}{AUC on validation dataset}
}
\seealso{
See \code{\link{PheCAP-package}} for code examples.
}
