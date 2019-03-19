\name{phecap_predict_phenotype}
\alias{phecap_predict_phenotype}
\title{
Predict Phenotype
}
\description{
Compute predicted probability of having the phenotype
for each patient in the dataset.
}
\usage{
phecap_predict_phenotype(data, model)
}
\arguments{
  \item{data}{
an object of class \code{PhecapData}, obtained by calling \code{PhecapData(...)}.
}
  \item{model}{
an object of class \code{PhecapModel}, probably returned from
\code{phecap_train_phenotyping_model}.
}
}
\value{
A \code{data.frame} with two columns:
  \item{patient_index}{patient identifier},
  \item{prediction}{predicted phenotype}.
}
\seealso{
See \code{\link{PheCAP-package}} for code examples.
}
