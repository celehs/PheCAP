\name{phecap_plot_roc_curves}
\alias{phecap_plot_roc_curves}
\title{
Plot ROC and Related Curves for Phenotyping Models
}
\description{
Plot ROC-like curves to illustrate phenotyping accuracy.
}
\usage{
phecap_plot_roc_curves(
  x, axis_x = "1 - spec", axis_y = "sen",
  what = c("training", "random-splits", "validation"),
  ggplot = TRUE, ...)
}
\arguments{
  \item{x}{
either a single object of class PhecapModel or PhecapValidation
(returned from \code{phecap_train_phenotyping_model} or
\code{phecap_validate_phenotyping_model}), or a named list of such objects
}
  \item{axis_x}{
an expression that leads to the \code{x} coordinate.
Recognized quantities include:
\code{cut} (probability cutoff),
\code{pct} (percent of predicted cases),
\code{acc} (accuracy),
\code{tpr} (true positive rate),
\code{fpr} (false positive rate),
\code{tnr} (true negative rate),
\code{ppv} (positive predictive value),
\code{fdr} (false discovery rate),
\code{npv} (negative predictive value),
\code{sen} (sensitivity),
\code{spec} (specificity),
\code{prec} (precision),
\code{rec} (recall),
\code{f1} (F1 score).
}
  \item{axis_y}{
an expression that leads to the \code{y} coordinate.
Recognized quantities are the same as those in \code{axis_x}.
}
  \item{what}{
The curves to be included in the figure.
}
  \item{ggplot}{
if TRUE and ggplot2 is installed, ggplot will be used for the figure.
Otherwise, the base R graphics functions will be used.
}
  \item{\dots}{
arguments to be ignored.
}
}
\seealso{
See \code{\link{PheCAP-package}} for code examples.
}
