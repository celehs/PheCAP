\name{PhecapSurrogate}
\alias{PhecapSurrogate}
\title{
Define a Surrogate Variable used in
Surrogate-Assisted Feature Extraction (SAFE)
}
\description{
Define a surrogate varible from existing features, and specify
associated lower and upper cutoffs.
}
\usage{
PhecapSurrogate(variable_names, lower_cutoff = 1L, upper_cutoff = 10L)
}
\arguments{
  \item{variable_names}{
a character scalar or vector consisting of variable names.
If a vector is given, the value of the surrogate is defined as the sum of
the values of each variable.
}
  \item{lower_cutoff}{
a numeric scalar. If the surrogate value of a patient is less than or equal to this cutoff,
then this patient is treated as a control in SAFE.
}
  \item{upper_cutoff}{
a numeric scalar. If the surrogate value of a patient is greater than or equal to this cutoff,
then this patient is treated as a case in SAFE.
}
}
\details{
This function only stores the definition. No calculation is done.
}
\value{
An object of class \code{PhecapSurrogate}.
}
\seealso{
See \code{\link{PheCAP-package}} for code examples.
}
