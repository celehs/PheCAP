\name{PhecapData}
\alias{PhecapData}
\title{
Define or Read Datasets for Phenotyping
}
\description{
Specify the data to be used for phenotyping.
}
\usage{
PhecapData(
  data, hu_feature, label, validation,
  patient_id = NULL, subject_weight = NULL,
  seed = 12300L, feature_transformation = log1p)
}
\arguments{
  \item{data}{
A data.frame consisting of all the variables needed for phenotyping,
or a character scalar of the path to the data,
or a list consisting of either character scalar or data.frame.
If a list is given, patient_id cannot be NULL.
All the datasets in the list will be joined into a single dataset
according to the columns specified by patient_id.
}
  \item{hu_feature}{
A character scalar or vector specifying the names of
one of more healthcare utilization (HU) variables.
There variables are always included in the phenotyping model.
}
  \item{label}{
A character scalar of the column name that gives the phenotype status
(1 or TRUE: present, 0 or FALSE: absent).
If label is not ready yet, just put a column filled with NA in data.
In such cases only the feature extraction step can be done.
}
  \item{validation}{
A character scalar, a real number strictly between 0 and 1,
or an integer not less than 2.
If a character scalar is used, it is treated as the column name
in the data that specifies whether this observation
belongs to the validation samples
(1 or TRUE: validation, 0 or FALSE: training).
If a real number strictly between 0 and 1 is used, it is treated as
the proportion of the validation samples. The actual validation samples
will be drawn from all labeled samples.
If an integer not less than 2 is used, it is treated as
the size of the validation samples. The actual validation samples
will be drawn from all labeled samples.
}
  \item{patient_id}{
A character vector for the column names, if any, that uniquely identifies
each patient.
Such variables must appear in the data.
patient_id can be NULL if such fields are not contained in the data.
}
  \item{subject_weight}{
An optional numeric vector of weights for observations.
}
  \item{seed}{
If validation samples need to be drawn from all labeled samples,
seed specifies the random seed for sampling.
}
  \item{feature_transformation}{
A function that will be applied to all the features.
Since count data are typically right-skewed,
by default \code{log1p} will be used.
feature_transformation can be NULL, in which case
no transformation will be done on any of the feature.
}
}
\value{
An object of class \code{PhecapData}.
}
\seealso{
See \code{\link{PheCAP-package}} for code examples.
}
