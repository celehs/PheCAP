\name{phecap_perform_majority_voting}
\alias{phecap_perform_majority_voting}
\title{
Perform Majority Voting on the CUIs from Multiple Knowledge Sources
}
\description{
Read parsed knowledge sources and identify CUIs.
Generate a list of CUIs that appear in at least half of the sources.
}
\usage{
phecap_perform_majority_voting(
  input_folder)
}
\arguments{
  \item{input_folder}{
a character scalar for the path to the folder that contains the parsed knowledge sources
}
}
\value{
A character vector consisting of CUIs that
pass the majority voting criterion.
}
