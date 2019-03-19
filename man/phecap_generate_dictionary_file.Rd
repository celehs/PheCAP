\name{phecap_generate_dictionary_file}
\alias{phecap_generate_dictionary_file}
\title{
Generate a Dictionary File for Note Parsing
}
\description{
Given a list of CUIs,
connect to the UMLS database stored in MySQL,
extract CUIs and associated terms,
and write a dictionary file for use in note parsing.
}
\usage{
phecap_generate_dictionary_file(
  cui_list, dict_file,
  user = "username", password = "password",
  host = "localhost", dbname = "umls", ...)
}
\arguments{
  \item{cui_list}{
a character vector consisting of CUIs of interest.
}
  \item{dict_file}{
a character scalar for the path to the dictionary file that will be generated.
}
  \item{user}{
a character scalar for the username for database connection; passed to \code{RMySQL::dbConnect} as it is.
}
  \item{password}{
a character scalar for the password for database connection; passed to \code{RMySQL::dbConnect} as it is.
}
  \item{host}{
a character scalar for the host (or URL) for database connection; passed to \code{RMySQL::dbConnect} as it is.
}
  \item{dbname}{
a character scalar for the database name for database connection; passed to \code{RMySQL::dbConnect} as it is.
}
  \item{...}{
Other arguments passed to \code{RMySQL::dbConnect} as they are.
}
}
\value{
The dictionary will be written to the location given by \code{dict_file}.
Return the dictionary invisibly.
}
