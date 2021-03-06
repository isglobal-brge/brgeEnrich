% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mapCpdbAcessionNumbers.R
\name{mapCpdbAccessionNumbers}
\alias{mapCpdbAccessionNumbers}
\title{Map accession numbers}
\usage{
mapCpdbAccessionNumbers(accType, accNumbers)
}
\arguments{
\item{accType, }{string with a valid accession number type (such as 'uniprot').}

\item{accNumbers, }{string with a list of accession numbers.}
}
\value{
cpdbIds, mapped accession numbers
}
\description{
Maps accession numbers of a valid type to CPDB entity IDs.
}
\examples{

library(PathData)

accNumbers <- c("MDHM_HUMAN", "MDHC_HUMAN", "DLDH_HUMAN")

mapCpdbAccessionNumbers("uniprot", accNumbers)
}
