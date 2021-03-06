% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cpdbOverRepresentationAnalysis.R
\name{cpdbOverRepresentationAnalysis}
\alias{cpdbOverRepresentationAnalysis}
\title{Get over-representation analysis}
\usage{
cpdbOverRepresentationAnalysis(
  entityType,
  fsetType,
  accNumbers,
  accType,
  cpdbIdsBg = NULL,
  pThreshold = 0.05
)
}
\arguments{
\item{entityType, }{should be either 'genes' or 'metabolites'.}

\item{fsetType}{is the type of the functional sets to be tested (as obtained with the getCpdbAvailableFsetTypes function; e.g. 'P' for pathways).}

\item{accNumbers}{is a list of interesting accNumbers (e.g. differentially expressed)}

\item{accType}{is a valid accession number type.}

\item{cpdbIdsBg}{is a list of CPDB entity IDs in the background. If empty, the default background is used (all different entities present in at least one functional set of the type fsetType and identifiable with accession numbers of type 'accType').}

\item{pThreshold}{is a p-value threshold, only sets with significant over-representation below or equal to this threshold will be provided.}
}
\value{

}
\description{
Performs over-representation analysis of functional sets with provided physical entities.
}
\examples{

library(PathData)

accType <- "uniprot"
accNumbers <- c("MDHM_HUMAN", "MDHC_HUMAN", "DLDH_HUMAN", "DHSA_HUMAN", "DHSB_HUMAN", "C560_HUMAN", "DHSD_HUMAN", "ODO2_HUMAN", "ODO1_HUMAN", "CISY_HUMAN", "ACON_HUMAN", "IDH3A_HUMAN", "IDH3B_HUMAN", "IDH3G_HUMAN", "SUCA_HUMAN", "SUCB1_HUMAN", "FUMH_HUMAN", "OGDHL_HUMAN", "ACOC_HUMAN", "DHTK1_HUMAN", "AMAC1_HUMAN")

overrep <- cpdbOverRepresentationAnalysis("genes", "C", accNumbers, accType)
}
