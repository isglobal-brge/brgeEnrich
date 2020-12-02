#' CPDB Enrichment Analysis
#'
#' Performs Wilcoxon enrichment analysis of functional sets with provided physical entities from CPDB.
#'
#' @param entityType, should be either 'genes' or 'metabolites'.
#' @param accNumbers is a list of interesting accNumbers (e.g. differentially expressed)
#' @param mv1 is a list of mesured value 1 or fold change if m2 is null
#' @param mv2 is a list of mesured value 2
#' @param accType is a valid accession number type.
#' @param fsetType, the type of the functional sets to be tested (as obtained with the getCpdbAvailableFsetTypes function; e.g. 'P' for pathways).
#' @param pThreshold, a p-value threshold, only sets with significant enrichment below or equal to this threshold will be provided.
#'
#' @return cpdbIds, mapped accession numbers
#'
#' @examples
#'
#' library(PathData)
#'
#' accType <- "uniprot"
#' accNumbers <- c("MDHM_HUMAN", "MDHC_HUMAN", "DLDH_HUMAN", "DHSA_HUMAN", "DHSB_HUMAN", "C560_HUMAN", "DHSD_HUMAN", "ODO2_HUMAN", "ODO1_HUMAN", "CISY_HUMAN", "ACON_HUMAN", "IDH3A_HUMAN", "IDH3B_HUMAN", "IDH3G_HUMAN", "SUCA_HUMAN", "SUCB1_HUMAN", "FUMH_HUMAN", "OGDHL_HUMAN", "ACOC_HUMAN", "DHTK1_HUMAN", "AMAC1_HUMAN")
#'
#' mv1 <- runif(length(accNumbers))
#' mv2 <- runif(length(accNumbers))
#'
#' eanalysis <- cpdbEnrichmentAnalysis("genes", accNumbers, mv1, mv2, "uniprot", "C", 1)
#' @export
cpdbEnrichmentAnalysis <- function(entityType, accNumbers, mv1, mv2 = NULL, accType, fsetType, pThreshold = 1) {
  if (!entityType %in% c("genes", "metabolites")) {
    stop("entityType should be 'genes' or 'metabolites'")
  }

  if (!fsetType %in% getCpdbAvailableFsetTypes(entityType)$ID) {
    stop("fsetType should be in Available FsetTypes for the current entityType")
  }

  if (!accType %in% getCpdbAccessionTypes(entityType)) {
    stop("fsetType should be in Available FsetTypes for the current entityType")
  }

  # Map accession numbers
  cpdbIds <- mapCpdbAccessionNumbers(accType, accNumbers)
  # Add values to accession numbers (values or fold change)
  if (!is.null(mv2)) {
    cpdbIds <- cpdbIds %>% mutate(
      m1 = mv1,
      m2 = mv2
    )
  } else {
    cpdbIds <- cpdbIds %>% mutate(
      m1 = mv1
    )
  }

  # Clean null accession numbers
  cpdbIds <- cpdbIds %>% filter(cpdbId != "")


  # Create data attending to xml description
  cpdbIds <- apply(cpdbIds[, 2:dim(cpdbIds)[2]], 1, function(dat) {
    paste0(
      "<cpd:cpdbIdsMeasurements>",
      paste0(dat, collapse = " "),
      "</cpd:cpdbIdsMeasurements>"
    )
  })


  body <- paste0('
   <soapenv:Envelope xmlns:soapenv="http://schemas.xmlsoap.org/soap/envelope/" xmlns:cpd="cpdbns">
   <soapenv:Header/>
   <soapenv:Body>
      <cpd:enrichmentAnalysis>
         <cpd:entityType>', entityType, "</cpd:entityType>
         <cpd:fsetType>", fsetType, "</cpd:fsetType>
            ", paste0(cpdbIds, collapse = " "), "
            <cpd:pThreshold>", pThreshold, "</cpd:pThreshold>
         </cpd:enrichmentAnalysis>
      </soapenv:Body>
   </soapenv:Envelope>", collapse = NULL)

  reader <- basicTextGatherer()

  curlPerform(
    url = "http://cpdb.molgen.mpg.de/ws2",
    httpheader = c(
      Accept = "text/xml", Accept = "multipart/*",
      SOAPAction = "urn:enrichmentAnalysis",
      "Content-Type" = "text/xml; charset=utf-8"
    ),
    postfields = body,
    writefunction = reader$update
  )

  data <- xml2::read_xml(reader$value())

  name <- xml_text(xml2::xml_find_all(data, "//ns1:name"))

  details <- xml_text(xml2::xml_find_all(data, "//ns1:details"))
  details <- as_tibble(data.frame(matrix(unlist(str_split(details, ";")), nrow = length(str_split(details, ";")), byrow = T)))
  details <- details %>% mutate(
    fsetId = str_remove(details$X1, "^fsetId:"),
    CPDBurl = str_remove(details$X3, "^CPDBurl:"),
    X1 = NULL,
    X2 = NULL,
    X3 = NULL
  )

  measuredEntitiesNum <- xml_text(xml2::xml_find_all(data, "//ns1:measuredEntitiesNum"))
  allEntitiesNum <- xml_text(xml2::xml_find_all(data, "//ns1:allEntitiesNum"))
  pValue <- xml_text(xml2::xml_find_all(data, "//ns1:pValue"))
  qValue <- xml_text(xml2::xml_find_all(data, "//ns1:qValue"))

  results <- bind_cols(
    "name" = name,
    "fsetId" = details$fsetId,
    "CPDBurl" = details$CPDBurl,
    "measuredEntitiesNum" = as.numeric(measuredEntitiesNum),
    "allEntitiesNum" = as.numeric(allEntitiesNum),
    "pValue" = as.numeric(pValue),
    "qValue" = as.numeric(qValue)
  )

  return(results)
}
