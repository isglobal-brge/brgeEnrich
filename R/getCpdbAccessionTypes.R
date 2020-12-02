#' Get CPDB Available Accession Types
#'
#' Provides a list of different types of accession numbers (e.g. 'uniprot') that are mappable to CPDB entity IDs.
#' Parameter 'entityType' should be either 'genes' or 'metabolites'.
#'
#' @param entityType, string with entity typme, should be 'genes' or metabolites
#'
#' @return array with accession types associated to genes or metabolites
#'
#' @examples
#'
#' library(PathData)
#'
#' getCpdbAccessionTypes(entityType = "genes")
#' @export
getCpdbAccessionTypes <- function(entityType) {
  if (!entityType %in% c("genes", "metabolites")) {
    stop("entityType should be 'genes' or 'metabolites'")
  }

  body <- paste0('
   <soapenv:Envelope xmlns:soapenv="http://schemas.xmlsoap.org/soap/envelope/" xmlns:cpd="cpdbns"> \
      <soapenv:Header/> \
      <soapenv:Body> \
         <cpd:getAvailableAccessionTypes> \
            <cpd:entityType>', entityType, "</cpd:entityType> \
         </cpd:getAvailableAccessionTypes> \
      </soapenv:Body> \
   </soapenv:Envelope>")

  reader <- basicTextGatherer()

  curlPerform(
    url = "http://cpdb.molgen.mpg.de/ws2",
    httpheader = c(
      Accept = "text/xml", Accept = "multipart/*",
      SOAPAction = "urn:getAvailableAccessionTypes",
      "Content-Type" = "text/xml; charset=utf-8"
    ),
    postfields = body,
    writefunction = reader$update
  )

  data <- xml2::read_xml(reader$value())

  results <- xml_text(xml2::xml_find_all(data, "//ns1:accType"))

  return(results)
}
