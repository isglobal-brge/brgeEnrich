#' Get all CPDB entity IDs in a functional set.
#'
#' Returns all CPDB entity IDs in a functional set. Note that functional sets of type 'N' (NESTs) are protected and cannot be retrieved.
#'
#' @param fsetId, is the ID of the functional set
#' @param fsetType, is its type (e.g. 'P').
#' @param entityType, entsetType
#'
#' @return array with CPDB entity IDs
#'
#' @examples
#'
#' library(PathData)
#'
#' getCpdbVersion()
#' @export
getCpdbIdsInFset <- function(fsetId, fsetType, entityType) {

  if (!entityType %in% c("genes", "metabolites")) {
    stop("entityType should be 'genes' or 'metabolites'")
  }

  body <- paste0('
   <soapenv:Envelope xmlns:soapenv="http://schemas.xmlsoap.org/soap/envelope/" xmlns:cpd="cpdbns">
      <soapenv:Header/>
      <soapenv:Body>
         <cpd:getCpdbIdsInFset>
            <cpd:fsetId>', fsetId, "</cpd:fsetId>
            <cpd:fsetType>", fsetType, "</cpd:fsetType>
            <cpd:entsetType>", entityType, "</cpd:entsetType>
         </cpd:getCpdbIdsInFset>
      </soapenv:Body>
   </soapenv:Envelope>")

  reader <- basicTextGatherer()

  curlPerform(
    url = "http://cpdb.molgen.mpg.de/ws2",
    httpheader = c(
      Accept = "text/xml", Accept = "multipart/*",
      SOAPAction = "urn:getCpdbIdsInFset",
      "Content-Type" = "text/xml; charset=utf-8"
    ),
    postfields = body,
    writefunction = reader$update
  )

  data <- xml2::read_xml(reader$value())
  results <- xml_text(xml2::xml_find_all(data, "//ns1:cpdbIds"))

  return(results)
}
