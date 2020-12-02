#' Get default Background Size
#'
#' Provides the default background size for over-representation analysis. Note that it is accession number-specific.
#'
#' @param fsetType, string with functional Set types ID associated to genes or metabolites
#' @param accType, string with accession type
#'
#' @return numeric with background size for over-representation analysis
#'
#' @examples
#'
#' library(PathData)
#'
#' getDefaultBackgroundSize("P", "hgnc-symbol")
#' @export
getDefaultBackgroundSize <- function(fsetType, accType) {
  if (!entityType %in% c("genes", "metabolites")) {
    stop("entityType should be 'genes' or 'metabolites'")
  }

  body <- paste0('
   <soapenv:Envelope xmlns:soapenv="http://schemas.xmlsoap.org/soap/envelope/" xmlns:cpd="cpdbns">
      <soapenv:Header/>
      <soapenv:Body>
         <cpd:getDefaultBackgroundSize>
            <cpd:fsetType>', fsetType, "</cpd:fsetType>
            <cpd:accType>", accType, "</cpd:accType>
         </cpd:getDefaultBackgroundSize>
      </soapenv:Body>
   </soapenv:Envelope>")

  reader <- basicTextGatherer()

  curlPerform(
    url = "http://cpdb.molgen.mpg.de/ws2",
    httpheader = c(
      Accept = "text/xml", Accept = "multipart/*",
      SOAPAction = "urn:getDefaultBackgroundSize",
      "Content-Type" = "text/xml; charset=utf-8"
    ),
    postfields = body,
    writefunction = reader$update
  )

  data <- xml2::read_xml(reader$value())
  result <- xml_text(xml2::xml_find_all(data, "//ns1:bgSize"))

  return(result)
}
