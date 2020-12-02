# Get CpdbVersion - Provides the version of ConsensusPathDB which is the source for functional sets.
#
#' Get CpdbVersion
#'
#' Provides the version of ConsensusPathDB which is the source for functional sets.
#'
#' @return string with ConsensusPathDB version
#'
#' @examples
#'
#' library(PathData)
#'
#' getCpdbVersion()
#' @export
getCpdbVersion <- function() {
  body <- '
   <soapenv:Envelope xmlns:soapenv="http://schemas.xmlsoap.org/soap/envelope/" xmlns:cpd="cpdbns"> \
      <soapenv:Header/> \
      <soapenv:Body> \
         <cpd:getCpdbVersion/> \
      </soapenv:Body> \
   </soapenv:Envelope>'

  reader <- basicTextGatherer()

  curlPerform(
    url = "http://cpdb.molgen.mpg.de/ws2",
    httpheader = c(
      Accept = "text/xml", Accept = "multipart/*",
      SOAPAction = "urn:getCpdbVersion",
      "Content-Type" = "text/xml; charset=utf-8"
    ),
    postfields = body,
    writefunction = reader$update
  )

  data <- xml2::read_xml(reader$value())
  results <- xml_text(xml2::xml_find_all(data, "//ns1:cpdbVersion"))

  return(results)
}
