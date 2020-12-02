#' Map accession numbers
#'
#' Maps accession numbers of a valid type to CPDB entity IDs.
#'
#' @param accType, string with a valid accession number type (such as 'uniprot').
#' @param accNumbers, string with a list of accession numbers.
#'
#' @return cpdbIds, mapped accession numbers
#'
#' @examples
#'
#' library(PathData)
#'
#' accNumbers <- c("MDHM_HUMAN", "MDHC_HUMAN", "DLDH_HUMAN")
#'
#' mapAccessionNumbers("uniprot", accNumbers)
#' @export
mapAccessionNumbers <- function(accType, accNumbers) {
  accNumbersQ <- sapply(accNumbers, function(dat) {
    paste0("<cpd:accNumbers>", dat, "</cpd:accNumbers>")
  })


  body <- paste0(
    '
   <soapenv:Envelope xmlns:soapenv="http://schemas.xmlsoap.org/soap/envelope/" xmlns:cpd="cpdbns">
      <soapenv:Header/>
      <soapenv:Body>
         <cpd:mapAccessionNumbers>
            <cpd:accType>', accType, "</cpd:accType>",
    paste0(accNumbersQ, collapse = "\n"), "
         </cpd:mapAccessionNumbers>
      </soapenv:Body>
   </soapenv:Envelope>"
  )

  reader <- basicTextGatherer()

  curlPerform(
    url = "http://cpdb.molgen.mpg.de/ws2",
    httpheader = c(
      Accept = "text/xml", Accept = "multipart/*",
      SOAPAction = "urn:mapAccessionNumbers",
      "Content-Type" = "text/xml; charset=utf-8"
    ),
    postfields = body,
    writefunction = reader$update
  )

  data <- xml2::read_xml(reader$value())

  accNumberR <- xml_text(xml2::xml_find_all(data, "//ns1:accNumber"))
  cpdbId <- xml_text(xml2::xml_find_all(data, "//ns1:cpdbId"))

  results <- bind_cols("accNumber" = accNumberR, "cpdbId" = cpdbId)

  return(results)
}
