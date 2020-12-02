#' Get available Functional Set Types
#'
#' Provides a list of available functional set types such as pathways, GO categories, NESTs, ...
#' Parameter 'entityType' should be either 'genes' or 'metabolites'.
#'
#' @param entityType, string with entity typme, should be 'genes' or metabolites
#'
#' @return Dataframe with Functional Set types ID and descriptions associated to genes or metabolites
#'
#' @examples
#'
#' library(PathData)
#'
#' getCpdbAvailableFsetTypes("genes")
#' @export
getCpdbAvailableFsetTypes <- function(entityType) {
  if (!entityType %in% c("genes", "metabolites")) {
    stop("entityType should be 'genes' or 'metabolites'")
  }

  body <- paste0('
   <soapenv:Envelope xmlns:soapenv="http://schemas.xmlsoap.org/soap/envelope/" xmlns:cpd="cpdbns"> \
      <soapenv:Header/> \
      <soapenv:Body> \
      <cpd:getAvailableFsetTypes> \
         <cpd:entityType>', entityType, "</cpd:entityType> \
      </cpd:getAvailableFsetTypes> \
      </soapenv:Body> \
   </soapenv:Envelope>")

  reader <- basicTextGatherer()

  curlPerform(
    url = "http://cpdb.molgen.mpg.de/ws2",
    httpheader = c(
      Accept = "text/xml", Accept = "multipart/*",
      SOAPAction = "urn:getAvailableFsetTypes",
      "Content-Type" = "text/xml; charset=utf-8"
    ),
    postfields = body,
    writefunction = reader$update
  )

  data <- xml2::read_xml(reader$value())

  ids <- xml_text(xml2::xml_find_all(data, "//ns1:fsetType"))
  desc <- xml_text(xml2::xml_find_all(data, "//ns1:description"))

  results <- bind_cols("ID" = ids, "description" = desc)

  return(results)
}
