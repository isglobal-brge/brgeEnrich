#' Get over-representation analysis
#'
#' Performs over-representation analysis of functional sets with provided physical entities.
#'
#' @param entityType, should be either 'genes' or 'metabolites'.
#' @param fsetType is the type of the functional sets to be tested (as obtained with the getCpdbAvailableFsetTypes function; e.g. 'P' for pathways).
#' @param accNumbers is a list of interesting accNumbers (e.g. differentially expressed)
#' @param accType is a valid accession number type.
#' @param cpdbIdsBg is a list of CPDB entity IDs in the background. If empty, the default background is used (all different entities present in at least one functional set of the type fsetType and identifiable with accession numbers of type 'accType').
#' @param pThreshold is a p-value threshold, only sets with significant over-representation below or equal to this threshold will be provided.
#'
#' @return
#'
#' @examples
#'
#' library(PathData)
#'
#' accType <- "uniprot"
#' accNumbers <- c("MDHM_HUMAN", "MDHC_HUMAN", "DLDH_HUMAN", "DHSA_HUMAN", "DHSB_HUMAN", "C560_HUMAN", "DHSD_HUMAN", "ODO2_HUMAN", "ODO1_HUMAN", "CISY_HUMAN", "ACON_HUMAN", "IDH3A_HUMAN", "IDH3B_HUMAN", "IDH3G_HUMAN", "SUCA_HUMAN", "SUCB1_HUMAN", "FUMH_HUMAN", "OGDHL_HUMAN", "ACOC_HUMAN", "DHTK1_HUMAN", "AMAC1_HUMAN")
#'
#' overrep <- cpdbOverRepresentationAnalysis("genes", "C", accNumbers, accType)
#' @export
cpdbOverRepresentationAnalysis <- function(entityType, fsetType, accNumbers, accType, cpdbIdsBg = NULL, pThreshold = 0.05) {
  if (!entityType %in% c("genes", "metabolites")) {
    stop("entityType should be 'genes' or 'metabolites'")
  }

  if (!fsetType %in% getCpdbAvailableFsetTypes(entityType)$ID) {
    stop("fsetType should be in Available FsetTypes for the current entityType")
  }

  if (is.null(cpdbIdsBg) & is.null(accType)) {
    stop("valid accession number type (accType). Should be specified if parameter 'cpdbIdsBg' is not set")
  }

  cpdbIds <- mapCpdbAccessionNumbers(accType, accNumbers)
  cpdbIds <- cpdbIds %>% filter(cpdbId != "")

  cpdbIdsFg <- apply(cpdbIds[, c(2)], 1, function(dat) {
    paste0(
      "<cpd:cpdbIdsFg>",
      paste0(dat, collapse = " "),
      "</cpd:cpdbIdsFg>"
    )
  })


  body <- paste0('
   <soapenv:Envelope xmlns:soapenv="http://schemas.xmlsoap.org/soap/envelope/" xmlns:cpd="cpdbns">
      <soapenv:Header/>
      <soapenv:Body>
         <cpd:overRepresentationAnalysis>
            <cpd:entityType>', entityType, "</cpd:entityType>
            <cpd:fsetType>", fsetType, "</cpd:fsetType>
            ", paste0(cpdbIdsFg, collapse = " "), collapse = NULL)

  if (!is.null(cpdbIdsBg)) {
    cpdbIdsBg <- apply(cpdbIds[, c(2)], 1, function(dat) {
      paste0(
        "<cpd:cpdbIdsBg>",
        paste0(dat, collapse = " "),
        "</cpd:cpdbIdsBg>"
      )
    })

    body <- paste0(body, cpdbIdsBg, collapse = NULL)
  }

  body <- paste0(body, "
            <cpd:accType>", accType, "</cpd:accType>
            <cpd:pThreshold>", pThreshold, "</cpd:pThreshold>
         </cpd:overRepresentationAnalysis>
      </soapenv:Body>
   </soapenv:Envelope>", collapse = NULL)

  reader <- basicTextGatherer()

  curlPerform(
    url = "http://cpdb.molgen.mpg.de/ws2",
    httpheader = c(
      Accept = "text/xml", Accept = "multipart/*",
      SOAPAction = "urn:overRepresentationAnalysis",
      "Content-Type" = "text/xml; charset=utf-8"
    ),
    postfields = body,
    writefunction = reader$update
  )

  data <- xml2::read_xml(reader$value())

  name <- xml_text(xml2::xml_find_all(data, "//ns1:name"))

  details <- xml_text(xml2::xml_find_all(data, "//ns1:details"))
  if(length(details)>0)
  {
    splited <- str_split(details, ";;")
    # splited2 <- lapply(splited, insert_pmids)

    if(fsetType=='P'){
      details <- as_tibble(data.frame(matrix(unlist(lapply(splited, insert_pmids)), nrow = length(str_split(details, ";;")), byrow = T)))
    }else if( grep("G?",fsetType)==1) {
      # details <- as_tibble(data.frame(matrix(unlist(lapply(splited, insert_pmids)), nrow = length(str_split(details, ";;")), byrow = T)))
      details <- as_tibble(data.frame(matrix(unlist(splited), nrow = length(str_split(details, ";;")), byrow = T)))
    }

    ## TODO :
    ##    Replace this code for an optimal one with a function that controls all fields in only on step

    details <- details %>% separate(X1, c("X1k", "X1v"), ":", extra = "merge") %>%
      separate(X2, c("X2k", "X2v"), ":", extra = "merge")
    xk1 <- unique(details$X1k)
    xk2 <- unique(details$X2k)
    details <- details %>%
      mutate(
        !!as.name((xk1)) := X1v,
        !!as.name((xk2)) := X2v,
        X1k = NULL,
        X1v = NULL,
        X2k = NULL,
        X2v = NULL
      )

    if(fsetType=='G2' | fsetType=='G3' | fsetType=='G4' |fsetType=='G5' | fsetType=='P')
    {
      details <- details %>% separate(X3, c("X3k", "X3v"), ":", extra = "merge")

      xk3 <- unique(details$X3k)
      details <- details %>%
        mutate(
          !!as.name((xk3)) := X3v,
          X3k = NULL,
          X3v = NULL,
        )

      if(fsetType=='P')
      {
        details <- details %>% separate(X4, c("X4k", "X4v"), ":", extra = "merge") %>%
          separate(X5, c("X5k", "X5v"), ":", extra = "merge")

        xk4 <- unique(details$X4k)
        xk5 <- unique(details$X5k)
        details <- details %>%
          mutate(
            !!as.name((xk4)) := X4v,
            !!as.name((xk5)) := X5v,
            X4k = NULL,
            X4v = NULL,
            X5k = NULL,
            X5v = NULL,
          )
      }
    }

    overlappingEntitiesNum <- xml_text(xml2::xml_find_all(data, "//ns1:overlappingEntitiesNum"))
    allEntitiesNum <- xml_text(xml2::xml_find_all(data, "//ns1:allEntitiesNum"))
    pValue <- xml_text(xml2::xml_find_all(data, "//ns1:pValue"))
    qValue <- xml_text(xml2::xml_find_all(data, "//ns1:qValue"))

    results <- bind_cols(
      "name" = name,
      details,
      "overlappingEntitiesNum" = overlappingEntitiesNum,
      "allEntitiesNum" = allEntitiesNum,
      "pValue" = as.numeric(pValue),
      "qValue" = as.numeric(qValue)
    )
  }else {
    results <- data.frame(name = character(),
                          fsetId = character(),
                          CPDBurl = character(),
                          overlappingEntitiesNum = character(),
                          allEntitiesNum = character(),
                          pValue = numeric(),
                          qValue = numeric(),
                          stringsAsFactors=FALSE)
  }

  return(results)
}



# Inserts element not availabled
insert_pmids <- function(x)
{
  val <- x
  if(length(grep('pmids',x[4]))==0) {
    val <- insert(x, ats=4,values = c('pmids:<NA>'))
  }
  return(val)
}
