
createMzQC = function(rprt_fns)
{
  cv_obo_file = system.file("./mzQC/qc-cv.obo", package="PTXQC")
  
  ##save required data of each metric in JSON-format
  cat("Creating metric json file ...", '\n')
  
  analysisSoftware <- list(name = "PTXQC",
                           version = as.character(packageVersion("PTXQC")),
                           uri = "github.com/cbielow/PTXQC")
  
  mzQC_version <- "0.1.1"
  creationDate <- as.character(Sys.time(), format="%Y-%m-%dT%H:%M:%S" )
  cv_infos = list(json_cvInfo("Proteomics Standards Initiative Quality Control Ontology",
                           "https://github.com/HUPO-PSI/qcML-development/blob/master/cv/v0_1_0/qc-cv.obo",
                           "0.1.0"),
               json_cvInfo("Proteomics Standards Initiative Mass Spectrometry Ontology",
                           "https://github.com/HUPO-PSI/psi-ms-CV/blob/master/psi-ms.obo",
                           "4.1.7"))
  
  jsonlite::prettify(jsonlite::toJSON(cv_infos, auto_unbox = TRUE))
  
  ## create output file
  createmzQC(version = mzQC_version, creationDate = creationDate, runQualities = runQuality, setQualities = setQuality,
             controlledVocabularies = controlledVocabularies, file = metrics_file, data_input = inputfiles_info,
             data_analysis = analysisSoftware, rawData_name = rawData_name)
  
  cat("done",'\n')
  cat(paste("mzQC file created at\n\n    ", rprt_fns$metrics_file, ".*\n\n", sep=""))
}

json_cvInfo = function(name, uri, version)
{
  return (list(name = name, uri = uri, version = version))
}


#'
#' Get the information of each CV term from an obo file.
#' 
#' @param cv_obo_file "xxx.obo"
#' @return A list containing cv term information
#' 
#' @export
#' 
parseOBO = function(cv_obo_file){
  ontology = ontologyIndex::get_ontology(cv_obo_file)
  obo = scan(file = cv_obo_file, what = "character")
  return(obo)
}

mzQC = setRefClass(
  'mzQC',
  fields = list(version = 'character',
                creationDate = 'MzQCDateTime',
                contactName = 'character',
                contactAddress = 'character',
                readMe = 'character',
                runQualities = 'runQuality',
                setQualities = 'setQuality',
                controlledVocabularies = 'controlledVocabulary'
                ),
  methods = list(
    toJSON = function(.self, ...)
    {
      return (jsonlite::toJSON(list(version = .self$version,
                                    creationDate = .self$creationDate
                                    ## todo, add the rest...
                                    ), ...))
    }
  )
)




asJSON <- jsonlite:::asJSON
setMethod('asJSON', 'BankAccount', function(x, ...) x$toJSON())
jsonlite::toJSON(bb)


