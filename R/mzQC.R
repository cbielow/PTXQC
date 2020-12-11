
createMzQC = function(rprt_fns)
{
  cv_obo_file = system.file("./mzQC/qc-cv.obo", package="PTXQC")
  
  ##save required data of each metric in JSON-format
  cat("Creating metric json file ...", '\n')
  ##load the cv Term
  metric_info <- read_cvTerm(cv_obo_file)
  ##copy cv Term table
  metric_info <- data.frame(matrix(unlist(metric_info), ncol = 4))
  colnames(metric_info) <- c("cvRef", " accession", "name", "unit")
  ##Read the value stored in each class
  metrics_file = rprt_fns$metrics_file
  metric_info = readValue(qcMetric = lst_qcMetrics, qc_cv = metric_info)
  
  input_rawData <- list(df_evd, df_msms, df_msmsScans, d_parAll, df_pg, d_smy)
  rawData_name <- c("df_evd", "df_msms", "df_msmsScans", "d_parAll", "df_pg", "d_smy")
  inputfiles_info <- get_inputInfo(input_rawData , rawData_name)
  
  analysisSoftware <- list(name = "PTXQC",
                           version = as.character(packageVersion("PTXQC")),
                           uri = "github.com/cbielow/PTXQC")
  
  mzQC_version <- "0.1.1"
  creationDate <- as.character(Sys.time(), format="%Y-%m-%dT%H:%M:%S" )
  runQuality = metric_info[which(metric_info$quality_type == "runQuality"), ]
  setQuality = metric_info[which(metric_info$quality_type == "setQuality"), ]
  controlledVocabularies <- list(list(name = "Quality control metrics generating software",
                                      uri = "github.com/HUPO-PSI/mzQC/blob/bulk-cvterms/cv/qc-cv.obo",
                                      version = "0.1.2" ))
  ## create output file
  createmzQC(version = mzQC_version, creationDate = creationDate, runQualities = runQuality, setQualities = setQuality,
             controlledVocabularies = controlledVocabularies, file = metrics_file, data_input = inputfiles_info,
             data_analysis = analysisSoftware, rawData_name = rawData_name)
  
  cat("done",'\n')
  cat(paste("mzQC file created at\n\n    ", rprt_fns$metrics_file, ".*\n\n", sep=""))
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
  ontology <- get_ontology(cv_obo_file)
  obo = scan(file = cv_obo_file, what = "character")
  return(obo)
}



