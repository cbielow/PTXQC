
#'
#' Get the information of each CV term from an obo file.
#' 
#' @param cv_obo_file A path to an .obo file
#' @return A list containing cv term information
#' 
#' @export
#' 
parseOBO = function(cv_obo_file){
  ontology = ontologyIndex:::get_ontology(cv_obo_file, extract_tags = "everything")
  obo = scan(file = cv_obo_file, what = "character")
  return(ontology)
}



#'
#' Parse the content of 'qc-cv.obo' and 'psi-ms.obo' from the 'PTXQC/cv/' folder and return their union as ontology
#' 
#' @return a list with 'id', 'name', 'def', 'parents', 'children' which contains the CV entries
#' 
#' @export
#' 
getCVDictionary = function()
{
  qc = parseOBO(system.file("./cv/qc-cv.obo", package="PTXQC"))
  ms = parseOBO(system.file("./cv/psi-ms.obo", package="PTXQC"))
  return(rbind(qc, ms))
}

# mzcv_dict = getCVDictionary()
# accession = "QC:4000020"




#' Fills a MzQCqualityMetric object with id(accession) and name.
#' The value (if any) and unit (if any) need to be set afterwards.
getQualityMetricTemplate = function(accession, mzcv_dict)
{
  idx = which(accession == mzcv_dict$id)
  if (length(idx) == 0) stop("Accession '", accession, "' is not a valid CV term in the current dictionary.")
  
  out = MzQCqualityMetric$new(accession, mzcv_dict$name[idx], mzcv_dict$def[idx])
  return(out)
}

#' Fills a MzQCcvParameter object with id(accession) and name.
#' The value (if any) needs to be set afterwards.
getCVTemplate = function(accession, mzcv_dict)
{
  idx = which(accession == mzcv_dict$id)
  if (length(idx) == 0) stop("Accession '", accession, "' is not a valid CV term in the current dictionary.")
  
  out = MzQCcvParameter$new(accession, mzcv_dict$name[idx])
  return(out)
}


#'
#' Collects all 'mzQC' members from each entry in lst_qcMetrics and stores them in an overall mzQC object, which can be written to disk or augmented otherwise
#'
#' @param lst_qcMetrics A list of qcMetric objects which have their mzQC member populated
#' 
#' @export
#' 
assembleMZQC = function(lst_qcMetrics)
{
  out = MzQCmzQC$new(version = "1.0.0", 
                     #creationDate = MzQCDateTime$new(), 
                     contactName = Sys.info()["user"], 
                     contactAddress = NA_character_, 
                     readMe = NA_character_,
                     runQualities = list(),
                     setQualities = list(), 
                     controlledVocabularies = list(
                                              MzQCcontrolledVocabulary$new(
                                                "Proteomics Standards Initiative Quality Control Ontology",
                                                "https://github.com/HUPO-PSI/mzQC/blob/master/cv/qc-cv.obo",
                                                "1.2.0"),
                                              MzQCcontrolledVocabulary$new(
                                                "Proteomics Standards Initiative Mass Spectrometry Ontology",
                                                "https://github.com/HUPO-PSI/psi-ms-CV/blob/master/psi-ms.obo",
                                                "4.1.7")))

  run_qualities = list(
    ## TEMP!!
    ## 
    ## 
    ## !!!
    #mzML = getCVTemplate("MS:1000584", mzcv_dict)
    ## 
    #                    MzQCrunQuality$new(MzQCmetadata$new("label",
    #                                                        list(MzQCinputFile$new("a.mzML", "c:\temp\a.mzML", MzQCcvParameter$new()))),
    #                                       list(MzQCqualityMetric$new("acc", "name", "desc", 1, "millisecond"))
  )                
                       
  set_qualities = list()
  for (metric in lst_qcMetrics)
  {
    mzqc_data = metric$mzQC
    if (is.null(mzqc_data)) next
    if (class(mzqc_data) != "list") stop("mzQC member of metric must be of class 'list()'")
    
    cl = lapply(mzqc_data, class)
    if (!all(cl %in% c("MzQCrunQuality", "MzQCsetQuality"))) stop("mzQC member must contain 'MzQCsetQuality' or 'MzQCsetQuality' elements in its list. Given: '", paste(cl, collapse = ","), "'.")
    
    run_qualities = append(run_qualities, mzqc_data[cl == 'MzQCrunQuality'])
    set_qualities = append(set_qualities, mzqc_data[cl == 'MzQCsetQuality'])
  }
  out$runQualities = run_qualities
  out$setQualities = set_qualities
  
  return(out)
}


#'
#' Checks if filepath ends in suffix. If suffix does not start with a '.' it is prepended automatically.
#' 
#' @return TRUE if yes, FALSE otherwise
#' 
#' @export
#' 
#' @examples 
#'   hasFileSuffix("bla.txt", "txt")    # TRUE
#'   hasFileSuffix("bla.txt", ".txt")   # TRUE
#'   hasFileSuffix("bla.txt", "doc")    # FALSE
#'   hasFileSuffix("bla.txt", ".doc")   # FALSE
#'   hasFileSuffix("fo", ".doc")        # FALSE
#'   hasFileSuffix("", ".doc")          # FALSE
#'   hasFileSuffix("foo", "")           # FALSE
#' 
hasFileSuffix = function(filepath, suffix)
{
  if (substr(suffix,1,1) != '.') suffix = paste0('.', suffix)
  
  filepath = tolower(filepath)
  suffix = tolower(suffix)
  
  return(suffix == substring(filepath, first = nchar(filepath) - nchar(suffix) + 1))
}

writeMZQC = function(filepath, mzqc_obj)
{
  
  if (!hasFileSuffix(filepath, ".mzQC")) warning("'", filepath, "' does not end in '.mzQC'. Please fix the output filename.")
  
  
  content = jsonlite::toJSON(mzqc_obj)
  
  cat(content, file = filepath)
  
}

