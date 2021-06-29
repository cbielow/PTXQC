
#'
#' Get the information of each CV term from an obo file.
#' 
#' @param cv_obo_file A path to an .obo file
#' @return A list containing cv term information
#' 
#' 
parseOBO = function(cv_obo_file){
  ontology = ontologyIndex:::get_ontology(cv_obo_file, extract_tags = "everything")
  #obo = scan(file = cv_obo_file, what = "character")
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
  both = list()
  for (name in names(qc))
  {
    both[[name]] = append(qc[[name]], ms[[name]])
  }
  return(both)
}

#' 
#' Define a Singleton class which can hold a CV dictionary (so we do not have to load the .obo files over and over again)
#' 
CVDictionarySingleton <- R6::R6Class("CVDictionarySingleton", inherit = R6P::Singleton, public = list(
  data = NULL
))


#' Fills a MzQCqualityMetric object with id(accession) and name.
#' The value (if any) and unit (if any) need to be set afterwards.
#' 
#' @param accession The ID (=accession) of the term in the CV
#' @param mzcv_dict A CV dictionary, as obtained by getCVDictionary(); defaults to a singleton, which needs to be filled manually beforehand
#' 
#' @return An instance of MzQCqualityMetric
#' 
getQualityMetricTemplate = function(accession, mzcv_dict = CVDictionarySingleton$new())
{
  idx = which(accession == mzcv_dict$id)
  if (length(idx) == 0) stop("Accession '", accession, "' is not a valid CV term in the current dictionary.")
  
  out = MzQCqualityMetric$new(accession, mzcv_dict$name[idx], mzcv_dict$def[idx])
  return(out)
}

#' Fills a MzQCcvParameter object with id(accession) and name.
#' The value (if any) needs to be set afterwards.
#' 
#' @param accession The ID (=accession) of the term in the CV
#' @param mzcv_dict A CV dictionary, as obtained by getCVDictionary(); defaults to a singleton, which needs to be filled manually beforehand
#' 
#' @return An instance of MzQCcvParameter
#' 
getCVTemplate = function(accession, mzcv_dict = CVDictionarySingleton$new())
{
  idx = which(accession == mzcv_dict$id)
  if (length(idx) == 0) stop("Accession '", accession, "' is not a valid CV term in the current dictionary.")
  
  out = MzQCcvParameter$new(accession, mzcv_dict$name[idx])
  return(out)
}


#'
#' Get an mzQC runQuality without actual metrics, but with full metadata
#' 
#' @param fc.raw.file For which run
#' @param raw_file_mapping A data.frame with cols 'from', to' and maybe 'best.effort' (if shorting was unsuccessful), as e.g. obtained by a FilenameMapper$raw_file_mapping
#' @return An MzQCrunQuality object
#'
getRunQualityTemplate = function(fc.raw.file, raw_file_mapping)
{
  
  idx = which(raw_file_mapping$to == fc.raw.file)
  if (length(idx) != 1) stop("fc.raw.file '", fc.raw.file, "' is not (or not unique) in mapping table.")
  
  raw_file = raw_file_mapping$from[id]

    ## todo: we're just guessing here...
  filename = paste0(raw_file, ".mzML"); 
  fullpath = paste0("???/", filename);
  
  mzML_format = getCVTemplate("MS:1000584")
  ptxqc_software = MzQCanalysisSoftware$new("MS:1003162", "PTX-QC", as.character(utils::packageVersion("PTXQC")), "https://github.com/cbielow/PTXQC/", "Proteomics (PTX) - QualityControl (QC) software for QC report generation and visualization.", "Proteomics Quality Control")
  
  out = MzQCrunQuality$new(MzQCmetadata$new(raw_file,  ## label
                                      list(MzQCinputFile$new(filename, fullpath, mzML_format)),
                                      list(ptxqc_software)),
                          list())
  
  return(out)
}

#'
#' Collects all 'mzQC' members from each entry in lst_qcMetrics and stores them in an overall mzQC object, which can be written to disk or augmented otherwise
#'
#' @param lst_qcMetrics A list of qcMetric objects which have their mzQC member populated with "MzQCrunQuality" and/or "MzQCsetQuality" objects
#' @param raw_file_mapping A data.frame with cols 'from', to' and maybe 'best.effort' (if shorting was unsuccessful), as e.g. obtained by a FilenameMapper$raw_file_mapping
#' @return An MzQCmzQC object
#' 
#' @export
#' 
assembleMZQC = function(lst_qcMetrics, raw_file_mapping)
{
  out = MzQCmzQC$new(version = "1.0.0", 
                     creationDate = MzQCDateTime$new(), 
                     contactName = Sys.info()["user"], 
                     contactAddress = NA_character_, 
                     description = NA_character_,
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

  run_qualities = list()
  set_qualities = list()
  
  for (metric in lst_qcMetrics)
  {
    mzqc_data = metric$mzQC
    if (is.null(mzqc_data)) next
    if (class(mzqc_data) != "list") stop("mzQC member of metric must be of class 'list()'")
    
    ## either an fc.raw.file or a concatenation of them
    for (name in names(mzqc_data))
    {
      if (name %in% raw_file_mapping$to)
      { ##
        ## runQuality
        ##
        rc = run_qualities[[name]];
        if (is.null(rc)) rc = getRunQualityTemplate(name, raw_file_mapping)
        ## append ...
        l = length(rc$qualityMetrics)
        rc$qualityMetrics[[l]] = mzqc_data[[name]]
        ## write back
        run_qualities[[name]] = rc
      } else {
        if (!grepl(";", name, fixed=TRUE)) stop("mzQC metric data must be an fc.raw.file or a contatenation of those using ';'. No ';' found!")
        ##
        ## setQuality
        ## 
        
      }
    }
    
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

#'
#' Writes a full mzQC object to disk
#' 
#' @param filepath A filename (with path) to write to. Should have '.mzQC' as suffix.
#' @param mzqc_obj An mzQC object, which is serialized to JSON and then written to disk
#'
writeMZQC = function(filepath, mzqc_obj)
{
  
  if (!hasFileSuffix(filepath, ".mzQC")) warning("'", filepath, "' does not end in '.mzQC'. Please fix the output filename.")
  
  
  content = jsonlite::toJSON(mzqc_obj)
  
  cat(content, file = filepath)
  
}




