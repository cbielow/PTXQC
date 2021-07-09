
#'
#' Get the information of each CV term from an obo file.
#' 
#' @param cv_obo_file A path to an .obo file
#' @return A list containing cv term information
#' 
#' 
parseOBO = function(cv_obo_file){
  ontology = ontologyIndex::get_ontology(cv_obo_file, extract_tags = "everything")
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
#' @export
#' 
CVDictionarySingleton <- R6::R6Class("CVDictionarySingleton", inherit = R6P::Singleton, public = list(
  #' @field data Stores the data of the singleton. Set the data once before using the singleton all over the place
  data = NULL
))



#' 
#' Define a Singleton class which holds the full raw filenames (+path) and their PSI-MS CV terms for usage in the mzQC metadata
#' 
#' @export
#' 
QCMetaFilenames <- R6::R6Class("QCMetaFilenames", inherit = R6P::Singleton, public = list(
  #' @field data Stores the data of the singleton. Set the data once before using the singleton all over the place
  data = NULL
))


getMetaFilenames = function(mqpar_file, base_folder)
{
  out = NA
  ## mqpar_file = "Z:\\projects\\QC\\PTXQC\\data\\ecoli_small\\mqpar.xml"
  ## base_folder = "Z:\\projects\\QC\\PTXQC\\data\\ecoli_small\\combined\\txt\\"
  xml_rawfiles = getMQPARValue(mqpar_file, "//string[parent::filePaths|parent::Filenames]", allow_multiple = TRUE)
  if (is.null(xml_rawfiles)) {
    ## try again using parent directory
    warning("No mqpar.xml found in '", mqpar_file, "'. Trying two folders up.")
    up_dir = paste0(base_folder, "/../../mqpar.xml")
    xml_rawfiles = getMQPARValue(up_dir, "//string[parent::filePaths|parent::Filenames]", allow_multiple = TRUE)
  }
  ## second try...
  if (is.null(xml_rawfiles)) {
    warning("No mqpar.xml found in '", up_dir, "'. Giving up to read full file paths The mqQC file will contain incomplete information and may not validate.")
  } else {
    ## cannot use 'eval(expr_fn_map)$raw_file_mapping' yet, since we do not have read any .txt files which fills the mapping
    out = data.frame(file = basename(xml_rawfiles), path = xml_rawfiles, file_no_suffix = removeSuffix(basename(xml_rawfiles)), CV = suffixToCV(xml_rawfiles))
  }
  return (out)
}


#' Fills a MzQCqualityMetric object with id(accession) and name.
#' The value (if any) and unit (if any) need to be set afterwards.
#' 
#' @param accession The ID (=accession) of the term in the CV
#' @param mzcv_dict A CV dictionary, as obtained by getCVDictionary(); defaults to a singleton, which needs to be filled manually beforehand
#' 
#' @return An instance of MzQCqualityMetric
#' 
getQualityMetricTemplate = function(accession, mzcv_dict = CVDictionarySingleton$new()$data)
{
  idx = which(accession == mzcv_dict$id)
  if (length(idx) == 0) stop("Accession '", accession, "' is not a valid CV term in the current dictionary (", length(mzcv_dict$id), " entries].")
  
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
getCVTemplate = function(accession, mzcv_dict = CVDictionarySingleton$new()$data)
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
#' @param raw_file_mapping A data.frame with cols 'from', 'to' and maybe 'best.effort' (if shorting was unsuccessful), as e.g. obtained by a FilenameMapper$raw_file_mapping
#' @return An MzQCrunQuality object
#'
getRunQualityTemplate = function(fc.raw.file, raw_file_mapping)
{
  
  idx = which(raw_file_mapping$to == fc.raw.file)
  if (length(idx) != 1) stop("fc.raw.file '", fc.raw.file, "' is not (or not unique) in mapping table.")
  
  raw_file = as.character(raw_file_mapping$from[idx])
  meta = QCMetaFilenames$new()$data
  if (is.null(meta) || is.na(meta) || sum(meta$file_no_suffix == raw_file) == 0) {
    ## we're just guessing here...
    warning("Cannot properly fill metadata of mzQC file, since full filenames are unknown. Using placeholders.")
    filename = paste0(raw_file, ".raw"); 
    fullpath = paste0("???/", filename);
    accession = suffixToCV(".raw")
  } else {
    idx_meta = which(meta$file_no_suffix == raw_file)
    filename = as.character(meta$file[idx_meta])
    fullpath = as.character(meta$path[idx_meta])
    accession = as.character(meta$CV[idx_meta])
  }
  
  file_format = getCVTemplate(accession = accession)
  ptxqc_software = MzQCanalysisSoftware$new("MS:1003162", "PTX-QC", 
                                            as.character(utils::packageVersion("PTXQC")), 
                                            "https://github.com/cbielow/PTXQC/",
                                            "Proteomics (PTX) - QualityControl (QC) software for QC report generation and visualization.",
                                            "Proteomics Quality Control")
  
  out = MzQCrunQuality$new(MzQCmetadata$new(raw_file,  ## label
                                            list(MzQCinputFile$new(filename, fullpath, file_format)),
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
    if (length(mzqc_data) == 0) next
    if (any(is.null(names(mzqc_data)))) stop("mzQC member of metric '", metric$qcName, "' is a list, but has no names (must be fc.raw.file names or a concatenation of those)")
    
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
        l = 1 + length(rc$qualityMetrics)
        rc$qualityMetrics[[l]] = mzqc_data[[name]]
        ## write back
        run_qualities[[name]] = rc
      } else {
        if (!grepl(";", name, fixed=TRUE)) stop("mzQC metric data must be an fc.raw.file or a contatenation of those using ';'. No ';' found!")
        ##
        ## setQuality
        ## 
        stop("setQuality not implemented yet")
      }
    }
    
  }
  ## remove the names from lists (to make it a JSON array; otherwise it would be an object)
  names(run_qualities) = NULL
  names(set_qualities) = NULL
  
  out$runQualities = run_qualities
  out$setQualities = set_qualities
  
  return(out)
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
  
  
  content = jsonlite::toJSON(mzqc_obj, pretty = TRUE, auto_unbox = TRUE)
  
  cat(content, file = filepath)
  
}


#'
#' For a given filename, check the suffix and translate it to an PSI-MS CV term, e.g. 'MS:1000584'
#' 
#' The following mapping is currently known:
#' .raw    : MS:1000563 ! Thermo RAW format
#' .mzML   : MS:1000584 ! mzML format
#' .mzData : MS:1000564 ! PSI mzData format
#' .wiff   : MS:1000562 ! ABI WIFF format
#' .pkl    : MS:1000565 ! Micromass PKL format
#' .mzXML  : MS:1000566 ! ISB mzXML format
#' .yep    : MS:1000567 ! Bruker/Agilent YEP format
#' .dta    : MS:1000613 ! Sequest DTA format
#' .mzMLb  : MS:1002838 ! mzMLb format
#' 
#' Falls back to 'MS:1000560 ! mass spectrometer file format' if no match could be found.
#' 
#' @param filepath A filename (with optional path)
#' @return A CV term accession as string, e.g. 'MS:1000584'
#' 
#' @examples 
#'   suffixToCV("test.mZmL")  # MS:1000584
#'   suffixToCV("test.raw")  # MS:1000563
#'   suffixToCV(c("test.raw", "bla.mzML"))
#'
#' @export
#' 
suffixToCV = function(filepath)
{
  if (length(filepath) > 1) return(sapply(filepath, suffixToCV))
  
  if (hasFileSuffix(filepath, ".raw")) return ("MS:1000563");
  if (hasFileSuffix(filepath, ".mzML")) return ("MS:1000584");
  if (hasFileSuffix(filepath, ".mzData")) return ("MS:1000564");
  if (hasFileSuffix(filepath, ".wiff")) return ("MS:1000562");
  if (hasFileSuffix(filepath, ".pkl")) return ("MS:1000565");
  if (hasFileSuffix(filepath, ".mzXML")) return ("MS:1000566");
  if (hasFileSuffix(filepath, ".yep")) return ("MS:1000567");
  if (hasFileSuffix(filepath, ".dta")) return ("MS:1000613");
  if (hasFileSuffix(filepath, ".mzMLb")) return ("MS:1002838");
  
  warning("File '", filepath, "' has an unknown suffix. Falling back to 'MS:1000560 ! mass spectrometer file format'.")
  return("MS:1000560")
}
  



