

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
  out = NULL
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
    out = data.frame(file = basename(xml_rawfiles), path = xml_rawfiles, file_no_suffix = removeFileSuffix(basename(xml_rawfiles)), CV = rmzqc::filenameToCV(xml_rawfiles))
  }
  return (out)
}



#'
#' Get an mzQC runQuality without actual metrics, but with full metadata
#' 
#' @param fc.raw.file For which run
#' @param raw_file_mapping A data.frame with cols 'from', 'to' and maybe 'best.effort' (if shorting was unsuccessful), as e.g. obtained by a FilenameMapper$raw_file_mapping
#' @return An MzQCrunQuality object
#' 
#' @import rmzqc
#'
getRunQualityTemplate = function(fc.raw.file, raw_file_mapping)
{
  
  idx = which(raw_file_mapping$to == fc.raw.file)
  if (length(idx) != 1) stop("fc.raw.file '", fc.raw.file, "' is not (or not unique) in mapping table.")
  
  raw_file = as.character(raw_file_mapping$from[idx])
  meta = QCMetaFilenames$new()$data
  if (is.null(meta) || sum(meta$file_no_suffix == raw_file) == 0) {
    ## we're just guessing here...
    warning("Cannot properly fill metadata of mzQC file, since full filenames are unknown. Using placeholders.")
    filename = paste0(raw_file, ".raw"); 
    fullpath = paste0("???/", filename);
    accession = rmzqc::filenameToCV(filename)
  } else {
    idx_meta = which(meta$file_no_suffix == raw_file)
    filename = as.character(meta$file[idx_meta])
    fullpath = as.character(meta$path[idx_meta])
    accession = as.character(meta$CV[idx_meta])
  }
  
  file_format = rmzqc::getCVTemplate(accession = accession)
  ptxqc_software = rmzqc::toAnalysisSoftware(id = "MS:1003162", version = as.character(utils::packageVersion("PTXQC")))
  
  out = rmzqc::MzQCrunQuality$new(rmzqc::MzQCmetadata$new(raw_file,  ## label
                                                          list(rmzqc::MzQCinputFile$new(filename, fullpath, file_format)),
                                                          list(ptxqc_software)),
                                  list())
  
  return(out)
}

#'
#' Collects all 'mzQC' members from each entry in lst_qcMetrics and stores them in an overall mzQC object, which can be written to disk (see writeMZQC()) or augmented otherwise
#'
#' @param lst_qcMetrics A list of qcMetric objects which have their mzQC member populated with "MzQCrunQuality" and/or "MzQCsetQuality" objects
#' @param raw_file_mapping A data.frame with cols 'from', to' and maybe 'best.effort' (if shorting was unsuccessful), as e.g. obtained by a FilenameMapper$raw_file_mapping
#' @return An MzQCmzQC object (root object of an mzQC document)
#' 
#' @export
#' 
assembleMZQC = function(lst_qcMetrics, raw_file_mapping)
{
  out = rmzqc::MzQCmzQC$new(version = "1.0.0", 
                            creationDate = MzQCDateTime$new(), 
                            contactName = Sys.info()["user"], 
                            contactAddress = NA_character_, 
                            description = NA_character_,
                            runQualities = list(),
                            setQualities = list(), 
                            controlledVocabularies = list(rmzqc::getDefaultCV()))

  run_qualities = list()
  set_qualities = list()
  
  for (metric in lst_qcMetrics)
  {
    mzqc_data = metric$mzQC
    if (is.null(mzqc_data)) next
    if (!inherits(mzqc_data, "list")) stop("mzQC member of metric must be of class 'list()'")
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





