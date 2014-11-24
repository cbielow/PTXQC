#'
#' Return list of raw file names which were used as MaxQuant alignment reference.
#' 
#' Usually there is only one reference -- except when fractions are used.
#' The reference can be identified by using the 'retention.time.calibration' column in evidence.txt.
#' The raw.file with very small values is usually the one (they are not exactly zero, even though they should be!)
#'
#' This function might return multiple, or no raw file names. In this case the result should be treated with
#' caution or (better) regarded as failure.
#' 
#' @param data The data.frame with columns 'retention.time.calibration' and 'raw.file'
#' @param threshold Maximum range (min to max) of values required to be recognized as reference
#' @return List of reference raw files (usually just one)
#' 
#' @importFrom plyr ddply
#'
findAlignReference = function(data, threshold = 1e-4)
{
  colnames(data) = tolower(colnames(data))
  if (!("retention.time.calibration" %in% colnames(data)))
  {
    stop("findAlignReference(): Error, could not find column 'retention.time.calibration' in data. Aborting!")  
  }
  if (!("raw.file" %in% colnames(data)))
  {
    stop("findAlignReference(): Error, could not find column 'raw.file' in data. Aborting!")  
  }
  fr = ddply(data, "raw.file", function(x) data.frame(range = diff(range(x$retention.time.calibration, na.rm=T))))
  ref = as.character(fr$raw.file[fr$range < threshold])
  return (ref)  
}


#'
#' Verify an alignment by checking the retention time differences of identical peptides across Raw files
#' 
#' The input is a data frame containing feature evidence with corrected retention times,
#' e.g. a 'calibrated.retention.time' column.
#'
#' Note that this function must be given real MS/MS identifications only (type "MULTI-MSMS")
#' in order to work correctly!
#' 
#' For each peptide sequence (and charge) in the reference Raw file, this function looks up the
#' already calibrated retention time difference of the same feature in all other files. For every comparison made,
#' we report the RT difference. If alignment worked perfectly, the differences are very small (<1 min).
#' 
#' An 'id' column must be present, to enable matching the result of this function to the original data frame.
#' 
#' A reference Raw file can be identified using 'findAlignReference()'.
#' 
#' @param data A data.frame with columns 'calibrated.retention.time', 'id', 'modified.sequence', 'charge' and 'raw.file'
#' @param referenceFile A raw file name as occuring in data$raw.file, serving as alignment reference
#' @return A data.frame containing the RT diff for each feature found in a Raw file and the reference.
#'
alignmentCheck = function(data, referenceFile) {
  colnames(data) = tolower(colnames(data))
  
  if (!all(c('calibrated.retention.time', 'id', 'modified.sequence', 'charge', 'raw.file') %in% colnames(data)))
  {
    stop("alignmentCheck(): columns missing!")  
  }
  print(dim(data))
  ## prefilter data: only peptides which occur in reference
  data = data[data$modified.sequence %in% data$modified.sequence[data$raw.file==referenceFile],]
  
  alignQ = ddply(data, c("modified.sequence", "charge"), function(x)
  {
    ## reference must be present
    if ((nrow(x)==1) | (sum(x$raw.file==referenceFile)==0) | (sum(duplicated(x$raw.file))>0)) {
      return(data.frame())
    }
    
    ## we could rescue this case, but for performance reasons we don't
    ## duplicates per raw file... ignore this raw file
    #if (sum(duplicated(x$raw.file))>0) {
    #  x$calibrated.retention.time[x$raw.file %in% x$raw.file[duplicated(x$raw.file)]] = NA
    #}
    
    m = x$calibrated.retention.time[x$raw.file == referenceFile]
    x$rtdiff = x$calibrated.retention.time - m;
    return(x[,c("raw.file", "id", "rtdiff", 'calibrated.retention.time')])
  })
  
  return(alignQ)
}


#'
#' Compute the fraction of features per Raw file which have an acceptable RT difference after alignment
#' 
#' Using the result from 'alignmentCheck()', score the features of every Raw file and see if they 
#' have been properly aligned.
#' Returned value is between 0 (bad) and 1 (all aligned).
#' 
#' @param data A data.frame with columns 'rtdiff' and 'raw.file'
#' @param allowed.deltaRT The allowed matching difference (1 minute by default)
#' @return A data.frame with one row for each raw.file and columns 'raw.file' and 'withinRT' (0-1)
#' 
#' 
ScoreInAlignWindow = function(data, allowed.deltaRT = 1)
{
  colnames(data) = tolower(colnames(data))
  
  if (!all(c('rtdiff', 'raw.file') %in% colnames(data)))
  {
    stop("alignmentCheck(): columns missing!")  
  }
  alignQC = ddply(data, "raw.file", function(x) {
    withinRT = sum(abs(x$rtdiff) < allowed.deltaRT, na.rm=T) / sum(!is.na(x$rtdiff))
    return(data.frame(withinRT = withinRT))
  })
  
  return(alignQC)
}
  
  
  
  