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
#' An 'id' column must be present, to enable mapping the result of this function to the original data frame.
#' 
#' A reference Raw file can be identified using 'findAlignReference()'. If Maxquants experimental design included
#' pre-fractionation, a column named 'fraction' should be given and 'referenceFile' should be empty. This function will
#' pick the one Raw file for each fraction (the first in order) to use as reference. Only the immediately neighbouring
#' fractions will be matched to this reference.
#' 
#' @param data A data.frame with columns 'calibrated.retention.time', 'id', 'modified.sequence', 'charge', 'raw.file' and 'fraction' (if present)
#' @param referenceFile A raw file name as occuring in data$raw.file, serving as alignment reference (when no fractions are used).
#' @return A data.frame containing the RT diff for each feature found in a Raw file and the reference.
#'
alignmentCheck = function(data, referenceFile) {
  colnames(data) = tolower(colnames(data))
  
  if (!all(c('calibrated.retention.time', 'id', 'modified.sequence', 'charge', 'raw.file') %in% colnames(data)))
  {
    stop("alignmentCheck(): columns missing!")  
  }
  if (is.na(referenceFile) & !('fraction' %in% colnames(data)))
  {
    stop("alignmentCheck(): Either a referenceFile or 'fraction' column must be provided!")  
  }
  if (!is.na(referenceFile) & ('fraction' %in% colnames(data)))
  {
    stop("alignmentCheck(): referenceFile and 'fraction' cannot be provided at the same time!")
  }
  
  if (!is.na(referenceFile)) {
    RefRounds = list( list("file_ref" = referenceFile, "file_clients" = unique(data$raw.file)) )
  } else {
    RefRounds = list()
    ## get Raw files and their fractions
    df.f = data.frame(raws = unique(data$raw.file))
    df.f$fraction = data$fraction[match(raws, data$raw.file)]
    df.f
    ## make a list of references:
    df.ref = df.f[!duplicated(df.f$fraction),]
    df.ref
    ## for each fraction: get list of references and the files that are allowed to match against them
    for (idx_ref in 1:nrow(df.ref))
    { # idx_ref = 1
      frac = df.ref$fraction[idx_ref]
      referenceFile = df.ref$raws[idx_ref]
      referenceFile
      clients = df.f$raws[df.f$fraction %in% (frac-1):(frac+1)]
      RefRounds[[idx_ref]] = list("file_ref" = referenceFile, "file_clients" = clients)
    }
    RefRounds
  }

  alignQ = data.frame()

  for (lRef in RefRounds)
  {
    file_ref = lRef[["file_ref"]]
    file_clients = lRef[["file_clients"]]
    
    ## prefilter data: only peptides which occur in reference
    ##            and are part of client list
    data_round = data[data$modified.sequence %in% data$modified.sequence[data$raw.file==file_ref] &
                      data$raw.file %in% file_clients,]
    
    alignQ_round = ddply(data_round, c("modified.sequence", "charge"), function(x)
    {
      ## reference must be present
      if ((nrow(x)==1) | (sum(x$raw.file == file_ref)==0) | (sum(duplicated(x$raw.file))>0)) {
        return (data.frame())
      }
      
      ## we could rescue this case, but for performance reasons we don't
      ## duplicates per raw file... ignore this raw file
      #if (sum(duplicated(x$raw.file))>0) {
      #  x$calibrated.retention.time[x$raw.file %in% x$raw.file[duplicated(x$raw.file)]] = NA
      #}
      
      ## RT of reference (we use 'mean' since there could be multiple ones when using fractions)
      m = x$calibrated.retention.time[x$raw.file == file_ref]
      
      x$rtdiff = x$calibrated.retention.time - m;
      return(x[,c("raw.file", "id", "rtdiff", 'calibrated.retention.time')])
    })
    alignQ = rbind(alignQ, alignQ_round)
  }
  
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
  
  
  
  