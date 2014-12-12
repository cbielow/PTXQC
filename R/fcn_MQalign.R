#'
#' Return list of raw file names which were reported by MaxQuant as reference point for alignment.
#' 
#' There is only one reference point which has '0' as corrected RT (in all MQ versions we've seen so far). 
#' This is also true for fractions.
#' The reference can be identified by using the 'retention.time.calibration' column in evidence.txt.
#' The raw.file with very small values is usually the one (they are not exactly zero, even though they should be!)
#' 
#' Note that MaxQuant uses a guide tree to align the Raw files, so the order of files does not influence the 
#' alignment. But the first file will always be used as reference point when reporting delta-RTs. And this file is also
#' used by PTXQC as reference file vs all other files to find the real calibration function (see alignmentCheck()).
#'
#' This function might return multiple, or no raw file names (if MQ decides to change its mind at some point in the future).
#' In this case the result should be treated with caution or (better) regarded as failure.
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
#' @importFrom plyr ddply
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
    df.f$fraction = data$fraction[match(df.f$raws, data$raw.file)]
    df.f
    ## make a list of references:
    df.ref = df.f[!duplicated(df.f$fraction),]
    df.ref
    ## for each fraction: get list of references and the files that are allowed to match against them
    for (idx_ref in 1:nrow(df.ref))
    { # idx_ref = 1
      frac = df.ref$fraction[idx_ref]
      rFile = df.ref$raws[idx_ref]
      rFile
      clients = df.f$raws[df.f$fraction %in% (frac-1):(frac+1)]
      RefRounds[[idx_ref]] = list("file_ref" = rFile, "file_clients" = clients)
    }
    RefRounds
  }

  alignQ = data.frame()

  for (lRef in RefRounds)
  { ## lRef = RefRounds[[1]]
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
    alignQ_round$refFile = file_ref
    ## remove matches to self (i.e. with rtDiff==0)
    alignQ_round = alignQ_round[ alignQ_round$raw.file != file_ref, ]
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
#' @importFrom plyr ddply
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
  

#'
#' Check how close transferred ID's after alignment are to their genuine IDs within one Raw file.
#' 
#' The input is a data frame containing feature evidence with corrected retention times,
#' e.g. a 'calibrated.retention.time' column.
#'
#' Note that this function must be given all MS/MS identifications only (type "MULTI-MSMS" and "MSMS-MATCH")
#' in order to work correctly!
#' 
#' We compare for each peptide sequence (and charge) the RT difference between genuine vs. tranferred pairs.
#' For every comparison made, we report the RT difference. If alignment worked perfectly, the differences are very small (<1 min),
#' i.e. the pairs are accidentally split 3D peaks.
#' 
#' @param data A data.frame with columns 'type', 'calibrated.retention.time', 'modified.sequence', 'charge', 'raw.file'
#' @return A data.frame containing the RT diff for each pair found in a Raw file.
#'
#' @importFrom plyr ddply
#' 
idTransferCheck = function(data) {
  colnames(data) = tolower(colnames(data))
  
  if (!all(c('type', 'calibrated.retention.time', 'modified.sequence', 'charge', 'fc.raw.file') %in% colnames(data)))
  {
    stop("idTransferCheck(): columns missing!")  
  }
  if (!all(c("MULTI-MSMS", "MULTI-MATCH") %in% unique(data$type)))
  {
    stop('idTransferCheck(): scan types missing! Required: "MULTI-MSMS" and "MULTI-MATCH".')  
  }
  #data = d_evd
  data$seq_charge = paste(factor(data$modified.sequence), data$charge, sep="_")
  alignQ = ddply(data[,c("fc.raw.file", "type", "calibrated.retention.time", "seq_charge")], "fc.raw.file", function(x) {
    ## retain only stuff which has the potential to be a pair
    seqs = intersect(x$seq_charge[x$type=="MULTI-MSMS"], x$seq_charge[x$type=="MULTI-MATCH"])
    #seqs = x$seq_charge[x$type=="MULTI-MSMS"]
    #seqs = x$seq_charge[x$type=="MULTI-MATCH"]
    x_sub = x[x$seq_charge %in% seqs,]
    if (nrow(x_sub)==0) return(data.frame())
    ## for all pairs...
    align = ddply(x_sub, "seq_charge", function(x2) {
      d = diff(range(x2$calibrated.retention.time))
      d_bg = diff(range(x2$calibrated.retention.time[x2$type=='MULTI-MSMS']))  ## track background from genuine 3D peaks as well
      return (data.frame(rtdiff = d, rtdiff_bg = d_bg))
    })  
    return (align)    
  })
  alignQ$rtdiff_bg[alignQ$rtdiff_bg==0] = NA ## there was only one genuine feature
  #head(alignQ)
  #hist(log(1+alignQ$rtdiff), 100)
  #hist(log(1+alignQ$rtdiff_bg), 100, add=T, col="red")
  
  return(alignQ)
}


#'
#' Compute the fraction of 3D-peak pairs per Raw file which have an acceptable RT difference after alignment
#' 
#' Using the result from 'idTransferCheck()', compute the fraction of pairs which are within a certain RT tolerance.
#' 
#' Returned value is between 0 (bad) and 1 (all within tolerance).
#' 
#' @param data A data.frame with columns 'fc.raw.file' and !colname (param)
#' @param colname Name of column which contains the RT difference of a pair
#' @param df.allowed.deltaRT The allowed matching difference for each Raw file (as data.frame(fc.rawfile, m))
#' @return A data.frame with one row for each raw.file and columns 'raw.file' and score 'withinRT' (0-1)
#' 
#' @importFrom plyr ddply
#' 
inMatchWindow = function(data, colname, df.allowed.deltaRT)
{
  colnames(data) = tolower(colnames(data))
  
  if (!all(c(colname, 'fc.raw.file') %in% colnames(data)))
  {
    stop("ScoreInMatchWindow(): 'data' has missing columns!")  
  }
  alignQC = ddply(data, "fc.raw.file", function(x) {
    allowed.deltaRT = df.allowed.deltaRT$m[match(x$fc.raw.file[1], df.allowed.deltaRT$fc.raw.file)]
    withinRT = sum(abs(x[,colname]) < allowed.deltaRT, na.rm=T) / sum(!is.na(x[,colname]))
    return(data.frame(withinRT = withinRT))
  })
  
  return(alignQC)
}


#'
#' Determine fraction of evidence which causes segmentation, i.e. sibling peaks at different RTs
#' confirmed either by genuine or transferred MS/MS.
#'
#' Sometimes, MQ split a feature into 2 or more if the chromatograpic conditions are not optimal and there
#' is a drop in RT intensity.
#' If both features contain successful MS/MS scans, we will find the same peptide twice (with slightly different RT)
#' in the same charge state. This constitutes a natively split peak and is rare (95% of all peaks are unique).
#' 
#' If Match-between-runs is used and the RT alignment is not perfect, then a peptide might be inferred at a wrong
#' RT position, even though this Raw file already contains MS/MS evidence of this peptide.
#' Usually the number of peak duplicates rises drastically (e.g. only 75% of peaks are unique after MBR was used).
#' In most cases, the RT is too far off to be a split peak. It's rather a lucky hit with accidentally the same mass-to-charge,
#' and thus the intensity is random.
#' To find by how much these peak pairs differ in RT, use idTransferCheck() and inMatchWindow().
#' 
#' Required columns are 'match.time.difference', 'fc.raw.file', 'modified.sequence', 'charge', 'type'.
#' 
#' @param d_evd A data.frame of evidences containing the above columns
#' @return A data.frame with one row per Raw file and two columns: natively split peaks (%) and % of split peaks using matching.
#'
#' @importFrom plyr ddply
#'
peakSegmentation = function(d_evd)
{
  if (!all(c("match.time.difference", "fc.raw.file", "modified.sequence", "charge", 'type') %in% colnames(d_evd)))
  {
    stop("qualMBR(): columns missing!")  
  }
  
  d_evd$hasMTD = !is.na(d_evd$match.time.difference)
  
  cols = c("hasMTD", "fc.raw.file", "modified.sequence", "charge")
  splitInferred = ddply(d_evd[d_evd$type!="MSMS", cols], cols[-1], function(x)
  {
    #ratio =  NA
    #if (nrow(x)==2 & sum(x$hasMTD)==1) ratio = x$Intensity[!x$hasMTD] / x$Intensity[x$hasMTD]
    #r = rep(NA, nrow(x))
    ## mixed state
    #if (sum(x$hasMTD) >0 & sum(x$hasMTD)<nrow(x))  {
    #  r[!x$hasMTD]
    #} 
    return(data.frame(nNative = sum(!x$hasMTD), nMatched = sum(x$hasMTD)))#, ratio = ratio))
  })
  
  mbr_score = ddply(splitInferred, "fc.raw.file", function(splitInferred)
  {
    ddt = table(splitInferred[, c("nMatched", "nNative")])
    ### ddt might look like this:
    #                 nNative
    #     nMatched    0    1   2  3  4 5 6 7 8 9 10 11 13 14 15 18 20 26 27 28 31 37 39
    #            0    0 6249 230 41 16 4 1 4 2 0  0  0  1  1  0  0  0  1  1  0  0  0  0
    #            1 1352  303  34 16  6 4 2 1 1 0  1  1  0  1  0  1  0  0  0  0  0  0  0
    #            2   49   29   4  2  1 0 0 1 0 0  1  
    
    ### scale ddt matrix from counts of events to counts of peaks,
    ### e.g. event: each 1native+2matched affects three peaks
    d1 = dimnames(ddt)$nMatched
    d2 = dimnames(ddt)$nNative
    #ddt_b = ddt
    for (i in 1:length(d1)) {
      for (j in 1:length(d2)) {
        ddt[i,j] = ddt[i,j] * (as.numeric(d1[i]) + as.numeric(d2[j]))
      }
    }
    ## add nMatched=1 row, if not present
    ddt = rbind(ddt, 0)

    ## segmentation genuine (ignore the nMatched, i.e. project everything onto first row)
    n.perf = sum(ddt[,colnames(ddt)==1]) ## uniquely matched
    n.all = sum(ddt[,!colnames(ddt)==0]) ## all natives (==0 is the matched-only column)
    corr.nat = n.perf/n.all ## e.g. 94%
    
    ## segmentation of matched-only
    i.perf = ddt[rownames(ddt)==1, colnames(ddt)==0] ## matched singlets (nMatched=1,nNative=0)
    i.all = sum(ddt[rownames(ddt)!=0, ])             ## all groups which contain a match (==without first row)
    corr.matched = max(0, i.perf) / i.all  ## e.g. 50%
    
    ## segmentation of all evidence (weighed average of above)
    comb.perf = ddt[rownames(ddt)==1, colnames(ddt)==0] + ## matched singlets (nMatched=1,nNative=0)
                ddt[rownames(ddt)==0, colnames(ddt)==1]
    comb.all = sum(ddt)             ## all groups
    corr.combined = max(0, comb.perf) / comb.all       ## e.g. 87%
    
    ## final result
    data.frame(corr.nat = corr.nat, corr.matched = corr.matched, corr.combined = corr.combined)
  })
  return (mbr_score = mbr_score)
}

#' Combine several data structs into a final picture for segmentation incurred by 'Match-between-runs'.
#' 
#' ...
#'
#' @param qMBR              A data.frame as computed by peakSegmentation()
#' @param qMBRSeg_Dist_r    A data.frame as computed by inMatchWindow() for using the 'rtdiff' (=all) column
#' @param qMBRSeg_Dist_r_bg A data.frame as computed by inMatchWindow() for using the 'rtdiff_bg' (=genuine only) column
#' @return A data.frame which details the distribution of singlets and pairs (inRT and outRT) for each Raw file and genuine vs. all
#' 
#' @importFrom plyr ddply
#' 
computeMatchRTFractions = function(qMBR, qMBRSeg_Dist_r, qMBRSeg_Dist_r_bg)
{
  ## data might look like this:
#     fc.raw.file  corr.nat corr.matched corr.combined
#   1  323_G.._01 0.9740162    0.3393057     0.9517918
#   2  521_G.._02 0.9813251    0.7689088     0.9433782
#   3  522_G.._01 0.9671893    0.6986601     0.9121866
# 
# qMBRSeg_Dist_r
#                      fc.raw.file   withinRT
# 1 20100730_Velos1_TaGe_SA_K562_1 0.17647059
# 2 20100730_Velos1_TaGe_SA_K562_2 0.28368794
# 3 20100730_Velos1_TaGe_SA_K564_3 0.38181818
# 4 20100730_Velos1_TaGe_SA_K565_4 0.37730871
# 5 20100730_Velos1_TaGe_SA_K565_5 0.33333333
# .. same for qMBRSeg_Dist_r_bg...
  
  ## compute percentage of outside dRT peaks in genuine. Same for matched evidence. And combined(=all)
  ## then calc the drop.
  f = ddply(qMBR, "fc.raw.file", function(x) {
    nat.inRT = qMBRSeg_Dist_r_bg$withinRT[qMBRSeg_Dist_r_bg$fc.raw.file==x$fc.raw.file]
    if (length(nat.inRT) == 0) nat.inRT=NA
    x$corr.nat.inRT = (1-x$corr.nat)*nat.inRT
    x$corr.nat.outRT = (1-x$corr.nat)*(1-nat.inRT)

    matched.inRT = qMBRSeg_Dist_r$withinRT[qMBRSeg_Dist_r$fc.raw.file==x$fc.raw.file]
    if (length(matched.inRT) == 0) matched.inRT=NA
    x$corr.matched.inRT = (1-x$corr.matched)*matched.inRT
    x$corr.matched.outRT = (1-x$corr.matched)*(1-matched.inRT)
    
    combined.inRT = matched.inRT ## same as matched (takes all groups)
    if (length(combined.inRT) == 0) combined.inRT=NA
    x$corr.combined.inRT = (1-x$corr.combined)*combined.inRT
    x$corr.combined.outRT = (1-x$corr.combined)*(1-combined.inRT)
    
    r = data.frame(fc.raw.file=x$fc.raw.file, 
                   single=c(x$corr.nat, x$corr.matched, x$corr.combined),
                   multi.inRT = c(x$corr.nat.inRT, x$corr.matched.inRT, x$corr.combined.inRT),
                   multi.outRT = c(x$corr.nat.outRT, x$corr.matched.outRT, x$corr.combined.outRT),
                   sample = factor(c("genuine", "matched", "all"), levels=c("genuine", "matched", "all")))
    return(r)    
  })
  return(f)
}
  