#'
#' Return list of raw file names which were reported by MaxQuant as reference point for alignment.
#' 
#' There is only one reference point which has '0' in 'retention.time.calibration' column in evidence.txt 
#' as corrected RT. This is true for most MaxQuant versions
#' and also true for fractions. However, some evidence.txt files show 0.03 as an averaged minimum per Raw file.
#' We use the raw.file with the smallest average as reference.
#' 
#' Note that MaxQuant uses a guide tree to align the Raw files, so the order of files does not influence the 
#' alignment. But the first file will always be used as reference point when reporting delta-RTs. And this file is also
#' used by PTXQC as reference file vs all other files to find the real calibration function (see alignmentCheck()).
#'
#' This function might return multiple raw file names (if MQ decides to change its mind at some point in the future).
#' In this case the result should be treated with caution or (better) regarded as failure.
#' 
#' @param data The data.frame with columns 'retention.time.calibration' and 'raw.file'
#' @return List of reference raw files (usually just one)
#'
findAlignReference = function(data)
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
  fr = plyr::ddply(data, "raw.file", function(x) data.frame(range = diff(range(x$retention.time.calibration, na.rm = TRUE))))
  ref = as.character(fr$raw.file[fr$range <= min(fr$range)])
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
    ## Fractions ...
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
    
    alignQ_round = plyr::ddply(data_round, c("modified.sequence", "charge"), function(x)
    {
      ## reference must be present
      if ((nrow(x)==1) | (sum(x$raw.file == file_ref)==0) | (sum(duplicated(x$raw.file))>0)) {
        return (NULL)
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
    ## alignQ_round might be an empty data.frame now, if
    ## there was only one fraction with no neighbours
    if (!plyr::empty(alignQ_round))
    {
      alignQ_round$refFile = file_ref
      ## remove matches to self (i.e. with rtDiff==0)
      alignQ_round = alignQ_round[ alignQ_round$raw.file != file_ref, ]
      alignQ = rbind(alignQ, alignQ_round)
    }
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
ScoreInAlignWindow = function(data, allowed.deltaRT = 1)
{
  colnames(data) = tolower(colnames(data))
  
  if (!all(c('rtdiff', 'raw.file') %in% colnames(data)))
  {
    stop("alignmentCheck(): columns missing!")  
  }
  alignQC = plyr::ddply(data, "raw.file", function(x) {
    withinRT = sum(abs(x$rtdiff) < allowed.deltaRT, na.rm = TRUE) / sum(!is.na(x$rtdiff))
    return(data.frame(withinRT = withinRT))
  })
  
  return(alignQC)
}
  

#'
#' Check how close transferred ID's after alignment are to their genuine IDs within one Raw file.
#' 
#' The input is a data.frame containing feature evidence with corrected retention times,
#' e.g. a 'calibrated.retention.time' column.
#'
#' Note that this function must be given MS/MS identifications of type "MULTI-MSMS" and "MSMS-MATCH".
#' It will stop() otherwise.
#'  
#' We compare for each peptide sequence (and charge) the RT difference within groups of either genuine as well as mixed pairs.
#' For every comparison made, we report the RT span If alignment worked perfectly, the span are very small (<1 min),
#' for the mixed group, i.e. the pairs are accidentally split 3D peaks. Alignment performance has no influence on the
#' genuine-only groups.
#' 
#' Note: We found early MaxQuant versions (e.g. 1.2.2.5) to have an empty 'modified.sequence' column for 'MULTI-MATCH' entries.
#' The sequence which SHOULD be present is equal to the immediate upper row. This is what we use to guess the sequence.
#' However, this relies on the data.frame not being subsetted before (we can sort using the 'id' column)!
#' 
#' @param df_evd_all A data.frame with columns 'type', 'calibrated.retention.time', 'modified.sequence', 'charge', 'raw.file'
#' @return A data.frame containing the RT diff for each ID-group found in a Raw file (bg = genuine).
#'
idTransferCheck = function(df_evd_all) {
  colnames(df_evd_all) = tolower(colnames(df_evd_all))
  
  if (!checkInput(c('id', 'type', 'calibrated.retention.time', 'modified.sequence', 'charge', 'fc.raw.file'), df_evd_all)) return()
 
  
  if (!all(c("MULTI-MSMS", "MULTI-MATCH") %in% unique(df_evd_all$type)))
  {
    stop('idTransferCheck(): scan types missing! Required: "MULTI-MSMS" and "MULTI-MATCH".')  
  }
  
  
  df_evd_all$seq_charge = paste(factor(df_evd_all$modified.sequence), df_evd_all$charge, sep="_")
  alignQ = plyr::ddply(df_evd_all[,c("fc.raw.file", "type", "calibrated.retention.time", "seq_charge")],
                       "fc.raw.file",
                       function(x) {
    # unique(df_evd_all$fc.raw.file)
    # x = df_evd_all[ df_evd_all$fc.raw.file == "file 01", ]
    
    ## genuine groups only (within this Raw file):
    x_genuine = x[x$type=="MULTI-MSMS",]
    rt_diffs_genuine = plyr::ddply(x_genuine, "seq_charge", 
                             function(x2) {
                               if (nrow(x2)==1) return(NULL) ## we do not want singlets
                               return (data.frame(rtdiff_genuine = diff(range(x2$calibrated.retention.time))))
                             })
    # rt_diffs_genuine might be empty, if no oversampling is seen
    if (nrow(rt_diffs_genuine) == 0) rt_diffs_genuine = data.frame(seq_charge = character(0), rtdiff_genuine = numeric(0)) ## remove columns from empty DF; avoid merge()ing NA columns
    
    ## mixed class:
    ## retain only IDs which have at least one transferred ID
    x_mixed = x[x$seq_charge %in% x$seq_charge[x$type=="MULTI-MATCH"], ]
    if (nrow(x_mixed)>0) {    
      rt_diffs_mixed = plyr::ddply(x_mixed, "seq_charge", 
                             function(x2) {
                               if (nrow(x2)==1) return(NULL) ## we do not want singlets
                               return (data.frame(rtdiff_mixed = diff(range(x2$calibrated.retention.time))))
                             })
      ## rtdiff_mixed might be empty, if only singlets were transferred
      ## only merge if non-empty (otherwise the whole merge is empty)
      if (nrow(rt_diffs_mixed) > 0) {
        rt_diffs_genuine = merge(rt_diffs_genuine, rt_diffs_mixed, all = TRUE)
      }
    }
    return (rt_diffs_genuine)    
  })
  #head(alignQ)
  #hist(log(1+alignQ$rtdiff), 100)
  #hist(log(1+alignQ$rtdiff_bg), 100, add=TRUE, col="red")
  
  return(alignQ)
}


#'
#' For grouped peaks: separate them into in-width vs. out-width class.
#' 
#' Looking at groups only: Compute the fraction of 3D-peak pair groups per Raw file which 
#' have an acceptable RT difference after alignment using the result from 'idTransferCheck()',
#' i.e. compute the fraction of groups which are within a certain RT tolerance.
#' 
#' Returned value is between 0 (bad) and 1 (all within tolerance).
#' 
#' @param data A data.frame with columns 'fc.raw.file', 'rtdiff_mixed', 'rtdiff_genuine'
#' @param df.allowed.deltaRT The allowed matching difference for each Raw file (as data.frame(fc.rawfile, m))
#' @return A data.frame with one row for each raw.file and columns 'raw.file' and score 'withinRT' (0-1)
#'
inMatchWindow = function(data, df.allowed.deltaRT)
{
  ## 'data' columns of interest : 
  cols_diffs = c("rtdiff_mixed", "rtdiff_genuine")
  
  colnames(data) = tolower(colnames(data))
  
  if (!all(c(cols_diffs, 'fc.raw.file') %in% colnames(data)))
  { ## be error tolerant and return an empty data frame.
    return (data.frame(fc.raw.file = NA, withinRT_genuine = NA, withinRT_mixed = NA, withinRT_all = NA)[numeric(0), ])
  }
  
  alignQC = plyr::ddply(data, "fc.raw.file", function(x) {
    # x=data[ data$fc.raw.file=="file 01",]
    allowed.deltaRT = df.allowed.deltaRT$m[match(x$fc.raw.file[1], df.allowed.deltaRT$fc.raw.file)]
    
    withinRT_genuine = sum(abs(x$rtdiff_genuine) < allowed.deltaRT, na.rm = TRUE) / sum(!is.na(x$rtdiff_genuine))
    withinRT_mixed   = sum(abs(x$rtdiff_mixed) < allowed.deltaRT, na.rm = TRUE) / sum(!is.na(x$rtdiff_mixed))
    ## compute the RT-span if all evidence is considered 
    ##  --> simply the max of mixed and genuine, since each evidence must be in either of these two
    x$rtdiff_all = pmax(x$rtdiff_mixed, x$rtdiff_genuine, na.rm = TRUE)
    withinRT_all = sum(abs(x$rtdiff_all) < allowed.deltaRT, na.rm = TRUE) / sum(!is.na(x$rtdiff_all))
    return(data.frame(withinRT_genuine, withinRT_mixed, withinRT_all))
  })
  
  return(alignQC)
}


#'
#' Determine fraction of evidence which causes segmentation, i.e. sibling peaks at different RTs
#' confirmed either by genuine or transferred MS/MS.
#'
#' Sometimes, MQ splits a feature into 2 or more if the chromatograpic conditions are not optimal and there
#' is a drop in RT intensity.
#' If both features contain successful MS/MS scans, we will find the same peptide twice (with slightly different RT)
#' in the same charge state. This constitutes a natively split peak and is rare (95% of all genuine peaks are unique).
#' 
#' If Match-between-runs is used and the RT alignment is not perfect, then a peptide might be inferred at a wrong
#' RT position, even though this Raw file already contains MS/MS evidence of this peptide.
#' Usually the number of peak duplicates rises drastically (e.g. only 75% of peaks are unique after MBR was used).
#' In most cases, the RT is too far off to be a split peak. It's rather a lucky hit with accidentally the same mass-to-charge,
#' and thus the intensity is random.
#' To find by how much these peak pairs differ in RT, use idTransferCheck() and inMatchWindow().
#' 
#' Required columns are 'is.transferred', 'fc.raw.file', 'modified.sequence', 'charge', 'type'.
#'
#' Note that this function must be given MS/MS identifications of type "MULTI-MSMS" and "MSMS-MATCH".
#' It will stop() otherwise.
#' 
#' @param df_evd_all A data.frame of evidences containing the above columns
#' @return A data.frame with one row per Raw file and 
#'         three columns: 
#'           1) % of native single peaks (ignoring transferred IDs)
#'           2) % of single peaks (group of size=1) using only groups which have one transferred evidence
#'           3) % of single peaks using all groups
#'
peakSegmentation = function(df_evd_all)
{
  if (!checkInput(c("is.transferred", "fc.raw.file", "modified.sequence", "charge", 'type'), df_evd_all)) return()

  if (!all(c("MULTI-MSMS", "MULTI-MATCH") %in% unique(df_evd_all$type)))
  {
    stop('peakSegmentation(): scan types missing! Required: "MULTI-MSMS" and "MULTI-MATCH".')  
  }
  
  fc.raw.files = unique(df_evd_all$fc.raw.file)
  
  ## just keep "MULTI-MATCH" and "MULTI-MSMS", to keep results comparable to idTransferCheck()
  df_evd_all = df_evd_all[df_evd_all$type %in% c("MULTI-MSMS", "MULTI-MATCH"), ]

  cols = c("is.transferred", "fc.raw.file", "modified.sequence", "charge")
  countSeqs = plyr::ddply(df_evd_all[, cols], cols[-1], function(x)
  {
    return(data.frame(nNative = sum(!x$is.transferred), nMatched = sum(x$is.transferred)))#, ratio = ratio))
  })
  
  mbr_score = plyr::ddply(countSeqs, "fc.raw.file", function(countSeqs_sub)
  {
    #unique(countSeqs$fc.raw.file)
    #countSeqs_sub = countSeqs[countSeqs$fc.raw.file == "file 02", ]
    
    ddt = table(countSeqs_sub[, c("nMatched", "nNative")])
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
    ## add nMatched=1 row -- if already present, we just insert an empty row at the end (no harm done)
    ddt = rbind(ddt, 0)

    ## segmentation genuine (ignore the nMatched, i.e. project everything onto first row)
    n.perf = sum(ddt[,colnames(ddt)==1]) ## uniquely matched
    n.all = sum(ddt[,!colnames(ddt)==0]) ## all natives (==0 is the matched-only column)
    single.nat = max(0, n.perf/n.all, na.rm=TRUE) ## e.g. 94%, but n.all could be 0, yielding NaN
    
    ## segmentation of matched-only
    i.perf = ddt[rownames(ddt)==1, colnames(ddt)==0] ## transferred singlets (nMatched=1,nNative=0)
    i.all = sum(ddt[rownames(ddt)!=0, ])             ## all groups which contain a match (==without first row)
    single.matched = max(0, i.perf) / i.all  ## e.g. 50%
    ## no ID's were transferred (div/0) -- this remains NaN
    
    ## segmentation of all evidence (weighed average of above)
    comb.perf = max(0, ddt[rownames(ddt)==1, colnames(ddt)==0]) + ## transferred singlets (nMatched=1,nNative=0) -- could be empty
                ddt[rownames(ddt)==0, colnames(ddt)==1]
    comb.all = sum(ddt)             ## all groups
    single.all = max(0, comb.perf) / comb.all       ## e.g. 87%
    
    ## final result
    data.frame(single.nat = single.nat, single.matched = single.matched, single.all = single.all)
  })
  
  ## some Raw files might be REALLY sparse and up to now, they would be missing
  empty_fc.raw.files = fc.raw.files[!(fc.raw.files %in% mbr_score$fc.raw.file)]
  if (length(empty_fc.raw.files) > 0)
  { ## ... add them with NA values
    mbr_score = merge(mbr_score, data.frame(fc.raw.file = empty_fc.raw.files), all = TRUE)   
  }
  
  return (mbr_score = mbr_score)
}

#'
#' Combine several data structs into a final picture for segmentation incurred by 'Match-between-runs'.
#' 
#' qMBRSeg_Dist_inGroup might be empty if there are only singlets (transferred and genuine), but then the scores will be pretty
#' boring as well (100%).
#'
#' @param qMBR                 A data.frame as computed by peakSegmentation()
#' @param qMBRSeg_Dist_inGroup A data.frame as computed by inMatchWindow()
#' @return A data.frame which details the distribution of singlets and pairs (inRT and outRT) for each Raw file and genuine vs. all
#' 
computeMatchRTFractions = function(qMBR, qMBRSeg_Dist_inGroup)
{
  ## data might look like this:
#     fc.raw.file single.nat  single.matched   single.all
#   1  file 01     0.9740162       0.3393057    0.9517918
#   2  file 02     0.9813251       0.7689088    0.9433782
#   3  file 03     0.9671893       0.6986601    0.9121866
# 
# qMBRSeg_Dist_inGroup
#      fc.raw.file withinRT_genuine withinRT_mixed withinRT_all
#   1      file 01        0.5000000            NaN    0.5000000
#   2      file 02        0.3589744      0.3076923    0.3750000
#   3      file 03        0.2068966      0.1764706    0.1851852
#   4      file 04        0.2500000      0.1034483    0.1818182
#   5      file 05        0.3125000      0.3333333    0.3333333
#   6      file 06        0.2500000      0.4090909    0.3529412
  
  ## compute percentage of outside dRT peaks in genuine, matched and combined(=all)
  ## then calc the drop.
  f = plyr::ddply(qMBR, "fc.raw.file", function(x) {
    #x = qMBR[3, , drop = FALSE]
    rr = qMBRSeg_Dist_inGroup$fc.raw.file==x$fc.raw.file
    
    nat.inRT = qMBRSeg_Dist_inGroup$withinRT_genuine[rr] ## e.g. 0.87
    if (length(nat.inRT) == 0) nat.inRT = NA
    x$single.nat.inRT = (1-x$single.nat)*nat.inRT
    x$single.nat.outRT = (1-x$single.nat)*(1-nat.inRT)

    matched.inRT = qMBRSeg_Dist_inGroup$withinRT_mixed[rr]
    if (length(matched.inRT) == 0) matched.inRT = NA
    x$single.matched.inRT = (1-x$single.matched)*matched.inRT
    x$single.matched.outRT = (1-x$single.matched)*(1-matched.inRT)
    
    combined.inRT = qMBRSeg_Dist_inGroup$withinRT_all[rr]
    if (length(combined.inRT) == 0) combined.inRT = NA
    x$single.all.inRT = (1-x$single.all)*combined.inRT
    x$single.all.outRT = (1-x$single.all)*(1-combined.inRT)
    
    r = data.frame(fc.raw.file=x$fc.raw.file, 
                   single=c(x$single.nat, x$single.matched, x$single.all),
                   multi.inRT = c(x$single.nat.inRT, x$single.matched.inRT, x$single.all.inRT),
                   multi.outRT = c(x$single.nat.outRT, x$single.matched.outRT, x$single.all.outRT),
                   sample = factor(c("genuine", "transferred", "all"), levels=c("genuine", "transferred", "all", ordered = TRUE)))
    return(r)    
  })
  return(f)
}

#'
#' Return a tree plot with a possible alignment tree.
#' 
#' This allows the user to judge which Raw files have similar corrected RT's (i.e. where aligned successfully).
#' If there are clear sub-clusters, it might be worth introducing artifical fractions into MaxQuant,
#' to avoid ID-transfer between these clusters (use the MBR-Align and MBR-ID-Transfer metrics to support the decision).
#' 
#' If the input contains fractions, leaf nodes will be colored accordingly.
#' Distinct sub-clusters should have their own color.
#' If not, MaxQuant's fraction settings should be optimized.
#' Note that introducing fractions in MaxQuant will naturally lead to a clustering here (it's somewhat circular).
#' 
#' @param df_evd  Evidence table containing calibrated retention times and sequence information.
#' @param col_fraction Empty vector or 1-values vector giving the name of the fraction column (if existing)
#' @return ggplot object containing the correlation tree
#' 
#' @import ggplot2
#' @export
#'
RTalignmentTree = function(df_evd, col_fraction = c())
{
  #df_evd$fc.raw.file=df_evd$raw.file
  head(df_evd)
  
  req_cols = c("calibrated.retention.time", "fc.raw.file", col_fraction, "modified.sequence", "charge")
  if (!all(req_cols %in% colnames(df_evd)))
  {
    stop("RTalignmentTree: Missing columns! Please fix the code: ", 
         setdiff(req_cols, colnames(df_evd)), "!")
  }
  
  d_cast = reshape2::dcast(df_evd, modified.sequence + charge ~ fc.raw.file, mean, value.var = "calibrated.retention.time")
  
  head(d_cast[,-(1:2)])
  d_cast.m = as.matrix(d_cast[,-(1:2)])
  head(d_cast.m)
  
  #Dissimilarity = 1 - Correlation
  #Dissimilarity = (1 - Correlation)/2
  #Dissimilarity = 1 - Abs(Correlation)
  #Dissimilarity = Sqrt(1 - Correlation2)
  dissimilarity = (1 - abs(cor(d_cast.m, use="pairwise.complete.obs")))
  ## if some samples have no overlap, their cell is NA --> set to 1 (max distance)
  dissimilarity[is.na(dissimilarity)] = 1
  #plot(hclust(as.dist(dissimilarity), method="ward.D"))
  ddata = ggdendro::dendro_data(hclust(as.dist(dissimilarity)), type = "rectangle")
  
  if (length(col_fraction))
  {
    idx_raw = match(ddata$labels$label, df_evd$fc.raw.file)
    ddata$labels$col = factor(df_evd[idx_raw, col_fraction])
  } else {
    ddata$labels$col = "black"
  }
  p = ggplot(ggdendro::segment(ddata)) + 
      geom_segment(aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend)) +
      scale_x_continuous(breaks = ddata$labels$x, labels = ddata$labels$label) +
      theme_blank() +
      theme(axis.text.y = element_text(colour = ddata$labels$col), ## is a vector not officially supported...
                     axis.text.x = element_blank()) +
      coord_flip() +
      ggtitle("[experimental] EVD: Clustering Tree of Raw files", "by Correlation of Corrected Retention Times")
  #p
  return(p)
}

