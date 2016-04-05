
qcMetric_EVD_UserContaminant =  setRefClass(
  "qcMetric_EVD_UserContaminant",
  contains = "qcMetric",
  methods = list(initialize=function() {  callSuper(  
    helpTextTemplate = 
      "User defined contaminant search (usually used for Mycoplasma detection, but can be used for an arbitrary (set of) proteins.

Two abundance measures are computed per Raw file:
  - fraction of intensity (used for score)
  - fraction of spectral counts (as comparison for user)

An additional plot with peptide score distributions will be shown if the threshold was reached (i.e. suspected contamination).
This allows to decide if the contamination is true, i.e. achieves good MS/MS Andromeda scores (requires a recent MQ version).

Heatmap score [EVD: Contaminant ...]: boolean, i.e. 0% (fail) if the threshold was reached. 100%% (pass) otherwise.
",
    workerFcn = function(.self, df_evd, df_pg, lst_contaminants)
    {
      ## completeness check
      stopifnot(c("id", "fasta.headers") %in% colnames(df_pg))
      ## "score" might not be present (e.g. missing in MQ 1.0.13.13)
      stopifnot(c("protein.group.ids", "type", "intensity", "fc.raw.file") %in% colnames(df_evd))

      local_qcScores = data.frame()
      
      lpl = list()
      
      ca_entry = lst_contaminants[[1]]
      for (ca_entry in lst_contaminants)
      {
        ca = ca_entry[1]
        ## 
        if (ca == FALSE) {
          cat("No special contaminants requested!\n")
          break;
        }
        
        ca_thresh = as.numeric(ca_entry[2])
        
        not_found = TRUE
        
        pg_id = df_pg$id[grep(ca, df_pg$fasta.headers, ignore.case = TRUE)]
        
        if (length(pg_id) > 0)
        {
          
          ## we might or might not have found something... we plot it anyways, so the user can be sure that we searched for it
          
          ## find peptides which only have one group (ignoring razor peptides where we cannot be sure)
          evd_uniqueGroup = !grepl(";", df_evd$protein.group.ids)
          ## do not trust MBR here. We want real evidence!
          evd_realMS = !grepl("MATCH", df_evd$type)
          ## for each Raw file: find unique peptides of our contaminant
          cont_data.l = dlply(df_evd[evd_uniqueGroup & evd_realMS, ], "fc.raw.file",
                              function(x) {
                                if (length(grep(";", x$protein.group.ids))) stop("more than one proteinGroup for supposedly unique peptide...")
                                
                                x$idx_cont = x$protein.group.ids %in% pg_id
                                
                                sc = sum(x$idx_cont) / nrow(x) * 100
                                int = sum(as.numeric(x$intensity[x$idx_cont]), na.rm = TRUE) / sum(as.numeric(x$intensity), na.rm = TRUE) * 100
                                
                                above.thresh = (sc > ca_thresh) | (int > ca_thresh)
                                cont_scoreECDF = NULL;
                                if ("score" %in% colnames(x)) {
                                  cont_scoreECDF = ddply(x, "idx_cont", function(xx) {
                                    if (length(unique(xx$score)) < 2) return(NULL) ## not enough data for ECDF
                                    r = getECDF(xx$score)
                                    r$condition = c("sample", "contaminant")[xx$idx_cont[1]+1]
                                    return(r)
                                  })
                                }
                                if (!any(x$idx_cont)){
                                  ks_p = NA
                                } else { ## no contaminant peptide 
                                  ks_p = suppressWarnings(  ## brags about '-value will be approximate in the presence of ties'
                                    ks.test(x$score[x$idx_cont], x$score[!x$idx_cont], alternative = "greater")$p.value
                                  )
                                }
                                return (list(cont_data = data.frame(spectralCount = sc, intensity = int,
                                                                    above.thresh = above.thresh, fc.raw.file = x$fc.raw.file[1],
                                                                    score_KS = ks_p),
                                             cont_scoreECDF = cont_scoreECDF))
                              })
          head(cont_data.l)
          
          ## melt
          cont_data = ldply(cont_data.l, function(l) { l$cont_data })
          cont_data.long = melt(cont_data, id.vars="fc.raw.file")
          
          not_found = all(cont_data.long$value[cont_data.long$variable == "above.thresh"] == FALSE)
        }
        
        if (not_found)
        { ## identifier was not found in any sample
          pl_cont = ggText("PG: Contaminants",
                           paste0("Contaminant '", ca, "' was not found in any sample.\n\nDid you use the correct database?"),
                           "red")
          lpl = append(lpl, list(pl_cont))
        } else {
          ## plot User-Contaminants
          lpl_i = byXflex(data = cont_data.long, indices = cont_data.long$fc.raw.file, subset_size = 120, 
                          FUN = plot_ContUser, sort_indices = TRUE, name_contaminant = ca, extra_limit = ca_thresh)
          lpl = append(lpl, lpl_i)
          
          ## plot Andromeda score distribution of contaminant vs. sample
          llply(cont_data.l, function(l)
          {
            if (l$cont_data$above.thresh == FALSE ||
                is.null(l$cont_scoreECDF))
            {
              return(NULL)
            } 
            p = plot_ContUserScore(l$cont_scoreECDF, l$cont_data$fc.raw.file, l$cont_data$score_KS)
            lpl = append(lpl, list(p))
            #print(p)
            return(NULL)
          })
          
          ## add heatmap column
          cname = sprintf(.self$qcName, ca)
          cont_data[,cname] = as.numeric(!cont_data$above.thresh) ## inverse (0 is 'bad')
          
          qcScore = cont_data[, c("fc.raw.file", cname)]
          if (ncol(local_qcScores) == 0){
            local_qcScores = qcScore
          } else {
            local_qcScores = merge(local_qcScores, qcScore)
          }
          
        }
      } ## contaminant loop
      
      return(list(plots = lpl, qcScores = local_qcScores))
    }, 
    qcCat = "Prep", 
    qcName = "EVD:Contaminant~(%s)", 
    orderNr = 0020
  )
    return(.self)
  })
)


#####################################################################

qcMetric_EVD_PeptideInt =  setRefClass(
  "qcMetric_EVD_PeptideInt",
  contains = "qcMetric",
  methods = list(initialize=function() {  callSuper(    
    helpTextTemplate = 
      "Peptide intensity per Raw file from evidence.txt.
Failing to reach the intensity threshold is usually due to unfavorable column conditions, inadequate 
column loading or ionization issues. If the study is not a dilution series or pulsed SILAC experiment, we 
would expect every condition to have about the same median log-intensity (of 2<sup>%1.1f</sup>).
The relative standard deviation (RSD) gives an indication about reproducibility across files and should be below 5%%.

Depending on your setup, your thresholds might vary from PTXQC's defaults.
Change the threshold using the YAML configuration file.

Heatmap score [EVD: Pep Intensity (>%1.1f)]: 
  Linear scale of the median intensity reaching the threshold, i.e. reaching 2<sup>21</sup> of 2<sup>23</sup> gives score 0.25.
",
    workerFcn = function(.self, df_evd, thresh_intensity)
    {
      ## completeness check
      stopifnot(c("fc.raw.file", "intensity") %in% colnames(df_evd))
      
      ## update helpText
      .self$helpText = sprintf(.self$helpTextTemplate, thresh_intensity, thresh_intensity)
      
      medians_pep = ddply(df_evd[ , c("fc.raw.file", "intensity")], "fc.raw.file",
                          function(x) data.frame(med = log2(quantile(x$intensity, probs=0.5, na.rm = TRUE))))
      
      int_dev_pep = RSD((medians_pep$med))
      int_dev.s = pastet("INT RSD [%]", round(int_dev_pep, 3))
      lpl = boxplotCompare(data = df_evd[, c("fc.raw.file", "intensity", "contaminant")],
                           log2 = TRUE, 
                           mainlab="EVD: peptide intensity distribution",
                           ylab = expression(log[2]*" intensity"),
                           sublab=paste0("RSD ", round(int_dev_pep, 1),"% (expected < 5%)\n"),
                           abline = thresh_intensity)
      #for (pl in lpl) print(pl)
      
      ## QC measure for peptide intensity
      qc_pepint = medians_pep
      cname = sprintf(.self$qcName, thresh_intensity)
      qc_pepint[,cname] = qualLinThresh(2^qc_pepint$med, 2^thresh_intensity) ## use non-log space 
      qcScore = qc_pepint[, c("fc.raw.file", cname)]
      
      return(list(plots = lpl, qcScores = qcScore))
    }, 
    qcCat = "prep", 
    qcName = "EVD:~Pep~Intensity~(\">%1.1f\")", 
    orderNr = 0030
  )
    return(.self)
  })
)



#####################################################################

qcMetric_EVD_ProteinCount =  setRefClass(
  "qcMetric_EVD_ProteinCount",
  contains = "qcMetric",
  methods = list(initialize=function() {  callSuper(  
    helpTextTemplate = 
      "Number of Protein groups (after FDR) per Raw file. A configurable target threshold is indicated as dashed line.

If MBR was enabled, three categories ('genuine (exclusive)', 'genuine + transferred', 'transferred (exclusive)'
are shown, so the user can judge the gain that MBR provides. If the gain is low and the MBR scores are bad (see other metrics),
MBR should be switched off for the Raw files which are affected (could be a few or all). 

Heatmap score [EVD: Prot Count (>%1.0f)]: Reaching the target threshold gives a score of 100%%.
",
    workerFcn = function(.self, df_evd, thresh_protCount)
    {
      ## completeness check
      stopifnot(c("fc.raw.file", "protein.group.ids", "match.time.difference") %in% colnames(df_evd))
      
      .self$helpText = sprintf(.self$helpTextTemplate, thresh_protCount)
      
      protC = getProteinCounts(df_evd[, c("fc.raw.file", "protein.group.ids", "match.time.difference")])
      protC$block = factor(assignBlocks(protC$fc.raw.file, 30))
      
      max_prot = max(unlist(dlply(protC, "fc.raw.file", function(x) sum(x$counts))))
      ## average gain in percent
      reportMTD = any(!is.na(df_evd$match.time.difference))
      gain_text = ifelse(reportMTD, sprintf("MBR gain: +%.0f%%", mean(protC$MBRgain, na.rm = TRUE)), "")
      
      lpl = dlply(protC, "block", .fun = function(x)
      {
        p = plot_CountData(data = x, 
                           y_max = max(thresh_protCount, max_prot)*1.1,
                           thresh_line = thresh_protCount,
                           title = c("EVD: ProteinGroups count", gain_text))
        #print(p)
        return (p)
      })
      
      ## QC measure for protein ID performance
      qc_protc = ddply(protC, "fc.raw.file", function(x){
        if (nrow(x) == 3 && length(grep("^genuine", x$category))!= 2){
          stop("expected two categories to start with 'genuine...'")
        }
        r = data.frame(genuineAll = sum(x$counts[grep("^genuine", x$category)]))
        return (r)
      })
      cname = sprintf(.self$qcName, thresh_protCount)
      qc_protc[,cname] = qualLinThresh(qc_protc$genuineAll, thresh_protCount)
      qcScore = qc_protc[, c("fc.raw.file", cname)]
      
      return(list(plots = lpl, qcScores = qcScore))
    }, 
    qcCat = 'general', 
    qcName = "EVD:~Prot~Count~(\">%1.0f\")", 
    orderNr = 0450
  )
    return(.self)
  })
)


#####################################################################

qcMetric_EVD_PeptideCount =  setRefClass(
  "qcMetric_EVD_PeptideCount",
  contains = "qcMetric",
  methods = list(initialize=function() {  callSuper(  
    helpTextTemplate = 
      "Number of peptides (after FDR) per Raw file. A configurable target threshold is indicated as dashed line.

If MBR was enabled, three categories ('genuine (exclusive)', 'genuine + transferred', 'transferred (exclusive)'
are shown, so the user can judge the gain that MBR provides. If the gain is low and the MBR scores are bad (see other metrics),
MBR should be switched off for the Raw files which are affected (could be a few or all). 

Heatmap score [EVD: Pep Count (>%1.0f)]: Reaching the target threshold gives a score of 100%%.
",
    workerFcn = function(.self, df_evd, thresh_pepCount)
    {
      ## completeness check
      stopifnot(c("fc.raw.file", "modified.sequence", "match.time.difference") %in% colnames(df_evd))
      
      .self$helpText = sprintf(.self$helpTextTemplate, thresh_pepCount)
      
      pepC = getPeptideCounts(df_evd[, c("fc.raw.file", "modified.sequence", "match.time.difference")])
      pepC$block = factor(assignBlocks(pepC$fc.raw.file, 30))
      
      max_pep = max(unlist(dlply(pepC, "fc.raw.file", function(x) sum(x$counts))))
      ## average gain in percent
      reportMTD = any(!is.na(df_evd$match.time.difference))
      gain_text = ifelse(reportMTD, sprintf("MBR gain: +%.0f%%", mean(pepC$MBRgain, na.rm = TRUE)), "")
      
      lpl = dlply(pepC, "block", .fun = function(x)
      {
        p = plot_CountData(data = x, 
                           y_max = max(thresh_pepCount, max_pep)*1.1,
                           thresh_line = thresh_pepCount,
                           title = c("EVD: Peptide ID count", gain_text))
        #print(p)
        return (p)
      })
      
      ## QC measure for peptide ID performance
      qc_pepc = ddply(pepC, "fc.raw.file", function(x){
        if (nrow(x) == 3 && length(grep("^genuine", x$category))!= 2){
          stop("expected two categories to start with 'genuine...'")
        }
        r = data.frame(genuineAll = sum(x$counts[grep("^genuine", x$category)]))
        return (r)
      })
      cname = sprintf(.self$qcName, thresh_pepCount)
      qc_pepc[,cname] = qualLinThresh(qc_pepc$genuineAll, thresh_pepCount)
      qcScore = qc_pepc[, c("fc.raw.file", cname)]
      
      return(list(plots = lpl, qcScores = qcScore))
    }, 
    qcCat = 'general', 
    qcName = "EVD:~Pep~Count~(\">%1.0f\")", 
    orderNr = 0400
  )
    return(.self)
  })
)  


#####################################################################

qcMetric_EVD_RTPeakWidth =  setRefClass(
  "qcMetric_EVD_RTPeakWidth",
  contains = "qcMetric",
  methods = list(initialize=function() {  callSuper(    
    helpTextTemplate = 
      "One parameter of optimal and reproducible chromatographic separation is the distribution of widths of 
peptide elution peaks, derived from the evidence table. Ideally, all Raw files show a similar 
distribution, e.g. to allow for equal conditions during dynamic precursor exclusion, RT alignment or 
peptide quantification. 
    
Heatmap score [EVD: RT Peak Width]: Scored using BestKS function, i.e. the D statistic of a Kolmogoriv-Smirnoff test.
",
    workerFcn = function(.self, df_evd)
    {
      ## completeness check
      stopifnot(c("retention.time", "retention.length", "fc.raw.file") %in% colnames(df_evd))
      
      ## compute some summary stats before passing data to ggplot (performance issue for large experiments) 
      df_evd.m.d = ddply(df_evd[,c("retention.time", "retention.length", "fc.raw.file")], "fc.raw.file", .fun = peakWidthOverTime)
      head(df_evd.m.d)
      ## median peak width
      df_evd.m.d_avg = ddply(df_evd[,c("retention.length","fc.raw.file")], "fc.raw.file", .fun = function(x) {
        #fcr = as.character(x$fc.raw.file[1])
        #cat(fcr)
        m = median(x$retention.length, na.rm = TRUE);
        return(data.frame(median = m))
      })
      df_evd.m.d_avg$fc.raw.file_aug = paste0(df_evd.m.d_avg$fc.raw.file, " (~", round(df_evd.m.d_avg$median, 1)," min)")
      .self$outData[["avg_peak_width"]] = df_evd.m.d_avg
      
      ## augment Raw filename with avg. RT peak width
      df_evd.m.d$fc.raw.file = mapvalues(df_evd.m.d$fc.raw.file, df_evd.m.d_avg$fc.raw.file, df_evd.m.d_avg$fc.raw.file_aug)
      df_evd.m.d$block = factor(assignBlocks(df_evd.m.d$fc.raw.file, 6)) ## color set is 9, so do not increase this (6*150%)
      ## identical limits for all plots
      df_evd.xlim = range(df_evd.m.d$RT, na.rm = TRUE)
      ## ignore top peaks, since they are usually early non-peptide eluents
      df_evd.ylim = c(0, quantile(df_evd.m.d$peakWidth, 0.99, na.rm = TRUE))
      
      ## plot peak width
      lpl = list()
      for (bl in unique(df_evd.m.d$block))
      { ## needs to be within a function, otherwise rep_data$add and print() somehow have delayed eval's which confused ggplot...
        lpl[[bl]] = plot_RTPeakWidth(data = df_evd.m.d[df_evd.m.d$block==bl,], x_lim = df_evd.xlim, y_lim = df_evd.ylim)
      }
      
      ## QC measure for reproducibility of peak shape
      ##.. create a list of distributions
      l_dists = dlply(df_evd[,c("retention.length", "fc.raw.file")], "fc.raw.file", function(x) return(x$retention.length))
      qc_evd_PeakShape = qualBestKS(l_dists)
      colnames(qc_evd_PeakShape) = c("fc.raw.file", .self$qcName)
      
      return(list(plots = lpl, qcScores = qc_evd_PeakShape))
    }, 
    qcCat = "LC", 
    qcName = "EVD:~RT~Peak~Width", 
    orderNr = 0170
  )
    return(.self)
  })
)  


#####################################################################

qcMetric_EVD_MBRAlign =  setRefClass(
  "qcMetric_EVD_MBRAlign",
  contains = "qcMetric",
  methods = list(initialize=function() {  callSuper(    
    helpTextTemplate = 
      "MBR Alignment: First of two steps (1=align, 2=transfer) during Match-between-runs.
This plot is based purely on real MS/MS ids. Ideally, RTs of identical peptides should be equal (i.e. very small residual RT delta)
across Raw files after alignment.

MaxQuants RT correction is shown in blue -- it should be well within the alignment search window (20min by default) set
during MaxQuant configuration.
The resulting residual RT delta after RT alignment (compared to a reference Raw file), is shown as green/red dots. One dot represents
one peptide (incl. charge). Every dot (peptide) outside an allowed residual delta RT (1min by default) is colored red.
All others are green.

If moving 'red' dots to the horizontal zero-line (to make them green) requires large RT shifts, then increasing the alignment
search window might help MaxQuant to find a better alignment.

Heatmap score [EVD: MBR Align]: fraction of 'green' vs. 'green+red' peptides.
",
    workerFcn = function(.self, df_evd, tolerance_matching, raw_file_mapping)
    {
      ## completeness check
      stopifnot(c("type", "calibrated.retention.time", "id", "raw.file", "modified.sequence", "charge") %in% colnames(df_evd))
      
      ## find reference
      if (('fraction' %in% colnames(df_evd)) && (length(unique(df_evd$fraction)) > 1)) {
        ## fractions: there must be more than one, otherwise MQ will treat the samples as unfractionated
        refRaw = NA
        col_fraction = "fraction"
        txt_subtitle = "fraction: neighbour comparison"
        evd_has_fractions = TRUE
        df_evd$fraction[is.na(df_evd$fraction)] = 32000
      } else {
        refRaw = findAlignReference(df_evd)
        col_fraction = c()
        txt_subtitle = paste("alignment reference:", refRaw)
        evd_has_fractions = FALSE
      }
      
      lpl = list()
      qcScore = .self$qcScores
      
      if (length(refRaw) != 1) {
        lpl[[1]] = ggText("EVD: Alignment check", paste0("Cannot find a unique reference Raw file (files: ", paste(refRaw, collapse=", "), ")"))
      } else {
        ## find RT curve based on genuine 3D peaks (should be flat)
        d_alignQ = alignmentCheck(df_evd[(df_evd$type %in% c("MULTI-MSMS")), 
                                         c("calibrated.retention.time", 
                                           "id", "raw.file", col_fraction, "modified.sequence", "charge")], 
                                  referenceFile = refRaw)
        ## augment more columns
        d_alignQ$retention.time.calibration = df_evd$retention.time.calibration[match(d_alignQ$id, df_evd$id)]
        
        if (nrow(d_alignQ)==0)
        { ## very unusual case: reference contains no evidence -- e.g. pull-down experiment
          lpl[[1]] = ggText("EVD: RT Distance of peptides from reference after alignment", "Alignment cannot be verfied -- no data.")
        } else {
          ## filter data (reduce PDF file size)
          evd_RT_t = thinOutBatch(d_alignQ,
                                  "calibrated.retention.time",
                                  "raw.file")
          
          evd_RT_t$fc.raw.file = renameFile(evd_RT_t$raw.file, raw_file_mapping)
          
          ## QC measure for alignment quality
          ## compute % of matches within matching boundary (1 min by default)
          qcAlign = ScoreInAlignWindow(d_alignQ, tolerance_matching)
          if (!is.na(refRaw)) { ## rescue reference file (it will not show up in fraction-less data, and would otherwise be scored 'red')
            qcAlign = rbind(qcAlign, data.frame(raw.file=refRaw, withinRT=1))
          }
          qcAlign[, .self$qcName] = qcAlign$withinRT
          qcScore = qcAlign[, c("raw.file", .self$qcName)]
          
          qcAlign$fc.raw.file = renameFile(qcAlign$raw.file, raw_file_mapping)
          qcAlign$newlabel = qcAlign$fc.raw.file
          if (evd_has_fractions)
          { ## amend fc.raw.file with fraction number
            qcAlign$fraction = df_evd$fraction[match(qcAlign$fc.raw.file, df_evd$fc.raw.file)]
            qcAlign$newlabel = paste0(qcAlign$fc.raw.file, " - frc", qcAlign$fraction)
          }
          ## amend fc.raw.file with % good ID pairs
          qcAlign$newlabel = paste0(qcAlign$newlabel, " (sc: ", round(qcAlign$withinRT*100), "%)")
          evd_RT_t$fc.raw.file_ext = mapvalues(evd_RT_t$fc.raw.file, qcAlign$fc.raw.file, qcAlign$newlabel)
          
          evd_RT_t$RTdiff_in = c("green", "red")[(abs(evd_RT_t$rtdiff) > tolerance_matching)+1]
          
          ## plot alignment result
          y_lim = quantile(c(evd_RT_t$rtdiff, evd_RT_t$retention.time.calibration), probs = c(0.01,0.99), na.rm = TRUE) * 1.1
          lpl =
            byX(evd_RT_t, evd_RT_t$fc.raw.file, 3*3, plot_MBRAlign, sort_indices = FALSE, 
                y_lim = y_lim, title_sub = txt_subtitle, match_tol = tolerance_matching)
          
        } ## no data
      } ## ambigous reference file
      
      
      return(list(plots = lpl, qcScores = qcScore))
    }, 
    qcCat = "LC", 
    qcName = "EVD:~MBR~Align", 
    orderNr = 0210
  )
    return(.self)
  })
)


#####################################################################

qcMetric_EVD_MBRIdTransfer =  setRefClass(
  "qcMetric_EVD_MBRIdTransfer",
  contains = "qcMetric",
  methods = list(initialize=function() {  callSuper(    
    helpTextTemplate = 
      "MBR Transfer: Last of two steps (1=align, 2=transfer) during Match-between-runs.
If MaxQuant only transfers peptide ID's which are not present in the target file, then
each Raw file should not have any duplicates of identical peptides (incl. charge).
Sometimes, a single or split 3D-peak gets annotated multiple times, that's ok. However, the same peptide should not
be annotated twice (or more) at vastly different points in RT.

This plot shows three columns:
 - left: the 'genuine' situation (pretending that no MBR was computed)
 - middle: looking only at transferred IDs
 - right: combined picture (a mixture of left+middle, usually)

Each peptide falls into three categories (the colors):
 - single (good, because it has either one genuine OR a transferred ID).
 - in-group (also good, because all ID's are very close in RT)
 - out-group (bad, spread across the RT gradient -- should not be possible; a false ID)

Heatmap score [EVD: MBR ID-Transfer]: The fraction of non-out-group peptides (i.e. good peptides) in the middle column.
This score is 'pessimistic' because if few ID's were transferred, but all of them are bad, the score is bad, even though
the majority of peptides is still ok (because they are genuine). However, in this case MBR
provides few (and wrong) additional information, and should be disabled.
",
    workerFcn = function(.self, df_evd, avg_peak_width)
    {
      ## completeness check
      #stopifnot(c("...") %in% colnames(df_evd))
      
      ## increase of segmentation by MBR:
      ## three values returned: single peaks(%) in genuine, transferred and all(combined)
      qMBR = peakSegmentation(df_evd)
      head(qMBR)
      ## for groups: get their RT-spans
      ## ... genuine ID's only (as 'rtdiff_genuine') 
      ##  or genuine+transferred (as 'rtdiff_mixed'))
      ## Could be empty (i.e. no groups, just singlets) if data is really sparse ..
      qMBRSeg_Dist = idTransferCheck(df_evd)
      #head(qMBRSeg_Dist)
      #head(qMBRSeg_Dist[qMBRSeg_Dist$fc.raw.file=="file 13",])
      
      
      ## Check which fraction of ID-pairs belong to the 'in-width' group.
      ## The allowed RT delta is given in 'avg_peak_width' (estimated from global peak width for each file)
      qMBRSeg_Dist_inGroup = inMatchWindow(qMBRSeg_Dist, df.allowed.deltaRT = avg_peak_width)
      ## puzzle together final picture
      scoreMBRMatch = computeMatchRTFractions(qMBR, qMBRSeg_Dist_inGroup)
      #head(scoreMBRMatch)
      #scoreMBRMatch[scoreMBRMatch$fc.raw.file=="file 3",]
      
      ## plot ID-transfer
      lpl =
        byX(scoreMBRMatch, scoreMBRMatch$fc.raw.file, 12, plot_MBRIDtransfer, sort_indices = FALSE)
      
      ##
      ##  Quality
      ##
      qualMBR.m = merge(scoreMBRMatch[scoreMBRMatch$sample=="genuine",], 
                        scoreMBRMatch[scoreMBRMatch$sample=="transferred",], by="fc.raw.file")
      qualMBR.m = merge(qualMBR.m, scoreMBRMatch[scoreMBRMatch$sample=="all",], by="fc.raw.file")
      cname = .self$qcName
      qualMBR.m[, cname] = 1 - qualMBR.m$multi.outRT.y # could be NaN if: no-transfer at all, or: no groups but only singlets transferred
      qualMBR.m[is.na(qualMBR.m$multi.outRT.y) & !is.na(qualMBR.m$single.y), cname] = 1 ## only singlets transferred, wow...
      qualMBR.m[is.na(qualMBR.m[, cname]), cname] = HEATMAP_NA_VALUE
      qcScore = qualMBR.m[, c("fc.raw.file", cname)]
      
      return(list(plots = lpl, qcScores = qcScore))
    }, 
    qcCat = "LC", 
    qcName = "EVD:~MBR~ID-Transfer", 
    orderNr = 0220
  )
    return(.self)
  })
)  


#####################################################################

qcMetric_EVD_MBRaux =  setRefClass(
  "qcMetric_EVD_MBRaux",
  contains = "qcMetric",
  methods = list(initialize=function() {  callSuper(    
    helpTextTemplate = 
      "Auxililiary plots -- experimental -- without scores.

Heatmap score: none.
",
    workerFcn = function(.self, df_evd)
    {
      ## completeness check
      stopifnot(c("type", "match.time.difference", "calibrated.retention.time", "fc.raw.file", "modified.sequence", "charge") %in% colnames(df_evd))
      
      if (('fraction' %in% colnames(df_evd)) && (length(unique(df_evd$fraction)) > 1)) {
        ## fractions: there must be more than one, otherwise MQ will treat the samples as unfractionated
        col_fraction = "fraction"
      } else {
        col_fraction = c()
      }
      
      lpl = list()
      lpl[["tree"]] =
        RTalignmentTree(df_evd[(df_evd$type %in% c("MULTI-MSMS")), 
                               c("calibrated.retention.time", "fc.raw.file", col_fraction, "modified.sequence", "charge")],
                        col_fraction = col_fraction)
      
      ## MBR: additional evidence by matching MS1 by AMT across files
      if (any(!is.na(df_evd$match.time.difference))) {
        ## gain for each raw file: absolute gain, and percent gain
        mtr.df = ddply(df_evd, "fc.raw.file", function(x) {
          match_count_abs = sum(!is.na(x$match.time.difference))
          match_count_pc  = round(100*match_count_abs/(nrow(x)-match_count_abs)) ## newIDs / oldIDs
          return (data.frame(abs = match_count_abs, pc = match_count_pc))
        })
        lpl[["gain"]] =
          plot_MBRgain(data = mtr.df, title_sub = "")
      }
      
      return(list(plots = lpl))
    }, 
    qcCat = "LC", 
    qcName = "EVD:~MBR~auxilliary", 
    orderNr = 0221
  )
    return(.self)
  })
)  



#####################################################################

qcMetric_EVD_Charge =  setRefClass(
  "qcMetric_EVD_Charge",
  contains = "qcMetric",
  methods = list(initialize=function() {  callSuper(    
    helpTextTemplate = 
      "Charge distribution per Raw file. Should be dominated by charge 2 
(one N-terminal and one at tryptic C-terminal R or K residue) and have the same fraction in each Raw file.
Consistent charge distribution is paramount for comparable 3D-peak intensities across samples.

Heatmap score [EVD: Charge]: Deviation of the charge 2 proportion from a representative Raw file ('qualMedianDist' function).
",
    workerFcn = function(.self, df_evd, int_cols, MAP_pg_groups)
    {
      ## completeness check
      stopifnot(c("hasMTD", "fc.raw.file", "charge") %in% colnames(df_evd))
      
      d_charge = mosaicize(df_evd[!df_evd$hasMTD, c("fc.raw.file", "charge")])
      lpl =
        byXflex(d_charge, d_charge$Var1, 30, plot_Charge, sort_indices = TRUE)
      
      ## QC measure for charge centeredness
      qc_charge = ddply(df_evd[!df_evd$hasMTD, c("charge",  "fc.raw.file")], "fc.raw.file", function(x) data.frame(c = (sum(x$charge==2)/nrow(x))))
      qc_charge[, .self$qcName] = qualMedianDist(qc_charge$c)
      
      return(list(plots = lpl, qcScores = qc_charge[, c("fc.raw.file", .self$qcName)]))
    }, 
    qcCat = "prep", 
    qcName = "EVD:~Charge", 
    orderNr = 0100
  )
    return(.self)
  })
)  


#####################################################################

qcMetric_EVD_IDoverRT =  setRefClass(
  "qcMetric_EVD_IDoverRT",
  contains = "qcMetric",
  methods = list(initialize=function() {  callSuper(  
    helpTextTemplate = 
      "Judge column occupancy over retention time. 
Ideally, the LC gradient is chosen such that the number of identifications (here, after FDR filtering) is 
uniform over time, to ensure consistent instrument duty cycles. Sharp peaks and uneven distribution of 
identifications over time indicate potential for LC gradient optimization. 
See [http://www.ncbi.nlm.nih.gov/pubmed/24700534](Moruz et al., GradientOptimizer: An open-source graphical environment for calculating optimized gradients in reversed-phase
liquid chromatography, Proteomics, 06/2014; 14) for details.

Heatmap score [EVD: ID rate over RT]: Scored using 'Uniform' scoring function, i.e. constant receives good score, extreme shapes are bad.
",
    workerFcn = function(.self, df_evd)
    {
      ## completeness check
      stopifnot(c("retention.time", "fc.raw.file") %in% colnames(df_evd))
      
      raws_perPlot = 6
      
      rt_range = range(df_evd$retention.time, na.rm = TRUE)
      df_idRT = ddply(df_evd, "fc.raw.file", function(x) {
        h = hist(x$retention.time, breaks=seq(from=rt_range[1]-3, to=rt_range[2]+3, by=3), plot = FALSE)
        return(data.frame(RT = h$mid, counts = h$counts))
      })
      lpl =
        byXflex(df_idRT, df_idRT$fc.raw.file, raws_perPlot, plot_IDsOverRT, sort_indices = TRUE)
      
      ## QC measure for uniform-ness
      qcScore = ddply(df_evd[, c("retention.time",  "fc.raw.file")], "fc.raw.file", 
                      function(x) data.frame(metric = qualUniform(x$retention.time)))
      colnames(qcScore)[colnames(qcScore)=="metric"] = .self$qcName
      
      return(list(plots = lpl, qcScores = qcScore))
    }, 
    qcCat = "LC", 
    qcName = "EVD:~ID~rate~over~RT", 
    orderNr = 0150
  )
    return(.self)
  })
)  



#####################################################################

qcMetric_EVD_PreCal =  setRefClass(
  "qcMetric_EVD_PreCal",
  contains = "qcMetric",
  methods = list(initialize=function() {  callSuper(  
    helpTextTemplate = 
      "Mass accurary before calibration. Outliers are marked as such ('out-of-search-tol') using ID rate and standard deviation as additional information (if available).
If any Raw file is flagged 'failed', increasing MaxQuant's first-search tolerance (20ppm by default, here: %1.1f ppm) might help
to enable successful recalibration.
A bug in MaxQuant sometimes leads to excessively high ppm mass errors (>10<sup>4</sup>) reported in the output 
data. However, this can sometimes be corrected for by re-computing the delta mass error from other data. If this is 
the case, a warning ('bugfix applied') will be shown.

Heatmap score [EVD: MS Cal Pre (%1.1f)]: the centeredness (function CenteredRef) of uncalibrated masses in relation to the search window size.
",
    workerFcn = function(.self, df_evd, df_idrate, tolerance_pc_ppm, tolerance_sd_PCoutOfCal)
    {
      ## completeness check
      #stopifnot(c("...") %in% colnames(df_pg))
      
      .self$helpText = sprintf(.self$helpTextTemplate, tolerance_pc_ppm, tolerance_pc_ppm)
      
      fix_cal = fixCalibration(df_evd, df_idrate, tolerance_sd_PCoutOfCal)
      
      ## some outliers can have ~5000ppm, blowing up the plot margins
      ## --> remove outliers 
      ylim_g = range(boxplot.stats(fix_cal$df_evd$uncalibrated.mass.error..ppm.)$stats[c(1, 5)], c(-tolerance_pc_ppm, tolerance_pc_ppm) * 1.05)
      ## PLOT
      lpl =
        byXflex(fix_cal$df_evd, fix_cal$df_evd$fc.raw.file, 20, plot_UncalibratedMSErr, sort_indices = TRUE, 
                MQBug_raw_files = fix_cal$affected_raw_files, 
                y_lim = ylim_g,
                stats = fix_cal$stats,
                extra_limit = tolerance_pc_ppm,
                title_sub = fix_cal$recal_message)
      
      ## scores
      qc_MS1deCal = ddply(fix_cal$df_evd, "fc.raw.file", 
                          function(x) {
                            xd = na.omit(x$uncalibrated.mass.error..ppm.)
                            if (length(xd)==0) {
                              r = HEATMAP_NA_VALUE ## if empty, give the Raw file an 'NA' score
                            } else if (fix_cal$stats$outOfCal[fix_cal$stats$fc.raw.file == x$fc.raw.file[1]]) {
                              r = 0 ## if we suspect out-of-calibration, give lowest score
                            } else {
                              r = qualCenteredRef(xd, tolerance_pc_ppm)
                            } 
                            return (data.frame(med_rat = r))
                          })
      
      cname = sprintf(.self$qcName, tolerance_pc_ppm)
      colnames(qc_MS1deCal) = c("fc.raw.file", cname)
      
      
      return(list(plots = lpl, qcScores = qc_MS1deCal))
    }, 
    qcCat = "MS", 
    qcName = "EVD:~MS~Cal-Pre~(%1.1f)", 
    orderNr = 0260
  )
    return(.self)
  })
)


#####################################################################

qcMetric_EVD_PostCal =  setRefClass(
  "qcMetric_EVD_PostCal",
  contains = "qcMetric",
  methods = list(initialize=function() {  callSuper(  
    helpTextTemplate = 
      "Precursor mass accuracy after calibration. Failed samples from precalibration data are still marked here.
Ppm errors should be centered on zero and their spread is expected to be significantly smaller than before calibration.

Heatmap score [EVD: MS Cal-Post]: The variance and centeredness around zero of the calibrated distribution (function GaussDev).
",
    workerFcn = function(.self, df_evd, df_idrate, tolerance_pc_ppm, tolerance_sd_PCoutOfCal, tol_ppm_mainSearch)
    {
      ## completeness check
      #stopifnot(c("...") %in% colnames(df_pg))
      
      fix_cal = fixCalibration(df_evd, df_idrate, tolerance_sd_PCoutOfCal)
      
      ylim_g = range(na.rm = TRUE, boxplot.stats(fix_cal$df_evd$mass.error..ppm.)$stats[c(1, 5)], c(-tol_ppm_mainSearch, tol_ppm_mainSearch) * 1.05)
      ## PLOT
      lpl =
        byXflex(fix_cal$df_evd, fix_cal$df_evd$fc.raw.file, 20, plot_CalibratedMSErr, sort_indices = TRUE,
                MQBug_raw_files = fix_cal$affected_raw_files,
                y_lim = ylim_g,
                stats = fix_cal$stats,
                extra_limit = tol_ppm_mainSearch,
                title_sub = fix_cal$recal_message_post)
      
      ## QC measure for post-calibration ppm error
      ## .. assume 0 centered and StdDev of observed data
      obs_par = ddply(fix_cal$df_evd[, c("mass.error..ppm.", "fc.raw.file")], "fc.raw.file", 
                      function(x) data.frame(mu = mean(x$mass.error..ppm., na.rm = TRUE), 
                                             sd = sd(x$mass.error..ppm., na.rm = TRUE)))
      qc_MS1Cal = data.frame(fc.raw.file = obs_par$fc.raw.file, 
                             val = sapply(1:nrow(obs_par), function(x) qualGaussDev(obs_par$mu[x], obs_par$sd[x])))
      ## if we suspect out-of-calibration, give lowest score
      qc_MS1Cal$val[qc_MS1Cal$fc.raw.file %in% fix_cal$stats$fc.raw.file[ fix_cal$stats$outOfCal ]] = 0 
      ## MQ mass bugfix will not work for postCalibration, since values are always too low
      qc_MS1Cal$val[qc_MS1Cal$fc.raw.file %in% fix_cal$stats$fc.raw.file[ fix_cal$stats$hasMassErrorBug ]] = HEATMAP_NA_VALUE
      colnames(qc_MS1Cal)[colnames(qc_MS1Cal) == "val"] = .self$qcName
      
      return(list(plots = lpl, qcScores = qc_MS1Cal))
    }, 
    qcCat = "MS", 
    qcName = "EVD:~MS~Cal-Post", 
    orderNr = 0270
  )
    return(.self)
  })
)  

#####################################################################

qcMetric_EVD_Top5Cont =  setRefClass(
  "qcMetric_EVD_Top5Cont",
  contains = "qcMetric",
  methods = list(initialize=function() {  callSuper(  
    helpTextTemplate = 
      "PTXQC will explicitly show the five most abundant external protein contaminants
(as detected via MaxQuant's contaminants FASTA file) by Raw file, and summarize the 
remaining contaminants as 'other'. This allows to track down which proteins exactly contaminate your sample.
Low contamination is obviously better.
Transparency indicates the average peptide intensity in each Raw file. It is
not unusual to see samples with low sample content to have higher contamination.
    
Heatmap score [EVD: Contaminants]: as fraction of summed intensity
 - 0 = sample full of contaminants
 - 1 = no contaminants

",
    workerFcn = function(.self, df_evd)
    {
      ## completeness check
      stopifnot(c("intensity", "contaminant", "fc.raw.file") %in% colnames(df_evd))
      
      ##
      ## elaborate contaminant fraction per Raw.file (this is not possible from PG, since raw files could be merged)
      ## find top 5 contaminants (globally)
      ##
      ## if possible, work on protein names (since MQ1.4), else use proteinIDs
      if ("protein.names" %in% colnames(df_evd))
      {
        evd_pname = "protein.names"        
      } else if ("proteins" %in% colnames(df_evd)) {
        evd_pname = "proteins" 
      } else {
        stop("Top5-Contaminants: Neither 'protein.names' nor 'proteins' column was found in data but is required.")
      }
      
      ## protein.names are sometimes not unique, e.g. if a contaminant is involved:
      ## "P02768;CON__P02768-1" and "P02768" will both give the same name (since contaminant name is empty)
      ## Thus, the distribution of bars will look slightly different (but summed percentages are identical)
      
      ## some protein.names are empty (usually the CON__ ones) ... so we substitute with ID
      df_evd$pname = df_evd[, evd_pname];
      df_evd$pname[df_evd$pname==""] = df_evd$proteins[df_evd$pname==""] ## a NOP if it already is 'proteins', but ok
      
      df_evd.totalInt = sum(as.numeric(df_evd$intensity), na.rm = TRUE)
      df_evd.cont.only = df_evd[df_evd$contaminant,]
      cont.top = by(df_evd.cont.only, df_evd.cont.only$pname, function(x) sum(as.numeric(x$intensity), na.rm = TRUE) / df_evd.totalInt*100)
      cont.top.sort = sort(cont.top, decreasing = TRUE)
      #head(cont.top.sort)
      cont.top5.names = names(cont.top.sort)[1:5]
      
      lpl = list()
      if (is.null(cont.top5.names))
      {
        lpl[["noCont"]] = ggText("EVD: Top5 Contaminant per Raw file",
                                 "No contaminants found in any sample.\n\nIncorporating contaminants during search is highly recommended!",
                                 "red")
      } else {
        lpl =
          byXflex(df_evd[, c("intensity", "pname", "fc.raw.file", "contaminant")], df_evd$fc.raw.file, 40, sort_indices = TRUE, 
                  plot_ContEVD, top5=cont.top5.names)
      }
      
      ## QC measure for contamination
      qc_cont = ddply(df_evd[, c("intensity", "contaminant", "fc.raw.file")], "fc.raw.file", 
                      function(x) {
                        val = ifelse(is.null(cont.top5.names), 
                                     HEATMAP_NA_VALUE, ## use NA in heatmap if there are no contaminants
                                     1-qualLinThresh(sum(as.numeric(x$intensity[x$contaminant]), na.rm = TRUE)/
                                                       sum(as.numeric(x$intensity), na.rm = TRUE)))
                        return(data.frame(val = val, check.names = FALSE))
                      }
      )
      colnames(qc_cont)[colnames(qc_cont) == "val"] = .self$qcName
      
      return(list(plots = lpl, qcScores = qc_cont))
    }, 
    qcCat = "Prep", 
    qcName = "EVD:~Contaminants", 
    orderNr = 0010
  )
    return(.self)
  })
)  

#####################################################################

qcMetric_EVD_MS2OverSampling =  setRefClass(
  "qcMetric_EVD_MS2OverSampling",
  contains = "qcMetric",
  methods = list(initialize=function() {  callSuper(  
    helpTextTemplate = 
      "An oversampled 3D-peak is defined as a peak whose peptide ion (same sequence and same charge 
state) was identified by at least two distinct MS<sup>2</sup> spectra in the same Raw file. 
    For high complexity samples, oversampling of individual 3D-peaks automatically leads to undersampling 
    or even omission of other 3D-peaks, reducing the number of identified peptides. Oversampling occurs in 
    low-complexity samples or long LC gradients, as well as undersized dynamic exclusion windows for data 
    independent acquisitions. 

Heatmap score [EVD: MS<sup>2</sup> Oversampling]: The percentage of non-oversampled 3D-peaks. 
",
    workerFcn = function(.self, df_evd)
    {
      ## completeness check
      stopifnot(c("fc.raw.file", "ms.ms.count") %in% colnames(df_evd))
      
      d_dups = ddply(df_evd, "fc.raw.file", function(x) {
        tt = as.data.frame(table(x$ms.ms.count), stringsAsFactors = FALSE)
        tt$Count = as.numeric(tt$Var1)
        ## remove "0", since this would be MBR-features
        tt = tt[tt$Count!=0,]
        ## summarize everything above 3 counts
        if (any(tt$Count >= 3)) {
          tt$Count[tt$Count >= 3] = "3+"
          tt = ddply(tt, "Count", function(x) data.frame(Freq=sum(x$Freq)))
        }
        ## make counts relative
        fraction = tt$Freq / sum(tt$Freq) * 100
        return (data.frame(n=as.character(tt$Count), fraction = fraction))
      })
      
      lpl =
        byXflex(d_dups, d_dups$fc.raw.file, 30, plot_MS2Oversampling, sort_indices = TRUE)
      
      ## QC measure for how many peaks were fragmented only once
      qc_evd_twin = d_dups[d_dups$n==1,]
      cname = .self$qcName
      qc_evd_twin[, cname] = qualLinThresh(qc_evd_twin$fraction/100)
      
      return(list(plots = lpl, qcScores = qc_evd_twin[, c("fc.raw.file", cname)]))
    }, 
    qcCat = "MS", 
    qcName = "EVD:~MS^2~Oversampling", 
    orderNr = 0250
  )
    return(.self)
  })
)  

#####################################################################

qcMetric_EVD_MissingValues =  setRefClass(
  "qcMetric_EVD_MissingValues",
  contains = "qcMetric",
  methods = list(initialize=function() {  callSuper(    
    helpTextTemplate = 
      "Missing peptide intensities per Raw file from evidence.txt.
This metric shows the fraction of missing peptides compared to all peptides seen in the whole experiment.
The more Raw files you have, the higher this fraction is going to be (because there is always going
to be some exotic [low intensity?] peptide which gets [falsely] identified in only a single Raw file).
A second plot shows how many peptides (Y-axis) are covered by at least X Raw files.
A third plot shows the density of the observed (line) and the missing (filled area) data.
To reconstruct the distribution of missing values, an imputation strategy is required, so the argument is somewhat
circular here. If all Raw files are (technical) replicates, i.e. we can expect that missing peptides are indeed
present and have an intensity similar to the peptides we do see, then the median is a good estimator.
This method performs a global normalization across Raw files (so their observed intensitiy distributions have the same mean),
before computing the imputed values. Afterwards, the distributions are de-normalized again (shifting them back to their)
original locations -- but this time with imputed peptides.

Peptides obtained via Match-between-run (MBR) are accounted for (i.e. are considered as present = non-missing).
Thus, make sure that MBR is working as intended (see MBR metrics).

<b>Warning:</b> this metric is meaningless for fractionated data!
<b>TODO:</b> compensate for lower scores in large studies (with many Raw files), since peptide FDR is accumulating!?

Heatmap score [EVD: Pep Missing]: Linear scale of the fraction of missing peptides.
",
    workerFcn = function(.self, df_evd)
    {
      ## completeness check
      stopifnot(c("fc.raw.file", "modified.sequence", "intensity") %in% colnames(df_evd))
      
      if (('fraction' %in% colnames(df_evd)) && (length(unique(df_evd$fraction)) > 1)) {
        lpl = list(ggText("Skipped", "Missing values calculation skipped. Fractionated data detected!"))
        return(list(plots = lpl))
      }
      
      ## make peptides unique per Raw file
      df_u = ddply(df_evd[ , c("fc.raw.file", "modified.sequence")], "fc.raw.file",
                   function(x) {
                     return(x[!duplicated(x$modified.sequence),])
                  })
            
      global_peps = unique(df_u$modified.sequence)
      global_peps_count = length(global_peps)
      
      ## percent identified in each Raw file
      pep_set = ddply(df_u[ , c("fc.raw.file", "modified.sequence")], "fc.raw.file",
                      function(x) {
                        score = 100*length(intersect(global_peps, x$modified.sequence)) / global_peps_count
                        return(data.frame(idFraction = score))
                      })
      
      lpl = byXflex(pep_set, pep_set$fc.raw.file, subset_size = 50, FUN = function(dx) {
        p = ggplot(dx) + 
          geom_bar(aes_string(x = "fc.raw.file", y = "idFraction"), stat = "identity") +
          addGGtitle("[experimental] EVD: Non-Missing Peptides", "compared to all peptides seen in experiment") +
          xlab("") +
          ylab("Fraction of total peptides [%]") +
          ylim(0, 100) +
          scale_x_discrete_reverse(dx$fc.raw.file) +
          coord_flip()
        return(p)
      })                      
      
      #for (pl in lpl) print(pl)
      
      tbl = table(df_u$modified.sequence)
      head(tbl)
      tbl_smry = as.data.frame(table(tbl))
      tbl_smry$FreqRel = tbl_smry$Freq / global_peps_count
      tbl_smry = tbl_smry[nrow(tbl_smry):1,] ## invert
      tbl_smry$FreqCum = cumsum(tbl_smry$FreqRel) * 100
      tbl_smry$x = as.numeric(tbl_smry$tbl)
      
      p = ggplot(tbl_smry, aes_string(x = "x", y = "FreqCum")) + 
        geom_line() +
        geom_point() +
        addGGtitle("[experimental] EVD: Non-missing by set", "") +
        xlab("Minimum # Raw files") +
        ylab("Fraction of total peptides [%]") +
        ylim(0, 100)
      lpl[["missingCul"]] = p
      
      ## intensity distribution of missing values
      df_evd$logInt = log2(df_evd$intensity)
      
      lpl_dens = byXflex(df_evd[, c("modified.sequence", "fc.raw.file", "logInt")], df_evd$fc.raw.file,
                         subset_size = 5, FUN = function(dx) {
             
        d_mat = dcast(dx, modified.sequence ~ fc.raw.file, fun.aggregate = mean, value.var = "logInt")
        
        ## ... normalization factors
        d_mat_mult = sapply(2:ncol(d_mat), function(x) {
          mult = mean(d_mat[, x] / d_mat[, 2], na.rm = TRUE)
          return(mult)
        })
        df_mult = data.frame(fc.raw.file = colnames(d_mat)[-1], mult = d_mat_mult)
        ## .. normalize data
        d_mat_n = d_mat
        d_mat_n[, -1] = sweep( d_mat_n[, -1], 2, d_mat_mult, '/')
        ## 
        head(d_mat_n)
        ## find impute value
        pep_mean = rowMeans(d_mat_n[, -1], na.rm = TRUE)
        df_missing = ddply(df_mult, "fc.raw.file", function(x) {
          ## get set of missing values
          values = pep_mean[is.na(d_mat_n[, as.character(x$fc.raw.file)])]
          ## de-normalize (back to old intensity range)
          values = values * x$mult
          return(data.frame(missingVals = values))
        })
        head(df_missing)
        
        pl = ggplot(df_missing) + 
          geom_freqpoly(data = df_evd, aes_string(x = "logInt", col="fc.raw.file"), binwidth = 0.5, size = 1.2) +
          geom_area(aes_string(x = "missingVals", col="fc.raw.file", fill = "fc.raw.file"), position = position_dodge(width=0), binwidth = 0.5, stat="bin", alpha=0.5) +
          xlab("Intensity [log2]") +
          ggtitle(" [experimental] EVD: Imputed Peptide Intensity Distribution of Missing Values") +
          scale_fill_manual(values = rep(brewer.pal(6,"Accent"), times=40), guide = guide_legend("")) +
          scale_colour_manual(values = rep(brewer.pal(6,"Accent"), times=40), guide = "none")
        
        return(pl)
      })
      
      lpl = append(lpl, lpl_dens)
      
      
      ## QC measure for fraction of missing values
      cname = .self$qcName
      pep_set[, cname] = qualLinThresh(pep_set$idFraction, 100) ## a no-op, just for clarity
      qcScore = pep_set[, c("fc.raw.file", cname)]
      
      return(list(plots = lpl, qcScores = qcScore))
    }, 
    qcCat = "prep", 
    qcName = "EVD:~Pep~Missing~Values", 
    orderNr = 0390  # just before peptide count
  )
    return(.self)
  })
)

#####################################################################


