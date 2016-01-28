
qcMetric_EVD_UserContaminant = qcMetric$new(
  helpText = 
"User defined contaminant search (usually used for Mycoplasma detection, but can be used for an arbitrary (set of) proteins.

Two abundance measures are computed per Raw file:
  - fraction of intensity
  - fraction of spectral counts
",
  workerFcn=function(.self, df_evd, df_pg, lst_contaminants)
  {
    ## completeness check
    stopifnot(c("id", "fasta.headers") %in% colnames(df_pg))
    stopifnot(c("protein.group.ids", "type", "score", "intensity", "fc.raw.file") %in% colnames(df_evd))
    
    
    local_qcScores = data.frame()
    
    idx_pl = 1
    plots = list()
    
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
                              cont_scoreECDF = ddply(x, "idx_cont", function(xx) {
                                if (length(unique(xx$score)) < 2) return(NULL) ## not enough data for ECDF
                                r = getECDF(xx$score)
                                r$condition = c("sample", "contaminant")[xx$idx_cont[1]+1]
                                return(r)
                              })
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
        plots[[idx_pl]] = pl_cont
        idx_pl = idx_pl + 1
      } else {
        ## plot User-Contaminants
        plots[[idx_pl]] = byXflex(data = cont_data.long, indices = cont_data.long$fc.raw.file, subset_size = 120, 
                                  FUN = plot_ContUser, sort_indices = FALSE, name_contaminant = ca, extra_limit = ca_thresh)
        idx_pl = idx_pl + 1
          
        ## plot Andromeda score distribution of contaminant vs. sample
        llply(cont_data.l, function(l)
        {
          if (l$cont_data$above.thresh == FALSE) return(NULL)
          p = plot_ContUserScore(l$cont_scoreECDF, l$cont_data$fc.raw.file, l$cont_data$score_KS)
          plots[[idx_pl]] = p
          idx_pl = idx_pl + 1
          #print(p)
          return(NULL)
        })
        
        ## add heatmap column
        cname = sprintf(qcName, ca)
        cont_data[,cname] = as.numeric(!cont_data$above.thresh) ## inverse (0 is 'bad')
        
        qcScore = cont_data[, c("fc.raw.file", cname)]
        if (ncol(local_qcScores) == 0){
          local_qcScores = qcScore
        } else {
          local_qcScores = merge(local_qcScores, qcScore)
        }

      }
    } ## contaminant loop
    
    return(list(plots = plots, qcScores = local_qcScores))
  }, 
  qcCat = "Prep", 
  qcName = "X002X_catPrep_EVD:Contaminant~(%s)", 
  heatmapOrder = 0020)


#####################################################################

qcMetric_EVD_PeptideInt = qcMetric$new(
  helpText = 
"Peptide intensity ...",
  workerFcn=function(.self, df_evd, thresh_intensity)
  {
    ## completeness check
    stopifnot(c("fc.raw.file", "intensity") %in% colnames(df_evd))
    
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
  qcName = "X003X_catPrep_EVD:~Pep~Intensity~(\">%1.1f\")", 
  heatmapOrder = 0030)



#####################################################################

qcMetric_EVD_ProteinCount = qcMetric$new(
  helpText = 
"Number of Protein groups (after FDR) per Raw file. If MBR was enabled, three categories ('genuine (exclusive)', 'genuine + transferred', 'transferred (exclusive)'
 are shown, so the user can judge the gain that MBR provides. If the gain is low and the MBR scores are bad,
MBR should be switched off for the Raw files which are affected (could be a few or all).",
  workerFcn=function(.self, df_evd, thresh_intensity)
  {
    ## completeness check
    stopifnot(c("fc.raw.file", "protein.group.ids", "match.time.difference") %in% colnames(df_evd))

    protC = getProteinCounts(df_evd[, c("fc.raw.file", "protein.group.ids", "match.time.difference")])
    protC$block = factor(assignBlocks(protC$fc.raw.file, 30))
    
    max_prot = max(unlist(dlply(protC, "fc.raw.file", function(x) sum(x$counts))))
    ## average gain in percent
    gain_text = ifelse(reportMTD, sprintf("MBR gain: +%.0f%%", mean(protC$MBRgain, na.rm = TRUE)), "")
    
    lpl = dlply(protC, "block", .fun = function(x)
    {
      p = plot_CountData(data = x, 
                         y_max = max(thresh_intensity, max_prot)*1.1,
                         thresh_line = thresh_intensity,
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
    cname = sprintf(.self$qcName, thresh_intensity)
    qc_protc[,cname] = qualLinThresh(qc_protc$genuineAll, thresh_intensity)
    qcScore = qc_protc[, c("fc.raw.file", cname)]

    return(list(plots = lpl, qcScores = qcScore))
  }, 
  qcCat = 'general', 
  qcName = "X045X_catGen_EVD:~Prot~Count~(\">%1.0f\")", 
  heatmapOrder = 0450)


#####################################################################

qcMetric_EVD_PeptideCount = qcMetric$new(
  helpText = 
    "Number of peptides (after FDR) per Raw file. If MBR was enabled, three categories ('genuine (exclusive)', 'genuine + transferred', 'transferred (exclusive)'
  are shown, so the user can judge the gain that MBR provides. If the gain is low and the MBR scores are bad,
  MBR should be switched off for the Raw files which are affected (could be a few or all).",
  workerFcn=function(.self, df_evd, thresh_intensity)
  {
    ## completeness check
    stopifnot(c("fc.raw.file", "modified.sequence", "match.time.difference") %in% colnames(df_evd))
    
    pepC = getPeptideCounts(df_evd[, c("fc.raw.file", "modified.sequence", "match.time.difference")])
    pepC$block = factor(assignBlocks(pepC$fc.raw.file, 30))
    
    max_pep = max(unlist(dlply(pepC, "fc.raw.file", function(x) sum(x$counts))))
    ## average gain in percent
    gain_text = ifelse(reportMTD, sprintf("MBR gain: +%.0f%%", mean(pepC$MBRgain, na.rm = TRUE)), "")
    
    lpl = dlply(pepC, "block", .fun = function(x)
    {
      p = plot_CountData(data = x, 
                         y_max = max(param_EV_pepThresh, max_pep)*1.1,
                         thresh_line = param_EV_pepThresh,
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
    cname = sprintf(.self$qcName, thresh_intensity)
    qc_pepc[,cname] = qualLinThresh(qc_pepc$genuineAll, thresh_intensity)
    qcScore = qc_pepc[, c("fc.raw.file", cname)]
    
    return(list(plots = lpl, qcScores = qcScore))
  }, 
  qcCat = 'general', 
  qcName = "X040X_catGen_EVD:~Pep~Count~(\">%1.0f\")", 
  heatmapOrder = 0400)


#####################################################################

qcMetric_EVD_RTPeakWidth = qcMetric$new(
  helpText = 
    "RT peak width distribution ...",
  workerFcn=function(.self, df_evd)
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
  qcName = "X017X_catLC_EVD:~RT~Peak~Width", 
  heatmapOrder = 0170)


#####################################################################

qcMetric_EVD_MBRAlign = qcMetric$new(
  helpText = 
    "Match-between-runs Alignment (step 1/2, 1=align, 2=transfer) ...",
  workerFcn=function(.self, df_evd, tolerance_matching, raw_file_mapping)
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

        ## save some memory
        if (!exists("DEBUG_PTXQC")) {
          rm("evd_RT_t")
          rm("proj_align_h")
        }
      } ## no data
    } ## ambigous reference file
    
    
    return(list(plots = lpl, qcScores = qcScore))
  }, 
  qcCat = "LC", 
  qcName = "X021X_catLC_EVD:~MBR~Align", 
  heatmapOrder = 0210)


#####################################################################

qcMetric_EVD_... = qcMetric$new(
  helpText = 
    "...",
  workerFcn=function(.self, df_pg, int_cols, MAP_pg_groups)
  {
    ## completeness check
    stopifnot(c(int_cols, "contaminant") %in% colnames(df_pg))
    
    
    return(list(plots = pg_plots_cont, qcScores = qcScore))
  }, 
  qcCat = NA_character_, 
  qcName = NA_character_, 
  heatmapOrder = NaN)




#####################################################################
