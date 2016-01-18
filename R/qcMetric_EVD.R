
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
        evd_uniqueGroup = !grepl(";", d_evd$protein.group.ids)
        ## do not trust MBR here. We want real evidence!
        evd_realMS = !grepl("MATCH", d_evd$type)
        ## for each Raw file: find unique peptides of our contaminant
        cont_data.l = dlply(d_evd[evd_uniqueGroup & evd_realMS, ], "fc.raw.file",
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

qcMetric_EVD_... = qcMetric$new(
  helpText = 
    "...",
  workerFcn=function(.self, df_pg, int_cols, MAP_pg_groups)
  {
    ## completeness check
    stopifnot(c(int_cols, "contaminant") %in% colnames(df_pg))
    
    
    return(list(plots = pg_plots_cont))
  }, 
  qcCat = NA_character_, 
  qcName = NA_character_, 
  heatmapOrder = NaN)




#####################################################################
