
qcMetric_param = qcMetric$new(helpText="MaxQuant parameters, extracted from parameters.txt.", 
   workerFcn=function(.self, df_mqpar)
   {
     ##todo: read in mqpar.xml to get group information and ppm tolerances for all groups (parameters.txt just gives Group1)
     
     line_break = "\n"; ## use space to make it work with table
     ## remove AIF stuff
     df_mqpar = df_mqpar[!grepl("^AIF ", df_mqpar$parameter),]
     df_mqpar$value = gsub(";", line_break, df_mqpar$value)
     ## seperate FASTA files (usually they destroy the layout)
     idx_fastafile = grepl("fasta file", df_mqpar$parameter, ignore.case = TRUE)
     d_par_file = df_mqpar[idx_fastafile, ]
     fasta_files = sapply(unlist(strsplit(d_par_file$value, "\n")), function(x) rev(strsplit(x,"\\", fixed = TRUE)[[1]])[1])
     d_par = df_mqpar[!idx_fastafile, ]
     ## remove duplicates
     d_par = d_par[!duplicated(d_par$parameter),]
     rownames(d_par) = d_par$parameter
     
     ## trim long param names (the user should know what they mean)
     d_par$parameter = sapply(d_par$parameter, function (s) {
       allowed_len = nchar("Min. score for unmodified .."); 
       if (nchar(s) > allowed_len) {
         s = paste(substring(s, 1, allowed_len), "..", collapse = "", sep="")
       }
       return (s)
     })
     ## break long values into multiple lines (to preserve table width)
     d_par$value = sapply(d_par$value, function (s) 
     {
       allowed_len = nchar("Use least modified peptide"); ## this is a typical entry -- everything which is longer gets split
       r = paste(sapply(unlist(strsplit(s, line_break, fixed = TRUE)), function(s1) {
         if (nchar(s1) > allowed_len) {
           s_beg = seq(1, nchar(s1) - 1, allowed_len)
           s1 = paste(unlist(substring(s1, s_beg, s_beg + allowed_len)), collapse = line_break)
         }
         return(s1)
       }), collapse = line_break)
       return (r)
     })
     
     ## sort by name
     d_par = d_par[order(d_par$parameter), ]
     
     ## two column layout
     if (nrow(d_par) %% 2 != 0) d_par = rbind(d_par, "") ## make even number of rows
     mid = nrow(d_par) / 2
     d_par$page = 1
     d_par$page[1:mid] = 0
     
     parC = c("parameter", "value")
     d_par2 = cbind(d_par[d_par$page==0, parC], d_par[d_par$page==1, parC])
     
     par_pl = plotTable(d_par2, title = "PAR: parameters", footer = fasta_files)
     
     return(list(plots = list(par_pl)))
   }, 
   qcCat = NA_character_, 
   qcName = NA_character_, 
   heatmapOrder = NaN)  ## should not show up in heatmap

#####
qcMetric_MSMSIdRate = qcMetric$new(helpText="MaxQuant parameters, extracted from parameters.txt.", 
  workerFcn=function(.self, df_summary, id_rate_bad, id_rate_great)
  {
    dms = df_summary$"ms.ms.identified...."
    dms[is.na(dms)] = 0  ## ID rate can be NaN for some raw files if NOTHING was acquired
    lab_IDd = c("red", 
                "blue", 
                "green")
    names(lab_IDd) = c("bad (<" %+% id_rate_bad %+% "%)", 
                       "ok (" %+% id_rate_bad %+% "-" %+% id_rate_great %+% "%)",
                       "great (>" %+% id_rate_great %+% "%)")
    df_summary$cat = factor(cut(dms, breaks=c(-1, id_rate_bad, id_rate_great, 100), labels=names(lab_IDd)))

    pl = plot_IDRate(df_summary, id_rate_bad, id_rate_great, lab_IDd)
    plots = list(pl)
    
    ## table of files with 'bad' MS/MS id rate
    bad_id_count = sum(df_summary$cat==names(lab_IDd[1]))
    if (bad_id_count>0)
    {
      sm_badID = df_summary[df_summary$cat==names(lab_IDd[1]), c("raw.file","ms.ms.identified....")]
      if (nrow(sm_badID) > 40)
      {
        sm_badID[40, "raw.file"] = paste(nrow(sm_badID) - 39, "more ...");
        sm_badID[40, "ms.ms.identified...."] = ""
        sm_badID = sm_badID[1:40, ]
      }
      p_tbl = plotTable(sm_badID, 
                      title = paste0("SM: Files with '", lab_IDd[1],
                                   "' ID rate (", round(bad_id_count*100/nrow(df_summary)), "% of samples)"),
                      col_names = c("Raw file", "% identified"))
      plots[["bad_id_table"]] =  p_tbl
    }
    
    ## QC measure for ID rate (threshold reached?)
    ## update QC name based on parameter values
    .self$qcName = sprintf(.self$qcName, id_rate_great)
    
    return(list(plots = plots))
  }, 
  qcCat = "MS", 
  qcName = "SM:~MS^2~ID~rate (\">%1.0f\")", 
  heatmapOrder = 300)  ## should not show up in heatmap


## list of qcMetric objects
lst_qcMetrics = ls(pattern="qcMetric_.*")
