
qcMetric_SM_MSMSIdRate =  setRefClass(
  "qcMetric_SM_MSMSIdRate",
  contains = "qcMetric",
  methods = list(initialize=function() {  callSuper(  
    helpTextTemplate = 
      "MS/MS identification rate per Raw file from summary.txt (SM). Each Raw file is colored
according to its ID rate and categorized into performance bins as 'bad', 'ok' and 'great'. 
Raw files below 'ok', are listed separately on the next page of the report for convenient follow-up.

The thresholds for the bins are

%s


Heatmap score [SM: MS<sup>2</sup> IDrate (>%1.0f)]: reaches 1 (=100%%) if the threshold for 'great' is reached or exceeded. ",
    workerFcn = function(.self, df_summary, id_rate_bad, id_rate_great)
    {
      stopifnot(.self$checkInput(c("fc.raw.file", "ms.ms.identified...."), colnames(df_summary)))
      
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
      lpl = list(pl)
      
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
        plot_title = paste0("SM: Files with '", lab_IDd[1], "' ID rate ")
        p_tbl = plotTable(sm_badID, 
                          title = plot_title,
                          footer = paste0(round(bad_id_count*100/nrow(df_summary)), "% of samples)"),
                          col_names = c("Raw file", "% identified"))
        lpl[["bad_id_table"]] =  p_tbl
      }
      
      ## update help text according to actual limits
      inText = paste(paste0("  - ", names(lab_IDd)), collapse="\n", sep="")
      
      .self$helpText = sprintf(.self$helpTextTemplate, inText, id_rate_great)
      
      ## QC measure for ID rate (threshold reached?)
      ## update QC name based on parameter values
      ## QC measure for ID rate (threshold reached?)
      qc_sm_id = df_summary[, c("fc.raw.file", "ms.ms.identified....")]
      cname = sprintf(.self$qcName, id_rate_great)
      qc_sm_id[, cname] = qualLinThresh(qc_sm_id$ms.ms.identified.... , id_rate_great)
      qcScore = qc_sm_id[,c("fc.raw.file", cname)]
      
      return(list(plots = lpl, qcScores = qcScore))
    }, 
    qcCat = "MS", 
    qcName = "SM:~MS^2~ID~rate (\">%1.0f\")", 
    orderNr = 300
  )
    return(.self)
  })
  
  
  qcMetric_SM_TIC =  setRefClass(
    "qcMetric_SM_TIC",
    contains = "qcMetric",
    methods = list(initialize=function() {  callSuper(
      helpTextTemplate = 
        "Total Ion Count. Reflects the summed intensity of all masses identified at any point of time.
",
      workerFcn = function(.self, d_smy)
      {
        ## completeness check
        stopifnot(.self$checkInput(c("fc.raw.file", "TIC"), colnames(d_smy)))
        
        df_long = ddply(d_smy, "fc.raw.file", function(x) {
        
          n = lenght(x$TIC)
          return(data.frame(intensity=x$TIC[seq(1,n,2)], RT=x$TIC[seq(2,n,2)]))
        })
        
        head(df_long)
        
        lpl =
          byXflex(df_long, df_long$fc.raw.file, 6, plot_TIC, sort_indices = FALSE)
        
        ## QC measure for smoothness of TopN over RT
        qc_TIC = ddply(df_long, "fc.raw.file", function(x) data.frame(val = qualUniform(x$intensity)))
        colnames(qc_TIC)[colnames(qc_TIC) == "val"] = .self$qcName
        
        return(list(plots = lpl, qcScores = qc_TIC))
      }, 
      qcCat = "LC", 
      qcName = "TIC", 
      orderNr = 0025
    )
      return(.self)
    })
  )
  
)

