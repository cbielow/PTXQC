
#####################################################################

qcMetric_MSMS_MSMSDecal =  setRefClass(
  "qcMetric_MSMS_MSMSDecal",
  contains = "qcMetric",
  methods = list(initialize=function() {  callSuper(
    helpTextTemplate = 
      "MS/MS decalibration metric. Use it to judge if your MS/MS tolerance has the correct
interval, i.e. can you narrow it more while still capturing most of the true (green) hits or does it need
widening because you are truncating the distribution (Gaussian) too much, thus loosing fragments?
The scoring function rewards centeredness around 0 ppm/Da.",
    workerFcn = function(.self, df_msms, fc_raw_files)
    {
      ## completeness check
      stopifnot(.self$checkInput(c("fc.raw.file", "fragmentation", "reverse", "mass.deviations..da."), colnames(df_msms)))
      ## older MQ versions do not have 'mass.analyzer' or 'mass.deviations..ppm.'
      ## , so we use fragmentation instead (this is a little risky, since you could to CID fragmentation and forward to Orbi, but hey...)
      if (!("mass.analyzer" %in% colnames(df_msms))) df_msms$mass.analyzer = df_msms$fragmentation
      
      
      ms2_decal = ddply(df_msms, c("fc.raw.file", "mass.analyzer"), .fun = function(x) {
        idx_nr = which(!x$reverse)
        ## select a representative subset, otherwise the number of datapoints is just too large
        idx_nr_subset = idx_nr[seq(1,length(idx_nr), by=ceiling(length(idx_nr)/1000))]
        df.ms = getFragmentErrors(x[idx_nr_subset,])
        df.ms$type="forward"
        
        if (any(x$reverse))
        {
          idx_nr = which(x$reverse)
          ## select a representative subset, otherwise the number of datapoints is just too large
          idx_nr_subset = idx_nr[seq(1,length(idx_nr), by=ceiling(length(idx_nr)/1000))]
          df.ms_r = getFragmentErrors(x[idx_nr_subset,])
          df.ms_r$type="decoy"
          
          ## only merge if we have hits (reverse hits might be few and $mass.deviations..da. might be empty)
          if (nrow(df.ms_r)) df.ms = rbind(df.ms, df.ms_r)
        }
        
        return (df.ms)
      })
      
      ms2_decal$msErr = as.numeric(as.character(ms2_decal$msErr))
      #ms2_range = diff(range(ms2_decal$msErr, na.rm = TRUE))
      #ms2_binwidth = ms2_range/20
      ## precision (plotting is just so much quicker, despite using a fixed binwidth)
      #ms2_decal$msErr = round(ms2_decal$msErr, digits=ceiling(-log10(ms2_binwidth)+1))
      ms2_decal$file = paste(ms2_decal$fc.raw.file, paste(ms2_decal$mass.analyzer, ms2_decal$unit), sep="\n")
      
      ## separate plots for each mass analyzer, since we want to keep 'fixed' scales for all raw.files (comparability)
      lpl = dlply(ms2_decal, "mass.analyzer", function(ms2_decal) {
        byXflex(ms2_decal, ms2_decal$fc.raw.file, 9, plot_MS2Decal, sort_indices = FALSE)
      })
      ## currently lpl is a list of lists() -- flatten
      lpl = Reduce(append, lpl)

      ##
      ## QC measure for centered-ness of MS2-calibration
      ##
      head(ms2_decal)
      qcScore = list()
      for (analyzer in unique(ms2_decal$mass.analyzer)) {
        qc_name = sprintf(.self$qcName, analyzer)
        qc_MS2_decal = ddply(ms2_decal[ms2_decal$mass.analyzer==analyzer, ], "fc.raw.file", 
                             function(x)
                             {
                               xx = na.omit(x$msErr);
                               return(data.frame(X1 = qualCentered(xx)))
                             })
        ## augment fragmentation methods with -Inf for missing raw files (otherwise they would become 'red'=fail)
        if (length( setdiff(fc_raw_files, qc_MS2_decal$fc.raw.file) )) {
          frag_missing = data.frame(fc.raw.file = setdiff(fc_raw_files, qc_MS2_decal$fc.raw.file), X1=-Inf)
          qc_MS2_decal = rbind(qc_MS2_decal, frag_missing)
        }
        colnames(qc_MS2_decal)[colnames(qc_MS2_decal)=="X1"] = qc_name
        qcScore[[analyzer]] = qc_MS2_decal[, c("fc.raw.file", qc_name)]
      }
      qcScore_all = Reduce(function(x,y) merge(x,y,all=TRUE), qcScore)
      
      return(list(plots = lpl, qcScores = qcScore_all))
    }, 
    qcCat = "MS", 
    qcName = "MSMS:~MS^2~Cal~(%s)", 
    orderNr = 0280
  )
    return(.self)
  })
)

#####################################################################

qcMetric_MSMS_MissedCleavages =  setRefClass(
  "qcMetric_MSMS_MissedCleavages",
  contains = "qcMetric",
  methods = list(initialize=function() {  callSuper(  
    helpTextTemplate = 
      "Metric for digestion efficiency. Fewer missed cleavages are better.
Also, the fraction of MC's across Raw files should be comparable (to ensure that we see
the same peptide species, enabling a more accurate comparison.",
    workerFcn = function(.self, df_msms, df_evd = NULL)
    {
      ## completeness check
      stopifnot(.self$checkInput(c("fc.raw.file"), colnames(df_msms)))
      if (!is.null(df_evd)) stopifnot(.self$checkInput(c("contaminant", "id"), colnames(df_evd)))
      
      max_mc = max(-Inf, df_msms$missed.cleavages, na.rm = TRUE) ## will be -Inf iff enzyme was not specified and columns is 100% NA
      if (!is.infinite(max_mc))
      { ## MC's require an enzyme to be set
        ## remove contaminants
        msg_cont_removed = "(includes contaminants -- no evidence.txt read)"
        if (!is.null(df_evd)) {
          msg_cont_removed = "(excludes contaminants)"
          df_msms$contaminant = df_evd$contaminant[match(df_msms$evidence.id, df_evd$id)]
          summary(df_msms$contaminant)
        }
        
        st_bin = ddply(df_msms[!df_msms$contaminant, c("missed.cleavages", "fc.raw.file")], "fc.raw.file", .fun = function(x) {
          t = table(x$missed.cleavages)/nrow(x)
          r = rep(0, max_mc + 1)
          names(r) = as.character(0:max_mc)
          r[names(t)] = t
          return (r)
        })
        lpl =
          byXflex(st_bin, st_bin$fc.raw.file, 25, plot_MissedCleavages, title_sub = msg_cont_removed, sort_indices = FALSE)
        
        ## QC measure for missed-cleavages variation
        qc_mc = data.frame(fc.raw.file = st_bin$fc.raw.file, valMC = st_bin[, "0"])
        qc_mc$valMCVar = qualMedianDist(qc_mc$valMC)
        colnames(qc_mc)[colnames(qc_mc) == "valMC"] = sprintf(.self$qcName, "MC")
        colnames(qc_mc)[colnames(qc_mc) == "valMCVar"] = sprintf(.self$qcName, "MC~Var")
        qc_score = qc_mc[, grep("valMC", colnames(qc_mc), invert=TRUE)]
      } ## end MC check
      
      return(list(plots = lpl, qcScores = qc_score))
    }, 
    qcCat = "Prep", 
    qcName = "MSMS:~%s", 
    orderNr = 0040
  )
    return(.self)
  })
)  



