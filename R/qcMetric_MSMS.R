
#####################################################################

qcMetric_MSMS_MSMSDecal =  setRefClass(
  "qcMetric_MSMS_MSMSDecal",
  contains = "qcMetric",
  methods = list(initialize=function() {  callSuper(
    helpTextTemplate = 
      "MS/MS decalibration metric. If most of the fragments are within tighter bounds, 
you can reduce the fragment mass tolerance to obtain more 
identifications under the same FDR. On the other hand, if the fragment mass errors are not centered on 
zero, a recalibration of the instrument should be performed.
If the (Gaussian-like) distribution is cut too severely on either side by the search tolerance window in MaxQuant,
you might be able to increase the number of identifications by allowing for a wider MS/MS search window when re-running MaxQuant.
However, the number of decoy identifications will increase as well, potentially offsetting any gain when FDR is applied.

Heatmap score [MSMS: MS<sup>2</sup> Cal (Analyzer)]: rewards centeredness around 0 ppm/Da (function Centered).
",
    workerFcn = function(.self, df_msms, fc_raw_files)
    {
      ## completeness check
      if(!(.self$checkInput(c("fc.raw.file", "fragmentation", "reverse", "mass.deviations..da."), colnames(df_msms)))){return(NULL)}
      ## older MQ versions do not have 'mass.analyzer' or 'mass.deviations..ppm.'
      ## , so we use fragmentation instead (this is a little risky, since you could do CID fragmentation and forward to Orbi, but hey...)
      if (!("mass.analyzer" %in% colnames(df_msms))) df_msms$mass.analyzer = df_msms$fragmentation
      
      sampleMax = function(x, max = 300) {
        ## sample at most 'max' items from 1:x
        ## if 'x' < max, return 1:x
        if (x < max) return(1:x)
        return (sort(sample.int(x, size = max)))
      }
      
      ms2_decal = ddply(df_msms, c("fc.raw.file", "mass.analyzer"), .fun = function(x) {
        df.ms = NULL
        ##
        ##  Forwards
        ##
        if (any(!x$reverse))
        {
          idx_nr = which(!x$reverse)
          ## select a representative subset, otherwise the number of datapoints is just too large
          idx_nr_subset = idx_nr[sampleMax(length(idx_nr))]
          df.ms = getFragmentErrors(x[idx_nr_subset, , drop=FALSE])
          if (!is.null(df.ms)) df.ms$type="forward"
        } 
        
        if (any(x$reverse))
        {
          idx_nr = which(x$reverse)
          ## select a representative subset, otherwise the number of datapoints is just too large
          idx_nr_subset = idx_nr[sampleMax(length(idx_nr))]
          df.ms_r = getFragmentErrors(x[idx_nr_subset, , drop=FALSE])
          if (!is.null(df.ms) & !is.null(df.ms_r)) {
            df.ms_r$type="decoy"
            ## only merge if we have hits
            df.ms = rbind(df.ms, df.ms_r)
          }
        }
        
        return (df.ms)
      })

      #ms2_range = diff(range(ms2_decal$msErr, na.rm = TRUE))
      #ms2_binwidth = ms2_range/20
      ## precision (plotting is just so much quicker, despite using a fixed binwidth)
      #ms2_decal$msErr = round(ms2_decal$msErr, digits=ceiling(-log10(ms2_binwidth)+1))
      
      ## separate plots for each mass analyzer, since we want to keep 'fixed' scales for all raw.files (comparability)
      lpl = dlply(ms2_decal, "mass.analyzer", function(ms2_decal) {
        ## create filename inside, since we need to retain the factor levels (i.e. ordering)
        ## and this only works if raw file + massanalyzer is unique
        ms2_decal$new_filename = paste(ms2_decal$fc.raw.file, paste(ms2_decal$mass.analyzer, ms2_decal$unit), sep="\n")
        ## i.e. change name without underlying value
        ms2_decal$file = ms2_decal$fc.raw.file
        levels(ms2_decal$file) = ms2_decal$new_filename[ match(levels(ms2_decal$file), ms2_decal$file) ]
        byXflex(ms2_decal, ms2_decal$fc.raw.file, 9, plot_MS2Decal, sort_indices = TRUE)
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
      "Under optimal digestion conditions (high enzyme grade etc.), only few missed cleavages (MC) are expected. In 
general, increased MC counts also increase the number of peptide signals, thus cluttering the available 
space and potentially provoking overlapping peptide signals, biasing peptide quantification.
Thus, low MC counts should be favored. Interestingly, it has been shown recently that 
incorporation of peptides with missed cleavages does not negatively influence protein quantification (see 
[http://pubs.acs.org/doi/abs/10.1021/pr500294d](Chiva, C., Ortega, M., and Sabido, E. Influence of the Digestion Technique, Protease, and Missed 
Cleavage Peptides in Protein Quantitation. J. Proteome Res. 2014, 13, 3979-86) ). 
However this is true only if all samples show the same degree of digestion. High missed cleavage values 
can indicate for example, either a) failed digestion, b) a high (post-digestion) protein contamination, or 
c) a sample with high amounts of unspecifically degraded peptides which are not digested by trypsin. 

If MC>=1 is high (>20%) you should increase the missed cleavages settings in MaxQuant and compare the number of peptides.
Usually high MC correlates with bad identification rates, since many spectra cannot be matched to the forward database.

In the rare case that 'no enzyme' was specified in MaxQuant, neither scores nor plots are shown.

Heatmap score [MSMS: MC]: the fraction (0% - 100%) of fully cleaved peptides per Raw file

Heatmap score [MSMS: MC Var]: each Raw file is scored for its deviation (score: MedianDist) from the 'average' digestion state of the 
current study. ",
    workerFcn = function(.self, df_msms, df_evd = NULL)
    {
      ## completeness check
      if(!(.self$checkInput(c("fc.raw.file"), colnames(df_msms)))){return(NULL)}
      if (!is.null(df_evd)) if(!.self$checkInput(c("contaminant", "id"), colnames(df_evd))){return(NULL)}
      
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
          byXflex(st_bin, st_bin$fc.raw.file, 25, plot_MissedCleavages, title_sub = msg_cont_removed, sort_indices = TRUE)
        
        ## QC measure for missed-cleavages variation
        qc_score = data.frame(fc.raw.file = st_bin$fc.raw.file, valMC = st_bin[, "0"])
        qc_score$valMCVar = qualMedianDist(qc_score$valMC)
        
      } else {
        lpl = list(ggText("MSMS: Missed cleavages per Raw file",
                          "No enzyme was specified.\nDigestion efficiency cannot be scored."))
        qc_score = data.frame(fc.raw.file = unique(df_msms$fc.raw.file),
                              valMC = HEATMAP_NA_VALUE,
                              valMCVar = HEATMAP_NA_VALUE)
      }## end enyzme check
      
      colnames(qc_score)[colnames(qc_score) == "valMC"] = sprintf(.self$qcName, "MC")
      colnames(qc_score)[colnames(qc_score) == "valMCVar"] = sprintf(.self$qcName, "MC~Var")

      return(list(plots = lpl, qcScores = qc_score))
    }, 
    qcCat = "Prep", 
    qcName = "MSMS:~%s", 
    orderNr = 0040
  )
    return(.self)
  })
)  



