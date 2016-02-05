
#####################################################################

qcMetric_MSMSScans_TopNoverRT =  setRefClass(
  "qcMetric_MSMSScans_TopNoverRT",
  contains = "qcMetric",
  methods = list(initialize=function() {  callSuper(
    helpTextTemplate = 
      "TopN over retention time. Similar to ID over RT, this metric reflects the complexity of the sample
at any point in time. Ideally complexity should be made roughly equal (constant) by choosing a proper (non-linear) LC gradient.
See [http://www.ncbi.nlm.nih.gov/pubmed/24700534](Moruz et al., GradientOptimizer: An open-source graphical environment for calculating optimized gradients in reversed-phase
liquid chromatography, Proteomics, 06/2014; 14) for details.
    
Heatmap score [MS<sup>2</sup> Scans: TopN over RT]: Rewards uniform (function Uniform) TopN events over time.
",
    workerFcn = function(.self, df_msmsScans)
    {
      ## completeness check
      stopifnot(.self$checkInput(c("fc.raw.file", "retention.time", "scan.event.number"), colnames(df_msmsScans)))
      
      ## maximum scan event over time
      DFmse = ddply(df_msmsScans, c("fc.raw.file"), function(x) {
        ## sort by RT
        if (is.unsorted(x$retention.time))
        { ## should not happen, but just to make sure
          x = x[, order(x$retention.time)]
        }
        ## only take the highest scan event (SE) in a series
        maxN = getMaxima(x$scan.event.number, thresh_rel = 0.0)
        df.max = x[maxN,]
        ## use median of that within one RT bin
        meanN = ddply(df.max, "rRT", function(xi) data.frame(topN = median(xi$scan.event.number)))
        return (meanN)    
      })
      
      lpl =
        byXflex(DFmse, DFmse$fc.raw.file, 6, plot_TopNoverRT, sort_indices = FALSE)
      
      ## QC measure for smoothness of TopN over RT
      qc_TopNRT = ddply(DFmse, "fc.raw.file", function(x) data.frame(val = qualUniform(x$topN)))
      colnames(qc_TopNRT)[colnames(qc_TopNRT) == "val"] = .self$qcName
      
      return(list(plots = lpl, qcScores = qc_TopNRT))
    }, 
    qcCat = "LC", 
    qcName = "MS^2*Scans:~TopN~over~RT", 
    orderNr = 0120
  )
    return(.self)
  })
)

#####################################################################

qcMetric_MSMSScans_IonInjTime =  setRefClass(
  "qcMetric_MSMSScans_IonInjTime",
  contains = "qcMetric",
  methods = list(initialize=function() {  callSuper(
    helpTextTemplate = 
      "Ion injection time score - should be as low as possible to allow fast cycles. Correlated with peptide intensity.
Note that this threshold needs customization depending on the instrument used (e.g., ITMS vs. FTMS).

Heatmap score [MS<sup>2</sup> Scans: Ion Inj Time]: Linear score as fraction of MS/MS below the threshold. 
",
    workerFcn = function(.self, df_msmsScans, threshold_iit)
    {
      ## completeness check
      stopifnot(.self$checkInput(c("fc.raw.file", "ion.injection.time"), colnames(df_msmsScans)))
      
      ## average injection time over RT
      DFmIIT = ddply(df_msmsScans, c("fc.raw.file"), function(x) {
        meanN = ddply(x, c("rRT"), function(x) data.frame(medIIT = median(x$ion.injection.time)))
        return (meanN) 
      })
      head(DFmIIT)
      ## average injection time overall
      DFmIITglob = ddply(df_msmsScans, c("fc.raw.file"), function(x) {
        return (data.frame(mean = mean(x$ion.injection.time)))
      })
      head(DFmIITglob)
      
      
      lpl =
        byXflex(DFmIIT, DFmIIT$fc.raw.file, 6, plot_IonInjectionTimeOverRT, sort_indices = FALSE,
                stats = DFmIITglob,
                extra_limit = threshold_iit)
      
      
      ## QC measure for injection times below expected threshold
      DFmIIT_belowThresh = ddply(df_msmsScans, c("fc.raw.file"), function(x) {
        return (data.frame(belowThresh_IIT = sum(x$ion.injection.time < threshold_iit, na.rm = TRUE) / nrow(x)))
      })
      head(DFmIIT_belowThresh)
      qc_IIT = ddply(DFmIIT_belowThresh, "fc.raw.file", 
                     function(x) data.frame(val = qualLinThresh(x$belowThresh_IIT, t = 1)))
      colnames(qc_IIT)[colnames(qc_IIT) == "val"] = .self$qcName
      
      return(list(plots = lpl, qcScores = qc_IIT))
    }, 
    qcCat = "MS", 
    qcName = "MS^2*Scans:~Ion~Inj~Time", 
    orderNr = 0240
  )
    return(.self)
  })
)

#####################################################################

qcMetric_MSMSScans_TopN =  setRefClass(
  "qcMetric_MSMSScans_TopN",
  contains = "qcMetric",
  methods = list(initialize=function() {  callSuper(
    helpTextTemplate = 
      "Reaching TopN on a regular basis indicates that all sections of the LC gradient
deliver a sufficient number of peptides to keep the instrument busy. This metric somewhat summarizes 'TopN over RT'.

Heatmap score [MS<sup>2</sup> Scans: TopN high]: rewards if TopN was reached on a regular basis (function qualHighest)
",
    workerFcn = function(.self, df_msmsScans)
    {
      ## completeness check
      stopifnot(.self$checkInput(c("fc.raw.file", "scan.event.number"), colnames(df_msmsScans)))
      
      ## check if scan.event.number requires fixing
      ## (e.g. when MS3 events are recorded between MS2 events, there are gaps in the numbering)
      ## we close the gaps by requiring consecutive scan event numbers in MS2
      scan.events = df_msmsScans[, c("scan.event.number", "fc.raw.file")]
      while (TRUE) { ## should be at most max(scan.even.number) iterations
        se_pos = 1 + which(diff(scan.events$scan.event.number) > 1) ## position of gaps>1
        if (length(se_pos) == 0) break;
        scan.events$scan.event.number[se_pos] = scan.events$scan.event.number[se_pos] - 1
      }
      DFc = ddply(scan.events, c("scan.event.number", "fc.raw.file"), function(x) data.frame(n = length(x$scan.event.number)))
      dfc.ratio = ddply(DFc, "fc.raw.file", function(x, maxn)
      {
        ## sort x by scan event
        event_count = x$n
        ## verify its monotonically increasing
        if (is.unsorted(rev(event_count))) {
          #print(x)
          stop("Scan event distribution is not monotonically increasing!")
        } 
        ## verify that there are no gaps
        if (max(x$scan.event.number) != nrow(x)) {
          #print(x)
          stop("Scan event distribution has unexpected holes...!")
        }
        
        event_pre = c(event_count[-1], 0)
        event_diff = event_count - event_pre
        
        ## build new DF of fixed length
        sn = x$scan.event.number
        if (max(sn) < maxn) 
        {
          event_diff = c(event_diff, rep(0, maxn-max(sn)))
          sn = c(sn, (max(sn)+1):maxn)
        }
        DF.new = data.frame(scan.event.number = sn, n = event_diff)
        return (DF.new)
      }, maxn = max(DFc$scan.event.number))
      head(dfc.ratio)
      
      lpl =
        byXflex(dfc.ratio, dfc.ratio$fc.raw.file, 9, plot_TopN, sort_indices = FALSE)
      
      ## QC measure for always reaching the maximum TopN
      maxTopN = max(dfc.ratio$scan.event.number)
      qc_TopN = ddply(dfc.ratio, "fc.raw.file", function(x) data.frame(val = qualHighest(x$n, maxTopN)))
      colnames(qc_TopN)[colnames(qc_TopN) == "val"] = .self$qcName
      
      return(list(plots = lpl, qcScores = qc_TopN))
    }, 
    qcCat = "MS", 
    qcName = "MS^2*Scans:~TopN~high", 
    orderNr = 0350
  )
    return(.self)
  })
)

#####################################################################

qcMetric_MSMSScans_TopNID =  setRefClass(
  "qcMetric_MSMSScans_TopNID",
  contains = "qcMetric",
  methods = list(initialize=function() {  callSuper(
    helpTextTemplate = 
      "Looking at the 
identification rates per scan event (i.e. the MS/MS scans after a survey scan) can give 
hints on how well scheduled precursor peaks could be fragmented and identified.
If performance drops for the later MS/MS scans, then the LC peaks are probably not wide enough to deliver
enough eluent or the intensity threshold to trigger the MS/MS event should be lowered (if LC peak is already over),
or increased (if LC peak is still to weak to collect enough ions).

Heatmap score [MS<sup>2</sup> Scans: TopN ID over N]: Rewards uniform identification performance across all scan events.
",
    workerFcn = function(.self, df_msmsScans)
    {
      ## completeness check
      stopifnot(.self$checkInput(c("fc.raw.file", "scan.event.number", "identified"), colnames(df_msmsScans)))
      
      DF = ddply(df_msmsScans, c("scan.event.number", "identified", "fc.raw.file"), summarise, n = length(scan.event.number))
      # try KS on underlying data instead of using qualUniform()
      #   DF2= ddply(df_msmsScans, "fc.raw.file", function(rf){
      #     cat(class(rf))
      #     cat(rf$fc.raw.file[1])
      #     idx_p = rf$identified=="+"
      #     cat(length(idx_p) %+% "\n")
      #     kk = ks.test(rf$scan.event.number[idx_p], rf$scan.event.number[-idx_p])
      #     cat(kk$p.value)
      #     kk$statistic  # = 'D' ,,  p.value is much smaller (~0)
      #   })
      #   --> fail, 'D' and p-values are too low
      df.ratio = ddply(DF, c("scan.event.number", "fc.raw.file"), function(x)
      {
        xp = xm = 0
        if ("+" %in% x$identified) xp = x$n[x$identified=="+"]
        if ("-" %in% x$identified) xm = x$n[x$identified=="-"]
        ratio = xp * 100 / sum(xp, xm)
        return (data.frame(ratio = ratio, count = sum(x$n)))
      })
      head(df.ratio)
      
      lpl = byXflex(df.ratio, df.ratio$fc.raw.file, 9, plot_ScanIDRate, sort_indices = FALSE)
      
      ## QC measure for constantly identifiying peptides, irrespective of scan event number
      ## -- we weight scan events by their number of occurence
      qc_TopN_ID = ddply(df.ratio, "fc.raw.file", function(x) data.frame(val = qualUniform(x$ratio, x$count)))
      colnames(qc_TopN_ID)[colnames(qc_TopN_ID) == "val"] = .self$qcName
      
      return(list(plots = lpl, qcScores = qc_TopN_ID))
    }, 
    qcCat = "MS", 
    qcName = "MS^2*Scans:~TopN~ID~over~N", 
    orderNr = 0380
  )
    return(.self)
  })
)

