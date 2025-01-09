#####################################################################

qcMetric_EVD_UserContaminant =  setRefClass(
  "qcMetric_EVD_UserContaminant",
  contains = "qcMetric",
  methods = list(initialize=function() {  callSuper(  
    helpTextTemplate = 
      "User defined contaminant plot based on peptide intensities and counts.
Usually used for Mycoplasma detection, but can be used for an arbitrary (set of) proteins.

All proteins (and their peptides) which contain the search string from the YAML file are considered contaminants. 
The contaminant's search string is searched in the full FASTA header in proteinGroups.txt.
If proteinGroups.txt is not available/found,
only protein identifiers can be considered. The search realm used is given in the plot subtitle.
You should choose the contaminant name to be distinctive.
Only peptides belonging to a single protein group are considered when computing the fractions (contaminant vs. all),
since peptides shared across multiple groups are potentially false positives.

Two abundance measures are computed per Raw file:

   - fraction of contaminant intensity (used for scoring of the metric)
   - fraction of contaminant spectral counts (as comparison; both should be similar)

If the intensity fraction exceeds the threshold (indicated by the dashed horizontal line) a contamination is assumed. 

For each Raw file exceeding the threshold an additional plot giving cumulative Andromeda peptide 
score distributions is shown.
This allows to decide if the contamination is true. Contaminant scores
should be equally high (or higher), i.e. to the right, compared to the sample scores.
Each graph's subtitle is augmented with a p-value of the Kologorov-Smirnoff test of this data
(Andromeda scores of contaminant peptides vs. sample peptides).
If the p-value is high, there is no score difference between the two peptide populations.
In particular, the contaminant peptides are not bad-scoring, random hits.
These p-values are also shown in the first figure for each Raw file. Note that the p-value is purely based
on Andromeda scores and is independent of intensity or spectral counts.
    

Heatmap score [EVD: Contaminant <name>]: boolean score, i.e. 0% (fail) if the intensity threshold was exceeded; otherwise 100% (pass).
",
    workerFcn = function(.self, df_evd, df_pg, lst_contaminants)
    {
      #lst_contaminants = yaml_contaminants
      ## completeness check
      ## PG is either missing, or has the correct data
      if (!is.null(df_pg) | !checkInput(c("id", "fasta.headers"), df_pg)) return()
      ## "score" might not be present (e.g. missing in MQ 1.0.13.13)
      if (!checkInput(c("protein.group.ids", "type", "intensity", "fc.raw.file"),df_evd)) return()

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
        
        if (is.null(df_pg)) {
          ## only search in protein IDs
          pg_id = df_evd$protein.group.ids[grep(ca, df_evd$proteins)] 
          ## this could be multiple PGs ("PG1; PG2") per cell, but we require unique peptides below, so its not a problem
          search_realm = "protein name only"
        } else {
          ## search in FASTA headers and protein IDs
          pg_id = df_pg$id[c(grep(ca, df_pg$fasta.headers, ignore.case = TRUE),
                             grep(ca, df_pg$protein.ids, ignore.case = TRUE))]
          search_realm = "full FASTA header"
        }
        
        
        if (length(pg_id) > 0)
        {
          ## we might or might not have found something... we plot it anyways, so the user can be sure that we searched for it
          not_found = FALSE
          
          ## find peptides which only have one group (ignoring razor peptides where we cannot be sure)
          evd_uniqueGroup = !grepl(";", df_evd$protein.group.ids)
          ## do not trust MBR here. We want real evidence!
          evd_realMS = !grepl("MATCH", df_evd$type)
          ## for each Raw file: find unique peptides of our contaminant
          cont_data.l = plyr::dlply(df_evd[evd_uniqueGroup & evd_realMS, ], "fc.raw.file",
                              function(x) {
                                if (length(grep(";", x$protein.group.ids))) stop("more than one proteinGroup for supposedly unique peptide...")
                                
                                x$idx_cont = x$protein.group.ids %in% pg_id
                                
                                sc = sum(x$idx_cont) / nrow(x) * 100
                                int = sum(as.numeric(x$intensity[x$idx_cont]), na.rm = TRUE) / sum(as.numeric(x$intensity), na.rm = TRUE) * 100
                                
                                above.thresh = (sc > ca_thresh) | (int > ca_thresh)
                                cont_scoreECDF = NULL;
                                if ("score" %in% colnames(x)) {
                                  cont_scoreECDF = plyr::ddply(x, "idx_cont", function(xx) {
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
          cont_data = plyr::ldply(cont_data.l, function(l) { l$cont_data })
          cont_data.long = reshape2::melt(cont_data, id.vars="fc.raw.file")
          
          # 
          # old: not_found = all(cont_data.long$value[cont_data.long$variable == "above.thresh"] == FALSE)
        }
        
        if (not_found)
        { ## identifier was not found in any sample
          pl_cont = ggText("EVD: Contaminants",
                           paste0("Contaminant '", ca, "' was not found in any sample.\n\nDid you use the correct database?"),
                           "red")
          lpl = append(lpl, list(pl_cont))
        } else {
          ## plot User-Contaminants
          lpl_i = byXflex(data = cont_data.long, indices = cont_data.long$fc.raw.file, subset_size = 120,
                          FUN = plot_ContUser, sort_indices = TRUE,
                          name_contaminant = ca, extra_limit = ca_thresh, subtitle = paste("search realm:", search_realm))
          lpl = append(lpl, lpl_i)
          
          ## plot Andromeda score distribution of contaminant vs. sample
          pl_andr = plyr::llply(cont_data.l, function(l)
          {
            if (l$cont_data$above.thresh == FALSE ||
                is.null(l$cont_scoreECDF))
            {
              return(NULL) ## some entries might be skipped (not significant)
            } 
            p = plot_ContUserScore(l$cont_scoreECDF, l$cont_data$fc.raw.file, l$cont_data$score_KS)
            #print(p)
            return(p)
          })
          pl_andr_nonNull = plyr::compact(pl_andr) ## remove 'NULL' entries from plot list
          lpl = append(lpl, pl_andr_nonNull)
          
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
    qcName = "EVD:User~Contaminant~(%s)",
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
      "Peptide precursor intensity per Raw file from evidence.txt WITHOUT match-between-runs evidence.
Low peptide intensity usually goes hand in hand with low MS/MS identifcation rates and unfavourable signal/noise ratios,
which makes signal detection harder. Also instrument acquisition time increases for trapping instruments.

Failing to reach the intensity threshold is usually due to unfavorable column conditions, inadequate 
column loading or ionization issues. If the study is not a dilution series or pulsed SILAC experiment, we 
would expect every condition to have about the same median log-intensity (of 2<sup>%1.1f</sup>).
The relative standard deviation (RSD) gives an indication about reproducibility across files and should be below 5%%.

Depending on your setup, your target thresholds might vary from PTXQC's defaults.
Change the threshold using the YAML configuration file.

Heatmap score [EVD: Pep Intensity (>%1.1f)]: 
  Linear scale of the median intensity reaching the threshold, i.e. reaching 2<sup>21</sup> of 2<sup>23</sup> gives score 0.25.
",
    workerFcn = function(.self, df_evd, thresh_intensity)
    {
      ## completeness check
      if (!checkInput(c("fc.raw.file", "intensity", "contaminant"), df_evd)) return()
      
      ## update helpText
      .self$helpText = sprintf(.self$helpTextTemplate, thresh_intensity, thresh_intensity)
      
      medians_pep = plyr::ddply(df_evd[ , c("fc.raw.file", "intensity")], "fc.raw.file",
                          function(x) data.frame(med = log2(quantile(x$intensity, probs=0.5, na.rm = TRUE))))
      
      int_dev_pep = RSD((medians_pep$med))
      int_dev.s = pastet("INT RSD [%]", round(int_dev_pep, 3))
      lpl = boxplotCompare(data = df_evd[, c("fc.raw.file", "intensity", "contaminant")],
                           log2 = TRUE, 
                           mainlab = "EVD: peptide intensity distribution",
                           ylab = expression(log[2]*" intensity"),
                           sublab = paste0("RSD ", round(int_dev_pep, 1),"% (expected < 5%)\n"),
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
    qcName = "EVD:~Peptide~Intensity~(\">%1.1f\")", 
    orderNr = 0030
  )
    return(.self)
  })
)


#####################################################################

qcMetric_EVD_ReporterInt =  setRefClass(
  "qcMetric_EVD_ReporterInt",
  contains = "qcMetric",
  methods = list(initialize=function() {  callSuper(
    helpText = 
      "ITRAQ/TMT reporter intensity violin plots of all PSMs for each channel and Raw file.
The second subplot shows labeling efficiency (LE), i.e the fraction of PSMs with non-zero abundance (100% = full labeling of all PSMs; 0% = no reporter ions at all). This is used for heatmap scoring. See below.

There is a similar 'Experimental Group' based metric/plot based on proteins.txt.

PTXQC uses isotope-corrected intensities (eliminating channel carry-over) to allow for detection of empty channels, e.g. due to mis-labeling.
If MaxQuant did no isotope correction (i.e. corrected and uncorrected channels are equal), 
the plot title will show a warning. The scores are too optimistic in this case (since carry-over will be mistaken for actual signal).

Note: global labelling efficiency can only be judged indirectly with this metric, since isobaric reporters where set as
      fixed modification. Thus, MaxQuant. will only identify labeled peptides in the first place.
      Observing only very few peptides (see peptide count metric), is a good indicator.
      However, if only the labeling of a few channels failed, this will be noticable here!

Labeling can still be poor, even though identification was successful.

A labeling efficiency (LE) is computed per Raw file AND channel as: the percentage of PSMs which have non-zero reporter intensity.
Ideally LE reaches 100 percent (all peptides have an intensity in the channel; biological missingness ignored).

Heatmap score: minimum labeling efficiency per Raw file across all channels.
I.e. for 4-plex ITRAQ and two Raw files, there will be 8 labeling efficiency (LE) values. 
Each Raw file is now scored by the minimum LE of all its 4 channels.
",
    workerFcn=function(.self, df_evd)
    {
      ## completeness check
      if (!checkInput(c("fc.raw.file"), df_evd)) return()
      ## check if reporter.intensity.0... is present
      cols_reporter = grepv("^reporter.intensity.corrected.[0-9]", colnames(df_evd));
      cols_reporter.nc = grepv("^reporter.intensity.[0-9]", colnames(df_evd));
      if(length(cols_reporter) <= 1 || length(cols_reporter.nc) <= 1) {
        message("Info: Two reporter.intensity and two reporter.intensity.corrected columns are needed for metric ReporterIntensity.")
        return()}
      ## check if correction was done at all
      if (all(df_evd[1:1000, cols_reporter] == df_evd[1:1000, cols_reporter.nc], na.rm = TRUE))
      {
        title_subtext = "Warning: MaxQuant did NO isotope correction";  
        title_color = "red"
      } else {
        title_subtext = "";  
        title_color = "black"
      }
      g_title = "EVD: Reporter Intensities"  
      ## use data.table for aggregation, its MUCH faster than ddply() and uses almost no extra memory
      df_reps = reshape2::melt(df_evd[, c("fc.raw.file", cols_reporter)], 
                     id.vars ="fc.raw.file", 
                     value.name = "intensity",
                     variable.name = "channel")
      head(df_reps)
      dt_reps = data.table::data.table(df_reps)

      ## do NOT remove -inf and NA's and 0's -- we need them to count labeling-efficiency (#entries with intensity > 0 vs. ALL)

      ## rename 'reporter.intensity.corrected.0' to '0'
      dt_reps$channel = substring(dt_reps$channel, nchar('reporter.intensity.corrected.') + 1)
      ## invert the channel order (so that channel 0 is highest, i.e. appears on top in plot)
      dt_reps$channel = factor(dt_reps$channel, levels = sort(unique(dt_reps$channel), decreasing = TRUE))
      head(dt_reps)
      
      ylims_minmax = range(dt_reps$intensity[dt_reps$intensity>0])
      if (is.na(ylims_minmax[1])) {ylims_minmax = range(1)} ## data is all 0! Make range at least 1, so log-range does not crash
      
      fcn_boxplot_internal = function(data, title_subtext = title_subtext, title_color = title_color) 
      {
        ### first subplot (distribution of intensities)
        data_noZero = data[data$intensity!=0,]
        pl = ggplot(data=data_noZero) +
          geom_violin(aes(x = .data$fc.raw.file, 
                          y = .data$intensity,  
                          color = .data$channel,
                          fill = .data$channel
                         )) +
          xlab("") + 
          ylab("reporter intensity (zeros removed)") +
          guides(#alpha = guide_legend(title="Label Eff"), 
                 fill = guide_legend(reverse = TRUE), ## inverse label order, so that channel 0 is on top
                 color = guide_none()) + 
          theme(axis.text.x = element_text(angle = 45, vjust = 0.5), 
                legend.position = "right",
                plot.title = element_text(color = title_color)) +
          ggtitle(g_title, title_subtext) + 
          PTXQC:::scale_x_discrete_reverse(unique(data$fc.raw.file)) +
          scale_y_log10(limits = ylims_minmax) +
          coord_flip() 
        #pl
        
        ylims = data[, { #limits = boxplot.stats(intensity, coef = 0.7)$stats;
                          list(labEff_PC = sum(intensity > 0, na.rm = TRUE) / (.N)) 
                       }, by = c("fc.raw.file", "channel")]
        
        
        ### second subplot (labeling efficiency)
        pl_eff = ggplot(data = ylims) + geom_bar(aes(x = .data$fc.raw.file,
                                                     y = .data$labEff_PC * 100,
                                                     fill = .data$channel), 
                                                 stat = "identity",
                                                 position = "dodge") + 
          xlab("") + 
          ylab("labelling efficiency (%)") +
          ylim(0, 100) +
          guides(fill = guide_legend(reverse = TRUE), ## inverse label order, so that channel 0 is on top
                 color = guide_none()) + 
          theme(legend.position = "right") +
          ggtitle("Fraction of Non-Zero Intensities", "") + 
          PTXQC:::scale_x_discrete_reverse(unique(ylims$fc.raw.file)) +
          coord_flip() 
        #pl_eff
        pl_both = gridExtra::grid.arrange(pl, pl_eff, ncol=2)
        #print(pl)
        return(pl_both)
      }
      channel_count = length(cols_reporter)
      lpl = byXflex(data = dt_reps, indices = dt_reps$fc.raw.file, subset_size = round(40 / channel_count), 
                    sort_indices = TRUE, FUN = fcn_boxplot_internal, title_subtext = title_subtext, title_color = title_color)
      lpl
      # heatmap scoring
      ## .. take min score over all channels
      ylims = dt_reps[, { #limits = boxplot.stats(intensity, coef = 0.7)$stats;
        list(labEff_PC = sum(intensity > 0, na.rm = TRUE) / (.N)) 
      }, by = c("fc.raw.file", "channel")]
      qcScore = ylims[, list(score_min = min(labEff_PC)), by=c("fc.raw.file")]
      colnames(qcScore) = c("fc.raw.file", .self$qcName)
  
      ## add manual title, since we return a grid.arrange() where automatic extraction is hard
      return(list(plots = lpl, qcScores = qcScore, title = rep(list(g_title), length(lpl))))
    }, 
    qcCat = "prep", 
    qcName = "EVD:~Reporter~intensity", 
    orderNr = 0031  ## should not show up in heatmap
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
are shown, so the user can judge the gain that MBR provides. Here, 'transferred (exclusive)' means that this protein group
has peptide evidence which originates only from transferred peptide IDs. The quantification is (of course) always from the
local Raw file. 
Proteins in the 'genuine + transferred' category have peptide evidence from within the Raw file by MS/MS, but at the same time
also peptide IDs transferred to this Raw file using MBR were used. It is not unusual to see the 'genuine + transferred' category be the 
rather large, since a protein group usually has peptide evidence from both sources.
To see of MBR worked, it is better to look at the two MBR-related metrics.

If MBR would be switched off, you can expect to see the number of protein groups corresponding to 'genuine (exclusive)' + 'genuine + transferred'.
In general, if the MBR gain is low and the MBR scores are bad (see the two MBR-related metrics),
MBR should be switched off for the Raw files which are affected (could be a few or all).

Heatmap score [EVD: Prot Count (>%1.0f)]: Linear scoring from zero. Reaching or exceeding the target threshold gives a score of 100%%.
",
    workerFcn = function(.self, df_evd, df_evd_tf, thresh_protCount)
    {
      ## completeness check

      req_cols = c("fc.raw.file", "protein.group.ids", "is.transferred")
      
      if (!checkInput(req_cols, df_evd)) return()

      .self$helpText = sprintf(.self$helpTextTemplate, thresh_protCount)
      
      protC = getProteinCounts(rbind(df_evd[,req_cols], df_evd_tf[, req_cols]))
      
      protC$block = factor(assignBlocks(protC$fc.raw.file, 30))
      
      max_prot = max(unlist(plyr::dlply(protC, "fc.raw.file", function(x) sum(x$counts))))
      ## average gain in percent
      reportMTD = nrow(df_evd_tf) > 0
      gain_text = ifelse(reportMTD, sprintf("MBR gain: +%.0f%%", mean(protC$MBRgain, na.rm = TRUE)), "")
      
      lpl = plyr::dlply(protC, "block", .fun = function(x)
      {
        p = plot_CountData(data = x, 
                           y_max = max(thresh_protCount, max_prot)*1.1,
                           thresh_line = thresh_protCount,
                           title = c("EVD: ProteinGroups count", gain_text))
        #print(p)
        return (p)
      })
      
      ## QC measure for protein ID performance
      qc_protc = plyr::ddply(protC, "fc.raw.file", function(x){
        if (nrow(x) == 3 && length(grep("^genuine", x$category))!= 2){
          stop("expected two categories to start with 'genuine...'")
        }
        r = data.frame(genuineAll = sum(x$counts[grep("^genuine", x$category)]))
        return (r)
      })
      cname = sprintf(.self$qcName, thresh_protCount)
      qc_protc[,cname] = qualLinThresh(qc_protc$genuineAll, thresh_protCount)
      qcScore = qc_protc[, c("fc.raw.file", cname)]
      
      ## add mzQC metric
      template_proteinCount = rmzqc::getQualityMetricTemplate("MS:1002406") # count of identified clusters
      mzqc = lapply(1:nrow(qc_protc), function(row){
          out = template_proteinCount$copy();
          out$value = qc_protc$genuineAll[row];
          #cat(row, " ", qc_protc$genuineAll[row], " ",out$value, "\n");
          return(out) })
      names(mzqc) = qc_protc$fc.raw.file
      ## done
      return(list(plots = lpl, qcScores = qcScore, mzQC = mzqc))
    }, 
    qcCat = 'general', 
    qcName = "EVD:~Protein~Count~(\">%1.0f\")", 
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
      "Number of unique (i.e. not counted twice) peptide sequences including modifications (after FDR) per Raw file. A configurable target threshold is indicated as dashed line.

If MBR was enabled, three categories ('genuine (exclusive)', 'genuine + transferred', 'transferred (exclusive)'
are shown, so the user can judge the gain that MBR provides.    
Peptides in the 'genuine + transferred' category were identified within the Raw file by MS/MS, but at the same time
also transferred to this Raw file using MBR. This ID transfer can be correct (e.g. in case of different charge states),
or incorrect -- see MBR-related metrics to tell the difference.
Ideally, the 'genuine + transferred' category should be rather small, the other two should be large.

If MBR would be switched off, you can expect to see the number of peptides corresponding to 'genuine (exclusive)' + 'genuine + transferred'.
In general, if the MBR gain is low and the MBR scores are bad (see the two MBR-related metrics),
MBR should be switched off for the Raw files which are affected (could be a few or all). 

Heatmap score [EVD: Pep Count (>%1.0f)]: Linear scoring from zero. Reaching or exceeding the target threshold gives a score of 100%%.
",
    workerFcn = function(.self, df_evd, df_evd_tf, thresh_pepCount)
    {
      ## completeness check

      req_cols = c("fc.raw.file", "modified.sequence", "is.transferred")
      if (!checkInput(req_cols, df_evd)) return()
      if (nrow(df_evd_tf)>0 & !checkInput(req_cols, df_evd_tf)) return()

      .self$helpText = sprintf(.self$helpTextTemplate, thresh_pepCount)
      
      pepC = getPeptideCounts(rbind(df_evd[, req_cols], df_evd_tf[, req_cols]))
      pepC$block = factor(assignBlocks(pepC$fc.raw.file, 30))
      
      max_pep = max(unlist(plyr::dlply(pepC, "fc.raw.file", function(x) sum(x$counts))))
      ## average gain in percent
      reportMTD = any(df_evd$is.transferred)
      gain_text = ifelse(reportMTD, sprintf("MBR gain: +%.0f%%", mean(pepC$MBRgain, na.rm = TRUE)), "")
      
      lpl = plyr::dlply(pepC, "block", .fun = function(x)
      {
        p = plot_CountData(data = x, 
                           y_max = max(thresh_pepCount, max_pep)*1.1,
                           thresh_line = thresh_pepCount,
                           title = c("EVD: Peptide ID count", gain_text))
        #print(p)
        return (p)
      })
      
      ## QC measure for peptide ID performance
      qc_pepc = plyr::ddply(pepC, "fc.raw.file", function(x){
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
    qcName = "EVD:~Peptide~Count~(\">%1.0f\")", 
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
      if (!checkInput(c("retention.time", "retention.length", "fc.raw.file"), df_evd)) return()
      
      ## compute some summary stats before passing data to ggplot (performance issue for large experiments) 
      df_evd.m.d = plyr::ddply(df_evd[,c("retention.time", "retention.length", "fc.raw.file")], "fc.raw.file", .fun = peakWidthOverTime)
      head(df_evd.m.d)
      ## median peak width
      df_evd.m.d_avg = plyr::ddply(df_evd[,c("retention.length","fc.raw.file")], "fc.raw.file", .fun = function(x) {
        #fcr = as.character(x$fc.raw.file[1])
        #cat(fcr)
        m = median(x$retention.length, na.rm = TRUE);
        return(data.frame(median = m))
      })
      df_evd.m.d_avg$fc.raw.file_aug = paste0(df_evd.m.d_avg$fc.raw.file, " (~", round(df_evd.m.d_avg$median, 1)," min)")
      .self$outData[["avg_peak_width"]] = df_evd.m.d_avg
      
      ## augment Raw filename with avg. RT peak width
      df_evd.m.d$fc.raw.file = plyr::mapvalues(df_evd.m.d$fc.raw.file, df_evd.m.d_avg$fc.raw.file, df_evd.m.d_avg$fc.raw.file_aug)
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
      l_dists = plyr::dlply(df_evd[,c("retention.length", "fc.raw.file")], "fc.raw.file", function(x) return(x$retention.length))
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
The ratio of 'green' vs. 'green+red' peptides is annotated using 'sc: ' (for 'score') in the plot subtitles. High values are better (green peptides dominate).

If moving 'red' dots to the horizontal zero-line (to make them green) requires large RT shifts, then increasing the alignment
search window might help MaxQuant to find a better alignment.


Heatmap score [EVD: MBR Align]: ratio of 'green' vs. 'green+red' peptides 
",
    workerFcn = function(.self, df_evd, tolerance_matching, raw_file_mapping)
    {
      ## completeness check
      if (!checkInput(c("type", "calibrated.retention.time", "retention.time.calibration", "id", "raw.file", "modified.sequence", "charge"), df_evd)) return()
      
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
        txt_subtitle = paste("alignment reference:", gsub("\\", "/", refRaw, fixed = TRUE)) ## subtitles in ggplot must not contain '\'
        evd_has_fractions = FALSE
      }
      
      lpl = list()
      qcScore = .self$qcScores
      
      if (!evd_has_fractions && length(refRaw) == 0) {
        lpl[[1]] = ggText("EVD: Alignment check", paste0("Cannot find a reference Raw file!\nPlease report this as a 'bug'!"))
      } else {
        if (!evd_has_fractions & (length(refRaw) != 1))
        {
          refRaw = refRaw[1] ## take the first
          warning(paste0("Cannot find a unique reference Raw file (files: ", paste(refRaw, collapse=", "), "). Picking the first."), immediate. = TRUE)
        }
        ## find RT curve based on genuine 3D peaks (should be flat)
        d_alignQ = alignmentCheck(df_evd[(df_evd$type %in% c("MULTI-MSMS")), 
                                         c("calibrated.retention.time", 
                                           "id", "raw.file", col_fraction, "modified.sequence", "charge")], 
                                  referenceFile = refRaw)
        ## augment more columns
        d_alignQ$retention.time.calibration = df_evd$retention.time.calibration[match(d_alignQ$id, df_evd$id)]
        
        if (diff(range(na.omit(d_alignQ$retention.time.calibration))) < 1e-5)
        {
          txt_subtitle = paste0(txt_subtitle, " || WARNING: MaxQuant did not correct RTs in any way!");
          warning("EVD MBRAlign: MaxQuant did not correct RTs in any way, despite MBR=on")
        }
        
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
          evd_RT_t$fc.raw.file_ext = plyr::mapvalues(evd_RT_t$fc.raw.file, qcAlign$fc.raw.file, qcAlign$newlabel)
          
          evd_RT_t$RTdiff_in = c("green", "red")[(abs(evd_RT_t$rtdiff) > tolerance_matching)+1]
          
          ## plot alignment result
          y_lim = quantile(c(evd_RT_t$rtdiff, evd_RT_t$retention.time.calibration), probs = c(0.01, 0.99), na.rm = TRUE) * 1.1
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
    workerFcn = function(.self, df_evd, df_evd_tf, avg_peak_width)
    {
      ## completeness check
      #stopifnot(c("...") %in% colnames(df_evd))
      if (!checkInput(c("modified.sequence"), df_evd)) return()
      
      df_evd_all = merge(df_evd, df_evd_tf, all = TRUE)
            
      ## increase of segmentation by MBR:
      ## three values returned: single peaks(%) in genuine, transferred and all(combined)
      qMBR = peakSegmentation(df_evd_all)
      head(qMBR)
      ## for groups: get their RT-spans
      ## ... genuine ID's only (as 'rtdiff_genuine') 
      ##  or genuine+transferred (as 'rtdiff_mixed'))
      ## Could be empty (i.e. no groups, just singlets) if data is really sparse ..
      qMBRSeg_Dist = idTransferCheck(df_evd_all)
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
  
Return a tree plot with a possible alignment tree.
This allows the user to judge which Raw files have similar corrected RT's (i.e. where aligned successfully).
If there are clear sub-clusters, it might be worth introducing artifical fractions into MaxQuant,
to avoid ID-transfer between these clusters (use the MBR-Align and MBR-ID-Transfer metrics to support the decision).
 
If the input contains fractions, leaf nodes will be colored accordingly.
Distinct sub-clusters should have their own color.
If not, MaxQuant's fraction settings should be optimized.
Note that introducing fractions in MaxQuant will naturally lead to a clustering here (it's somewhat circular).

Heatmap score: none.
",
    workerFcn = function(.self, df_evd)
    {
      ## completeness check

      if (!checkInput(c("type", "is.transferred", "calibrated.retention.time", "fc.raw.file", "modified.sequence", "charge"), df_evd)) return()
    
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
      if (any(df_evd$is.transferred)) {
        ## gain for each raw file: absolute gain, and percent gain
        mtr.df = plyr::ddply(df_evd, "fc.raw.file", function(x) {
          match_count_abs = sum(x$is.transferred)
          ## if only matched IDs are present, this would be 'Inf' -- we limit that to 1e4
          match_count_pc  = min(1e4, round(100*match_count_abs/(nrow(x)-match_count_abs))) ## newIDs / oldIDs
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
      "Charge distribution per Raw file. For typtic digests, peptides of charge 2 
(one N-terminal and one at tryptic C-terminal R or K residue) should be dominant.
Ionization issues (voltage?), in-source fragmentation, missed cleavages and buffer irregularities can 
cause a shift (see Bittremieux 2017, DOI: 10.1002/mas.21544).
The charge distribution should be similar across Raw files.
Consistent charge distribution is paramount for comparable 3D-peak intensities across samples.

Heatmap score [EVD: Charge]: Deviation of the charge 2 proportion from a representative Raw file ('qualMedianDist' function).
",
    workerFcn = function(.self, df_evd)
    {
      ## completeness check
      if (!checkInput(c("is.transferred", "fc.raw.file", "charge"), df_evd)) return()
      
      d_charge = mosaicize(df_evd[!df_evd$is.transferred, c("fc.raw.file", "charge")])
      lpl =
        byXflex(d_charge, d_charge$Var1, 30, plot_Charge, sort_indices = TRUE)
      
      ## QC measure for charge centeredness
      qc_charge = plyr::ddply(df_evd[!df_evd$is.transferred, c("charge",  "fc.raw.file")], "fc.raw.file", function(x) data.frame(c = (sum(x$charge==2)/nrow(x))))
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
See [Moruz 2014, DOI: 10.1002/pmic.201400036](https://pubmed.ncbi.nlm.nih.gov/24700534/) for details.

Heatmap score [EVD: ID rate over RT]: Scored using 'Uniform' scoring function, i.e. constant receives good score, extreme shapes are bad.
",
    workerFcn = function(.self, df_evd)
    {
      ## completeness check
      if (!checkInput(c("retention.time", "fc.raw.file"), df_evd)) return()
      
      raws_perPlot = 6
      
      rt_range = range(df_evd$retention.time, na.rm = TRUE)
      df_idRT = plyr::ddply(df_evd, "fc.raw.file", function(x) {
        h = hist(x$retention.time, breaks=seq(from=rt_range[1]-3, to=rt_range[2]+3, by=3), plot = FALSE)
        return(data.frame(RT = h$mid, counts = h$counts))
      })
      lpl =
        byXflex(df_idRT, df_idRT$fc.raw.file, raws_perPlot, plot_IDsOverRT, sort_indices = TRUE)
      
      ## QC measure for uniform-ness
      qcScore = plyr::ddply(df_evd[, c("retention.time",  "fc.raw.file")], "fc.raw.file", 
                      function(x) data.frame(metric = qualUniform(na.omit(x$retention.time))))
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
      
      if (!checkInput(c("fc.raw.file", "uncalibrated.mass.error..ppm."), df_evd)) return()
      
      ## for some mzTab (not recalibrated) 'mass.error..ppm.' is not there... but we only need a dummy
      if (!("mass.error..ppm." %in% colnames(df_evd))) df_evd$mass.error..ppm. = 0
      
      fix_cal = fixCalibration(df_evd, df_idrate, tolerance_sd_PCoutOfCal)
      if (is.null(fix_cal)) {
        warning("Internal error. Data missing. Skipping metric!", immediate. = TRUE)
        return()
      }
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
      qc_MS1deCal = plyr::ddply(fix_cal$df_evd, "fc.raw.file", 
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
      if (!checkInput(c("uncalibrated.mass.error..ppm.", "mass", "mass.error..ppm."), df_evd)) return()
      
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
      obs_par = plyr::ddply(fix_cal$df_evd[, c("mass.error..ppm.", "fc.raw.file")], "fc.raw.file", 
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
The 'Abundance class' models the average peptide intensity in each Raw file and is visualized using varying degrees of
transparency. It is not unusual to see samples with low sample content to have higher contamination.
If you see only one abundance class ('mid'), this means all your Raw files have roughly
the same peptide intensity distribution.

If you see less than 5 contaminants, it either means there are actually less, or that one (or more) of the shortened contaminant names
subsume multiple of the top5 contaminants (since they start with the same prefix).
    
Heatmap score [EVD: Contaminants]: as fraction of summed intensity with 0 = sample full of contaminants; 1 = no contaminants

",
    workerFcn = function(.self, df_evd)
    {
      ## completeness check
      if (!checkInput(c("intensity", "contaminant", "fc.raw.file", "proteins"), df_evd)) return()
      
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
      df_evd.cont.only = df_evd[df_evd$contaminant > 0,]

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
      qc_cont = plyr::ddply(df_evd[, c("intensity", "contaminant", "fc.raw.file")], "fc.raw.file", 
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
    
    If DIA-Data: this metric is skipped

Heatmap score [EVD: MS<sup>2</sup> Oversampling]: The percentage of non-oversampled 3D-peaks. 
",
    workerFcn = function(.self, df_evd)
    {
      ## completeness check
      if (!checkInput(c("fc.raw.file", "ms.ms.count"), df_evd)) return()
      
      ## this is DIA data -- omit metric
      if (all(range(df_evd$ms.ms.count) == c(0,0))) return()
      
      d_dups = plyr::ddply(df_evd, "fc.raw.file", function(x) {
        tt = as.data.frame(table(x$ms.ms.count), stringsAsFactors = FALSE)
        tt$Count = as.numeric(tt$Var1)
        ## remove "0", since this would be MBR-features
        tt = tt[tt$Count!=0,]
        ## summarize everything above 3 counts
        if (any(tt$Count >= 3)) {
          tt$Count[tt$Count >= 3] = "3+"
          tt = plyr::ddply(tt, "Count", function(x) data.frame(Freq=sum(x$Freq)))
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

qcMetric_EVD_CarryOver =  setRefClass(
  "qcMetric_EVD_CarryOver",
  contains = "qcMetric",
  methods = list(initialize=function() {  callSuper(  
    helpTextTemplate = 
      "Sample carryover ... not ready yet... 

Heatmap score [EVD: Carryover]: The percentage of peptide identifications whose sequence gives rise to a large 'retention span'. 
",
    workerFcn = function(.self, df_evd)
    {
      ## completeness check
      if (!checkInput(c("fc.raw.file", "ms.ms.count", "retention.length"), df_evd)) return()
      
      raws_perPlot = 6
      
      aboveThresh_fn = function(data) {
        #t = quantile(data, probs=0.75, na.rm=TRUE) + 3*IQR(data, na.rm = TRUE)
        t = median(data, na.rm=TRUE) * 10
        return(t)
      }
      
      rt_range = range(df_evd$retention.time, na.rm = TRUE)
      df_carry = plyr::ddply(df_evd, "fc.raw.file", function(x) {
        thresh = round(aboveThresh_fn(x$retention.length), 1)
        weirds = x[x$retention.length > thresh,]
        ## remove "0", since this would be MBR-features
        weirds = weirds[weirds$ms.ms.count!=0,]
        if (nrow(weirds) == 0) return(data.frame())
        weirds$fc.raw.file = paste0(weirds$fc.raw.file, " (>", thresh, " min)")
        h = hist(weirds$retention.time, breaks=seq(from=rt_range[1]-3, to=rt_range[2]+3, by=3), plot = FALSE)
        result = data.frame(RT = h$mid, counts = h$counts, fn = weirds$fc.raw.file[1])
        return(result)
      })
      #df_carry = result
      df_carry$fc.raw.file = df_carry$fn
      lpl =
        byXflex(df_carry, df_carry$fc.raw.file, raws_perPlot, plot_DataOverRT, sort_indices = TRUE, 
                title = "EVD: Peptides with wide RT span", y_lab = "# of Peptide Sequences")
      lpl
      
      
      qc_evd_carry = plyr::ddply(df_evd, "fc.raw.file", function(x) {
        thresh = aboveThresh_fn(x$retention.length);
        pc = sum(x$retention.length > thresh, na.rm = TRUE) / nrow(x)
        return (data.frame(larger_pc=pc))
      })
      
      ## QC measure for how many IDs are part of a large span
      cname = .self$qcName
      qc_evd_carry[, cname] = qualLinThresh(1 - qc_evd_carry$larger_pc)
      
      return(list(plots = lpl, qcScores = qc_evd_carry[, c("fc.raw.file", cname)]))
    }, 
    qcCat = "MS", 
    qcName = "EVD:~CarryOver", 
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
      if (!checkInput(c("fc.raw.file", "modified.sequence", "intensity"), df_evd)) return()
      
      if (('fraction' %in% colnames(df_evd)) && (length(unique(df_evd$fraction)) > 1)) {
        lpl = list(ggText("Missing Values Skipped", "Missing values calculation skipped. Fractionated data detected!"))
        return(list(plots = lpl))
      }
      
      if (length(unique(df_evd$fc.raw.file)) < 2) {
        lpl = list(ggText("Missing Values Skipped", "Need more than one Raw file!"))
        return(list(plots = lpl))
      }
      
      ## make peptides unique per Raw file
      df_u = plyr::ddply(df_evd[ , c("fc.raw.file", "modified.sequence")], "fc.raw.file",
                   function(x) {
                     return(x[!duplicated(x$modified.sequence),])
                  })
            
      global_peps = unique(df_u$modified.sequence)
      global_peps_count = length(global_peps)
      
      ## percent identified in each Raw file
      pep_set = plyr::ddply(df_u[ , c("fc.raw.file", "modified.sequence")], "fc.raw.file",
                      function(x) {
                        score = 100*length(intersect(global_peps, x$modified.sequence)) / global_peps_count
                        return(data.frame(idFraction = score))
                      })
      
      lpl = byXflex(pep_set, pep_set$fc.raw.file, subset_size = 50, FUN = function(dx) {
        p = ggplot(dx) + 
          geom_bar(aes(x = .data$fc.raw.file, y = .data$idFraction), stat = "identity") +
          ggtitle("[experimental] EVD: Non-Missing Peptides", "compared to all peptides seen in experiment") +
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
      
      p = ggplot(tbl_smry, aes(x = .data$x, y = .data$FreqCum)) + 
        geom_line() +
        geom_point() +
        ggtitle("[experimental] EVD: Non-missing by set") +
        xlab("Minimum # Raw files") +
        ylab("Fraction of total peptides [%]") +
        ylim(0, 100)
      lpl[["missingCul"]] = p
      
      ## intensity distribution of missing values
      df_evd$logInt = log2(df_evd$intensity)
      
      lpl_dens = byXflex(df_evd[, c("modified.sequence", "fc.raw.file", "logInt")], df_evd$fc.raw.file,
                         subset_size = 5, FUN = function(dx) {
        d_mat = reshape2::dcast(dx, modified.sequence ~ fc.raw.file, fun.aggregate = mean, value.var = "logInt")
        
        ## ... normalization factors
        d_mat_mult = sapply(2:ncol(d_mat), function(x) {
          mult = mean(d_mat[, x] / d_mat[, 2], na.rm = TRUE)
          return(mult)
        })
        df_mult = data.frame(fc.raw.file = colnames(d_mat)[-1], mult = d_mat_mult)
        ## .. normalize data
        d_mat_n = d_mat
        d_mat_n[, -1] = sweep( d_mat_n[, -1, drop=FALSE], 2, d_mat_mult, '/')
        ## 
        head(d_mat_n)
        ## find impute value
        pep_mean = rowMeans(d_mat_n[, -1, drop=FALSE], na.rm = TRUE)
        df_missing = plyr::ddply(df_mult, "fc.raw.file", function(x) {
          ## get set of missing values
          values = pep_mean[is.na(d_mat_n[, as.character(x$fc.raw.file)])]
          ## de-normalize (back to old intensity range)
          values = values * x$mult
          return(data.frame(missingVals = values))
        })
        head(df_missing)
        
        pl = ggplot(df_missing, aes(x = .data$missingVals, col = .data$fc.raw.file, fill = .data$fc.raw.file)) + 
          geom_area(position = position_dodge(width=0), binwidth = 0.5, stat="bin", alpha=0.5) +
          geom_freqpoly(data = dx, aes(x = .data$logInt, col = .data$fc.raw.file), binwidth = 0.5, linewidth = 1.2) +
          xlab("Intensity [log2]") +
          ggtitle(" [experimental] EVD: Imputed Peptide Intensity Distribution of Missing Values") +
          scale_fill_manual(values = rep(RColorBrewer::brewer.pal(6,"Accent"), times=40), guide = guide_legend("")) +
          scale_colour_manual(values = rep(RColorBrewer::brewer.pal(6,"Accent"), times=40), guide = "none")
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

qcMetric_EVD_UpSet =  setRefClass(
  "qcMetric_EVD_UpSet",
  contains = "qcMetric",
  methods = list(initialize=function() {  callSuper(    
    helpTextTemplate = 
      'The metric shows an upSet plot based on the number of modified peptide sequences per Raw file, intersected or merged with other Raw files (see below for details).<br>

If the number of Raw files is >=6, only the `distinct` plot is generated (the other two are skipped for performance reasons).
<a href="https://raw.githubusercontent.com/cbielow/PTXQC/master/inst/reportTemplate/modes_UpSet.png" target="_blank" rel="noopener"><span>See here for an example plot showing how the set size is computed</span> </a>.

Definition: An `active set` is the set of black dots in a column of the plot -- as opposed to the grey dots (you will understand when you see it).

<p>
<b>distinct:</b> shows the number of sequences that are present in ALL active sets. For three Raw files and active sets A and B, this would mean all sequences which occur in A and B (intersect), but not in C (setdiff).<br>
<b>intersection:</b> shows the number of sequences that occurs in all active sets (intersection).<br>
<b>union:</b> shows the number of sequences that occurs in total. For two files that are all sequences that occurs either in A or in B (union).<br>
<p>
Heatmap score [EVD: UpSet]: The proportion of sequences that the file has in common with all other files.
',
    workerFcn = function(.self, df_evd)
    {
      if (!checkInput(c("modified.sequence", "fc.raw.file"), df_evd)) return()
      
      getOutputWithMod = function(dl, mode){
        unlist(sapply(1:length(dl), function(numElem){
          comb = combn(names(dl),numElem)
          sapply(1:ncol(comb), function(x){
              sets = comb[,x]
              exp = as.expression(paste(sets, collapse = "&"))
              value = length(Reduce(mode, dl[sets]))
              names(value) = exp
              return(value)
          })
        }))
      }
      
      lf = tapply(df_evd$modified.sequence, df_evd$fc.raw.file, function(x){return(list(unique(x)))})
      # get rid of rawfiles without any PepIDs
      lf = Filter(function(l) length(l)>0 && any(!is.na(l)), lf)
      if (length(lf) <= 1)
      {
        lpl = list(ggText("UpSetR", "Only single Raw file detected. Cannot compute unions/intersections."))
        return(list(plots = lpl, title = list("EVD: UpSet")))
      }
      
      
      lpl = list(UpSetR::upset(UpSetR::fromList(lf), nsets = min(20, length(lf)), keep.order = TRUE, mainbar.y.label = "distinct size"))
      if (length(lf) < 6)
      { ## performance for enumerating all supersets forbids doing it on larger sets until we make this code smarter...
        lpl[[2]] = UpSetR::upset(UpSetR::fromExpression(getOutputWithMod(lf, intersect)), mainbar.y.label = "intersection size")
        lpl[[3]] = UpSetR::upset(UpSetR::fromExpression(getOutputWithMod(lf, union)), mainbar.y.label = "union size")
      }
      titles = list("EVD: UpSet distinct", 
                    "EVD: UpSet intersect",
                    "EVD: UpSet union")[1:length(lpl)]
      
      score = sapply(1:length(names(lf)), function(x){
        union = unique(unlist(lf[-x]))
        inters = intersect(lf[[x]], union)
        score = length(inters)/length(union)
        return(score)
      })
      
      qcScore = data.frame(fc.raw.file = names(lf), score = score)
      colnames(qcScore)[2] = .self$qcName
      
      return(list(plots = lpl, title = titles, qcScores = qcScore))
    }, 
    qcCat = "LC",
    qcName = "EVD:~UpSet", 
    orderNr = 0500  # just before peptide count
  )
    return(.self)
  })
)
