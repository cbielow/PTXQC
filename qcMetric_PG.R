
qcMetric_PG_Cont =  setRefClass(
  "qcMetric_PG_Cont",
  contains = "qcMetric",
  methods = list(initialize=function() {  callSuper(
    helpText = 
      "External protein contamination should be controlled for, therefore MaxQuant ships with a 
comprehensive, yet customizable protein contamination database, which is searched by MaxQuant by default. PTXQC 
generates a contamination plot derived from the proteinGroups (PG) table showing the fraction of total 
protein intensity attributable to contaminants. The plot employs transparency to discern differences in
the group-wise summed protein abundance. This allows to delineate a high contamination in high complexity samples from a high 
contamination in low complexity samples (e.g. from in-gel digestion). If you see only one 
abundance class ('mid'), this means all your groups have roughly
the same summed protein intensity.
Note that this plot is based on experimental groups, and therefore may not correspond 1:1 to Raw files. 
    
Heatmap score: none (since data source proteinGroups.txt is not related 1:1 to Raw files)
",
    workerFcn=function(.self, df_pg, int_cols, MAP_pg_groups)
    {
      ## completeness check
      if(!(.self$checkInput(c(int_cols, "contaminant"),colnames(df_pg)))){return(NULL)}
      
      df.con_stats = adply(int_cols, .margins=1, function(group) {
        #cat(group)
        total_int = sum(as.numeric(df_pg[, group]), na.rm = TRUE)
        return(data.frame(group_long = as.character(group),
                          log10_int  = log10(total_int),
                          cont_pc    = sum(as.numeric(df_pg[df_pg$contaminant, group]), na.rm = TRUE) /
                            total_int * 100
        ))
      })
      df.con_stats$group = MAP_pg_groups$short[match(df.con_stats$group_long, MAP_pg_groups$long)]
      df.con_stats$logAbdClass = getAbundanceClass(df.con_stats$log10_int)
      df.con_stats
      
      pg_plots_cont = byXflex(df.con_stats, 1:nrow(df.con_stats), 90, plot_ContsPG, sort_indices = FALSE)
      
      return(list(plots = pg_plots_cont))
    }, 
    qcCat = "Prep", 
    qcName = "PG:~Contaminants", 
    orderNr = 0110  ## should not show up in heatmap
  )
    return(.self)
  })
)

#####################################################################

qcMetric_PG_RawInt =  setRefClass(
  "qcMetric_PG_RawInt",
  contains = "qcMetric",
  methods = list(initialize=function() {  callSuper(
    helpText = 
      "Intensity boxplots by experimental groups. Groups are user-defined during MaxQuant configuration.
This plot displays a (customizable) threshold line for the desired mean intensity of proteins. Groups
which underperform here, are likely to also suffer from a worse MS/MS id rate and higher contamination due to
the lack of total protein loaded/detected. If possible, all groups should show a high and consistent amount
of total protein.
The height of the bar correlates to the number of proteins with non-zero abundance.

Contaminants are shown as overlayed yellow boxes, whose height corresponds to the number of contaminant proteins.
The position of the box gives the intensity distribution of the contaminants.

Heatmap score: none (since data source proteinGroups.txt is not related 1:1 to Raw files)
",
    workerFcn=function(.self, df_pg, int_cols, MAP_pg_groups, thresh_intensity)
    {
      ## completeness check
      if(!(.self$checkInput(c(int_cols, "contaminant"), colnames(df_pg)))){return(NULL)}
      
      
      ## some stats (for plot title)
      medians_no0 = sort(apply(log2(df_pg[, int_cols, drop = FALSE]), 2, function(x) quantile(x[x>0], probs=0.5, na.rm = TRUE))) # + c(0,0,0,0,0,0))
      int_dev_no0 = RSD(medians_no0)
      ## do not remove zeros (but add +1 since RSD is 'NA' when 'inf' is included in log-data)
      medians = sort(apply(log2(df_pg[, int_cols, drop = FALSE]+1), 2, quantile, na.rm = TRUE, probs=0.5)) # + c(0,0,0,0,0,0))
      int_dev = RSD(medians)
      int_dev.s = pastet("INT RSD [%]", round(int_dev, 3))
      lpl = boxplotCompare(data = melt(df_pg[, c(int_cols, "contaminant"), drop = FALSE], id.vars=c("contaminant"))[,c(2,3,1)],
                           log2 = TRUE, 
                           mainlab = "PG: intensity distribution",
                           ylab = expression(log[2]*" intensity"),
                           sublab = paste0("RSD ", round(int_dev_no0, 1),"% (w/o zero int.; expected < 5%)\n",
                                           "RSD ", round(int_dev, 1),"% [high RSD --> few peptides])"),
                           abline = thresh_intensity,
                           names = MAP_pg_groups)
      
      return(list(plots = lpl))
    }, 
    qcCat = "prep", 
    qcName = "PG:~raw~intensity", 
    orderNr = 0032
  )
    return(.self)
  })
)

#####################################################################

qcMetric_PG_LFQInt =  setRefClass(
  "qcMetric_PG_LFQInt",
  contains = "qcMetric",
  methods = list(initialize=function() {  callSuper(
    helpText = 
      "Label-free quantification (LFQ) intensity boxplots by experimental groups. Groups are user-defined during MaxQuant configuration.
This plot displays a (customizable) threshold line for the desired mean of LFQ intensity of proteins. Raw files
which underperform in *Raw* intensity, are likely to show an *increased* mean here, since
only high-abundance proteins are recovered and quantifyable by MaxQuant in this Raw file. The remaining proteins
are likely to receive an LFQ value of 0 (i.e. do not contribute to the distribution).
The height of the bar correlates to the number of proteins with non-zero abundance.

Contaminants are shown as overlayed yellow boxes, whose height corresponds to the number of contaminant proteins.
The position of the box gives the intensity distribution of the contaminants.

Heatmap score: none (since data source proteinGroups.txt is not related 1:1 to Raw files)
",
    workerFcn=function(.self, df_pg, int_cols, MAP_pg_groups, thresh_intensity)
    {
      ## completeness check
      if(!(.self$checkInput(c(int_cols, "contaminant"), colnames(df_pg)))){return(NULL)}
      
      
      ## some stats (for plot title)
      medians_no0 = sort(apply(log2(df_pg[, int_cols, drop = FALSE]), 2, function(x) quantile(x[x>0], probs=0.5, na.rm = TRUE))) # + c(0,0,0,0,0,0))
      int_dev_no0 = RSD(medians_no0)
      ## do not remove zeros (but add +1 since RSD is 'NA' when 'inf' is included in log-data)
      medians = sort(apply(log2(df_pg[, int_cols, drop = FALSE]+1), 2, quantile, na.rm = TRUE, probs=0.5)) # + c(0,0,0,0,0,0))
      int_dev = RSD(medians)
      lpl = boxplotCompare(data = melt(df_pg[, c(int_cols, "contaminant"), drop = FALSE], id.vars=c("contaminant"))[,c(2,3,1)],
                           log2 = TRUE, 
                           mainlab = "PG: LFQ intensity distribution",
                           ylab = expression(log[2]*" intensity"),
                           sublab = paste0("RSD ", round(int_dev_no0, 1),"% (w/o zero int.; expected < 5%)\n",
                                           "RSD ", round(int_dev, 1),"% [high RSD --> few peptides])"),
                           abline = thresh_intensity,
                           names = MAP_pg_groups)
      
      return(list(plots = lpl))
    }, 
    qcCat = "prep",
    qcName = "PG:~LFQ~intensity", 
    orderNr = 0033  ## should not show up in heatmap
  )
    return(.self)
  })
)

#####################################################################

qcMetric_PG_ReporterInt =  setRefClass(
  "qcMetric_PG_ReporterInt",
  contains = "qcMetric",
  methods = list(initialize=function() {  callSuper(
    helpText = 
      "ITRAQ/TMT reporter intensity boxplots by experimental groups. Groups are user-defined during MaxQuant configuration.
This plot displays a (customizable) threshold line for the desired mean of reporter ion intensity of proteins.
The height of the bar correlates to the number of proteins with non-zero abundance.

Contaminants are shown as overlayed yellow boxes, whose height corresponds to the number of contaminant proteins.
The position of the box gives the intensity distribution of the contaminants.
Contaminants should be lower compared to label-free samples, since all contaminants introduced after the labeling 
should not be identified by Andromeda (since they lack the isobaric tag).
    
There is a similar 'Raw file' based metric/plot based on evidence.txt.

Heatmap score: none (since data source proteinGroups.txt is not related 1:1 to Raw files)
",
    workerFcn=function(.self, df_pg, int_cols, MAP_pg_groups, thresh_intensity)
    {
      ## completeness check
      if(!(.self$checkInput(c(int_cols, "contaminant"), colnames(df_pg)))){return(NULL)}
      
      
      ## some stats (for plot title)
      medians_no0 = sort(apply(log2(df_pg[, int_cols, drop = FALSE]), 2, function(x) quantile(x[x>0], probs=0.5, na.rm = TRUE))) # + c(0,0,0,0,0,0))
      reprt_dev_no0 = RSD(medians_no0)
      ## do not remove zeros (but add +1 since RSD is 'NA' when 'inf' is included in log-data)
      medians = sort(apply(log2(df_pg[, int_cols, drop = FALSE]+1), 2, quantile, probs=0.5, na.rm = TRUE)) # + c(0,0,0,0,0,0))
      reprt_dev = RSD(medians)
      lpl = boxplotCompare(   data = melt(df_pg[, c(int_cols, "contaminant"), drop = FALSE], id.vars=c("contaminant"))[,c(2,3,1)],
                              log2 = TRUE, 
                              ylab = expression(log[2]*" reporter intensity"),
                              mainlab = "PG: reporter intensity distribution",
                              sublab = paste0("RSD ", round(reprt_dev_no0, 1),"% (w/o zero int.; expected < 5%)\n",
                                              "RSD ", round(reprt_dev, 1),"% [high RSD --> few peptides])"),
                              abline = thresh_intensity,
                              names = MAP_pg_groups)
      
      return(list(plots = lpl))
    }, 
    qcCat = "prep", 
    qcName = "PG:~Reporter~intensity", 
    orderNr = 0034  ## should not show up in heatmap
  )
    return(.self)
  })
)

#####################################################################

qcMetric_PG_PCA =  setRefClass(
  "qcMetric_PG_PCA",
  contains = "qcMetric",
  methods = list(initialize=function() {  callSuper(
    helpText = 
      "Principal components plots of experimental groups (as defined during MaxQuant configuration).

This plot is shown only if more than one experimental group was defined.
If LFQ was activated in MaxQuant, an additional PCA plot for LFQ intensities is shown. Similarly, if iTRAQ/TMT
reporter intensities are detected.

Since experimental groups and Raw files do not necessarily correspond 1:1, this plot cannot use the abbreviated
Raw file names, but instead must rely on automatic shortening of group names.

Heatmap score: none (since data source proteinGroups.txt is not related 1:1 to Raw files)
",
    workerFcn=function(.self, df_pg, lst_cols, MAP_pg_groups)
    {
      ## completeness check
      if(!(.self$checkInput(c(unlist(lst_cols), "contaminant"), colnames(df_pg)))){return(NULL)}
      
      lpl = list()
      for (cond in names(lst_cols))
      {
        #cond = names(lst_cols)[3]
        #print(lst_cols[cond])
        if (length(lst_cols[[cond]]) <= 1) 
        { ## only one condition.. PCA does not make sense (and will not work)
          next;
        }
        ## remove contaminants
        data = t(df_pg[!df_pg$contaminant, unlist(lst_cols[cond]), drop = FALSE])
        ## remove constant/zero columns (== dimensions == proteins)
        data = data[, colSums(data, na.rm = TRUE) > 0, drop = FALSE]
        rownames(data) = MAP_pg_groups$short[match(rownames(data), MAP_pg_groups$long)]
        pl = try(getPCA(data = data,
                        gg_layer = addGGtitle(paste0("PG: PCA of '", sub(".", " ", cond, fixed = TRUE), "'"), "(excludes contaminants)")
                       )[["plots"]])
        #print(pl)
        if (!inherits(pl, "try-error")) lpl = append(lpl, pl);
      }
      
      return(list(plots = lpl))
    }, 
    qcCat = "General", 
    qcName = "PG:~Principal~Component", 
    orderNr = 0003
  )
    return(.self)
  })
)

#####################################################################

qcMetric_PG_Ratio =  setRefClass(
  "qcMetric_PG_Ratio",
  contains = "qcMetric",
  methods = list(initialize=function() {  callSuper(
    helpText = 
      "This plot shows log<sub>2</sub> ratios for SILAC-like experiments (whenever MaxQuant reports a set of 'ratio.*' columns).
Useful to spot unequal channel mixing during sample preparation. If equal mixing is expected, the 
distribution should be unimodal and its mode close to 1 (i.e., a 1:1 ratio), as indicated by a visual 
guidance line. Multimodal distributions are flagged as such automatically. If PTXQC detects ratios 
deviating strongly from 1:1 (parameterized by default beyond the range between 1:4 and 4:1), PTXQC 
automatically assumes a pulsed experiment and reports the label incorporation in percent for all groups. 

Heatmap score: none (since data source proteinGroups.txt is not related 1:1 to Raw files)
",
    workerFcn=function(.self, df_pg, ratio_cols, thresh_LabelIncorp, GL_name_min_length)
    {
      ## completeness check
      if(!(.self$checkInput(c(ratio_cols, "contaminant", "reverse"), colnames(df_pg)))){return(NULL)}
      
      
      ## remove reverse and contaminants (might skew the picture)
      idx_row = !df_pg$contaminant & !df_pg$reverse
      d_sub = log2(df_pg[idx_row, ratio_cols, drop = FALSE])
      ## rename "ratio.h.l" to "h.l" (same for m.l in tripleSILAC)
      idx_globalRatio = grep("ratio\\.[hm]\\.l$", colnames(d_sub))
      if (length(idx_globalRatio)) colnames(d_sub)[idx_globalRatio] = gsub("^ratio\\.", "", colnames(d_sub)[idx_globalRatio])
      ## simplify the rest
      if (ncol(d_sub) > length(idx_globalRatio))
      {
        idx_other = setdiff(1:ncol(d_sub), idx_globalRatio)
        colnames(d_sub)[idx_other] = shortenStrings(simplifyNames(delLCP(colnames(d_sub)[idx_other], 
                                                                         min_out_length = GL_name_min_length, 
                                                                         add_dots = TRUE), 
                                                                  min_out_length = GL_name_min_length))
      }
      #summary(d_sub)
      # 
      # plot(density(d_sub[,1], bw = "SJ", adjust=1, na.rm = TRUE, n=128))
      # plot(d_sub[,1])
      # h = density(d_sub[,1], bw = "SJ", adjust=1, na.rm = TRUE, n=128)
      
      
      ## get ranges to fix breaks for density intervals
      # breaks = seq(min(d_sub, na.rm = TRUE), max(d_sub, na.rm = TRUE), length.out=(max(dd, na.rm = TRUE)-min(dd, na.rm = TRUE))/0.5)
      # mid = hist(d_sub[, 1], breaks = breaks)$mids
      ratio.densities = do.call(rbind, (lapply(1:ncol(d_sub), function(x) {
        name = colnames(d_sub)[x]
        ## density estimation can fail if not enough data
        h = try(density(na.omit(d_sub[ ,x]), bw = "SJ", adjust=2, na.rm = TRUE), silent = TRUE)
        if (inherits(h, "try-error")) return (NULL)
        count = sum(getMaxima(h$y))
        if (count > 1) name = paste(name, "*")
        df = data.frame(x = h$x, y = h$y, col = name, multimodal = (count>1))
        return (df)
      })))
      if (is.null(ratio.densities)) {
        pl_cont = ggText("PG: ratio density",
                         paste0("No data for plotting!"),
                         "red")
        return(list(plots = list(pl_cont)))
      } 
      
      ratio.densities$alpha = c(0.8, 1)[ratio.densities$multimodal+1]
      ratio.densities$ltype = c("dotted", "solid")[ratio.densities$multimodal+1]
      #head((ratio.densities))
      
      
      ## compute label incorporation?
      ratio.mode = ddply(ratio.densities, "col", .fun = function(x) {
        mode = x$x[which.max(x$y)]
        return (data.frame(mode = mode))
      })
      
      ## on more than ratio 1:4 or 4:1 ratio, report label incorporation
      # thresh_LabelIncorp == 4 by default
      if (max(abs(ratio.mode$mode)) > abs(log2(thresh_LabelIncorp))) {
        cat(paste0("Maximum ratio (log2) was ", round(max(abs(ratio.mode$mode)),1) , ", reaching the threshold of ", abs(log2(thresh_LabelIncorp)), " for label-incorporation computation.\nComputing ratios ...\n"))
        ## back to normal scale
        ratio.norm = 2^ratio.mode$mode
        ## compute incorporation
        ratio.inc =  ratio.norm / (ratio.norm+1) * 100
        ## round
        ratio.mode$li = round(ratio.inc)
        ## new label
        ratio.mode$col_new = ratio.mode$col %+% " (" %+% ratio.mode$li %+% "%)"
        ## replace column names in data.frame
        ratio.densities$col = ratio.mode$col_new[match(ratio.densities$col, ratio.mode$col)]
        ## notify user via legend title
        legend_title = "group\n(with label inc)"
      } else
      {
        cat(paste0("Maximum ratio (log2) was ", max(abs(ratio.mode$mode)), ", NOT reaching the threshold of ", abs(log2(thresh_LabelIncorp)), " for label-incorporation computation.\nSkipping ratios ...\n"))
        legend_title = "group"
      }
      
      main_title = "PG: ratio density\n(w/o contaminants)"
      main_col = "black"
      if (any(ratio.densities$multimodal))
      {
        main_title = paste0(main_title, "\nWarning: multimodal densities detected")
        main_col = "red"
      }
      ##
      ## plot ratios
      ##
      lpl = 
        byXflex(ratio.densities, ratio.densities$col, 5, plot_RatiosPG, sort_indices = FALSE, d_range = range(d_sub, na.rm=TRUE), main_title, main_col, legend_title)
      
      return(list(plots = lpl))
    },
    qcCat = "prep", 
    qcName = "PG:~Ratio", 
    orderNr = 0019
  )
    return(.self)
  })
)

#####################################################################

