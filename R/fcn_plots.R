##
## All plotting functions for computeReport() [main function]
##

#'
#' Plot contaminants from proteinGroups.txt
#' 
#' @param data A data.frame with columns 'group', 'cont_pc', 'logAbdClass'
#' @return GGplot object
#' 
#' @import ggplot2
#'
#' @examples 
#' 
#'  data = data.frame( 'group' = letters[1:10], 'cont_pc' = 2:11, 'logAbdClass' = c("low","high"))
#'  plot_ContsPG(data)
#' 
plot_ContsPG = function(data)
{
  data$section = as.integer(seq(0, nrow(data)/correctSetSize(nrow(data),30)*0.999, length.out=nrow(data)))
  p = ggplot(data=data, aes_string(x = "group", y = "cont_pc", alpha="logAbdClass")) +
        scale_alpha_discrete(range = c(c(0.3, 1)[(length(unique(data$logAbdClass))==1) + 1], 1.0), ## ordering of range is critical!
                             name = "Abundance\nclass") +
        geom_bar(stat="identity") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        xlab("")  +
        ggtitle("PG: Contaminant per condition") +
        ylab("contaminant (% intensity)") +
        geom_hline(aes_string(yintercept = "5"), linetype = 'dashed')
  if (length(unique(data$section))>1) p = p + facet_wrap(~ section, ncol = 1, scales="free_x")
  #print(p)
  return (p)
}

#'
#' Plot ratios of labeled data (e.g. SILAC) from proteinGroups.txt
#'
#' The 'x' values are expected to be log2() transformed already.
#' 
#' @param df_ratios A data.frame with columns 'x', 'y', 'col', 'ltype'
#' @param d_range X-axis range of plot
#' @param main_title Plot title
#' @param main_col Color of title
#' @param legend_title Legend text
#' @return GGplot object
#' 
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal

#' @examples 
#' 
#'  x1 = seq(-3, 3, by = 0.1)
#'  y1 = dnorm(x1)
#'  x2 = seq(-5, 1, by = 0.1)
#'  y2 = dnorm(x2, mean = -1)
#'  data = data.frame( x = c(x1,x2),
#'                     y = c(y1,y2), 
#'                     col = c(rep("ok", length(x1)), rep("shifted", length(x2))), 
#'                     ltype = "dotted")
#'  plot_RatiosPG(data, range(data$x), "Ratio plot", "red", "group")
#' 
plot_RatiosPG = function(df_ratios, d_range, main_title, main_col, legend_title)
{
  br = c(2, 5, 10, 20);
  
  p =
    ggplot(data = df_ratios, aes_string(x = "x", y = "y", colour = "col")) + 
      facet_grid(col ~ ., scales = "free_y") +
      geom_line(size = 1.2) +
      geom_area(aes_string(alpha = "ltype", fill = "col")) +
      xlab("ratio")  +
      ylab("density")  +
      scale_fill_manual(values = rep(brewer.pal(6,"Accent"), times=40), guide_legend(legend_title)) + 
      scale_colour_manual(values = rep(brewer.pal(6,"Accent"), times=40)) +
      scale_alpha_discrete(range = c(1, 0.2), 
                           labels=c("dotted"="unimodal", "solid"="multimodal"),
                           guide_legend("shape")
      ) +
      scale_x_continuous(limits = d_range, trans = "identity", breaks = c(-br, 0, br), labels=c(paste0("1/",2^(br)), 1, 2^br)) +
      guides(colour = FALSE) +
      theme(plot.title = element_text(colour = main_col)) +
      theme_bw() +
      geom_vline(alpha = 0.5, xintercept = 0, colour = "green", linetype = "dashed", size = 1.5) +
      ggtitle(main_title)
  #print(p)
  return (p)
}


#'
#' Plot user-defined contaminants from evidence.txt
#' 
#' @param data A data.frame with columns 'fc.raw.file', 'variable', 'value'
#' @param name_contaminant Name of the contaminant shown in title
#' @param extra_limit Horizontal position where a h-line is plotted (for visual guidance)
#' @import ggplot2
#' @return GGplot object
#' 
#' @examples 
#' 
#'  data = data.frame(fc.raw.file = letters[1:3], 
#'                    variable = c(rep("spectralCount", 3),
#'                                  rep("intensity", 3),
#'                                  rep("above.thresh", 3),
#'                                  rep("score_KS", 3)),
#'                    value = c(10, 20, 15, 9, 21, 14, 0, 1, 1, 0.3, 0.01, 0.04))
#'  plot_ContUser(data, "myco", 5)
#' 
plot_ContUser = function(data, name_contaminant, extra_limit) {
  datav = subset(data, data$variable %in% c('spectralCount', "intensity"))
  dataAT = subset(data, data$variable %in% c('above.thresh'))
  contRaws = dataAT$fc.raw.file[ dataAT$value == TRUE]
  dataKS = subset(data, data$variable == 'score_KS' & (data$fc.raw.file %in% contRaws))
  dataKS$value = paste0("p = ", round(dataKS$value,2))
  #cat(paste0("CA entry is ", extra_limit, "\n"))
  maxY = max(datav$value, extra_limit)
  datav$section = as.integer(seq(0, nrow(datav)/40, length.out = nrow(datav)))
  p = ggplot(datav, aes_string(x = "fc.raw.file", y = "value")) +
        geom_bar(stat="identity", aes_string(fill = "variable"), position = "dodge", width=.7) +
        ggtitle(paste0("EVD: Contaminant '", name_contaminant, "'")) +
        xlab("")  +
        ylab("abundance (%)") +
        ylim(c(0, maxY * 1.1)) +
        theme(plot.title = element_text(colour = "red"),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        scale_fill_discrete(name = "Method") +
        geom_hline(yintercept = extra_limit, linetype = 'dashed') +
        geom_text(data = dataKS, aes_string(label = "value", y = maxY * 1.05)) +
        facet_wrap(~ section, ncol = 1, scales = "free_x")
  #print(p)
  return(p)
}



#'
#' Plot Andromeda score distribution of contaminant peptide vs. matrix peptides.
#' 
#' The data is expected to be an ECDF already, x being the Andromeda score, y being the culmulative probability.
#' The Score is the probability of a Kolm.-Smirnoff test that the contaminant scores are larger (i.e.
#' large p-values indicate true contamination).
#' You will only see this plot if the %-threshold (YAML config) was reached. This is a saveguard against false-positive,
#' but high-scoring contaminant peptides, which would erroneously give you a large p-value and make you believe
#' your sample is contaminated although that's not the case.
#' 
#' @param data A data.frame with columns 'x', 'y', 'condition'
#' @param raw.file Name of Raw file for which the data is displayed (will become part of the plot title)
#' @param score Score of how distinct the distributions are (will become part of the title)
#' @return GGplot object
#' 
#' @import ggplot2
#'
#' @examples 
#' 
#'  data = data.frame(x = 10:60,
#'                    y = c(seq(0,1,length=51), seq(0.1, 1, length=51)), 
#'                    condition = rep(c("sample","contaminant"), each=51))
#'  plot_ContUserScore(data, 'test file', 0.96)
#' 
plot_ContUserScore = function(data, raw.file, score) {
  p = ggplot(data) + 
        geom_line(aes_string(x = "x", y = "y", col = "condition")) + 
        ggtitle(paste0("Empirical CDF of '", raw.file, "'\np = ", round(score, 2))) + 
        ylab("Pr") +
        xlab("Andromeda score")
  return(p)
}


#'
#' Plot Protein groups per Raw file
#' 
#' The input is a data.frame with protein/peptide counts, where 'category' designates
#' the origin of information (genuine ID, transferred ID, or both).
#' 
#' @param data A data.frame with columns 'fc.raw.file', 'counts', 'category'
#' @param y_max Plot limit of y-axis
#' @param thresh_line Position of a threshold line, indicating the usual target value
#' @param title Main title, and optional subtitle (if vector of length 2 is provided)
#' @return GGplot object
#'
#' @import ggplot2
#'
#' @examples 
#' 
#'  data = data.frame(fc.raw.file = rep(c("file A", "file B"), each=3),
#'                    counts = c(3674, 593, 1120, 2300, 400, 600), 
#'                    category = c("genuine","genuine+transferred","transferred"))
#'  plot_CountData(data, 6000, 4000, c("EVD: Protein Groups count", "gain: 23%"))
#' 
plot_CountData = function(data, y_max, thresh_line, title)
{
  title_main = title[1]
  title_sub = ifelse(length(title) > 1,  title[2], "")
  p = ggplot(data, aes_string(x = "fc.raw.file", y = "counts", fill = "category")) +
        geom_bar(stat = "identity", position = "stack") +
        xlab("") +
        ylab("count") +
        scale_x_discrete_reverse(data$fc.raw.file) +
        ylim(0, y_max) +
        scale_fill_manual(values = c("green", "#BEAED4", "red")) +
        addGGtitle(title_main, title_sub) + 
        geom_abline(alpha = 0.5, intercept = thresh_line, slope = 0, colour = "black", linetype = "dashed", size = 1.5) +
        coord_flip()
  return(p)
}


#'
#' Plot RT peak width over time
#' 
#' The input is a data.frame with already averaged counts over binned RT-slices.
#' 
#' @param data A data.frame with columns 'fc.raw.file', 'RT', 'peakWidth'
#' @param x_lim Plot range of x-axis
#' @param y_lim Plot range of y-axis
#' @return GGplot object
#' @importFrom RColorBrewer brewer.pal
#'
#' @import ggplot2
#'
#' @examples 
#' 
#'  data = data.frame(fc.raw.file = rep(c("file A", "file B", "file C"), each=81),
#'                    RT = c(20:100), 
#'                    peakWidth = c(rnorm(81, mean=20), rnorm(81, mean=10), rnorm(81, mean=30)))
#'  plot_RTPeakWidth(data, c(10, 100), c(0, 40))
#' 
plot_RTPeakWidth = function(data, x_lim, y_lim)
{
  p = ggplot(data) +
    geom_line(aes_string(x = "RT", y = "peakWidth", colour = "fc.raw.file"), size=1, alpha=0.7) +
    scale_color_manual(values = brewer.pal(length(unique(data$fc.raw.file)), "Set1")) +
    guides(color = guide_legend(title = "Raw file\n(avg. peak width)")) +
    xlab("retention time [min]") +
    ylab("peak width [min]") +
    coord_cartesian(xlim = x_lim, ylim = y_lim) + ## zoom in y -- do not cut data (preserve lines)
    ggtitle("EVD: Peak width over RT") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  #print(p)
  return(p)
}


#'
#' Plot MaxQuant Match-between-runs alignment performance.
#' 
#' The plots shows the correction function applied by MaxQuant, and the
#' residual RT (ideally 0) of each peptide to its reference. Uncalibrated peptides
#' are shown in red, calibrated ones in green.
#' The MaxQuant RT correction which was applied prior is shown in blue. The range of this function
#' can give hints if the allowed RT search window (20min by default) is sufficient or if
#' MaxQuant should be re-run with more tolerant settings.
#' 
#' The input is a data.frame with columns
#'   'calibrated.retention.time' - resulting (hopefully) calibratated RT after MQ-recal (the X-axis of the plot)
#'   'retention.time.calibration' - delta applied by MaxQuant
#'   'rtdiff' - remaining RT diff to reference peptide of the same sequence
#'   'RTdiff_in' - is the feature aligned (within 'match_tol')?
#'   'fc.raw.file_ext' - raw file
#' where each row represents one peptide whose RT was corrected by MaxQuant.
#' 
#' @param data A data.frame with columns as described above
#' @param y_lim Plot range of y-axis
#' @param title_sub Subtitle
#' @param match_tol Maximal residual RT delta to reference (usually ~1 min)
#' @return GGplot object
#'
#' @import ggplot2
#'
#' @examples 
#' 
#'  data = data.frame(fc.raw.file_ext = "file A", ## more than one would be possible
#'                    calibrated.retention.time = c(20:100), 
#'                    retention.time.calibration = 6 + sin((20:100)/10))
#'  data$rtdiff = rnorm(nrow(data))
#'  data$RTdiff_in = c("green", "red")[1 + (abs(data$rtdiff) > 0.7)]
#'  
#'  plot_MBRAlign(data, c(-10, 10), "fancy subtitle", 0.7)
#' 
plot_MBRAlign = function(data, y_lim, title_sub, match_tol)
{
  #data = evd_RT_t[ evd_RT_t$fc.raw.file == "file 13",]
  p = ggplot(data, aes_string(x = "calibrated.retention.time", y = "retention.time.calibration")) + 
        ## the MaxQuant correction (plot real data, no spline, since it can be very irregular)
        geom_line(aes(alpha = 0.7), color = "blue") +
        ## PTXQC correction
        geom_point(aes_string(x = "calibrated.retention.time", y = "rtdiff", color = "RTdiff_in"), alpha = 0.5) + 
        scale_alpha(name = 'Alignment function', 
                    labels = list(expression("MaxQuant" ~ Delta*"RT")),
                    range = c(1,1)) + 
        scale_colour_manual(name = expression(bold("ID pairs ("*Delta*"RT to Ref)")), 
                            values = c("green" = "green", "red" = "red"),
                            labels=c("green" = paste0("good (<", match_tol, "min)"), 
                                     "red" = paste0("bad (>", match_tol, "min)"))) +
        guides(colour = guide_legend(order = 2), 
               alpha = guide_legend(order = 1)) +   ## alpha-legend on top, color below
        ylim(y_lim) +
        xlab("corrected RT [min]") +
        ylab(expression(Delta*"RT [min]")) +
        facet_wrap(~ fc.raw.file_ext) +
        addGGtitle("EVD: MBR - alignment", title_sub)  
  #print(p)
  return(p)
}


#'
#' Plot MaxQuant Match-between-runs id transfer performance.
#' 
#' The plots shows the different categories of peak classes
#' 
#' The input is a data.frame with columns
#'   'fc.raw.file' - raw file name
#'   'single' - fraction of peptides with are represent only once
#'   'multi.inRT' - fraction of peptides with are represent multiple times, 
#'                  but within a certain RT peak width
#'   'multi.outRT' - fraction of peptides with are represent multiple times,
#'                   with large RT distance
#'   'sample' - raw file
#' where each row represents one peptide sequence.
#' 
#' @param data A data.frame with columns as described above
#' @return GGplot object
#'
#' @import ggplot2
#'
#' @examples 
#'  data = data.frame(fc.raw.file = rep(c("file A", "file B"), each = 3),
#'                    single = c(0.9853628, 0.8323160, 0.9438375, 0.9825538, 0.8003763, 0.9329961), 
#'                    multi.inRT = c(0.002927445, 0.055101018, 0.017593087, 0.005636457, 0.099640044, 0.031870056),
#'                    multi.outRT = c(0.01170978, 0.11258294, 0.03856946, 0.01180972, 0.09998363, 0.03513386),
#'                    sample = rep(c("genuine", "transferred", "all"), 2))
#'  plot_MBRIDtransfer(data)
#' 
plot_MBRIDtransfer = function(data)
{
  data.m = melt(data, id.vars=c("fc.raw.file", "sample"))
  data.m$value = data.m$value * 100 ## used to be scores in [0-1]
  p = ggplot(data.m) + 
        geom_bar(aes_string(x="fc.raw.file", y="value", fill="variable"), stat="identity", position="stack") + 
        scale_fill_manual("peak class", 
                          values = c("single"="green", "multi.inRT"="lightgreen", "multi.outRT"="red"),
                          labels=c("single", "group (in width)", "group (out width)")) +
        ylim(0, 100.1) + ## ggplot might not show the last (red) group upon 100.0
        xlab("") +
        ylab("fraction of 3D-peaks [%]") +
        coord_flip() + 
        scale_x_discrete_reverse(factor(data$fc.raw.file)) +
        ggtitle("EVD: MBR - ID Transfer") + 
        facet_wrap(~sample)
  #print(p)
  return(p)
}


#'
#' Plot MaxQuant Match-between-runs id transfer performance.
#' 
#' The plots shows the different categories of peak classes
#' 
#' The input is a data.frame with columns
#'   'fc.raw.file' - raw file name
#'   'single' - fraction of peptides with are represent only once
#'   'multi.inRT' - fraction of peptides with are represent multiple times, 
#'                  but within a certain RT peak width
#'   'multi.outRT' - fraction of peptides with are represent multiple times,
#'                   with large RT distance
#'   'sample' - raw file
#' where each row represents one peptide sequence.
#' 
#' @param data A data.frame with columns as described above
#' @return GGplot object
#'
#' @import ggplot2
#' @import directlabels
#'
#' @examples
#'  data = data.frame(fc.raw.file = paste("file", letters[1:4]),
#'                    abs = c(5461, 5312, 3618, 502), 
#'                    pc = c(34, 32, 22, 2))
#'  plot_MBRgain(data, "MBR gain: 18%")
#' 
plot_MBRgain = function(data, title_sub)
{
  p = ggplot(data = data, aes_string(x = "abs", y = "pc", col = "fc.raw.file")) + 
    geom_point(size=2) + 
    addGGtitle("EVD: Peptides inferred by MBR", title_sub) +
    xlab("number of transferred ID's") +
    ylab("gain on top of genuine IDs [%]") +
    xlim(0, max(mtr.df$abs, na.rm = TRUE)*1.1) + ## accounting for labels near the border
    ylim(0, max(mtr.df$pc, na.rm = TRUE)*1.1)
  p = direct.label(p, list(cex=0.5, "smart.grid"))
  #print(p)
  return(p)
}


#'
#' Plot MaxQuant Match-between-runs id transfer performance.
#' 
#' The plots shows the charge distribution per Raw file.
#' The output of 'mosaicize()' can be used directly.
#' 
#' The input is a data.frame with columns
#'   'Var1' - name of the Raw file
#'   'Var2' - charge (used as fill color)
#'   'Var1_center' - contains X-position of the Raw file
#'   'Var2_height' - relative frequency of the charge
#'   'Margin_var1' - 
#' where each row represents one peptide sequence.
#' 
#' @param data A data.frame with columns as described above
#' @return GGplot object
#'
#' @import ggplot2
#'
#' @examples
#'  data = data.frame(raw.file = c(rep('file A', 100), rep('file B', 40)),
#'                      charge = c(rep(2, 60), rep(3, 30), rep(4, 10),
#'                                 rep(2, 30), rep(3, 7), rep(4, 3)))
#'  plot_Charge(mosaicize(data))
#' 
plot_Charge = function(data)
{
  p = ggplot(data, aes_string(x = "Var1_center", y = "Var2_height")) +
        geom_bar(stat = "identity", aes_string(width = "Margin_var1", fill = "Var2"), color = "black")  +
        geom_text(aes_string(label = "as.character(Var1)", x = "Var1_center", y = 1.05)) +
        xlab("Raw file") +
        ylab("fraction [%]") +
        guides(fill = guide_legend(title="charge"), color = FALSE) + # avoid black line in legend
        scale_x_reverse() +
        coord_flip() +
        theme(axis.text.y = element_blank(), axis.ticks = element_blank()) +
        ggtitle("EVD: charge distribution")
  #print(p)
  return(p)
}


#'
#' Plot IDs over time for each Raw file.
#' 
#' The plots shows the charge distribution per Raw file.
#' The output of 'mosaicize()' can be used directly.
#' 
#' The input is a data.frame with columns
#'   'RT' - RT in seconds, representing one bin
#'   'counts' - number of IDs at this bin
#'   'fc.raw.file' - name of the Raw file
#' where each row represents one bin in RT.
#' 
#' @param data A data.frame with columns as described above
#' @return GGplot object
#'
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' 
#' @examples
#'  data = data.frame(fc.raw.file = rep(paste('file', letters[1:3]), each=30),
#'                             RT = seq(20, 120, length.out = 30),
#'                         counts = c(rnorm(30, 400, 20), rnorm(30, 250, 15), rnorm(30, 50, 15)))
#'  plot_IDsOverRT(data)
#' 
plot_IDsOverRT = function(data, x_lim = range(data$RT), y_max = max(data$counts))
{
  nrOfRaws = length(unique(data$fc.raw.file))
  p = ggplot(data = data) +
    geom_line(aes_string(x = "RT", y = "counts", colour = "fc.raw.file", linetype = "fc.raw.file")) +
    xlim(x_lim) +
    xlab("RT [min]") + 
    ylim(from = 0, to = y_max) +
    ylab("ID count") +
    ggtitle("EVD: IDs over RT") +
    guides(colour = guide_legend(title="Raw file"), linetype = FALSE) +
    scale_linetype_manual(values = rep_len(c("solid", "dashed"), nrOfRaws)) +
    scale_color_manual(values = brewer.pal(nrOfRaws, "Set1")) 
  #print(p)
  return(p)
}

#'
#' Plot percent of identified MS/MS for each Raw file.
#' 
#' Useful for a first overall impression of the data.
#' 
#' The input is a data.frame with columns
#'   'fc.raw.file' - name of the Raw file
#'   'ms.ms.identified....' - fraction of identified MS/MS spectra in percent
#'   'cat' - identification category as arbitrary string
#' where each row represents one Raw file.
#' 
#' @param data A data.frame with columns as described above
#' @param id_rate_bad Number below which the ID rate is considered bad
#' @param id_rate_great Number above which the ID rate is considered great
#' @param label_ID Named vector with colors for the categories given in data$cat
#' @return GGplot object
#'
#' @import ggplot2
#'
#' @examples
#'  id_rate_bad = 20; id_rate_great = 35;
#'  label_ID = c("bad (<20%)" = "red", "ok (...)" = "blue", "great (>35%)" = "green")
#'  data = data.frame(fc.raw.file = paste('file', letters[1:3]),
#'                    ms.ms.identified.... = rnorm(3, 25, 15))
#'  data$cat = factor(cut(data$ms.ms.identified...., breaks=c(-1, id_rate_bad, id_rate_great, 100), labels=names(label_ID)))                  
#'  plot_IDRate(data, id_rate_bad, id_rate_great, label_ID)
#'
plot_IDRate = function(data, id_rate_bad, id_rate_great, label_ID)
{
  p = ggplot(data, aes_string(y = "fc.raw.file", x = "ms.ms.identified....")) +
        geom_point(aes_string(colour = "cat")) +
        geom_vline(xintercept = id_rate_bad, color=(label_ID)[1]) +
        geom_vline(xintercept = id_rate_great, color=(label_ID)[3]) +
        ylab("") + 
        xlab("MS/MS identified [%]") +
        scale_colour_manual(values=label_ID) + 
        ggtitle("SM: MS/MS identified per Raw file") + 
        xlim(0, max(data$ms.ms.identified...., id_rate_great)*1.1) + 
        guides(color=guide_legend(title="ID class")) +
        scale_y_discrete_reverse(data$fc.raw.file, breaks = ggAxisLabels)
  #print(p)
  return(p)
}


#'
#' Plot a table of Raw files which failed ID rate requirements.
#' 
#' Just to have their names explicitly in the report.
#' 
#' The input is a data.frame with two columns
#'   1st column - name of the Raw file
#'   2nd column - fraction of identified MS/MS spectra in percent
#' where each row represents one Raw file.
#' 
#' @param data A data.frame with columns as described above
#' @param title Table title
#' @return gTree object with class 'PTXQC_table'
#'
#' @import grid
#' @import gridExtra
#' @import gtable
#'
#' @examples
#'   data = data.frame(raw.file = letters[1:4],
#'                     id.rate = 3:6)
#'   plot_IDRateBad(data, "Bad files")
#' 
plot_IDRateBad = function(data, title)
{
  table = tableGrob(data, rows = NULL, cols = c("Raw file", "% identified"))
  title = textGrob(title, gp = gpar(fontsize = 20))
  padding = unit(0.5, "line")
  table = gtable_add_rows(table, 
                          heights = grobHeight(title) + padding,
                          pos = 0)
  table = gtable_add_grob(table, list(title), t = 1, l = 1, r = ncol(table))
  
  ## neat trick to enable calling print(g), to mimic ggplot-behaviour on this object
  ## in combination with print.PTXQC_table() -- see below
  p = gTree(children = gList(table), cl = c("PTXQC_table"))
  print(p)
  return(p)
}

#' helper S3 class, enabling print(some-plot_Table-object)
#' @import grid
#' @param x Some Grid object to plot
print.PTXQC_table = function(x) {grid.newpage(); grid.draw(x)}





