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
#' @export
#' 
#' @examples 
#' 
#'  data = data.frame( 'group' = letters[1:10], 'cont_pc' = 2:11, 'logAbdClass' = c("low","high"))
#'  plot_ContsPG(data)
#' 
plot_ContsPG = function(data)
{
  data$section = as.integer(seq(0, nrow(data)/correctSetSize(nrow(data),30)*0.999, length.out=nrow(data)))
  p = ggplot(data=data, aes(x = .data$group, y = .data$cont_pc, alpha = .data$logAbdClass)) +
        suppressWarnings(## suppresses 'Using alpha for a discrete variable is not advised'
        scale_alpha_discrete(range = c(c(0.3, 1)[(length(unique(data$logAbdClass))==1) + 1], 1.0), ## ordering of range is critical!
                             name = "Abundance\nclass")) +
        geom_col() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        xlab("")  +
        ggtitle("PG: Contaminant per condition") +
        ylab("contaminant (% intensity)") +
        geom_hline(aes(yintercept = 5), linetype = 'dashed')
  if (length(unique(data$section))>1) p = p + facet_wrap(~ section, ncol = 1, scales="free_x")
  #print(p)
  return (p)
}



#'
#' Plot user-defined contaminants from evidence.txt
#' 
#' Kolmogorov-Smirnoff p-values are plotted on top of each group.
#' High p-values indicate that Andromeda scores for contaminant peptides
#' are equal or higher compared to sample peptide scores, i.e. the probability that
#' sample peptides scores are NOT greater than contaminant peptide scores.
#' 
#' @param data A data.frame with columns 'fc.raw.file', 'variable', 'value'
#' @param name_contaminant Name of the contaminant shown in title
#' @param extra_limit Position where a h-line is plotted (for visual guidance)
#' @param subtitle Optional subtitle for plot
#' @return GGplot object
#' 
#' @import ggplot2
#' @export
#' 
#' @examples 
#' 
#'  data = data.frame(fc.raw.file = letters[1:3], 
#'                    variable = c(rep("spectralCount", 3),
#'                                  rep("intensity", 3),
#'                                  rep("above.thresh", 3),
#'                                  rep("score_KS", 3)),
#'                    value = c(10, 20, 15, 9, 21, 14, 0, 1, 1, 0.3, 0.01, 0.04))
#'  plot_ContUser(data, "myco", 5, "subtitle")
#' 
plot_ContUser = function(data, name_contaminant, extra_limit, subtitle = NULL)
{
  datav = subset(data, data$variable %in% c('spectralCount', "intensity"))
  datav$section = assignBlocks(datav$fc.raw.file, set_size = 40, sort_values = TRUE)
  dataAT = subset(data, data$variable %in% c('above.thresh'))
  ## contRaws might be empty
  contRaws = dataAT$fc.raw.file[ dataAT$value == TRUE]
  dataKS = subset(data, data$variable == 'score_KS' & (data$fc.raw.file %in% contRaws))
  if (nrow(dataKS)>0) {
    dataKS$value = paste0("p = ", round(dataKS$value,2))
    ## use the same section, so ggplot knows how to subset the data
    dataKS$section = datav$section[match(dataKS$fc.raw.file, datav$fc.raw.file)]
  } 
  #cat(paste0("CA entry is ", extra_limit, "\n"))
  maxY = max(datav$value, extra_limit)
  p = ggplot(datav, aes(x = .data$fc.raw.file, y = .data$value)) +
        geom_col(aes(fill = .data$variable), position = "dodge", width=.7) +
        ggtitle(paste0("EVD: Contaminant '", name_contaminant, "'"), subtitle) +
        xlab("")  +
        ylab("abundance fraction (%)") +
        ylim(c(0, maxY * 1.1)) +
        theme(plot.title = element_text(colour = "red"),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        scale_fill_discrete(name = "Method") +
        geom_hline(yintercept = extra_limit, linetype = 'dashed')
  ## group(NULL) seems important in geom_text()
  if (nrow(dataKS)>0) p = p + geom_text(data = dataKS, aes(label = .data$value, y = maxY * 1.05, group = NULL))
  p = p + facet_wrap(~ section, ncol = 1, scales = "free_x")
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
#' @export
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
    geom_line(aes(x = .data$x, y = .data$y, col = .data$condition)) + 
    ggtitle(paste0("Empirical CDF of '", raw.file, "'\np = ", round(score, 2))) + 
    ylab("Pr") +
    xlab("Andromeda score")
  return(p)
}


#'
#' Plot contaminants from evidence.txt, broken down into top5-proteins.
#' 
#' @param data A data.frame with columns 'fc.raw.file', 'contaminant', 'pname', 'intensity'
#' @param top5 Name of the Top-5 Proteins (by relative intensity or whatever seems relevant)
#' @return GGplot object
#' 
#' @import ggplot2
#' @export
#' 
#' @examples 
#' 
#'  data = data.frame(intensity = 1:12, 
#'                    pname = rep(letters[1:3], 4), 
#'                    fc.raw.file = rep(paste("f", 1:4), each=3),
#'                    contaminant = TRUE)
#'  ## providing more proteins than present... d,e will be ignored
#'  plot_ContEVD(data, top5 = letters[1:5])
#'  ## classify 'c' as 'other'
#'  plot_ContEVD(data, top5 = letters[1:2])
#' 
plot_ContEVD = function(data, top5) 
{ 
  #top5 = cont.top5.names
  if (is.null(top5)) stop("Function plot_ContEVD() called with invalid argument. Please report this bug.")
  if (length(top5) > 5) stop("Top5 protein list is longer than 5 (which is the maximum allowed).")
  
  intensity = NULL ## to make R CHECK happy...
  data.sub = data[data$contaminant > 0,]
  ## rewrite prot names, and subsume 6th and below as 'other'
  data.sub$pname = as.character(data.sub$pname)
  data.sub[!(data.sub$pname %in% top5), "pname"] = 'other'
  ## aggregate identical proteins
  ##  use sum(as.numeric(.)) to prevent overflow
  d_sum = plyr::ddply(data.sub[, c("intensity", "pname", "fc.raw.file")], c("pname", "fc.raw.file"), 
                function(x) plyr::summarise(x, s.intensity=sum(as.numeric(intensity), na.rm = TRUE)))
  ## normalize by total intensity of raw file
  d_norm = plyr::ddply(data[, c("intensity", "fc.raw.file")],  "fc.raw.file", 
                 function(x) plyr::summarise(x, total.intensity=sum(as.numeric(intensity), na.rm = TRUE)))
  
  d_sum$total.intensity = d_norm$total.intensity[match(d_sum$fc.raw.file, d_norm$fc.raw.file)]
  d_sum$Log10Diff = getAbundanceClass(log10(d_sum$total.intensity))
  d_sum$s.intensity = d_sum$s.intensity / d_sum$total.intensity * 100
  ## shorten protein-groups (at most two protein names)
  d_sum$pname = sapply(d_sum$pname, function(x) {
    p.split = unlist(strsplit(x, split=";"))
    ## shorten entries as well (at most 15 characters)
    p.split_s = sapply(p.split[1:(min(2, length(p.split)))], function(x) ifelse(nchar(x)>15, paste0(substr(x, start=1, stop=13), ".."), x))
    r = paste(p.split_s, sep="", collapse=";")
    if (length(p.split)>2) r=paste0(r, ";..")
    return(r)
  })
  ## order of pname determines order of bars    
  d_sum = rbind(d_sum[d_sum$pname!="other",], d_sum[d_sum$pname=="other",])

  ## value of factors determines order in the legend
  ## --> make proteins a factor, with 'other' being the first
  d_sum$Protein = factor(d_sum$pname, levels = unique(c("other", d_sum$pname)), ordered = TRUE)
  head(d_sum)
  
  ## plot
  p = ggplot(d_sum, aes(  x = .data$fc.raw.file,
                          y = .data$s.intensity, 
                       fill = .data$Protein)) +
        geom_col(aes(alpha = .data$Log10Diff)) +
        suppressWarnings(## suppresses 'Using alpha for a discrete variable is not advised'
          scale_alpha_discrete(range = c(c(0.3, 1)[(length(unique(d_sum$Log10Diff))==1) + 1], 1.0),
                               name = "Abundance\nclass")) +
        xlab("")  +
        theme_bw() +
        ggtitle("EVD: Top5 Contaminants per Raw file") +
        ylab("contaminant (% intensity)") +
        scale_fill_manual(values = RColorBrewer::brewer.pal(6,"Accent")) + 
        scale_colour_manual(values = RColorBrewer::brewer.pal(6,"Accent")) +
        geom_hline(aes(yintercept = 5), linetype='dashed') +
        #guides(alpha=NULL, fill = guide_legend(nrow = 2, ncol = 3, byrow = TRUE, reverse = TRUE)) +
        #theme(legend.position="top", legend.title=element_blank()) +
        coord_flip() +
        scale_x_discrete_reverse(d_sum$fc.raw.file)
  
  #print(p)
  return(p)  
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
#' @export
#' 
#' @examples 
#' 
#'  x1 = seq(-3, 3, by = 0.1)
#'  y1 = dnorm(x1)
#'  x2 = seq(-5, 1, by = 0.1)
#'  y2 = dnorm(x2, mean = -1)
#'  data = data.frame( x = c(x1,x2),
#'                     y = c(y1,y2), 
#'                     col = c(rep("ok", length(x1)), rep("shifted", length(x2))), 
#'                     ltype = c(rep("solid", length(x1)), rep("dotted", length(x2))))
#'  plot_RatiosPG(data, range(data$x), "Ratio plot", "red", "group")
#' 
plot_RatiosPG = function(df_ratios, d_range, main_title, main_col, legend_title)
{
  br = c(2, 5, 10, 20);
  
  p =
    ggplot(data = df_ratios, aes(x = .data$x, y = .data$y, colour = .data$col)) + 
    facet_grid(col ~ ., scales = "free_y") +
    geom_line(linewidth = 1.2) +
    geom_area(aes(alpha = .data$ltype, fill = .data$col)) +
    xlab("ratio")  +
    ylab("density")  +
    scale_fill_manual(values = rep(RColorBrewer::brewer.pal(6,"Accent"), times=40), guide = guide_legend(legend_title)) + 
    scale_colour_manual(values = rep(RColorBrewer::brewer.pal(6,"Accent"), times=40)) +
    suppressWarnings(scale_alpha_discrete(range = c(1, 0.2), 
                         labels=c("dotted"="unimodal", "solid"="multimodal"),
                         guide = guide_legend("shape")
    )) +
    scale_x_continuous(limits = d_range, trans = "identity", breaks = c(-br, 0, br), labels=c(paste0("1/",2^(br)), 1, 2^br)) +
    guides(color = "none") +
    theme(plot.title = element_text(colour = main_col)) +
    theme_bw() +
    geom_vline(alpha = 0.5, xintercept = 0, colour = "green", linetype = "dashed", linewidth = 1.5) +
    ggtitle(main_title)
  #print(p)
  return (p)
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
#' @export
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
  p = ggplot(data, aes(x = .data$fc.raw.file, y = .data$counts, fill = .data$category)) +
        geom_col(position = position_stack(reverse = TRUE)) +
        xlab("") +
        ylab("count") +
        scale_x_discrete_reverse(data$fc.raw.file) +
        ylim(0, y_max) +
        scale_fill_manual(values = c("green", "#BEAED4", "blue")) +
        ggtitle(title_main, title_sub) + 
        geom_abline(alpha = 0.5, intercept = thresh_line, slope = 0, colour = "black", linetype = "dashed", linewidth = 1.5) +
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
#'
#' @import ggplot2
#' @export
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
    geom_line(aes(x = .data$RT, y = .data$peakWidth, colour = .data$fc.raw.file), linewidth = 1, alpha = 0.7) +
    scale_color_manual(values = brewer.pal.Safe(length(unique(data$fc.raw.file)), "Set1")) +
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
#' @export
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
  p = ggplot(data, aes(x = .data$calibrated.retention.time, y = .data$retention.time.calibration)) + 
        ## the MaxQuant correction (plot real data, no spline, since it can be very irregular)
        geom_line(aes(alpha = 0.7), color = "blue") +
        ## PTXQC correction
        geom_point(aes(x = .data$calibrated.retention.time, y = .data$rtdiff, color = .data$RTdiff_in), alpha = 0.5) + 
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
        ggtitle("EVD: MBR - alignment", title_sub)  
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
#' @export
#' 
#' @examples 
#'  data = data.frame(fc.raw.file = rep(c("file A", "file B"), each = 3),
#'                    single = c(0.9853628, 0.8323160, 0.9438375,
#'                               0.9825538, 0.8003763, 0.9329961), 
#'                    multi.inRT = c(0.002927445, 0.055101018, 0.017593087,
#'                                   0.005636457, 0.099640044, 0.031870056),
#'                    multi.outRT = c(0.01170978, 0.11258294, 0.03856946,
#'                                    0.01180972, 0.09998363, 0.03513386),
#'                    sample = rep(c("genuine", "transferred", "all"), 2))
#'  plot_MBRIDtransfer(data)
#' 
plot_MBRIDtransfer = function(data)
{
  data.m = reshape2::melt(data, id.vars=c("fc.raw.file", "sample"))
  data.m$value = data.m$value * 100 ## used to be scores in [0-1]
  if (all(is.na(data.m$value)))
  {# the slice of Raw file we are looking at could have no MBR data -- and ggplot needs something to plot...
    data.m$value = 0
  }
  p = ggplot(data.m) + 
        geom_col(aes(x = .data$fc.raw.file, y = .data$value, fill = .data$variable), position = position_stack(reverse = TRUE)) + 
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
#' Plot MaxQuant Match-between-runs id transfer performance as a scatterplot.
#' 
#' Per Raw file, the absolute number of transferred IDs as well as the relative gain in percent.
#' 
#' The input is a data.frame with columns
#'   'fc.raw.file' - raw file name
#'   'abs' - absolute number of transferred ID's
#'   'pc' - gain on top of genuine IDs [%]
#' where each row represents one rawfile.
#' 
#' @param data A data.frame with columns as described above
#' @param title_sub Subtitle text
#' @return GGplot object
#'
#' @import ggplot2
#' @export
#' 
#' @examples
#'  data = data.frame(fc.raw.file = paste("file", letters[1:4]),
#'                    abs = c(5461, 5312, 3618, 502), 
#'                    pc = c(34, 32, 22, 2))
#'  plot_MBRgain(data, "MBR gain: 18%")
#' 
plot_MBRgain = function(data, title_sub = "")
{
  p = ggplot(data = data, aes(x = .data$abs, y = .data$pc, col = .data$fc.raw.file)) + 
        geom_point(size = 2) + 
        ggtitle("EVD: Peptides inferred by MBR", title_sub) +
        xlab("number of transferred ID's") +
        ylab("gain on top of genuine IDs [%]") +
        xlim(0, max(data$abs, na.rm = TRUE)*1.1) + ## accounting for labels near the border
        ylim(0, max(data$pc, na.rm = TRUE)*1.1) +
        guides(color = "none") +
        geom_text(aes(label = .data$fc.raw.file), hjust = -0.2, show.legend = FALSE, check_overlap = TRUE)
  #print(p)
  return(p)
}


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
#' @param d_charge A data.frame with columns as described above
#' @return GGplot object
#'
#' @import ggplot2
#' @export
#' 
#' @examples
#'  data = data.frame(raw.file = c(rep('file A', 100), rep('file B', 40)),
#'                        data = c(rep(2, 60), rep(3, 30), rep(4, 10),
#'                                 rep(2, 30), rep(3, 7), rep(4, 3)))
#'  plot_Charge(mosaicize(data))
#' 
plot_Charge = function(d_charge)
{
  p = ggplot(d_charge, aes(x = .data$Var1_center, y = .data$Var2_height, width = .data$Margin_var1)) +
        geom_col(aes(fill = .data$Var2), color = "black", position = position_stack(reverse = TRUE))  +
        geom_text(aes(label = .data$Var1, x = .data$Var1_center, y = 1.05)) +
        xlab("Raw file") +
        ylab("fraction [%]") +
        guides(fill = guide_legend(title = "charge"),
                                   color = "none") + # avoid black line in legend
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
#' Uses plot_DataOverRT() internally.
#' 
#' @param data A data.frame with columns as described above
#' @param x_lim Limits of the x-axis (2-tuple)
#' @param y_max Maximum of the y-axis (single value)
#' @return GGplot object
#'
#' @import ggplot2
#' @export
#' 
#' @examples
#'  data = data.frame(fc.raw.file = rep(paste('file', letters[1:3]), each=30),
#'                             RT = seq(20, 120, length.out = 30),
#'                         counts = c(rnorm(30, 400, 20), rnorm(30, 250, 15), rnorm(30, 50, 15)))
#'  plot_IDsOverRT(data)
#' 
plot_IDsOverRT = function(data, x_lim = range(data$RT), y_max = max(data$counts))
{
  return(plot_DataOverRT(data, "EVD: IDs over RT", "ID count", x_lim, y_max))
}


#'
#' Plot some count data over time for each Raw file.
#' 
#' The input is a data.frame with columns
#'   'RT' - RT in seconds, representing one bin
#'   'counts' - number of counts at this bin
#'   'fc.raw.file' - name of the Raw file
#' where each row represents one bin in RT.
#' 
#' At most nine(!) Raw files can be plotted. If more are given,
#' an error is thrown.
#' 
#' 
#' @param data A data.frame with columns as described above
#' @param title The plot title
#' @param y_lab Label of y-axis
#' @param x_lim Limits of the x-axis (2-tuple)
#' @param y_max Maximum of the y-axis (single value)
#' @return GGplot object
#'
#' @import ggplot2
#' @export
#' 
#' @examples
#'  data = data.frame(fc.raw.file = rep(paste('file', letters[1:3]), each=30),
#'                             RT = seq(20, 120, length.out = 30),
#'                         counts = c(rnorm(30, 400, 20), rnorm(30, 250, 15), rnorm(30, 50, 15)))
#'  plot_DataOverRT(data, "some title", "count data")
#' 
plot_DataOverRT = function(data, title, y_lab, x_lim = range(data$RT), y_max = max(data$counts))
{
  nrOfRaws = length(unique(data$fc.raw.file))
  p = ggplot(data = data) +
    geom_line(aes(x = .data$RT, y = .data$counts, colour = .data$fc.raw.file, linetype = .data$fc.raw.file)) +
    xlim(x_lim) +
    xlab("RT [min]") + 
    ylim(from = 0, to = y_max) +
    ylab(y_lab) +
    ggtitle(title) +
    guides(colour = guide_legend(title="Raw file"), linetype = "none") +
    scale_linetype_manual(values = rep_len(c("solid", "dashed"), nrOfRaws)) +
    scale_color_manual(values = brewer.pal.Safe(nrOfRaws, "Set1")) 
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
#' @export
#' 
#' @examples
#'  id_rate_bad = 20; id_rate_great = 35;
#'  label_ID = c("bad (<20%)" = "red", "ok (...)" = "blue", "great (>35%)" = "green")
#'  data = data.frame(fc.raw.file = paste('file', letters[1:3]),
#'                    ms.ms.identified.... = rnorm(3, 25, 15))
#'  data$cat = factor(cut(data$ms.ms.identified....,
#'                        breaks=c(-1, id_rate_bad, id_rate_great, 100),
#'                        labels=names(label_ID)))                  
#'  plot_IDRate(data, id_rate_bad, id_rate_great, label_ID)
#'
plot_IDRate = function(data, id_rate_bad, id_rate_great, label_ID)
{
    p = ggplot(data, aes(y = .data$fc.raw.file, x = .data$ms.ms.identified....)) +
        geom_point(aes(colour = .data$cat)) +
        geom_vline(xintercept = id_rate_bad, color=(label_ID)[1]) +
        geom_vline(xintercept = id_rate_great, color=(label_ID)[3]) +
        ylab("") + 
        xlab("MS/MS identified [%]") +
        scale_colour_manual(values=label_ID) + 
        ggtitle("SM: MS/MS identified per Raw file") + 
        xlim(0, max(data$ms.ms.identified...., id_rate_great)*1.1) + 
        guides(color = guide_legend(title="ID class")) +
        scale_y_discrete_reverse(data$fc.raw.file, breaks = ggAxisLabels)
  #print(p)
  return(p)
}



#'
#' Colored table plot.
#' 
#' Code taken from http://stackoverflow.com/questions/23819209/change-text-color-for-cells-using-tablegrob-in-r
#' 
#' @param data Table as Data.frame
#' @param colours Single or set of colours (col-wise)
#' @param fill Cell fill (row-wise)
#' @param just (ignored)
#' @return gTable
#'
plotTableRaw = function(data, colours="black", fill=NA, just="centre")
{
  
  label_matrix = as.matrix(data)
  
  nc = ncol(label_matrix)
  nr = nrow(label_matrix)
  n = nc*nr
  
  colours <- rep(colours, length.out = n)
  fill <- rep(fill, length.out = n)
  justs <- rep(just, length.out = n)
  
  ## text for each cell
  labels <- lapply(seq_len(n), function(ii)
    grid::textGrob(as.character(label_matrix[ii]), gp = grid::gpar(fontsize=8, col=colours[ii]), just="left", x = grid::unit(0.05, "npc")))
  label_grobs <- matrix(labels, ncol=nc)
  
  ## define the fill background of cells
  fill <- lapply(seq_len(n), function(ii) 
    grid::rectGrob(gp = grid::gpar(fill=fill[ii])))
  
  ## some calculations of cell sizes
  row_heights <- function(m){
    do.call(grid::unit.c, apply(m, 1, function(l)
      max(do.call(grid::unit.c, lapply(l, grid::grobHeight)))))
  }
  col_widths <- function(m){
    do.call(grid::unit.c, apply(m, 2, function(l)
      max(do.call(grid::unit.c, lapply(l, grid::grobWidth)))))
  }
  
  ## place labels in a gtable
  g <- gtable::gtable_matrix("table", grobs = label_grobs, 
                             widths = col_widths(label_grobs) + grid::unit(2,"mm"), 
                             heights = row_heights(label_grobs) + grid::unit(2,"mm"))
  
  ## add the background
  xt <- rep(seq_len(nr), each=nc)
  xl <- rep(seq_len(nc), times=nr)
  g <- gtable::gtable_add_grob(g, fill, t=xt, l=xl, z=0, name="fill")
  
  return(g)
}

#'
#' Create an HTML table with an extra header row
#' 
#' 
#' @param data A data.frame which serves as table
#' @param caption A set of headlines, e.g. c("top line", "bottom line")
#' @return table as html character string for cat()'ing into an html document
#'
#' @import htmlTable
#' @import magrittr
#'
#' @export 
#' 
#' @examples
#'   data = data.frame(raw.file = letters[1:4],
#'                     id.rate = 3:6)
#'   getHTMLTable(data, 
#'                caption = "some header line")
#' 
getHTMLTable = function(data, caption = NA)
{
  
  tbl = htmlTable::addHtmlTableStyle(data,
                                     align = 'l',  ## align columns left
                                     col.rgroup = c("none", "#F7F7F7"))
  tbl = htmlTable::htmlTable(tbl, rnames = FALSE,    ## no row names
                             caption = caption) 

  return(tbl)
}

#'
#' Plot a table with row names and title
#' 
#' Restriction: currently, the footer will be cropped at the table width.
#' 
#' @param data A data.frame with columns as described above
#' @param title Table title
#' @param footer Footer text
#' @param col_names Column names for Table
#' @param fill Fill pattern (by row)
#' @param col Text color (by column)
#' @param just (ignored)
#' @return gTree object with class 'PTXQC_table'
#'
#' @export 
#' 
#' @examples
#'   data = data.frame(raw.file = letters[1:4],
#'                     id.rate = 3:6)
#'   plotTable(data, 
#'             title = "Bad files",
#'             footer = "bottom", 
#'             col_names = c("first col", "second col"),
#'             col=c("red", "green"))
#' 
plotTable = function(data, title = "", footer = "", col_names = colnames(data), fill = c("grey90", "grey70"), col = "black", just="centre")
{
  ## add column names
  data2 = rbind(col_names, apply(data, 2, function(x) as.character(x)))
  ## create table
  n = nrow(data2)*ncol(data2)
  nd = nrow(data)*ncol(data)
  table = plotTableRaw(data2, 
                       fill = c(rep("grey50", ncol(data)), rep(fill, each=ncol(data), length.out=nd)), ## row-wise
                       colours = unlist(lapply(col, function(cc) c("black", rep(cc, nrow(data))))), ## col-wise
                       just = c(rep("centre", ncol(data)), rep(just, each=nrow(data), length.out=nd))) 
  
  colhead = lapply(col_names, function(ii) grid::textGrob(ii, gp = grid::gpar(fontsize=12, col="black", fontface="bold", fill="grey")))
  ## replace column names
  table = gtable::gtable_add_grob(table, colhead, t = 1, l = 1:ncol(data))

  if (nchar(title[1]) > 0)
  {
    gtitle = grid::textGrob(title, gp = grid::gpar(fontsize = 14))
    padding = grid::unit(1.5, "line")
    ## add heading (white space)
    table = gtable::gtable_add_rows(table, heights = grid::grobHeight(gtitle) + padding, pos = 0)
    ## add heading (text as overlay)
    table = gtable::gtable_add_grob(table, list(gtitle), t = 1, l = 1, r = ncol(table), clip = "off")
  }
  if (nchar(footer[1]) > 0)
  {
    gfooter = grid::textGrob(footer, gp = grid::gpar(fontsize = 10))
    padding = grid::unit(1.5, "line")
    ## add heading (white space)
    table = gtable::gtable_add_rows(table, heights = grid::grobHeight(gfooter) + padding, pos = -1) ## bottom
    ## add heading (text as overlay)
    table = gtable::gtable_add_grob(table, list(gfooter), t = nrow(table), l = 1, r = ncol(table), clip = "off")
  }
  
  
  ## neat trick to enable calling print(g), to mimic ggplot-behaviour on this object
  ## in combination with print.PTXQC_table() -- see below
  p = grid::gTree(children = grid::gList(table), cl = c("PTXQC_table"))
  
  ## hide the table name inside (for qcMetric::getTitles())
  p$labels$title = title
    
  #print(p)
  return(p)
}

#' helper S3 class, enabling print(some-plot_Table-object)
#' @param x Some Grid object to plot
#' @param ... Further arguments (not used, but required for consistency with other print methods)
#' @return NULL
#' 
#' @export
#' 
print.PTXQC_table = function(x, ...) {
  grid::grid.newpage();
  grid::grid.draw(x)
  return(NULL)
}

#'
#' A boxplot of uncalibrated mass errors for each Raw file.
#'
#' Boxes are optionally colored to indicate that a MQ bug was detected or 
#' if PTXQC detected a too narrow search window.
#' 
#' @param data A data.frame with columns 'fc.raw.file', 'uncalibrated.mass.error..ppm.'
#' @param MQBug_raw_files List of Raw files with invalid calibration values
#' @param stats A data.frame with columns 'fc.raw.file', 'sd', 'outOfCal'
#' @param y_lim Range of y-axis
#' @param extra_limit Position where a v-line is plotted (for visual guidance) 
#' @param title_sub Subtitle
#' @return GGplot object
#' 
#' @import ggplot2
#' @export
#' 
#' @examples
#'   n = c(150, 1000, 1000, 1000)
#'   data = data.frame(fc.raw.file = repEach(letters[4:1], n),
#'                     uncalibrated.mass.error..ppm. = c(rnorm(n[1], 13, 2.4),
#'                                                       rnorm(n[2], 1, 0.5),
#'                                                       rnorm(n[3], 3, 0.7),
#'                                                       rnorm(n[4], 4.5, 0.8)))
#'   stats = data.frame(fc.raw.file = letters[4:1],
#'                      sd_uncal = c(2.4, 0.5, 0.7, 0.8),
#'                      outOfCal = c(TRUE, FALSE, FALSE, FALSE))           
#'   plot_UncalibratedMSErr(data, MQBug_raw_files = letters[1],
#'                          stats, y_lim = c(-20,20), 15, "subtitle")
#' 
plot_UncalibratedMSErr = function(data, MQBug_raw_files, stats, y_lim, extra_limit, title_sub)
{
  
  data$col = "default"
  if (length(MQBug_raw_files) > 0)
  {
    data$col[data$fc.raw.file %in% MQBug_raw_files] = "MQ bug"
  }
  ## add 'out-of-calibration' Raw files:
  data$col[data$fc.raw.file %in% stats$fc.raw.file[stats$outOfCal]] = "out-of-search-tol"
  ## only show legend if special things happen  
  showColLegend = ifelse(length(setdiff(data$col, "default")) > 0, "legend", "none")
  ## amend SD to fc.raw.file
  stats$fcr_new_lvl = paste0(stats$fc.raw.file, " (sd = ", stats$sd_uncal, "ppm)")
  
  ## use augmented name
  data$fc.raw.file_ext = stats$fcr_new_lvl[ match(data$fc.raw.file, stats$fc.raw.file) ]

  cols_sub = c("default"="black", "MQ bug"="red", "out-of-search-tol"="orange")
  cols_sub = cols_sub[names(cols_sub) %in% data$col]
  
  p = ggplot(data, col=data$col) +
        geom_boxplot(aes(x = .data$fc.raw.file_ext, y = .data$uncalibrated.mass.error..ppm., col = .data$col), varwidth = TRUE, outlier.shape = NA) +
        scale_colour_manual("", values = cols_sub, guide = showColLegend) +
        ylab(expression(Delta~"mass [ppm]")) +
        xlab("") +
        ylim(y_lim) +
        #scale_x_discrete_reverse(data$fc.raw.file_ext) +
        scale_x_discrete(limits=rev) +
        geom_hline(yintercept = c(-extra_limit, extra_limit), 
                   colour="red",
                   linetype = "longdash") +  ## == vline for coord_flip
        coord_flip() +
        ggtitle("EVD: Uncalibrated mass error", title_sub)

  #print(p)
  return(p)
}

#'
#' Plot bargraph of uncalibrated mass errors for each Raw file.
#'
#' Boxes are optionally colored to indicate that a MQ bug was detected or 
#' if PTXQC detected a too narrow search window.
#' 
#' @param data A data.frame with columns 'fc.raw.file', 'mass.error..ppm.'
#' @param MQBug_raw_files List of Raw files with invalid calibration values
#' @param stats A data.frame with columns 'fc.raw.file', 'outOfCal'
#' @param y_lim Range of y-axis
#' @param extra_limit Position where a v-line is plotted (for visual guidance) 
#' @param title_sub Subtitle
#' @return GGplot object
#' 
#' @import ggplot2
#' @export
#' 
#' @examples
#'   n = c(150, 1000, 1000, 1000)
#'   data = data.frame(fc.raw.file = repEach(letters[4:1], n),
#'                     mass.error..ppm. = c(rnorm(n[1], 1, 2.4),
#'                                          rnorm(n[2], 0.5, 0.5),
#'                                          rnorm(n[3], 0.1, 0.7),
#'                                          rnorm(n[4], 0.3, 0.8)))
#'   stats = data.frame(fc.raw.file = letters[4:1],
#'                      sd = c(2.4, 0.5, 0.7, 0.8),
#'                      outOfCal = c(TRUE, FALSE, FALSE, FALSE))           
#'   plot_CalibratedMSErr(data, MQBug_raw_files = letters[1], stats, y_lim = c(-20,20), 15, "subtitle")
#'
plot_CalibratedMSErr = function(data, MQBug_raw_files, stats, y_lim, extra_limit = NA, title_sub = "")
{
  data$col = "default"
  if (length(MQBug_raw_files) > 0) {
    data$col = c("default", "MQ bug")[(data$fc.raw.file %in% MQBug_raw_files) + 1]
    data$mass.error..ppm.[data$fc.raw.file %in% MQBug_raw_files] = 0
    if (all(data$mass.error..ppm.==0)) data$mass.error..ppm. = rnorm(nrow(data), sd=0.0001, mean=mean(y_lim))
  }
  ## add 'out-of-calibration' Raw files:
  data$col[data$fc.raw.file %in% stats$fc.raw.file[stats$outOfCal]] = "out-of-search-tol"
  ## only show legend if special things happen  
  showColLegend = ifelse(length(setdiff(data$col, "default")) > 0, "legend", "none")
  
  cols_sub = c("default"="black", "MQ bug"="red", "out-of-search-tol"="orange")
  cols_sub = cols_sub[names(cols_sub) %in% data$col]
  
  ## plot
  p = ggplot(data, col = data$col) +
    geom_boxplot(aes(x = .data$fc.raw.file, y = .data$mass.error..ppm., col = .data$col), varwidth = TRUE, outlier.shape = NA) +
    scale_colour_manual("", values = cols_sub, guide = showColLegend) +
    ylab(expression(Delta~"mass [ppm]")) +
    xlab("") +
    ylim(y_lim) +
    #scale_x_discrete_reverse(data$fc.raw.file) +
    scale_x_discrete(limits=rev) +
    coord_flip() +
    ggtitle("EVD: Calibrated mass error", title_sub)
  if (!is.na(extra_limit)) {
    p = p + geom_hline(yintercept = c(-extra_limit, extra_limit), colour="red", linetype = "longdash")  ## == vline for coord_flip
  }
  
  #print(p)
  return (p)
}


#'
#' Plot bargraph of oversampled 3D-peaks.
#'
#' Per Raw file, at most three n's must be given, i.e.
#' the fraction of 3D-peaks for n=1, n=2 and n=3(or more).
#' The fractions must sum to 1 (=100%).
#' 
#' 
#' @param data A data.frame with columns 'fc.raw.file', 'n', 'fraction'
#' @return GGplot object
#' 
#' @import ggplot2
#' @export
#' 
#' @examples
#'   data = data.frame(fc.raw.file = rep(letters[1:3], each=3),
#'                     n = 1:3,
#'                     fraction = c(0.8, 0.1, 0.1, 0.6, 0.3, 0.1, 0.7, 0.25, 0.05))
#'   plot_MS2Oversampling(data)
#'
plot_MS2Oversampling = function(data)
{
  stopifnot(length(unique(data$n)) <= 3) ## at most three -- to match color vector below
  #data = d_dups
  ## reorder factor, such that '10+' is last
  data$n = as.character(data$n)
  n_unique = sort(unique(data$n)) ## sort as character vector!
  data$n = factor(data$n, levels=n_unique[order(nchar(n_unique))], ordered = TRUE)
  
  p = ggplot(data) + 
        geom_col(position = position_stack(reverse = TRUE), aes(x = .data$fc.raw.file, y = .data$fraction, fill = .data$n)) +
        scale_fill_manual("MS/MS\ncounts", values =c("green", "blue", "red")) +
        scale_x_discrete_reverse(data$fc.raw.file) +
        xlab("") +
        ylab("MS/MS counts per 3D-peak [%]") +
        ggtitle(paste0("EVD: Oversampling (MS/MS counts per 3D-peak)")) + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        coord_flip()
  
  #print(p)
  return(p)
}


#'
#' Plot bargraph of oversampled 3D-peaks.
#'
#' Per Raw file, at most three n's must be given, i.e.
#' the fraction of 3D-peaks for n=1, n=2 and n=3(or more).
#' The fractions must sum to 1 (=100%).
#' 
#' 
#' @param data A data.frame with columns 'file', 'msErr', 'type'
#' @return GGplot object
#' 
#' @import ggplot2
#' @export
#' 
#' @examples
#'   n = c(100, 130, 50)
#'   data = data.frame(file = repEach(paste(letters[1:3],"\nLTQ [Da]"), n),
#'                     msErr = c(rnorm(n[1], 0.5), rnorm(n[2], 0.0), rnorm(n[3], -0.5)),
#'                     type = c("forward", "decoy")[1+(runif(sum(n))>0.95)])
#'   plot_MS2Decal(data)
#'
plot_MS2Decal = function(data)
{
  ## trim down the data to 2-98 percentiles (to avoid outliers far off)
  data2 = plyr::ddply(data, "file", function(x) {
    qnt = quantile(x$msErr, probs = c(0.02, 0.98), na.rm = TRUE)
    return (x[qnt[1] < x$msErr & x$msErr < qnt[2], ])
  })
  p = ggplot(data2, aes(x = .data$msErr, fill = .data$type)) + 
    geom_histogram(bins = 30) +
    xlab("fragment mass delta") +  
    ylab("count") + 
    scale_fill_manual(values = c(forward = "#99d594", decoy = "#ff0000")) +
    ggtitle("MSMS: Fragment mass errors per Raw file") +
    facet_wrap(~file, scales = "fixed")
  
  #print(p)
  return(p)
}

#'
#' Plot bargraph of missed cleavages.
#'
#' Per Raw file, an arbitrary number of missed cleavage classes (one per column) can be given.
#' The total fraction of 3D-peaks must sum to 1 (=100%).
#' Columns are ordered by name.
#' 
#' A visual threshold line is drawn at 75% (expected MC0 count).
#' 
#' @param data A data.frame with columns 'fc.raw.file', '...' (missed cleavage classes)
#' @param title_sub Plot's subtitle
#' @return GGplot object
#' 
#' @import ggplot2
#' @export
#' 
#' @examples
#'   data = data.frame(fc.raw.file = letters[1:5],
#'                     MC0 = c(0.8, 0.5, 0.85, 0.2, 0.9),
#'                     MC1 = c(0.1, 0.4, 0.05, 0.7, 0.0),
#'                     "MS2+" =  c(0.1, 0.1, 0.1, 0.1, 0.1),
#'                     check.names = FALSE)
#'   plot_MissedCleavages(data, "contaminant inclusion unknown")
#'
plot_MissedCleavages = function(data, title_sub = "")
{
  st_bin.m = reshape2::melt(data, id.vars = c("fc.raw.file"))
  p =
    ggplot(data = st_bin.m, aes(x = factor(.data$fc.raw.file), y = .data$value, fill = .data$variable)) + 
        geom_col(position = position_stack(reverse = TRUE)) +
        xlab("Raw file") +  
        ylab("missed cleavages [%]") + 
        theme(legend.title = element_blank()) +
        scale_fill_manual(values = rep(c("#99d594", "#ffffbf", "#fc8d59", "#ff0000", "#800080", "#000000"), 10)) +
        geom_abline(alpha = 0.5, intercept = 0.75, slope = 0, colour = "black", linetype = "dashed", linewidth = 1.5) +
        coord_flip() +
        scale_x_discrete_reverse(st_bin.m$fc.raw.file) +
        ggtitle("MSMS: Missed cleavages per Raw file", title_sub)
  
  #print(p)
  return(p)
}

#'
#' Plot line graph of TopN over Retention time.
#'
#' Number of Raw files must be 6 at most. Function will stop otherwise.
#' 
#' @param data A data.frame with columns 'fc.raw.file', 'rRT', 'topN'
#' @return GGplot object
#' 
#' @import ggplot2
#' @export
#' 
#' @examples
#'   data = data.frame(fc.raw.file = rep(letters[1:3], each=100),
#'                     rRT = seq(20, 120, length.out = 100),
#'                     topN = c(round(runif(100, min=3, max=5)),
#'                              round(runif(100, min=5, max=8)),
#'                              round(runif(100, min=1, max=3)))
#'                     )
#'   plot_TopNoverRT(data)
#'
plot_TopNoverRT = function(data)
{
  nrOfRaws = length(unique(data$fc.raw.file))
  p = ggplot(data, aes(x = .data$rRT, y = .data$topN, col = .data$fc.raw.file)) +
        geom_line() +
        scale_color_manual(values = brewer.pal.Safe(nrOfRaws, "Set1")) +
        xlab("retention time [min]") +
        ylab("highest N [median per RT bin]") +
        #stat_smooth(method = "loess", formula = y ~ x, se = FALSE, span = 0.1) +
        guides(color=guide_legend(title="")) +
        ggtitle("MSMSscans: TopN over RT")
    
  #print(p)
  return (p)
}

#'
#' Plot line graph of TopN over Retention time.
#'
#' Number of Raw files must be 6 at most. Function will stop otherwise.
#' 
#' @param data A data.frame with columns 'fc.raw.file', 'rRT', 'medIIT'
#' @param stats A data.frame with columns 'fc.raw.file', 'mean'
#' @param extra_limit Visual guidance line (maximum acceptable IIT)
#' @return GGplot object
#' 
#' @import ggplot2
#' @export
#' 
#' @examples
#'   data = data.frame(fc.raw.file = rep(c("d","a","x"), each=100),
#'                     rRT = seq(20, 120, length.out = 100),
#'                     medIIT = c(round(runif(100, min=3, max=5)),
#'                                round(runif(100, min=5, max=8)),
#'                                round(runif(100, min=1, max=3)))
#'                     )
#'   stats = data.frame(fc.raw.file = c("d","a","x"),
#'                      mean = c(4, 6.5, 2))
#'   plot_IonInjectionTimeOverRT(data, stats, 10)
#'
plot_IonInjectionTimeOverRT = function(data, stats, extra_limit)
{
  data$fc.raw.file = data$fc.raw.file[,drop = TRUE] ## drop unused factor levels
  nrOfRaws = length(unique(data$fc.raw.file))
  stats_sub = stats[stats$fc.raw.file %in% data$fc.raw.file, , drop = FALSE]
  ## augment legend with average II-time[ms]
  data$fc.raw.file = paste0(data$fc.raw.file, " (~", 
                            round(stats_sub$mean[match(data$fc.raw.file, stats_sub$fc.raw.file)]),
                            " ms)")
  ## manually convert to factor to keep old ordering (otherwise ggplot will sort it, since its a string)
  data$fc.raw.file = factor(data$fc.raw.file, levels = unique(data$fc.raw.file), ordered = TRUE)
  stats_sub$fc.raw.file = paste0(stats_sub$fc.raw.file, " (~", 
                                 round(stats_sub$mean[match(stats_sub$fc.raw.file, stats_sub$fc.raw.file)]),
                                 " ms)")
  p = ggplot(data) +
        geom_line(aes(x = .data$rRT, y = .data$medIIT, col = .data$fc.raw.file)) +
        scale_color_manual(values = brewer.pal.Safe(nrOfRaws, "Set1")) +
        xlab("retention time [min]") +
        ylab("ion injection time [ms]") +
        geom_hline(yintercept = extra_limit, linetype = 'dashed') +
        guides(color=guide_legend(title="Raw file with\naverage inj. time")) +
        ggtitle("MSMSscans: Ion Injection Time over RT") +
        pointsPutX(x_range = range(data$rRT), x_section = c(0.03, 0.08), y = stats_sub$mean, col = stats_sub$fc.raw.file[,drop = TRUE])
  
  #print(p)
  return(p)
}

#'
#' Plot line graph of TopN over Retention time.
#'
#' Number of Raw files must be 6 at most. Function will stop otherwise.
#' 
#' @param data A data.frame with columns 'fc.raw.file', 'scan.event.number', 'n'
#' @return GGplot object
#' 
#' @import ggplot2
#' @export
#' 
#' @examples
#'   data = data.frame(fc.raw.file = rep(c("d","a","x"), each=10),
#'                     scan.event.number = 1:10,
#'                     n = 11:20)
#'   plot_TopN(data)
#'
plot_TopN = function(data)
{
  
  p = ggplot(data, aes(x = .data$scan.event.number, y = .data$n)) +
        geom_col() +
        xlab("highest scan event") +
        ylab("count") +
        facet_wrap(~ fc.raw.file, scales = "free_y") +
        ggtitle(paste0("MSMSscans: TopN"))
  
  #print(p)
  return(p)
}

#'
#' Plot line graph of TopN over Retention time.
#'
#' Number of Raw files must be 6 at most. Function will stop otherwise.
#' 
#' @param data A data.frame with columns 'fc.raw.file', 'scan.event.number', 'ratio', 'count'
#' @return GGplot object
#' 
#' @import ggplot2
#' @export
#' 
#' @examples
#'   data = data.frame(fc.raw.file = factor(rep(c("d","a","x"), each=10), levels = c("d","a","x")),
#'                     scan.event.number = 1:10,
#'                     ratio = seq(40, 20, length.out=10),
#'                     count = seq(400, 200, length.out=10))
#'   plot_ScanIDRate(data)
#'
plot_ScanIDRate = function(data)
{
  
  p = ggplot(data, aes(x = .data$scan.event.number, y = .data$ratio, alpha = .data$count)) +
        geom_col() +
        xlab("scan event") +
        ylab("percent identified") +
        facet_wrap(~ fc.raw.file) +
        ggtitle(paste0("MSMSscans: TopN % identified over N"))
  return (p)
}


#'
#' Plot Total Ion Count over time
#' 
#' The input is a data.frame with already averaged counts over binned RT-slices.
#' 
#' @param data A data.frame with columns 'fc.raw.file', 'RT', 'intensity'
#' @param x_lim Plot range of x-axis
#' @param y_lim Plot range of y-axis
#' @return GGplot object
#'
#' @import ggplot2
#' @export
#' 
#' @examples 
#' 
#'  data = data.frame(fc.raw.file = rep(c("file A", "file B", "file C"), each=81),
#'                    RT = c(20:100), 
#'                    intensity = c(rnorm(81, mean=20), rnorm(81, mean=10), rnorm(81, mean=30)))
#'  plot_TIC(data, c(10, 100), c(0, 40))
#' 
plot_TIC = function(data, x_lim, y_lim)
{
  p = ggplot(data) +
    geom_line(aes(x = .data$RT, y = .data$intensity, colour = .data$fc.raw.file), linewidth = 1, alpha = 0.7) +
    scale_color_manual(values = brewer.pal.Safe(length(unique(data$fc.raw.file)), "Set1")) +
    guides(color = guide_legend(title = "Raw file\n(avg. peak width)")) +
    xlab("retention time [min]") +
    ylab("intensity") +
    coord_cartesian(xlim = x_lim, ylim = y_lim) + ## zoom in y -- do not cut data (preserve lines)
    ggtitle("SM: Total Ion Count") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  #print(p)
  return(p)
}