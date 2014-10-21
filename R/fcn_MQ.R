###
### Author: Chris Bielow
###
###

#install.packages("ggplot2")
#source("http://bioconductor.org/biocLite.R")
#biocLite("ggplot2")
#biocLite("plyr")
#require(plyr)


#' Boxplots - one for each condition (=column) in a data frame.
#' 
#' Given a data.frame with two columns in long format (name, value; in that order), each group (given from 1st column)
#' is plotted as a bar.
#' 
#' Boxes are shaded: many NA or Inf lead to more transparency. Allows to easily spot sparse groups
#' 
#' 
#' @param data    Data frame in long format with numerical expression data
#' @param log2    Apply log2 to the data (yes/no)
#' @param ylab    Label on Y-axis
#' @param mainlab Main title
#' @param sublab  Sub title
#' @param boxes_per_page  Maximum number of boxplots per plot. Yields multiple plots if more groups are given.
#' @param abline          Draw a horziontal green line at the specified y-position (e.g. to indicate target median values)
#' 
#' @return List of ggplot objects
#' 
#' @import ggplot2
#' @importFrom plyr ddply
#' 
#' @export
#' 
boxplotCompare <- function(data, 
                           log2 = T,
                           ylab = "intensity",
                           mainlab = ylab,
                           sublab = "",
                           boxes_per_page = 30,
                           abline = NA)
{
 
  colnames(data) = c("group", "value")
  if (log2) {
    data$value = log2(data$value)
  }
  

  ## actual number of entries in each column (e.g. LFQ often has 0)
  ncol.stat = ddply(data, colnames(data)[1], function(x){ notNA = sum(!is.infinite(x$value) & !is.na(x$value));
                                                          data.frame(n = nrow(x), 
                                                                    notNA = notNA, 
                                                                    newname = paste0(x$group[1], " (n=", notNA, ")"))})
  ## compute alpha value for plotting
  ncol.stat$alpha = ncol.stat$notNA / max(ncol.stat$n)
  
  ## rename (augment with 'n')
  data$group = ncol.stat$newname[match(data$group, ncol.stat$group)]
  
  ## remote -inf and NA's
  data = data[!is.infinite(data$value) & !is.na(data$value), ]

  groups = unique(data$group);
  ## add color for H vs L (if SILAC)
  cat = factor(c("light", "medium", "heavy"), levels=c("light", "medium", "heavy"))
  data$cat = cat[1]
  if (sum(grepl("^[^HLM]", groups )) == 0) { ## all start with either L, M or H
    data$cat[grep("^M", data$group)] = cat[2]
    data$cat[grep("^H", data$group)] = cat[3]
  }
  cols = c("black", "blue", "red")[unique(data$cat)]
  
  ## augment with alpha values
  data$alphav = ncol.stat$alpha[match(data$group, ncol.stat$group)] 
  
  #ex: datar$section = as.integer(as.numeric(datar$section)/boxes_per_page)
  
  ## compute global y-limits (so we can fix it across plots)
  ylims = boxplot.stats(data$value)$stats[c(1, 5)]
  ## make sure to inlude abline (if existing)
  if (!is.na(abline))
  {
    ylims = c(ylims[1], max(ylims[2], abline))
  }
  fcn_boxplot_internal = function(data, abline = NA) 
  {
    #require(ggplot2)
    pl = ggplot(data=data, aes_string(x = "group", y = "value")) +
      geom_boxplot(aes_string(fill = "cat", alpha = "alphav")) + 
      xlab("data set") + 
      ylab(ylab) +
      ylim(ylims) +
      scale_fill_manual(values=cols) + 
      theme(axis.text.x = element_text(angle=90, vjust = 0.5)) +
      theme(legend.position=ifelse(length(cols)==1, "none", "right")) +
      scale_alpha(range=c(min(ncol.stat$alpha), max(ncol.stat$alpha)))
    pl = addGGtitle(pl, mainlab, sublab)
    if (!is.na(abline))
    {
      pl = pl + geom_abline(alpha = 0.5, intercept = abline, slope = 0, colour = "green")
    }
    return(pl)
  }
  #ex: fcn_boxplot_internal(datar[datar$section<2,])
  lpl = byXflex(data = data, indices = data$group, subset_size = boxes_per_page, sort_indices = F, FUN = fcn_boxplot_internal, abline)
  return (lpl)
}



