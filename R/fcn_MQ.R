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
#' Given two data.frames, each with matching rows and columns and just containing expression data.
#' We plot the column ratios and annotate with correlation.
#' If just one data.frame is given, we assume that it contains the data 
#' (e.g. ratios, or LFQ) already (no correlation will be reported).
#' 
#' Boxes are shaded: many NA or Inf lead to more transparency. Allows to easily spot sparse columns.
#' 
#' 
#' @param data1 Data frame with numerical expression data
#' @param data2 Optional second data frame (with matching rows and columns)
#' @param log2_ratios Apply log2 to the data (yes/no)
#' @param ylab    Label on Y-axis
#' @param mainlab Main title
#' @param sublab  Sub title
#' @param boxes_per_page  Maximum number of boxplots per plot. Yields multiple plots if more columns are given.
#' @param abline Draw a horziontal green line at the specified y-position (e.g. to indicate target median values)
#' 
#' @return List of printed ggplots
#' 
#' @import ggplot2
#' @importFrom reshape2 melt
#' 
#' @export
#' 
boxplotCompare <- function(data1, data2 = NA, 
                           log2_ratios = T,
                           ylab = "intensity",
                           mainlab = ylab,
                           sublab = "",
                           boxes_per_page = 30,
                           abline = NA)
{
  if (is.na(data2)) {
    ratios = data1
  } else {
    ratios = data1/data2
    ylab = paste(ylab, "ratio")
  } 
  
  if (log2_ratios) {
    ratios = log2(ratios)
    ylab = paste("log2", ylab)
  }
  
  #lcp = nchar(lcPrefix( colnames(ratios) ))   # shorten name (remove common prefix)
  #colnames(ratios) = sapply(colnames(ratios), substr, lcp+1, 100000)
  colnames(ratios) = delLCP(colnames(ratios))
  
  ## maximum number of possible entries
  nmax = nrow(ratios)
  ## actual number of entries in each column (e.g. LFQ often has 0)
  ncol.stat = apply(ratios, 2, function(x) sum(!is.infinite(x) & !is.na(x)))
  
  if (!is.na(data2))
  {
    cors = sapply(1:ncol(data1), function(i) {
      x = data1[, i]
      y = data2[, i]
      ok <- is.finite(x) & is.finite(y)
      r <- (cor(x[ok], y[ok], use="pairwise.complete.obs"))
      txt <- format(c(r, 0.123456789), digits=2)[1]
      txt
    })
    colnames(ratios) = paste(colnames(ratios), " (n=", ncol.stat, ", c=", cors, ")", sep="")
  } else {
    colnames(ratios) = paste(colnames(ratios), " (n=", ncol.stat, ")", sep="")
  }
  
  ## compute alpha value for plotting
  ncol.stat.alpha.df = data.frame(variable = colnames(ratios), alphav = ncol.stat/nmax)
  
  ## long table
  #require(reshape2)
  datar = melt(ratios)
  head(datar)
  ## remote -inf and NA's
  datar = datar[!is.infinite(datar$value) & !is.na(datar$value), ]
    
  ## add color for H vs L (if SILAC)
  cat = factor(c("light", "medium", "heavy"), levels=c("light", "medium", "heavy"))
  datar$cat = cat[1]
  if (sum(grepl("^[^HLM]", colnames(ratios) )) == 0) { ## all start with either L, M or H
    datar$cat[grep("^M", datar$variable)] = cat[2]
    datar$cat[grep("^H", datar$variable)] = cat[3]
  }
  cols = c("black", "blue", "red")[unique(datar$cat)]
  
  ## augment with alpha values
  datar$alphav = ncol.stat.alpha.df$alphav[match(datar$variable, ncol.stat.alpha.df$variable)] 
  
  #ex: datar$section = as.integer(as.numeric(datar$section)/boxes_per_page)
  
  ## compute global y-limits (so we can fix it across plots)
  ylims = boxplot.stats(datar$value)$stats[c(1, 5)]
  fcn_boxplot_internal = function(data, abline = NA) 
  {
    #require(ggplot2)
    pl = ggplot(data=data, aes_string(x = "variable", y = "value")) +
      geom_boxplot(aes_string(fill = "cat", alpha = "alphav")) + 
      xlab("data set") + 
      ylab(ylab) +
      ylim(ylims) +
      scale_fill_manual(values=cols) + 
      theme(axis.text.x = element_text(angle=90)) +
      theme(legend.position=ifelse(length(cols)==1, "none", "right")) +
      scale_alpha(range=c(min(ncol.stat.alpha.df$alphav), max(ncol.stat.alpha.df$alphav)))
    pl = addGGtitle(pl, mainlab, sublab)
    if (!is.na(abline))
    {
      pl = pl + geom_abline(alpha = 0.5, intercept = abline, slope = 0, colour = "green")
    }
    print(pl)
  }
  #ex: fcn_boxplot_internal(datar[datar$section<2,])
  byXflex(data = datar, indices = datar$variable, subset_size = boxes_per_page, sort_indices = F, FUN = fcn_boxplot_internal, abline)
  
}



