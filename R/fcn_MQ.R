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
#' Given a data.frame with two/three columns in long format (name, value, [contaminant]; in that order), each group (given from 1st column)
#' is plotted as a bar.
#' Contaminants (if given) are seperated and plotted as yellow bars.
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
#' @param coord_flip      Exchange Y and X-axis for better readability
#' @param names           An optional data.frame(long=.., short=..), giving a renaming scheme for the 'group' column
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
                           abline = NA,
                           coord_flip = T,
                           names = NA)
{
 
  if (ncol(data) == 2) {
    data$contaminant = FALSE ## add a third column, if missing
  }
  colnames(data) = c("group", "value", "contaminant")
  
  if (log2) {
    data$value = log2(data$value)
  }
  
  ## shorten group names
  if (class(names)=='data.frame')
  {
    stopifnot(sort(colnames(names)) == c("long", "short"))
    data$group = names$short[match(data$group, names$long)]
    if (any(is.na(data$group))) 
    {
      print(names)
      stop("Group renaming is incomplete! Aborting...")
    }
  }
  
  ## actual number of entries in each column (e.g. LFQ often has 0)
  ncol.stat = ddply(data, colnames(data)[1], function(x){ notNA = sum(!is.infinite(x$value) & !is.na(x$value));
                                                          data.frame(n = nrow(x), 
                                                                     notNA = notNA, 
                                                                     newname = paste0(x$group[1], " (n=", notNA, ")"))})
  
  ## rename (augment with '#n')
  data$group = ncol.stat$newname[match(data$group, ncol.stat$group)]
  
  ## remove -inf and NA's
  data = data[!is.infinite(data$value) & !is.na(data$value), ]

  groups = unique(data$group);
  ## add color for H vs L (if SILAC)
  cols = c("sample" = "black", 
           "sample (light)" = "black", 
           "sample (medium)" = "blue", 
           "sample (heavy)" = "green",
           "contaminant" = "yellow")
  cat_names = names(cols)
  cat = factor(cat_names, levels=cat_names)
  data$cat = cat[1]
  if (sum(grepl("^[^HLM]", groups )) == 0 || sum(grepl("^intensity\\.[hlm]\\.", groups )) > 0) { ## all start with either L, M or H
    data$cat = cat[2]
    data$cat[grepl("^M", data$group) | grepl("^intensity\\.m\\.", data$group)] = cat[3]
    data$cat[grepl("^H", data$group) | grepl("^intensity\\.h\\.", data$group)] = cat[4]
  }
  data$cat[data$contaminant] = cat[5]
  
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
    pl = ggplot(data=data, aes_string(x = "group", y = "value", fill = "cat")) + ## do not use col="cat", since this will dodge bars and loose scaling
          geom_boxplot(varwidth=T) +
          xlab("") + 
          ylab(ylab) +
          ylim(ylims) +
          scale_alpha(guide=FALSE) +
          scale_fill_manual(values=cols, name = "Category") + 
          scale_color_manual(values=cols, name = "Category") + 
          theme(axis.text.x = element_text(angle=90, vjust = 0.5)) +
          theme(legend.position=ifelse(length(cols)==1, "none", "right")) +
          addGGtitle(mainlab, sublab) + 
          scale_x_discrete_reverse(unique(data$group))
    
    if (!is.na(abline))
    {
      pl = pl + geom_abline(alpha = 0.5, intercept = abline, slope = 0, colour = "green")
    }
    if (coord_flip == T)
    {
      pl = pl + coord_flip()
    }
    #print(pl)
    return(pl)
  }
  lpl = byXflex(data = data, indices = data$group, subset_size = boxes_per_page, sort_indices = F, FUN = fcn_boxplot_internal, abline)
  return (lpl)
}


#'
#' Extract fragment mass deviation errors from a data.frame from msms.txt
#' 
#' Given a data.frame as obtainable from a msms.txt with 
#'  - a 'mass.analyzer' column which contains only a single value for the whole column
#'  - a 'mass.deviations..da.' and (if available) 'mass.deviations..ppm.'
#'  - a 'masses' column (only required if 'mass.deviations..ppm.' is unavailable and the mass.analyzer
#'    indicates hig-res data)
#' 
#' Mass deviations are extracted from the columns, e.g. each cell containing values separated by
#' semicolons is split into single values. The appropriate unit is chosen (Da or ppm, depending on
#' ITMS or FTMS data). Also the fragmentation type can be used: CID indicates ITMS, HCD to FTMS.
#' This is not 100% safe, but older MQ versions do not report the mass analyzer properly.
#' 
#' If ppm mass deviations are not available, errors in Da will be converted to ppm using the corresponding mass values.
#' 
#' @param x Data frame in long format with numerical expression data
#' @return  Data frame with mass errors ('msErr') and their 'unit' (Da or ppm)
#' 
#' @export
#' 
getFragmentErrors = function(x)
{
  ## require only one mass analyzer type:
  stopifnot(length(unique(x$mass.analyzer))==1)
  stopifnot(all(c("mass.analyzer", "mass.deviations..da.") %in% colnames(x)))
  
  convert_Da2PPM = FALSE
  if (grepl("ITMS|TOF|CID", x$mass.analyzer[1]) & ("mass.deviations..da." %in% colnames(x)))
  {
    ms2_unit = "[Da]"; ms2_col = "mass.deviations..da."
  } else
    if (grepl("FTMS|HCD", x$mass.analyzer[1]) & ("mass.deviations..ppm." %in% colnames(x)))
    {
      ms2_unit = "[ppm]"; ms2_col = "mass.deviations..ppm."
    } else
      if (grepl("FTMS|HCD", x$mass.analyzer[1]) & ("mass.deviations..da." %in% colnames(x)))
      {
        ## we know its high resolution, but this MQ version only gave us Dalton mass deviations
        ## --> convert back to ppm
        ms2_unit = "[ppm]"; ms2_col = "mass.deviations..da."
        convert_Da2PPM = TRUE
      }
  err = unlist(strsplit(paste(x[, ms2_col], sep="", collapse=";"), split=";", fixed=T))
  if (convert_Da2PPM) {
    stopifnot("masses" %in% colnames(x))
    mass = unlist(strsplit(paste(x$masses, sep="", collapse=";"), split=";", fixed=T))
    err = as.numeric(err) / as.numeric(mass) * 1e6
  }
  ## return as character, otherwise it will get converted to factor by ddply?
  return(data.frame(msErr = as.character(err), unit = ms2_unit))
}

