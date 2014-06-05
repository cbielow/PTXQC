###
### Author: Chris Bielow
###
###

#install.packages("ggplot2")
#source("http://bioconductor.org/biocLite.R")
#biocLite("ggplot2")
#biocLite("plyr")
#require(plyr)

#' Wrapper to read a MQ txt file (e.g. proteinGroups.txt).
#' 
#' Since MaxQuant changes capitalization and sometimes even column names, it seemed convenient
#' to have a function which just reads a txt file and returns unified column names, irrespective of the MQ version.
#' So, it unifies access to columns (e.g. by using lower case for ALL columns) and ensures columns are
#' identically named across MQ versions. We only of one case: "protease" == "enzyme".
#'
#' If the file is empty, this function stops with an error.
#'
#' @param file   (relative) path to a MQ txt file ()
#' @param filter Searched for "C" and "R". If present, [c]ontaminants and [r]everse hits are removed
#' @param type   Allowed values are:
#'               "pg" (proteinGroups) [default], adds abundance index columns (*AbInd*, replacing 'intensity')
#'               "sm" (summary), splits into three row subsets (raw.file, condition, total)
#'               Any other value will not add any special columns
#' @param col_subset  A vector of column names as read by read.delim(), e.g., spaces are replaced by dot already.
#'                   If given, only columns with these names (ignoring lower/uppercase) will be returned (regex allowed)
#'                   E.g. col_subset=c("^lfq.intensity.", "protein.name")
#' @param add_fs_col If TRUE and a column 'raw.file' is present, an additional column 'fc.raw.file' will be added with 
#'                   common prefix AND common substrings removed (\code{\link{simplifyNames}})
#'                           E.g. two rawfiles named 'OrbiXL_2014_Hek293_Control', 'OrbiXL_2014_Hek293_Treated' will give
#'                                                   'Control', 'Treated'
#' @param LFQ_action [for type=='pg' only] An additional custom LFQ column ('cLFQ...') is created where
#'               zero values in LFQ columns are replaced by the following method IFF(!) the corresponding raw intensity is >0 (indicating that LFQ is erroneusly 0)
#'               "toNA": replace by NA
#'               "impute": replace by lowest LFQ value >0 (simulating 'noise')
#' @param ... Additional parameters passed on to read.delim()             
#' @return data frame

#' @import utils
#' @import graphics
#' 
#' @export
readMQ <- function(file, filter="C+R", type="pg", col_subset=NA, add_fs_col=T, LFQ_action=FALSE, ...)
{
  cat(paste("Reading file", file,"...\n"))
  msg_parse_error = paste0("\n\nParsing the file '", file, "' failed. See message above why. If the file is not usable but other files are ok, disable the corresponding section in the YAML config.")
  
  ## resolve set of columns which we want to keep
  #col_subset = c("Multi.*", "^Peaks$")
  if (sum(sapply(col_subset, function(x) !is.na(x))) > 0)
  { ## just read a tiny bit to get column names
    data_header = try(read.delim(file, na.strings=c("NA", "n. def."), stringsAsFactors=F, comment.char="#", nrows=2))
    if (inherits(data_header, 'try-error')) stop(msg_parse_error, call.=F);
    
    colnames(data_header) = tolower(colnames(data_header))
    idx_keep = rep(FALSE, ncol(data_header))    
    for (valid in col_subset)
    {
      idx_new = grepl(valid, colnames(data_header), ignore.case = T)
      if (sum(idx_new) == 0) cat(paste0("WARNING: Could not find column regex '", valid, "' using case-INsensitive matching.\n"))
      idx_keep = idx_keep | idx_new
    }
    #summary(idx_keep)
    col_subset = colnames(data_header)
    col_subset[ idx_keep] = NA     ## default action for selected columns
    col_subset[!idx_keep] = "NULL" ## skip over unselected columns during reading
    cat(paste("Keeping", sum(idx_keep),"of",length(idx_keep),"columns!\n"))
    #print (colnames(data_header))
  }
  
  data = try(read.delim(file, na.strings=c("NA", "n. def."), stringsAsFactors=F, comment.char="#", colClasses = col_subset, ...))
  if (inherits(data, 'try-error')) stop(msg_parse_error, call.=F);
  
  #colnames(data)
  
  cat(paste0("Read ", nrow(data), " entries from ", file,".\n"))
  
  ### just make everything lower.case (MQ versions keep changing it and we want it to be reproducible)
  colnames(data) = tolower(colnames(data))
  ## rename some columns since MQ 1.2 vs. 1.3 differ....
  colnames(data)[colnames(data)=="protease"] = "enzyme"
  
  if (add_fs_col & "raw.file" %in% colnames(data))
  {
    rf_name = unique(data$raw.file)
    ## remove prefix
    rf_name_s = delLCP(rf_name)
    ## remove infix (2 iterations)
    rf_name_s = simplifyNames(rf_name_s, 2, min_LCS_length=7)
    
    ## check if shorter filenames are still unique (they should be.. if not we have a problem!!)
    if (length(rf_name) != length(unique(rf_name_s)))
    {
      cat("Original names:\n")
      cat(rf_name)
      cat("Short names:\n")
      cat(rf_name_s)
      stop("While loading MQ data: shortened raw filenames are not unique! This should not happen. Please contact chris.bielow@mdc-berlin.de and provide the above names!")
    }
    
    data$fc.raw.file = as.factor(rf_name_s[match(data$raw.file, rf_name)])
  }
  
  ## proteingroups.txt special treatment
  if (type=="pg")
  {
    if (grepl("C", filter)) data = data[data$contaminant != "+",]
    if (grepl("R", filter)) data = data[data$reverse != "+",]

    stats = data.frame(n=NA, v=NA)
    if (LFQ_action!=FALSE)
    { ## replace erroneous zero LFQ values with something else
      cat("Starting LFQ action.\nReplacing ...\n")
      lfq_cols = grepv("^lfq", colnames(data))
      for (cc in lfq_cols)
      {
        ## get corresponding raw intensity column
        rawint_col = sub("^lfq\\.", "", cc)
        if (!(rawint_col %in% colnames(data))) {stop(paste0("Could not find column '", rawint_col, "' in dataframe with columns: ", paste(colnames(data), collapse=",")), "\n")}
        vals = data[, cc]
        ## affected rows
        bad_rows = (data[, rawint_col]>0 & data[, cc]==0)
        if (sum(bad_rows, na.rm=T)==0) {next;}
        ## take action
        if (LFQ_action=="toNA" | LFQ_action=="impute") {
          ## set to NA
          impVal = NA;
          vals[bad_rows] = impVal;
          if (LFQ_action=="impute") {
            impVal = min(vals[vals>0], na.rm=T)
            ## replace with minimum noise value (>0!)
            vals[bad_rows] = impVal;
          }
          cat(paste0("   '", cc, "' ", sum(bad_rows, na.rm=T), ' entries (', sum(bad_rows, na.rm=T)/nrow(data)*100,'%) with ', impVal, '\n'))
          ## add column
          data[, paste0("c", cc)] = vals;
          ##
          stats = rbind(stats, c(sum(bad_rows, na.rm=T), impVal))
          
        }
        else {
          stop(paste0("Unknown action '", LFQ_action, "' for LFQ_action parameter! Aborting!"));
        }
      }
      if (LFQ_action=="impute") {
        plot(stats, log="y", xlab="number of replaced zeros", ylab='imputed intensity', main='Imputation statistics')
      }
    }
    
    ### add abundance index columns (for both, intensity and lfq.intensity)
    int_cols= colnames(data)[grep("intensity", colnames(data))]
    data[, sub("intensity", "AbInd", int_cols)] = apply(data[,int_cols], 2, function(x)
                                                  {
                                                    x / data[,"mol..weight..kda."]
                                                  })
    
  } else  ## summary.txt special treatment
  if (type=="sm")
  { ## split
    idx_group = which(data$enzyme=="" | is.na(data$enzyme))[1]
    raw.files = data[1:(idx_group-1), ]
    groups = data[1:(idx_group-1), ]
    total = data[nrow(data), ]
    return (list(raw.files, groups, total))
  }
  
  return (data);
}


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



