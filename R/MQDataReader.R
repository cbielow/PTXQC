
## A proto class for handling consistent Raw file names while loading
## multiple MQ result files.
## If the names are too long, an alias name (eg 'f1', 'f2', ...) is used instead.
##
##
## [Occasional rage: since S4 is so very inadequate for basically everything which is important in OOP and the syntax is
##                   even more horrible, we use the 'proto' package to at least abstract away much of this S4 non-sense.
##            Read http://cran.r-project.org/doc/contrib/Genolini-S4tutorialV0-5en.pdf::10:2 Method to modify a field and you'll see
##            what I mean]
##

#install.packages("proto")
#require(proto)

## CLASS 'MQDataReader'
MQDataReader <- proto()

#' Constructor for class 'MQDataReader'.
#'
#' This class is used to read MQ data tables using readMQ() while holding
#' the internal raw file --> short raw file name mapping (stored in a member called 
#' 'raw_file_mapping') and updating/using it every time readMQ() is called.
#' 
#' @name MQDataReader$new
#' @import proto
#' 
MQDataReader$new <- function(.)
{
  proto(., raw_file_mapping = data.frame())
}

##
## Functions
##

#' Wrapper to read a MQ txt file (e.g. proteinGroups.txt).
#' 
#' Since MaxQuant changes capitalization and sometimes even column names, it seemed convenient
#' to have a function which just reads a txt file and returns unified column names, irrespective of the MQ version.
#' So, it unifies access to columns (e.g. by using lower case for ALL columns) and ensures columns are
#' identically named across MQ versions. We only of one case: "protease" == "enzyme".
#'
#' If the file is empty, this function stops with an error.
#'
#' @param .      A 'this' pointer (in C++ terms). Use it to refer/change internal members. It's implicitly added, thus not required too call the function!
#' @param file   (Relative) path to a MQ txt file ()
#' @param filter Searched for "C" and "R". If present, [c]ontaminants and [r]everse hits are removed
#' @param type   Allowed values are:
#'               "pg" (proteinGroups) [default], adds abundance index columns (*AbInd*, replacing 'intensity')
#'               "sm" (summary), splits into three row subsets (raw.file, condition, total)
#'               Any other value will not add any special columns
#' @param col_subset A vector of column names as read by read.delim(), e.g., spaces are replaced by dot already.
#'                   If given, only columns with these names (ignoring lower/uppercase) will be returned (regex allowed)
#'                   E.g. col_subset=c("^lfq.intensity.", "protein.name")
#' @param add_fs_col If TRUE and a column 'raw.file' is present, an additional column 'fc.raw.file' will be added with 
#'                   common prefix AND common substrings removed (\code{\link{simplifyNames}})
#'                           E.g. two rawfiles named 'OrbiXL_2014_Hek293_Control', 'OrbiXL_2014_Hek293_Treated' will give
#'                                                   'Control', 'Treated'
#'                   If 'add_fs_col' is a number AND the longest short-name is still longer, the names are discarded and replaced by
#'                   a running ID of the form 'f<x>', where <x> is a number from 1 to N.
#'                   If the function is called again and a mapping already exists, this mapping is used.
#'                   Should some raw.files be unknown (ie the mapping from the previous file is incomplete), an error is thrown.
#' @param LFQ_action [For type=='pg' only] An additional custom LFQ column ('cLFQ...') is created where
#'               zero values in LFQ columns are replaced by the following method IFF(!) the corresponding raw intensity is >0 (indicating that LFQ is erroneusly 0)
#'               "toNA": replace by NA
#'               "impute": replace by lowest LFQ value >0 (simulating 'noise')
#' @param ... Additional parameters passed on to read.delim()             
#' @return A data.frame of the respective file
#' 
#' @name MQDataReader$readMQ
#' @import utils
#' @import graphics
#' 
MQDataReader$readMQ <- function(., file, filter="C+R", type="pg", col_subset=NA, add_fs_col=10, LFQ_action=FALSE, ...)
{
  cat(paste("Reading file", file,"...\n"))
  ## error message if failure should occur below
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
    int_cols = grepv("intensity", colnames(data))
    data[, sub("intensity", "AbInd", int_cols)] = apply(data[,int_cols, drop=F], 2, function(x)
    {
      x / data[,"mol..weight..kda."]
    })
    
  } else if (type=="sm")
    ## summary.txt special treatment
  { ## split
    idx_group = which(data$enzyme=="" | is.na(data$enzyme))[1]
    raw.files = data[1:(idx_group-1), ]
    groups = data[idx_group:(nrow(data)-1), ]
    total = data
    data = raw.files ## temporary, until we have assigned the fc.raw.files
  }
  
  
  if (add_fs_col & "raw.file" %in% colnames(data))
  {
    ## check if we already have a mapping
    if (length(.$raw_file_mapping) > 0)
    {
     
    } else {
      
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
      .$raw_file_mapping = data.frame(from = rf_name, to = rf_name_s, stringsAsFactors = FALSE)
      ## check if the minimal length was reached
      add_fs_col = 10
      if (is.numeric(add_fs_col) & max(nchar((.$raw_file_mapping$to))) > add_fs_col)
      { ## resort to short naming convention
        .$raw_file_mapping[, "best effort"] = .$raw_file_mapping$to
        cat("Filenames are longer than the maximal allowed size of '" %+% add_fs_col %+% "'. Resorting to short versions 'f...'.\n\n")
        maxl = length(unique(data$raw.file))
        .$raw_file_mapping$to = paste("f", sprintf(paste0("%0", nchar(maxl), "d"), 1:maxl)) ## with leading 0's if required
        ## plot a mapping table ... (this might not be the best place to do it though..)
        .$plotNameMapping()
      }
      
    }
    ## do the mapping    
    data$fc.raw.file = as.factor(.$raw_file_mapping$to[match(data$raw.file, .$raw_file_mapping$from)])
    ## check for NA's
    if (any(is.na(data$fc.raw.file)))
    {
      missing = unique(data$raw.file[is.na(data$fc.raw.file)])
      stop("Generation of short Raw file names failed due to missing mapping entries:" %+% paste(missing, collapse=", ", sep=""))
    }
  }

  if (type=="sm") { ## post processing for summary
    data = list(raw = raw.files, groups = groups, total = total)
  }
  
  return (data);
} ## end readMQ()


#' Plots the current mapping of Raw file names to their shortened version.
#'
#' Convenience function to plot the mapping (e.g. to a PDF device for reporting).
#' The data frame can be accessed directly via '.$raw_file_mapping'.
#' If no mapping exists, the function prints a warning to console and omits the plot.
#'
#' @return Returns 'TRUE' if mapping is present. 'FALSE' otherwise (no plot is generated).
#'
#' @name MQDataReader$plotNameMapping
#' @import plotrix
#' 
MQDataReader$plotNameMapping <- function(.)
{
  if (length(.$raw_file_mapping) > 0)
  {
    extra = ""
    if ("best effort" %in% colnames(.$raw_file_mapping))
    {
      extra = "\n(automatic shortening of names was not sufficiently short - see 'best effort'"
    }
    plot(0:100, 0:100, type="n", axes=F, xlab="", ylab="")
    #tbl_cex = ifelse(bad_id_count<25, 1, 25/bad_id_count)
    addtable2plot("topleft", y=NULL, 
                  .$raw_file_mapping,
                  xjust=0, yjust=1,
                  xpad=0, ypad=0.1,
                  title="Info: mapping of raw files to their short names" %+% extra) 
  } else {
    cat("No mapping found. Omitting plot.")
    return (FALSE);
  }
    
  return (TRUE);  
}
