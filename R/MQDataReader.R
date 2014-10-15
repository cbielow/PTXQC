
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

#'
#' Convenience wrapper for MQDataReader when only a single MQ file should be read
#' and file mapping need not be stored
#' 
#' For params, see \code{\link{MQDataReader$readMQ}}.
#' 
#' #' @param file   (Relative) path to a MQ txt file ()
#' @param filter see \code{\link{MQDataReader$readMQ}}
#' @param type   see \code{\link{MQDataReader$readMQ}}
#' @param col_subset see \code{\link{MQDataReader$readMQ}}
#' @param add_fs_col see \code{\link{MQDataReader$readMQ}}
#' @param LFQ_action see \code{\link{MQDataReader$readMQ}}
#' @param ... see \code{\link{MQDataReader$readMQ}}
#' @return see \code{\link{MQDataReader$readMQ}}
#'
#' @export
#' 
read.MQ <- function(file, filter="", type="pg", col_subset=NA, add_fs_col=10, LFQ_action=FALSE, ...)
{
  mq = MQDataReader$new()
  mq$readMQ(file, filter, type, col_subset, add_fs_col, LFQ_action, ...)
}


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
  proto(., raw_file_mapping = NULL, mq.data = NULL, mapping.created = FALSE)
}


##
## Functions
##

#' Wrapper to read a MQ txt file (e.g. proteinGroups.txt).
#' 
#' Since MaxQuant changes capitalization and sometimes even column names, it seemed convenient
#' to have a function which just reads a txt file and returns unified column names, irrespective of the MQ version.
#' So, it unifies access to columns (e.g. by using lower case for ALL columns) and ensures columns are
#' identically named across MQ versions:
#' \preformatted{
#'  old term                  new term
#'  -----------------------------------------
#'  protease                  enzyme
#'  protein.descriptions      fasta.headers
#' }
#' 
#' 
#' Example of usage:
#' \preformatted{
#'   mq = MQDataReader$new()
#'   d_evd = mq$readMQ("evidence.txt", type="ev", filter="R", col_subset=c("proteins", "Retention.Length", "retention.time.calibration")) 
#' }
#' 
#' If the file is empty, this function stops with an error.
#'
#' @param .      A 'this' pointer. Use it to refer/change internal members. It's implicitly added, thus not required too call the function!
#' @param file   (Relative) path to a MQ txt file ()
#' @param filter Searched for "C" and "R". If present, [c]ontaminants and [r]everse hits are removed if the respective columns are present.
#'               E.g. to filter both, \code{filter = "C+R"}
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
#'                   If \code{add_fs_col} is a number AND the longest short-name is still longer, the names are discarded and replaced by
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
# (not exported!)
MQDataReader$readMQ <- function(., file, filter="", type="pg", col_subset=NA, add_fs_col=10, LFQ_action=FALSE, ...)
{
  ## it's either present already or will be created
  .$mapping.created = FALSE
  
  # . = MQDataReader$new() ## debug
  # ... = NULL
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
  
  .$mq.data = try(read.delim(file, na.strings=c("NA", "n. def."), stringsAsFactors=F, comment.char="#", colClasses = col_subset, ...))
  if (inherits(.$mq.data, 'try-error')) stop(msg_parse_error, call.=F);
  
  #colnames(.$mq.data)
  
  cat(paste0("Read ", nrow(.$mq.data), " entries from ", file,".\n"))
  
  ### checking for invalid rows
  if (type != "sm") ## summary.txt has irregular structure
  {
    inv_lines = .$getInvalidLines();
    if (length(inv_lines) > 0)
    {
      stop(paste0("\n\nError: file '", file, "' seems to have been edited in Microsoft Excel and",
                                             " has artificial line-breaks which destroy the data at lines (roughly):\n",
                                             paste(inv_lines, collapse="\n"), "\nPlease fix (e.g. try LibreOffice 4.0.x or above)!"))
    }
  }
  
  ### just make everything lower.case (MQ versions keep changing it and we want it to be reproducible)
  colnames(.$mq.data) = tolower(colnames(.$mq.data))
  ## rename some columns since MQ 1.2 vs. 1.3 differ....
  colnames(.$mq.data)[colnames(.$mq.data)=="protease"] = "enzyme"
  colnames(.$mq.data)[colnames(.$mq.data)=="protein.descriptions"] = "fasta.headers"
        
  ## work in-place on 'contaminant' column
  .$substitute("contaminant");
  .$substitute("reverse");
  if (grepl("C", filter) & ("contaminant" %in% colnames(.$mq.data))) .$mq.data = .$mq.data[!(.$mq.data$contaminant),]
  if (grepl("R", filter) & ("reverse" %in% colnames(.$mq.data))) .$mq.data = .$mq.data[!(.$mq.data$reverse),]
  
  ## proteingroups.txt special treatment
  if (type=="pg") {
    stats = data.frame(n=NA, v=NA)
    if (LFQ_action!=FALSE)
    { ## replace erroneous zero LFQ values with something else
      cat("Starting LFQ action.\nReplacing ...\n")
      lfq_cols = grepv("^lfq", colnames(.$mq.data))
      for (cc in lfq_cols)
      {
        ## get corresponding raw intensity column
        rawint_col = sub("^lfq\\.", "", cc)
        if (!(rawint_col %in% colnames(.$mq.data))) {stop(paste0("Could not find column '", rawint_col, "' in dataframe with columns: ", paste(colnames(.$mq.data), collapse=",")), "\n")}
        vals = .$mq.data[, cc]
        ## affected rows
        bad_rows = (.$mq.data[, rawint_col]>0 & .$mq.data[, cc]==0)
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
          cat(paste0("   '", cc, "' ", sum(bad_rows, na.rm=T), ' entries (', sum(bad_rows, na.rm=T)/nrow(.$mq.data)*100,'%) with ', impVal, '\n'))
          ## add column
          .$mq.data[, paste0("c", cc)] = vals;
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
    int_cols = grepv("intensity", colnames(.$mq.data))
    .$mq.data[, sub("intensity", "AbInd", int_cols)] = apply(.$mq.data[,int_cols, drop=F], 2, function(x)
    {
      x / .$mq.data[,"mol..weight..kda."]
    })
    
  } else if (type=="sm") {
    ## summary.txt special treatment
    ## split by a column which has values for raw.files but not for groups (e.g. multiplicity, NOT instrument (missing in 1.2), NOT enzyme since it can be all NA if disabled)
    idx_group = which(.$mq.data$multiplicity=="" | is.na(.$mq.data$multiplicity))[1]
    raw.files = .$mq.data[1:(idx_group-1), ]
    groups = .$mq.data[idx_group:(nrow(.$mq.data)-1), ]
    total = .$mq.data
    .$mq.data = raw.files ## temporary, until we have assigned the fc.raw.files
  }
  
  
  if (add_fs_col & "raw.file" %in% colnames(.$mq.data))
  {
    ## check if we already have a mapping
    if (is.null(.$raw_file_mapping))
    {
      rf_name = unique(.$mq.data$raw.file)
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
        maxl = length(unique(.$mq.data$raw.file))
        .$raw_file_mapping$to = paste("f", sprintf(paste0("%0", nchar(maxl), "d"), 1:maxl)) ## with leading 0's if required
        ## indicate to outside that a new table is ready
        .$mapping.created = FALSE
      }
      
    }
    ## do the mapping    
    .$mq.data$fc.raw.file = as.factor(.$raw_file_mapping$to[match(.$mq.data$raw.file, .$raw_file_mapping$from)])
    ## check for NA's
    if (any(is.na(.$mq.data$fc.raw.file)))
    {
      missing = unique(.$mq.data$raw.file[is.na(.$mq.data$fc.raw.file)])
      stop("Generation of short Raw file names failed due to missing mapping entries:" %+% paste(missing, collapse=", ", sep=""))
    }
  }

  if (type=="sm") { ## post processing for summary
    .$mq.data = list(raw = raw.files, groups = groups, total = total)
  }
  
  return (.$mq.data);
} ## end readMQ()


#' Plots the current mapping of Raw file names to their shortened version.
#'
#' Convenience function to plot the mapping (e.g. to a PDF device for reporting).
#' The data frame can be accessed directly via \code{.$raw_file_mapping}.
#' If no mapping exists, the function prints a warning to console and returns NULL (which is safe to use in print(NULL)).
#'
#' @return Returns a mapping plot if mapping is available, 'NULL' otherwise.
#'
#' @name MQDataReader$plotNameMapping
#' 
MQDataReader$plotNameMapping <- function(.)
{
  if (!is.null(.$raw_file_mapping))
  {
    extra = ""
    if ("best effort" %in% colnames(.$raw_file_mapping))
    {
      extra = "\n(automatic shortening of names was not sufficiently short - see 'best effort')"
    }
    
    #mq_mapping = mq$raw_file_mapping
    mq_mapping = .$raw_file_mapping
    mq_mapping$ypos = -(1:nrow(mq_mapping))
    head(mq_mapping)
    mq_mapping.long = melt(mq_mapping, id.vars = c("ypos"), measure.vars=c("from","to","best effort"))
    head(mq_mapping.long)
    mq_mapping.long$variable = as.character(mq_mapping.long$variable)
    mq_mapping.long$col = "#000000";
    mq_mapping.long$col[mq_mapping.long$variable=="to"] = "#5F0000"
    mq_mapping.long$variable[mq_mapping.long$variable=="to"] = 9
    mq_mapping.long$variable[mq_mapping.long$variable=="from"] = 8
    mq_mapping.long$variable[mq_mapping.long$variable=="best effort"] = 12
    mq_mapping.long$variable = as.numeric(mq_mapping.long$variable)
    mq_mapping.long$size = 2;
    df.header = data.frame(ypos = 0, variable = c(8,9,12), value = c("original", "short\nname", "best_effort"), col = "#000000", size=3)
    mq_mapping.long2 = rbind(mq_mapping.long, df.header)
    mqmap_pl = ggplot(mq_mapping.long2, aes_string(x = "variable", y = "ypos"))  +
      geom_text(aes_string(label="value"), color = mq_mapping.long2$col, hjust=1, size=mq_mapping.long2$size) +
      coord_cartesian(xlim=c(0,15)) +
      theme_bw() +
      theme(plot.margin = unit(c(1,1,1,1), "cm"), line = element_blank(), axis.title = element_blank(), panel.border = element_blank(),
            axis.text = element_blank(), strip.text = element_blank(), legend.position = "none") +
      ggtitle("Info: mapping of raw files to their short names" %+% extra)
    return(mqmap_pl)
  } else {
    cat("No mapping found. Omitting plot.")
    return (NULL);
  }
    
  return (NULL);  
}

#' Replaces values in the mq.data member with (binary) values.
#'
#' Most MQ tables contain columns like 'contaminants' or 'reverse', whose values are either empty strings
#' or "+", which is inconvenient and can be much better represented as TRUE/FALSE.
#' The params \code{valid_entries} and \code{replacements} contain the matched pairs, which determine what is replaced with what.
#' 
#' @param colname       Name of the column (e.g. "contaminants") in the mq.data table
#' @param valid_entries Vector of values to be replaced (must contain all values expected in the column -- fails otherwise)
#' @param replacements  Vector of values inserted with the same length as \code{valid_entries}.
#' @return Returns \code{TRUE} if successful.
#'
#' @name MQDataReader$substitute
#' 
MQDataReader$substitute <- function(., colname, valid_entries = c(NA, "","+"), replacements = c(FALSE, FALSE, TRUE))
{
  if (length(valid_entries) == 0)
  {
    stop("Entries given to $substitute() must not be empty.")
  }
  if (length(valid_entries) != length(replacements))
  {
    stop("In function $substitute(): 'valid_entries' and 'replacements' to not have the same length!")
  }
  if (colname %in% colnames(.$mq.data))
  {
    ## verify that there are only known entries (usually c("","+") )
    setD_c = setdiff(.$mq.data[, colname], valid_entries)
    if (length(setD_c) > 0) stop(paste0("'", colname, "' column contains unknown entry (", paste(setD_c, collapse=",", sep="") ,")."))
    ## replace with TRUE/FALSE
    .$mq.data[, colname] = replacements[ match(.$mq.data[, colname], valid_entries) ];
  }
  return (TRUE);
}


#' Detect broken lines (e.g. due to Excel import+export)
#'
#' When editing a MQ txt file in Microsoft Excel, saving the file can cause it to be corrupted,
#' since Excel has a single cell content limit of 32k characters 
#' (see http://office.microsoft.com/en-001/excel-help/excel-specifications-and-limits-HP010342495.aspx)
#' while MQ can easily reach 60k (e.g. in oxidation sites column).
#' Thus, affected cells will trigger a line break, effectively splitting one line into two (or more).
#' 
#' We detect this by counting the number of NA's per row and finding outliers.
#' The line break then must be in this line (plus the preceeding or following one). Depending on where
#' the break happened we can also detect both lines right away (if both have more NA's than expected).
#'
#' Currently, we have no good strategy to fix this problem since columns are not aligned any longer, which
#' leads to columns not having the class (e.g. numeric) they should have.
#' (thus one would need to un-do the linebreak and read the whole file again)
#' 
#' [Solution to the problem: try LibreOffice 4.0.x or above -- seems not to have this limitation]
#' 
#' @return Returns a vector of indices of broken (i.e. invalid) lines
#'
#' @name MQDataReader$getInvalidLines
#' 
MQDataReader$getInvalidLines <- function(.)
{
  if (!inherits(.$mq.data, 'data.frame'))
  {
    stop("In 'MQDataReader$getInvalidLines': function called before data was loaded. Internal error. Exiting.", call.=F);
  }
  
  counts = apply(.$mq.data, 1, function(x) sum(is.na(x)));
  ## NA counts should be roughly equal across rows
  expected_count = quantile(counts, probs = 0.75)
  broken_rows = which(counts > (expected_count * 3 + 10))
  if (length(broken_rows) > 0)
  {
    print("Table:")
    print(table(counts))
    print(paste0("NAn count limit: 3*", expected_count, " + 10 = ", expected_count * 3 + 10))    
  }
  
  return (broken_rows);
}
