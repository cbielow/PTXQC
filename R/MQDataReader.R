
## A proto class for handling consistent Raw file names while loading
## multiple MQ result files.
## If the names are too long, an alias name (eg 'file 1', 'file 2', ...) is used instead.
##
##
## [Occasional rage: since S4 is so very inadequate for basically everything which is important in OOP and the syntax is
##                   even more horrible, we use the 'proto' package to at least abstract away much of this S4 non-sense.
##            Read http://cran.r-project.org/doc/contrib/Genolini-S4tutorialV0-5en.pdf::10:2 Method to modify a field and you'll see
##            what I mean]
##

#'
#' Convenience wrapper for MQDataReader when only a single MQ file should be read
#' and file mapping need not be stored.
#' 
#' For params, see \code{\link{MQDataReader$readMQ}}.
#' 
#' @param file   see \code{\link{MQDataReader$readMQ}}
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
read.MQ <- function(file, filter = "", type = "pg", col_subset = NA, add_fs_col = 10, LFQ_action = FALSE, ...)
{
  mq = MQDataReader$new()
  mq$readMQ(file, filter, type, col_subset, add_fs_col, LFQ_action, ...)
}


## CLASS 'MQDataReader'
MQDataReader <- proto::proto()

#' Constructor for class 'MQDataReader'.
#'
#' This class is used to read MQ data tables using readMQ() while holding
#' the internal raw file --> short raw file name mapping (stored in a member called 
#' 'fn_map') and updating/using it every time readMQ() is called.
#' 
#' @name MQDataReader$new
#' @import proto
#' 
MQDataReader$new <- function(.)
{
  proto::proto(., mq.data = NULL, fn_map = FilenameMapper$new())
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
#'  alternative term          new term
#'  -----------------------------------------
#'  protease                  enzyme
#'  protein.descriptions      fasta.headers
#'  potential.contaminant     contaminant
#'  mass.deviations           mass.deviations..da.
#'  basepeak.intensity        base.peak.intensity
#' }
#' 
#' We also correct 'reporter.intensity.*' naming issues to MQ 1.6 convention, when 'reporter.intensity.not.corrected' is present.
#' MQ 1.5 uses: reporter.intensity.X and reporter.intensity.not.corrected.X
#' MQ 1.6 uses: reporter.intensity.X and reporter.intensity.corrected.X       
#'
#' Note: you must find a regex which matches both versions, or explicitly add both terms if you are requesting only a subset
#'       of columns!
#' 
#' Example of usage:
#' \preformatted{
#'   mq = MQDataReader$new()
#'   d_evd = mq$readMQ("evidence.txt", type="ev", filter="R", col_subset=c("proteins", "Retention.Length", "retention.time.calibration")) 
#' }
#' 
#' If the file is empty, this function shows a warning and returns NULL.
#' If the file is present but cannot be read, the program will stop.
#'
#' @param .      A 'this' pointer. Use it to refer/change internal members. It's implicitly added, thus not required too call the function!
#' @param file   (Relative) path to a MQ txt file.
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
#'                   a running ID of the form 'file <x>', where <x> is a number from 1 to N.
#'                   If the function is called again and a mapping already exists, this mapping is used.
#'                   Should some raw.files be unknown (ie the mapping from the previous file is incomplete), they will be augmented
#' @param check_invalid_lines After reading the data, check for unusual number of NA's to detect if file was corrupted by Excel or alike                 
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
#'
# (not exported!)
MQDataReader$readMQ <- function(., file, filter="", type="pg", col_subset=NA, add_fs_col=10, check_invalid_lines = TRUE, LFQ_action=FALSE, ...)
{
  # . = MQDataReader$new() ## debug
  # ... = NULL

  if (!file.exists(file)) {
    cat(paste0("MaxQuant file ", file, " was not found. Reading skipped.\n"))
    return (NULL);
  }
  cat(paste("Reading file", file,"...\n"))
  ## error message if failure should occur below
  msg_parse_error = paste0("\n\nParsing the file '", file, "' failed. See message above why. If the file is not usable but other files are ok, disable the corresponding section in the YAML config. You might also be running a foreign locale (e.g. Chinese) - switch to an English locale and make sure that txt files are encoded in ASCII (Latin-1)!")
  
  ## resolve set of columns which we want to keep
  #example: col_subset = c("Multi.*", "^Peaks$")
  colClasses = NA ## read.table: keep all columns by default
  if (sum(!is.na(col_subset)) > 0)
  { ## just read a tiny bit to get column names
    ## do not use data.table::fread for this, since it will read the WHOLE file and takes ages...
    data_header = try(read.delim(file, comment.char="", nrows=2))
    if (inherits(data_header, 'try-error')) stop(msg_parse_error, call. = FALSE);
    
    colnames(data_header) = tolower(colnames(data_header))
    colClasses = rep("NULL", ncol(data_header)) ## either "NULL" (discard), NA (auto-detect), or specific, i.e. "numeric", ...
    for (idx in 1:length(col_subset))
    {
      idx_new = grepl(col_subset[idx], colnames(data_header), ignore.case = TRUE)
      if (sum(idx_new) == 0)
      {
        cat(paste0("WARNING: Could not find column regex '", col_subset[idx], "' using case-INsensitive matching.\n"))
      } else {
        colClasses[idx_new] = NA  ## auto
        if (nchar(names(col_subset)[idx]) > 0) {
          colClasses[idx_new] = names(col_subset)[idx]
          cat(paste0("Requiring column(s) '", paste(colnames(data_header)[idx_new], sep="", collapse="', '"), "' to be of type '", names(col_subset)[idx], "'!\n"))
        } 
      }
      
    }
    ## keep the 'id' column if available (for checking data integrity: invalid line-breaks)
    if ("id" %in% colnames(data_header) & check_invalid_lines == TRUE) colClasses[which("id"==colnames(data_header))] = NA
    cat(paste("Keeping", sum(colClasses != "NULL", na.rm=TRUE), "of", ncol(data_header), "columns!\n"))
    if (sum(colClasses != "NULL", na.rm=TRUE) == 0) {
      ## can happen for very old MQ files without header, or if the user just gave the wrong colClasses
      .$mq.data = data.frame()
      return (.$mq.data)
    }
  }
  
  ## higher memory consumption during load (due to memory mapped files) compared to read.delim... but about 5x faster
  ## , but also different numerical results when parsing numbers!!!
  #.$mq.data = try(
  #  fread(file, header = TRUE, sep='\t', na.strings=c("NA", "n. def."), verbose = TRUE, select = idx_keep, data.table = FALSE, ...)
  #)
  #colnames(.$mq.data) = make.names(colnames(.$mq.data), unique = TRUE)

  ## comment.char should be "", since lines will be TRUNCATED starting at the comment char.. and a protein identifier might contain just anything...
  ## na.strings:
  ##  - use BOTH 'n. def.' and 'n.def.' even though only the former is present in the data
  ##    However, when the colClass is 'numeric', whitespaces are stripped, and only AFTERWARDS the string
  ##    is checked against na.strings
  ##  - the '\u975E\u6570\u5B57' na-string is the chinese UTF-8 representation of "NA"
  .$mq.data = try(read.delim(file, na.strings=c("NA", "n. def.", "n.def.", "\u975E\u6570\u5B57"), encoding="UTF-8", comment.char="", stringsAsFactors = FALSE, colClasses = colClasses, ...))
  if (inherits(.$mq.data, 'try-error')) stop(msg_parse_error, call. = FALSE);
  
  #colnames(.$mq.data)
  
  cat(paste0("Read ", nrow(.$mq.data), " entries from ", file,".\n"))

  ### checking for invalid rows
  if (check_invalid_lines == TRUE & type != "sm") ## summary.txt has irregular structure
  {
    inv_lines = .$getInvalidLines();
    if (length(inv_lines) > 0)
    {
      stop(paste0("\n\nError: file '", file, "' seems to have been edited in Microsoft Excel and",
                                             " has artificial line-breaks which destroy the data at lines (roughly):\n",
                                             paste(inv_lines, collapse="\n"), "\nPlease fix (e.g. try LibreOffice 4.0.x or above)!",
                                             "You might also have unknown UTF-8 characters (Asian?) in your file which PTXQC does not recognize. Try running MaxQuant on an english locale!"))
    }
  }
  
  cat(paste0("Updating colnames\n"))
  cn = colnames(.$mq.data)
  ### just make everything lower.case (MQ versions keep changing it and we want it to be reproducible)
  cn = tolower(cn)
  ## rename some columns since MQ 1.2 vs. 1.3 differ....
  cn[cn=="protease"] = "enzyme"
  cn[cn=="protein.descriptions"] = "fasta.headers"
  cn[cn=="basepeak.intensity"] = "base.peak.intensity"
  ## MQ 1.0.13 and 1.1.x have 'mass.deviations', later versions have 'mass.deviations..da.'
  cn[cn=="mass.deviations"] = "mass.deviations..da."
  ## MQ 1.5 uses 'potential.contaminant' instead of 'contaminant'
  cn[cn=="potential.contaminant"] = "contaminant"
  ## MQ 1.5 uses     : reporter.intensity.X and reporter.intensity.not.corrected.X
  ## whereas 1.6 has : reporter.intensity.X and reporter.intensity.corrected.X
  ## we rename to 1.6 convention:
  if (sum(grepl("reporter.intensity.not.corrected", cn)) > 0)
  {
    idx_nc = grep("reporter.intensity.not.corrected", cn);
    idx_c = grep("reporter.intensity.[0-9]", cn);
    cn[idx_nc] = gsub(".not.corrected", "", cn[idx_nc]);
    cn[idx_c] = gsub(".intensity", ".intensity.corrected", cn[idx_c]);
  }
  
  colnames(.$mq.data) = cn
  
  
  
  ## work in-place on 'contaminant' column
  cat(paste0("Simplifying contaminants\n"))
  .$substitute("contaminant");
  cat(paste0("Simplifying reverse\n"))
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
        if (sum(bad_rows, na.rm = TRUE)==0) {next;}
        ## take action
        if (LFQ_action=="toNA" | LFQ_action=="impute") {
          ## set to NA
          impVal = NA;
          vals[bad_rows] = impVal;
          if (LFQ_action=="impute") {
            impVal = min(vals[vals>0], na.rm = TRUE)
            ## replace with minimum noise value (>0!)
            vals[bad_rows] = impVal;
          }
          cat(paste0("   '", cc, "' ", sum(bad_rows, na.rm = TRUE), ' entries (', sum(bad_rows, na.rm = TRUE)/nrow(.$mq.data)*100,'%) with ', impVal, '\n'))
          ## add column
          .$mq.data[, paste0("c", cc)] = vals;
          ##
          stats = rbind(stats, c(sum(bad_rows, na.rm = TRUE), impVal))
          
        }
        else {
          stop(paste0("Unknown action '", LFQ_action, "' for LFQ_action parameter! Aborting!"));
        }
      }
      if (LFQ_action=="impute") {
        plot(stats, log="y", xlab="number of replaced zeros", ylab='imputed intensity', main='Imputation statistics')
      }
    }
    
    int_cols = grepv("intensity", colnames(.$mq.data))
    
    ##
    ## apply potential fix (MQ 1.5 writes numbers in scientific notation with ',' -- which parses as String :( )
    ##
    int_cols_nn = (apply(.$mq.data[,int_cols, drop = FALSE], 2, class) != "numeric")
    if (any(int_cols_nn))
    {
      .$mq.data[, int_cols[int_cols_nn]] = sapply(int_cols[int_cols_nn], function(x_name)
      {
        x = .$mq.data[, x_name]
        if (class(x) == "integer")
        {
          x = as.numeric(x)
        }
        if (class(x) != "numeric")
        {
          cat(paste(x[grep(",", x)], collapse = ",", sep=", "))
          cat(paste(unique(x), collapse = ",", sep=", "))
          warning(paste0("Column '", x_name, "' in file '", file, "' is expected to contain numbers, but was classified as '",
                         class(x), "'.",
                         " A bug in MaxQuant 1.5.2.8 (and related?) versions causes usage of comma instead of dot for large numbers in scientific notation,",
                         " e.g. '1,73E+011' instead of '1.73E+011'. We will try to fix this now ..."),
                  call. = FALSE, immediate. = TRUE)
          x = as.numeric(gsub(",",".", x))
        }
        return(x)
      })
    }
    
    ##
    ## add Abundance index
    ##
    if ("mol..weight..kda." %in% colnames(.$mq.data)){
      ### add abundance index columns (for both, intensity and lfq.intensity)
      .$mq.data[, sub("intensity", "AbInd", int_cols)] = apply(.$mq.data[,int_cols, drop = FALSE], 2, function(x)
      {
        return (x / .$mq.data[,"mol..weight..kda."])
      })
    } else {
      stop("MQDataReader$readMQ(): Cannot add abundance index since 'mol..weight..kda.' was not loaded from file. Did you use the correct 'type' or forgot to add the column in 'col_subset'?")
    }
    
  } else if (type=="sm") {
    ## summary.txt special treatment
    ## find the first row, which lists Groups (after Raw files): it has two non-zero entries only 
    ##                                                           (or even less if the group name is empty)
    ##dx <<- .$mq.data;
    idx_group = which(apply(.$mq.data, 1, function(x) sum(x!="", na.rm = TRUE))<=2)[1]
    ## summary.txt will not contain groups, if none where specified during MQ-configuration
    if (is.na(idx_group)) {
      idx_group = nrow(.$mq.data)
      groups= NA
    } else {
      groups = .$mq.data[idx_group:(nrow(.$mq.data)-1), ]
    }
    raw.files = .$mq.data[1:(idx_group-1), ]
    total = .$mq.data
    .$mq.data = raw.files ## temporary, until we have assigned the fc.raw.files
  }
  
  
  if (add_fs_col & "raw.file" %in% colnames(.$mq.data))
  {
    .$mq.data$fc.raw.file = .$fn_map$getShortNames(.$mq.data$raw.file, add_fs_col)
    
  }

  if (type=="sm") { ## post processing for summary
    ## .$mq.data is basically "raw.files#, but with fc.raw.files augmented
    .$mq.data = list(raw = .$mq.data, groups = groups, total = total)
  }
  
  return (.$mq.data);
} ## end readMQ()


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
#' If the table has an 'id' column, we can simply check the numbers are consecutive. If no 'id' column is available,
#' we detect line-breaks by counting the number of NA's per row and finding outliers.
#' The line break then must be in this line (plus the preceeding or following one). Depending on where
#' the break happened we can also detect both lines right away (if both have more NA's than expected).
#'
#' Currently, we have no good strategy to fix the problem since columns are not aligned any longer, which
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
    stop("In 'MQDataReader$getInvalidLines': function called before data was loaded. Internal error. Exiting.", call. = FALSE);
  }
  broken_rows = c()
  if ("id" %in% colnames(.$mq.data))
  {
    last_id = as.numeric(as.character(.$mq.data$id[nrow(.$mq.data)]))
    if (is.na(last_id) || (last_id+1)!=nrow(.$mq.data))
    {
      print(paste0("While checking ID column: last ID was '", last_id, "', while table has '", nrow(.$mq.data), "' rows."))
      broken_rows = which(!is.numeric(as.character(.$mq.data$id)))
    }
  } else
  {
    cols = !grepl("ratio", colnames(.$mq.data)) ## exclude ratio columns, since these can have regular NA's in unpredictable frequency
    counts = apply(.$mq.data[, cols], 1, function(x) sum(is.na(x)));
    ## NA counts should be roughly equal across rows
    expected_count = quantile(counts, probs = 0.75)
    broken_rows = which(counts > (expected_count * 3 + 10))
    if (length(broken_rows) > 0)
    {
      print("Table:")
      print(table(counts))
      print(paste0("NAn count limit: 3*", expected_count, " + 10 = ", expected_count * 3 + 10))    
    }
  }
  
  return (broken_rows);
}
