#'
#'  Make sure to call $readMappingFile(some_file) if you want to support a user-defined file mapping.
#'  Otherwise, calls to $getShortNames() will create/augment the mapping for filenames.
#' 
#'
#'
#' @field raw_file_mapping Data.frame with columns 'from', 'to' and maybe 'best.effort' (if shorting was unsuccessful)
#' @field mapping.creation how the current mapping was obtained (user or auto)
#' @field external.mapping.file Filename of user-defined mapping file; only defined if readMappingFile() was called
#'
#' @import methods
#' @import ggplot2
#' 
#' @exportClass FilenameMapper
#' @export FilenameMapper
#' 
#' @examples 
#' a = FilenameMapper$new()
#' a$readMappingFile('filenamemapping.txt') 
#' 
FilenameMapper = setRefClass("FilenameMapper",
                       
                       fields = list(raw_file_mapping = "data.frame", ## with cols 'from', to' and maybe 'best.effort' (if shorting was unsuccessful)
                                     mapping.creation = "character",  ## how the current mapping was obtained (user or auto)
                                     external.mapping.file = "character" ## filename of user-defined mapping file; only defined if readMappingFile() was called
                       ),
                       methods = list(

initialize=function() {
 .self$raw_file_mapping = data.frame()
 .self$mapping.creation = NA_character_
 .self$external.mapping.file = NA_character_
 return(.self)
},

mapRunsToShort = function(.self, ms_runs)
{
  "Given a vector of ms_runs, return a data.frame of identical length with columns 'raw.file' and 'fc.raw.file'."
  if (!"ms.run" %in% colnames(.self$raw_file_mapping)) stop("Mapping is missing 'ms.run' from mzTab!")
  
  res = .self$raw_file_mapping[ match(ms_runs, .self$raw_file_mapping$ms.run), c("from", "to")]
  colnames(res) = c("raw.file", "fc.raw.file")
  return (res)
},

getShortNames = function(.self, raw.filenames, max_length = 10, ms_runs = NULL)
{
  "Uses the internal mapping (or augments it if possible) and maps the input raw names to shorter output names. 
    Returns a vector of the same length."
  #rf <<- raw.filenames
  #raw.filenames = rf

  if (!is.null(ms_runs) && length(ms_runs) != length(raw.filenames)) stop("raw.filenames and ms_runs do not have the same length!")
  
  cat(paste0("Adding fc.raw.file column ..."))
  ## check if we already have a mapping
  if (nrow(.self$raw_file_mapping) == 0)
  {
    rfm = .self$getShortNamesStatic(unique(raw.filenames), max_length)
    if (!is.null(ms_runs)) {
      rfm$ms.run = ms_runs[ match(rfm$from, raw.filenames)  ]
    }
    .self$raw_file_mapping = rfm
    ## indicate to outside that a new table is ready
    .self$mapping.creation = .self$getMappingCreation()['auto']
  }
  ## do the mapping    
  v.result = as.factor(.self$raw_file_mapping$to[match(raw.filenames, .self$raw_file_mapping$from)])
  
  ## check for NA's
  if (any(is.na(v.result)))
  { ## if mapping is incomplete
    missing = unique(raw.filenames[is.na(v.result)])
    if (.self$mapping.creation != .self$getMappingCreation()['user'])
    {
      stop("Hithero unknown Raw files: " %+% paste(missing, collapse=", ", sep="") %+% " occurred in file '" %+% file %+% "' which were not present in previous txt files.")
    }
    ## the user has re-run the pipeline with more Raw files,
    ## but the old _filename_sort.txt file was used to read the (now incomplete mapping)
    warning("Incomplete mapping file '", .self$external.mapping.file, "'.\nAugmenting shortened Raw files:\n  " %+% 
              paste(missing, collapse="\n  ", sep="") %+% ".\nEdit the table if necessary and re-run PTXQC.")
    ## augment
    addon = .self$getShortNamesStatic(missing, max_length, nrow(.self$raw_file_mapping) + 1)
    if (!is.null(ms_runs)) {
      missing_msrun = ms_runs [ match( missing, raw.filenames ) ]
      addon$ms.run = missing_msrun [ match( addon$from, missing) ]
    }
    .self$raw_file_mapping = rbind(.self$raw_file_mapping, addon)
    ## redo mapping with full map
    v.result = as.factor(.self$raw_file_mapping$to[match(raw.filenames, .self$raw_file_mapping$from)])
    
  } 
  cat(paste0(" done\n"))
  return (v.result)
}, 
                         
getShortNamesStatic = function(raw.files, max_len, fallbackStartNr = 1)
{
  "Static method: Shorten a set of Raw file names and return a data frame with the mappings.
   Mapping will have: $from, $to and optionally $best.effort (if shorting was unsuccessful and numbers had to be used)
   \\itemize{
     \\item{\\verb{raw.files}  Vector of Raw files.}
     \\item{\\verb{max_len} Maximal length of shortening results, before resorting to canonical names (file 1,...).}
     \\item{\\verb{fallbackStartNr} Starting index for canonical names.}
   }
   \\subsection{Return Value}{ data.frame with mapping.}
  "
  rf_name = raw.files
  ## remove prefix
  rf_name_s = delLCP(rf_name, 
                     min_out_length = 8,
                     add_dots = TRUE)
  ## remove infix (2 iterations)
  rf_name_s = simplifyNames(rf_name_s, 
                            2, 
                            min_LCS_length = 7,
                            min_out_length = 8)
  
  ## check if shorter filenames are still unique (they should be.. if not we have a problem!!)
  if (length(rf_name) != length(unique(rf_name_s)))
  {
    cat("Original names:\n")
    cat(rf_name)
    cat("Short names:\n")
    cat(rf_name_s)
    stop("While loading MQ data: shortened raw filenames are not unique! This should not happen. Please contact the developers and provide the above names!")
  }
  df.mapping = data.frame(from = rf_name, to = rf_name_s, stringsAsFactors = FALSE)
  
  ## always include 'best.effort' column
  df.mapping[, "best.effort"] = df.mapping$to
  
  ## check if the minimal length was reached
  if (max(nchar(df.mapping$to)) > max_len)
  { ## resort to short naming convention
    cat("Filenames are longer than the maximal allowed size of '" %+% max_len %+% "'. Resorting to short versions 'file X'.\n\n")
    maxl = length(raw.files) - 1 + fallbackStartNr
    df.mapping$to = paste("file", sprintf(paste0("%0", nchar(maxl), "d"), fallbackStartNr:maxl)) ## with leading 0's if required
  }
  return(df.mapping)
},



plotNameMapping = function(.self)
{
  "Plots the current mapping of Raw file names to their shortened version.
  
   Convenience function to plot the mapping (e.g. to a PDF device for reporting).
   The data frame can be accessed directly via \\verb{.self$raw_file_mapping}.
   If no mapping exists, the function prints a warning to console and returns NULL (which is safe to use in print(NULL)).
  
   @return if mapping is available, returns a list of plots 'plots' and a Html table string 'htmlTable' ; 'NULL' otherwise.
  
  "
  if (nrow(.self$raw_file_mapping) == 0)
  {
    cat("No mapping found. Omitting plot.")
    return (NULL);
  }
  
  table_header = c("original", "short\nname")
  xpos = c(9, 11)
  extra = ""
  has_best_effort = FALSE
  if ("best.effort" %in% colnames(.self$raw_file_mapping))
  {
    has_best_effort = TRUE
    table_header = c(table_header, "best\neffort")
    xpos = c(9, 11, 13)
    if (all(.self$raw_file_mapping$to != .self$raw_file_mapping$best.effort)) {
      extra = "\n(automatic shortening of names was not sufficient - see 'best effort')"
    }
    
  }
  
  #mq_mapping = mq$raw_file_mapping
  mq_mapping = .self$raw_file_mapping
  pl_title = "Mapping of Raw files to their short names\nMapping source: " %+% .self$mapping.creation %+% extra;
  
  mappingChunk = function(mq_mapping)
  {
    mq_mapping$ypos = -(1:nrow(mq_mapping))
    head(mq_mapping)
    mq_mapping.long = melt(mq_mapping, id.vars = c("ypos"))
    head(mq_mapping.long)
    mq_mapping.long$variable = as.character(mq_mapping.long$variable)
    mq_mapping.long$col = "#000000";
    mq_mapping.long$col[mq_mapping.long$variable=="to"] = "#5F0000"
    mq_mapping.long$variable[mq_mapping.long$variable=="from"] = xpos[1]
    mq_mapping.long$variable[mq_mapping.long$variable=="to"] = xpos[2]
    if (nchar(extra)) mq_mapping.long$variable[mq_mapping.long$variable=="best.effort"] = xpos[3]
    mq_mapping.long$variable = as.numeric(mq_mapping.long$variable)
    mq_mapping.long$size = 2;
    
    df.header = data.frame(ypos = 1, variable = xpos, value = table_header, col = "#000000", size=3)
    mq_mapping.long2 = rbind(mq_mapping.long, df.header)
    mq_mapping.long2$hpos = 0 ## left aligned,  1=right aligned
    mq_mapping.long2$hpos[mq_mapping.long2$variable==xpos[1]] = 1
    mq_mapping.long2$hpos[mq_mapping.long2$variable==xpos[2]] = 0
    if (nchar(extra)) mq_mapping.long2$hpos[mq_mapping.long2$variable==xpos[3]] = 0
    
    mqmap_pl = ggplot(mq_mapping.long2, aes_string(x = "variable", y = "ypos"))  +
      geom_text(aes_string(label="value"), color = mq_mapping.long2$col, hjust=mq_mapping.long2$hpos, size=mq_mapping.long2$size) +
      coord_cartesian(xlim=c(0,20)) +
      theme_bw() +
      theme(plot.margin = grid::unit(c(1,1,1,1), "cm"), line = element_blank(), 
            axis.title = element_blank(), panel.border = element_blank(),
            axis.text = element_blank(), strip.text = element_blank(), legend.position = "none") +
      ggtitle(pl_title)
    return(mqmap_pl)
  }
  l_plots = byXflex(mq_mapping, 1:nrow(mq_mapping), 20, mappingChunk, sort_indices = FALSE);
  return (list(plots = l_plots, htmlTable = getHTMLTable(.self$raw_file_mapping, pl_title)))
  
}, 

getRawfm = function(.self)
{
  "Wrapper function for member 'raw_file_mapping', ensuring that $to is a factor"
  
  tmp = .self$raw_file_mapping
  tmp$to = factor(tmp$to)
  return(tmp)
},


readMappingFile = function(.self, filename)
{
  "Reads a mapping table of full Raw file names to shortened names.

  The internal structure \\verb{raw_file_mapping} is created using this file.
  If the file is missing, nothing is done.
 
  The file must have two columns named: 'orig.Name' and 'new.Name' and use Tab as separator.
  This file can be used to manually substitute Raw file names within the report.
  The ordering of Raw files in the report can be changed by re-arranging the rows.
  I.e.
  \\preformatted{
  orig.Name  new.Name
  2011_05_30_ALH_OT_21_VIL_TMT_FR01   myfile A
  2011_05_30_ALH_OT_22_VIL_TMT_FR02   another B
  }
 
  @param filename  Source filename to read.
  @return Returns \\verb{TRUE} if file was read, \\verb{FALSE} if it does not exist.
"
  
  if (file.exists(filename))
  {
    message(paste0("Reading mapping file '", filename, "'\n"))
    dfs = read.delim(filename, comment.char="#", stringsAsFactors = FALSE)
    colnames(dfs) = gsub("_", ".", colnames(dfs)) ## legacy support for old "best_effort" column (now "best.effort")
    req_cols = c(from = "orig.Name", to = "new.Name")
    if (!all(req_cols %in% colnames(dfs)))
    {
      stop("Input file '", filename, "' does not contain the columns '", paste(req_cols, collapse="' and '"), "'.",
           " Please fix and re-run PTXQC!")
    }
    req_cols = c(req_cols, best.effort = "best.effort", ms.run = "ms.run") ## augment 
    colnames(dfs) = names(req_cols)[match(colnames(dfs), req_cols)]
    
    if (any(duplicated(dfs$from)) | any(duplicated(dfs$to)))
    {
      dups = c(dfs$from[duplicated(dfs$from)], dfs$to[duplicated(dfs$to)])
      stop("Input file '", filename_sorting, "' has duplicate entries ('", paste(dups, collapse=", "), ")'!",
           " Please fix and re-run PTXQC!")
    }
    dfs
    dfs$to = factor(dfs$to, levels = unique(dfs$to), ordered = TRUE) ## keep the order
    dfs$from = factor(dfs$from, levels = unique(dfs$from), ordered = TRUE) ## keep the order
    ## set internal mapping
    .self$raw_file_mapping = dfs
    ## set who defined it
    .self$mapping.creation = 'file (user-defined)'
    .self$external.mapping.file = filename; ## remember filename for later error messages
    return (TRUE)
  }
  return (FALSE)
},


writeMappingFile = function(.self, filename)
{
  "Writes a mapping table of full Raw file names to shortened names.
  
  The internal structure \\verb{raw_file_mapping} is written to the
  file specified.
  File is only created if mapping exists (in .self$raw_file_mapping).
  
  @param filename  Target filename to create.
  @return Returns NULL.
  "
  if (nrow(.self$raw_file_mapping) == 0)
  {
    cat("No mapping found. Writing not possible!")
    return (FALSE)
  }
  
  dfs = data.frame(orig.Name = .self$raw_file_mapping$from, new.Name = .self$raw_file_mapping$to)
  if (nrow(dfs) == 0) return(NULL)
  
  if ("best.effort" %in% colnames(.self$raw_file_mapping)) {
    dfs$best.effort = .self$raw_file_mapping[, "best.effort"]
  }
  
  if ("ms.run" %in% colnames(.self$raw_file_mapping)) {
    dfs$ms.run = .self$raw_file_mapping[,"ms.run"]
  }
  
  ## use a file handle to avoid warning from write.table() when appending
  ## a table with column names 'Warning(): appending column names to file'
  FH = file(filename, "w")
  cat(file = FH,
      "# This file can be used to manually substitute Raw file names within the report.",
      "# The ordering of Raw files in the report can be changed by re-arranging the rows.",
      sep = "\n")
  write.table(x = dfs, file = FH, quote = FALSE, sep="\t", row.names = FALSE)
  close(FH) ## flush
  return (TRUE)
},



getMappingCreation = function(.self)
{
  "A static function"
  return(c(user = 'file (user-defined)', auto = 'automatic'))
}





) ## end methods list
) ## end RefClass
