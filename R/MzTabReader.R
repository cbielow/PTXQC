#'
#' Class to read an mzTab file and store the tables internally.
#' 
#' The 'sections' field is initialized after $readMzTab was called.
#' The 'fn_map' fields should be initialized via ...$fn_map$readMappingFile(...) manually if user-defined filename mappings are desired
#' and is automatically updated/queried when $readMzTab is called.
#' 
#' @field sections MzTab sections as list. Valid list entries are: "MTD", "PRT", "PEP", "PSM", "SML", "filename" and "comments"
#' @field fn_map FilenameMapper which can translate raw filenames into something shorter
#'
#'
#'
MzTabReader = setRefClass("MzTabReader",
                       
                       fields = list(sections = "list",
                                     fn_map = "FilenameMapper"
                       ),
                       methods = list(
                         initialize=function(sections = NA_character_,
                                             workerFcn = function(){},
                                             qcCat = NA_character_,
                                             qcName = NA_character_,
                                             orderNr = NaN) {
                           .self$sections = list();
                           .self$fn_map = NULL;
                           
                           return(.self)
                         },
                         #'
readMzTab = function(file) {
  "Read a mzTab file into a list of 5 data.frames (one df per mzTab section).
   Data.frames in the resulting list are named as follows:
     'MTD', 'PRT', 'PEP', 'PSM', 'SML',.
   Additionally, 'filename' and 'comments' are valid list elements."

  
  ## this implementation is derived from with minor modifications
  ## https://github.com/lgatto/MSnbase/blob/master/R/MzTab.R
  
  lines = readLines(file)
  # remove empty lines
  lines = lines[ nzchar(lines) ]
  
  ## Split on the first two characters (so headers stay in
  ## the same group as table content rows)
  lineType = substring(lines, 1, 2)
  
  ## Could be stricter in the type checking to check that all
  ## three of the first characters match the 10 allowed types
  ## but since it doesn't affect parsing, I don't think it's
  ## worth bothering.
  allowed_types = c("CO", "MT", "PR", "PE", "PS", "SM")
  stopifnot(all(lineType %in% allowed_types))
  linesByType = split(lines, lineType)
  
  ## Comments are easy: just strip the first four characters
  ## from each line.
  comments = substring(linesByType[["CO"]], 5)
  
  ## Parse the other five blocks in a loop, then fix up
  ## metadata afterwards
  res = setNames(
    lapply(
      linesByType[c("MT", "PR", "PE", "PS", "SM")],
      function(x) {
        if (length(x) == 0) return(data.frame())
        return(read.delim(text = x,
                          na.strings = c("", "null"),
                          stringsAsFactors = FALSE)[,-1])
      }),
    c("MTD", "PRT", "PEP", "PSM", "SML"))
  
  ## rewrite MetaData as named vector
  res[["MTD"]] = setNames(res[["MTD"]][,2], res[["MTD"]][, 1])
  
  res[["filename"]] = file
  res[["comments"]] = comments
  .self$sections = res
  return (NULL)
}

) # methods
) # class
