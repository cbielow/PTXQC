#'
#' Read a mzTab file into a list of 5 data.frames (one df per mzTab section).
#'
#' Data.frames in the resulting list are named as follows:
#' "MTD", "PRT", "PEP", "PSM", "SML"
#' 
#' Additionally, "filename" and "comments" are valid list elements.
#' 
#' @param file A path to an mzTab file to read in.
#'
readMzTab = function(file) {
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
  
  return (res)
}
