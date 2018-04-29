
#' Retrieve a parameter value from a mqpar.xml file
#' 
#' If the file has the param, then return it as string.
#' If the file is missing, warning is shown and NULL is returned.
#' If the param (i.e. XML tag) is unknown or cannot be extracted, the program will quit (since this is a hard error).
#' When multiple occurrences of the param are found (usually due to parameter groups), we test if the values are all identical.
#' If so, the value is returned. If the values are different, a warning is emitted and NULL is returned.
#' 
#' E.g. calling getMQPARValue("mqpar.xml", "firstSearchTol")
#' will look up the line
#' <firstSearchTol>20</firstSearchTol>
#' and return "20" (string!).
#' 
#' 
#' @param mqpar_filename Filename (incl. absolute or relative path) to the mqpar.xml file
#' @param param_name     XML tag name, e.g. 'firstSearchTol' from which to read the value
#' 
#' @return The stored value as string(!)
#'
#' @export
#'
getMQPARValue = function(mqpar_filename, param_name)
{
  #param_name = "firstSearchTol"
  #mqpar_filename = txt_files$mqpar
  ## TODO: at some point we might use a real XML parser, but for now, this would be overkill and
  ##       also add another dependency library
  if (!file.exists(pattern=mqpar_filename)) {
    warning("The file '", mqpar_filename, "' was not found. MaxQuant parameters could not be extracted. Will fall back to default value, which however is only an approximation.",
            " Please either: a) copy the mqpar.xml which was used for this MQ analysis into your TXT folder or,",
            " b) make sure that you configure all YAML parameters whose name starts with 'MQpar_' correctly.", immediate. = TRUE)
    return (NULL)
  }
  
  lines = readLines(con = mqpar_filename, warn = FALSE)
  idx = grep(param_name, lines)
  ## is the tag present multiple times? (if yes, we found parameter groups)
  results = gsub(paste0("[ ]*<", param_name, ">(.*)</", param_name, ">[ ]*"), "\\1", lines[idx])
  
  ## if regex did not work, the whole line will be returned, including the tag
  if (length(grep(param_name, results)) > 0)
  {
    stop("getMQPARValue(): The parameter '", param_name, "' was found but could not be extracted from the line(s)\n  ", paste(lines[idx], collapse="\n  "), "\n  Please contact the package support.", call. = FALSE)
  }
  
  if (length(unique(results)) > 1) {
    warning("getMQPARValue(): The parameter '", param_name, "' was found more than once in the file '", mqpar_filename, "' with different values (probably due to usage of parameter groups).",
            " PTXQC currently cannot deal with that -- the YAML param is going to be used. Sorry.", immediate. = TRUE);
    return (NULL);
  } else if (length(results) == 0) {
    stop("getMQPARValue(): The parameter '", param_name, "' was not found in the file '", mqpar_filename, "'. Please contact the package support.", call. = FALSE);
  }
  
  ## all tests passed, return the unique result
  return (results[1])
}