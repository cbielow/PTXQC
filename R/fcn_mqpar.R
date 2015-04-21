
#' Retrieve a parameter value from a mqpar.xml file
#' 
#' If the file has the param, then return it as string.
#' If the file is missing, warning is shown and NULL is returned.
#' If the param (i.e. XML tag) is unknown, ambiguous or cannot be extracted, the program will quit (since this is a hard error)
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
  #mqpar_filename = "mqpar_MQ15_MBRon.xml"
  ## TODO: at some point we might use a real XML parser, but for now, this would be overkill and
  ##       also add another dependency library
  if (!file.exists(pattern=mqpar_filename)) {
    warning("The file '", mqpar_filename, "' was not found. MaxQuant parameters could not be extracted. Will fall back to default value, which however is only an approximation.",
            " Please either: a) copy the mqpar.xml which was used for this MQ analysis into your TXT folder or,",
            " b) make sure that you configure all YAML parameters whose name contains '_mqpar' correctly.", immediate. = T)
    return (NULL)
  }
  
  lines = readLines(con = mqpar_filename, warn = F)
  idx = grep(param_name, lines)
  if (length(idx) != 1) {
    stop("getMQPARValue(): The parameter '", param_name, "' was found more than once in the file '", mqpar_filename, "'. Cannot resolve ambiguity. Please contact the package support.", immediate. = T)
  }
  result = gsub(paste0("[ ]*<", param_name, ">(.*)</", param_name, ">[ ]*"), "\\1", lines[idx])
  if (length(grep(param_name, result)) > 0)
  {
    stop("getMQPARValue(): The parameter '", param_name, "' was found but could not be extracted from the line '", lines[idx], "'. Please contact the package support.", immediate. = T)
  }
  
  return (result)
}