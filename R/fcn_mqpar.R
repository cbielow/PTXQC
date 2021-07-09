
#' Retrieve a parameter value from a mqpar.xml file
#' 
#' If the file has the param, then return it as string.
#' If the file is missing, warning is shown and NULL is returned.
#' If the param (i.e. XML tag) is unknown or cannot be extracted, the program will quit (since this is a hard error).
#' When multiple occurrences of the param are found (usually due to parameter groups), we test if the values are all identical.
#' If so, the value is returned. If the values are different, a warning is emitted and NULL is returned unless 'allow_multiple = TRUE'
#' 
#' E.g. calling getMQPARValue("mqpar.xml", "//firstSearchTol")
#' will look up the line
#' <firstSearchTol>20</firstSearchTol>
#' and return "20" (string!).
#' 
#' 
#' @param mqpar_filename Filename (incl. absolute or relative path) to the mqpar.xml file
#' @param xpath          An XPath to extract the content of XML tag(s), e.g. '//firstSearchTol'
#' @param allow_multiple If the XPath expression returns more than one value, all values must be identical (not allowing multiple different values) or 'stop()' is called
#' 
#' @return The stored value as string(!)
#'
#' @export
#'
getMQPARValue = function(mqpar_filename, xpath, allow_multiple = FALSE)
{
  #xpath = "//firstSearchTol"
  #mqpar_filename = txt_files$mqpar
  if (!file.exists(pattern=mqpar_filename)) {
    message("Info: The file '", mqpar_filename, "' was not found. MaxQuant parameters could not be extracted. Will fall back to default value, which however is only an approximation.",
            " Please either: a) copy the mqpar.xml which was used for this MQ analysis into your TXT folder or,",
            " b) make sure that you configure all YAML parameters whose name starts with 'MQpar_' correctly.", immediate. = TRUE)
    return (NULL)
  }
  
  lines = xml2::read_xml(mqpar_filename)
  ## is the tag present multiple times? (if yes, we found parameter groups)
  results = xml_text(xml_find_all(lines, xpath))

  if (length(results) == 0) {
    stop("getMQPARValue(): The XPath '", xpath, "' was not found in the file '", mqpar_filename, "'. Please contact the package support.", call. = FALSE);
  }
  
  if ((allow_multiple == FALSE)) {
    
    if (length(unique(results)) == 1) {
      return (results[1]);
    }
    
    warning("getMQPARValue(): The XPath '", xpath, "' was found more than once in the file '", mqpar_filename, "' with different values (probably due to usage of parameter groups).",
            " PTXQC currently cannot deal with that -- the YAML param is going to be used. Sorry.", immediate. = TRUE);
    return (NULL);
  }
  
  ## all tests passed, return the result(s)
  return (results)
}