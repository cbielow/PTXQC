# Hilfsfunktion für Mapping
getShortNames = function(raw.files, max_len, fallbackStartNr = 1)
{
  ##
  ## mapping will have: $from, $to and optionally $best.effort (if shorting was unsuccessful and numbers had to be used)
  ##
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
}
