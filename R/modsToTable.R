#' 
#' Convert list of (mixed)modifications to a frequency table
#' 
#' @param mod_list A vector with modifications, each for a specific peptide. Multiple mods per entry are allowed, each separated by comma.
#' @return A data.frame with 'modification_names' and 'Freq' (0-100)
#'
#' @export
#' 
#' @examples
#' modsToTable(c("Ox (M)", "Unmodified", "Ox (M),Acetyl (Prot N-term)", "2 Ox (M)", "Unmodified", "Unmodified"))
#' 
#'
modsToTable = function(mod_list)
{
  modification_names = unlist(strsplit(mod_list, ",", fixed=TRUE))
  tt = data.frame(table(modification_names) / length(mod_list) * 100)
  if (! ("modification_names" %in% colnames(tt))) stop("table has wrong column names!")
  return(tt)
}
