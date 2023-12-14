#' 
#' Convert list of (mixed)modifications to a frequency table
#' 
#' @param df_evd data.frame with 'fc.raw.file' and a 'modifications' column, which contains the modifications for each peptide.
#' @param name_unmod String in 'modifications' which represents an unmodified peptide
#' @param name_unmod_inverse If non-empty, then inverse the frequencies of the 'name_unmod' modifications (i.e. 100-x) IFF they are >=50\% on average (across Raw files) and rename them to this string
#' @return A data.table with 'fc.raw.file', 'modification_names' (factor), and 'Freq' (0-100)
#'
#' @export
#' 
#' @examples
#' data = data.frame(fc.raw.file = rep(c("file A", "file B"), each = 3),
#'                   modifications = c("Oxidation (M)", "Unmodified", "Oxidation (M),Acetyl (Protein N-term)", "2 Oxidation (M)", "Unmodified", "Unmodified"))
#' modsToTableByRaw(data)
#' 
#'
modsToTableByRaw = function(df_evd, name_unmod = "Unmodified", name_unmod_inverse = "Modified (total)")
{
  dt_evd = data.table::data.table(df_evd)
  mods_tbl = dt_evd[, modsToTable(.SD$modifications), by = fc.raw.file]
  
  mods_tbl$modification_names = as.factor(mods_tbl$modification_names) ## ensure all mods are known when subsetting the data later
  
  ## inverse frequencies for 'modified' if they are >=50% (to make the plot more compact)
  if (nchar(name_unmod_inverse)){
    if (mean(mods_tbl[modification_names==name_unmod, Freq]) > 50)
    {
      mods_tbl[modification_names==name_unmod, Freq := 100 - Freq]
      mods_tbl[modification_names==name_unmod, modification_names := name_unmod_inverse]
    }
    mods_tbl = droplevels(mods_tbl) ## drop 'Unmodified' level
  }
  
  
  return(as.data.frame(mods_tbl))
}
