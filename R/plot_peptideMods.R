
#'
#' Plot peptide modification frequencies
#' 
#' The input is a data.frame, as obtained from modsToTableByRaw().
#' 
#' @param tbl A data.frame with 'fc.raw.file', 'modification_names', and 'Freq' (0-100)
#' @param y_max The upper limit of the y-axis's (==Freq); useful for multiple plots with identical limits
#' @return GGplot object
#'
#' @import ggplot2
#' @export
#' 
#' @examples 
#' 
#'  data = data.frame(fc.raw.file = rep(c("file A", "file B"), each=3),
#'                    modifications = c("Oxidation (M)", "Unmodified", "Oxidation (M), Acetyl (Protein N-term)", "2 Oxidation (M)", "Unmodified", "Unmodified"))
#'  tbl = modsToTableByRaw(data)
#'  plot_peptideMods(tbl)
#' 
plot_peptideMods = function(tbl, y_max = NA)
{
  if (is.na(y_max)) y_max = max(tbl$Freq, na.rm = TRUE)
  
  p = ggplot(tbl, aes_string(x = "fc.raw.file", y = "Freq", fill = "modification_names")) +
    geom_col(position = "dodge") +
    xlab("") +
    ylab("Occurence [%]") +
    scale_x_discrete_reverse(tbl$fc.raw.file) +
    ylim(0, y_max) +
    ggtitle("Variable modifications per Raw file") + 
    coord_flip()
  return(p)
}
