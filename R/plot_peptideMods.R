
#'
#' Plot peptide modification frequencies
#' 
#' The input is a data.frame, as obtained from modsToTableByRaw().
#' 
#' 
#' 
#' @param tbl A data.frame with 'fc.raw.file', 'modification_names' (can be a factor), and 'Freq' (0-100)
#' @param y_max The upper limit of the y-axis's (==Freq); useful for multiple plots with identical limits; if 'NA' the limit is computed from the given 'tbl'
#' @param show_missing_modification_levels If 'tbl$modification_names' is a factor and has more (but missing) levels than actually used, should missing values be dropped or assumed as '0' frequency?
#' @return GGplot object
#'
#' @import ggplot2
#' @export
#' 
#' @examples 
#' 
#'  data = data.frame(fc.raw.file = rep(c("file A", "file B"), each=3),
#'                    modifications = c("Oxidation (M)",
#'                    "Unmodified",
#'                    "Oxidation (M), Acetyl (Protein N-term)",
#'                    "2 Oxidation (M)", 
#'                    "Unmodified", 
#'                    "Unmodified"))
#'  tbl = modsToTableByRaw(data)
#'  plot_peptideMods(tbl,show_missing_modification_levels = TRUE)
#' 
plot_peptideMods = function(tbl, y_max = NA, show_missing_modification_levels = TRUE)
{
  if (is.na(y_max)){
    y_max = max(tbl$Freq, na.rm = TRUE)
  } 

  ## augment '0'-frequency for missing factors
  if (show_missing_modification_levels) {
    tbl$fc.raw.file = droplevels(as.factor(tbl$fc.raw.file))  ## we need fc.raw.file to be a factor, but we do not want its unused levels
    all_combinations = expand.grid(lapply(tbl[c("fc.raw.file","modification_names")], levels))
    new_tbl = merge(tbl, all_combinations , all.y=TRUE)
    new_tbl$Freq[is.na(new_tbl$Freq)] = 0  ## replace new 'NA's with 0
    tbl = new_tbl
  }
    
  p = ggplot(tbl, aes(x = .data$fc.raw.file, y = .data$Freq, fill = .data$modification_names, colour = .data$modification_names)) + ## use 'colour' for outline of '0'-frequency features (invisible otherwise)
    geom_col(position = "dodge") +
    xlab("") +
    ylab("Occurence [%]") +
    scale_x_discrete_reverse(tbl$fc.raw.file) +
    guides(fill = guide_legend(reverse = TRUE), color = guide_none()) + ## reverse order of mods in legend, to match order in plot
    ylim(0, y_max) +
    ggtitle("EVD: variable modifications per Raw file") + 
    coord_flip()
  return(p)
}
