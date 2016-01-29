## default value for NA's in heatmap
HEATMAP_NA_VALUE = -Inf

#'
#' Generate a Heatmap from a list of QC measurements.
#' 
#' Each list entry is a data.frame with two columns. 
#' The first one contains the Raw file name (or the short version).
#' and should be named 'raw.file' (or 'fc.raw.file'). 
#' The second column's name must be an expression (see ?plotmath)
#' and contains quality values in the range [0,1]. If values are outside this range, 
#' a warning is issued and values are cut to the nearest allowed value (e.g. '1.2' becomes '1').
#' List entries are merged and columns are ordered by name.
#'   
#' All substrings enclosed by 'X[0-9]*X.' will be removed (can be used for sorting columns).
#' The resulting string is evaluated as an expression. 
#' E.g. parse(text = <colname>)
#' 
#' To judge the overall quality of each raw file a summary column is added, 
#' values being the mean of all other columns per row.
#'        
#' @param QCM List of data.frames, each having a 'raw.file' or 'fc.raw.file' column and a metric column with numeric values
#' @param raw_file_mapping Data.frame with 'from' and 'to' columns for name mapping to unify names from list entries
#' @return A ggplot object for printing
#'
#' @importFrom reshape2 melt
#' @importFrom reshape2 dcast
#' @importFrom plyr empty
#'
#' @export
#'
#' @examples
#'   mapping = data.frame(from=c("A.raw","B.raw"), to=c("A","B"))
#'   QC_data = list(
#'                somedata  = data.frame(raw.file=c("A.raw", "B.raw"),
#'                                       "X005X.EVD:~Some~later~MS^2" = 1:2),
#'                someother = data.frame(raw.file=c("A.raw", "B.raw"),
#'                                       "X002X.EVD:~First~col" = 5:6),
#'                   middle = data.frame(raw.file=c("A.raw", "B.raw"),
#'                                       "X003X.EVD:~Middle_col" = 3:4))
#'   getQCHeatMap(QC_data, mapping)
#' 
#' 
#' 
getQCHeatMap = function(QCM, raw_file_mapping)
{
  QCM_shortNames = lapply(QCM, function(x) {
    if (empty(x)) return(NULL) ## if metric was not computed, default DF is empty
    if ("raw.file" %in% colnames(x)) {
      x$fc.raw.file = renameFile(x$raw.file, raw_file_mapping)  ## create short name column
      x = x[, !(colnames(x) %in% "raw.file")]  ## remove raw.file column
    }
    if (!("fc.raw.file" %in% colnames(x))) {
      cat(paste("columns:", paste0(colnames(x), collapse=", ", sep="")))
      stop("Internal error in getQCHeatMap(): 'fc.raw.file' column missing from QC measure.")
    } 
    ## check if fc.raw.filenames are known (e.g. when column was named fc.raw.file but values are from raw.file)
    if (!(all(x$fc.raw.file %in% raw_file_mapping$to))) {
      stop("Internal error in getQCHeatMap(): 'fc.raw.file' column has invalid entries for '", grepv("^X", colnames(x)), "'!")
    }
  
    return(x)
  })
  
  ## final heat map of QC metrics
  QCM_final = Reduce(function(a,b) merge(a,b,all = TRUE), QCM_shortNames)
  ## add summary column
  QCM_final$"X999X_catGen_Average~Overall~Quality" = apply(QCM_final[,!grepl("fc.raw.file", colnames(QCM_final)), drop = FALSE], 1, function(row) {
    row[is.infinite(row)] = NA  ## mask explicitly missing values, since it will bias the mean otherwise
    return(mean(row, na.rm = TRUE))
  })
  
  ## reorder file names
  QCM_final = QCM_final[match(raw_file_mapping$to, QCM_final$fc.raw.file), ]
  ## ... fix factor levels
  QCM_final$fc.raw.file = factor(QCM_final$fc.raw.file, levels = QCM_final$fc.raw.file)
  
  ## reorder columns
  QCM_final = QCM_final[, order(colnames(QCM_final))]
  ## add column numbering (ignore first column, which is 'fc.raw.file')
  idx = 2:(ncol(QCM_final)-1)
  colnames(QCM_final)[idx] = paste0(colnames(QCM_final)[idx], "~\"[", idx-1, "]\"")
  
  QCM_final.m = melt(QCM_final, id.vars="fc.raw.file")
  QCM_final.m$variable = factor(QCM_final.m$variable, ordered = TRUE)
  
  ## some files might not be in the original list (will receive 'bad' score in table)
  QCM_final.m$value[is.na(QCM_final.m$value)] = 0
  ## some other files might be missing on purpose
  QCM_final.m$value[is.infinite(QCM_final.m$value)] = NA
  
  if (any(QCM_final.m$value > 1, na.rm = TRUE))
  {
    warning("getQCHeatMap() received quality data values larger than one! This can be corrected here, but should be done upstream.")
    QCM_final.m$value[QCM_final.m$value > 1] = 1
  }
  if (any(QCM_final.m$value < 0, na.rm = TRUE))
  {
    warning("getQCHeatMap() received quality data values smaller than zero! This can be corrected here, but should be done upstream.")
    QCM_final.m$value[QCM_final.m$value < 0] = 0
  }
  
  ##
  ## rename metrics (remove underscores and dots)
  ##
  ## sorting prefix
  QCM_final.m$variable = gsub("X[0-9]*X", "", as.character(QCM_final.m$variable))  ## prefix sorting
  ## axis color
  col_axis = c("_catPrep_" = "#606060", "_catLC_" = "black", "_catMS_" = "#606060", "_catGen_" = "black")
  QCM_final.m$axisCat = gsub("^(_cat[A-Za-z]*_).*", "\\1", QCM_final.m$variable)  ## assign axis color
  QCM_final.m$axisCol = col_axis[match(QCM_final.m$axisCat, names(col_axis))]
  QCM_final.m$variable = gsub("_cat[A-Za-z]*_", "", QCM_final.m$variable)            ## remove color category from name
  ## trim
  QCM_final.m$variable = gsub("^\\s+|\\s+$", "", QCM_final.m$variable)

  QCM_final.m$variable2 = factor(QCM_final.m$variable, levels = unique(QCM_final.m$variable))
  if (any(is.na(QCM_final.m$value))) QCM_final.m$dummy_col = "NA" ## use color legend for missing values
  
  ## replace the x-axis labels by expressions
  QCM_final.m$variable3 = sapply(as.character(QCM_final.m$variable2), function(x) {
    parse(text=x)
  })
  
  p = ggplot(QCM_final.m, aes_string(y="fc.raw.file", x="variable2"))
  if (any(is.na(QCM_final.m$value))) {
    p = p + geom_tile(aes_string(fill = "value", colour = "dummy_col")) +
            scale_colour_manual(name="Missing", values=c("NA" = "grey50")) +
            guides(colour = guide_legend(override.aes = list(fill = 'grey50')))
  } else {
    p = p + geom_tile(aes_string(fill = "value"))
  }  
  p = p + scale_y_discrete_reverse(QCM_final.m$fc.raw.file, breaks = ggAxisLabels) + 
          scale_x_discrete(labels = QCM_final.m$variable3) +
          scale_fill_gradientn(colours = c("red", "black", "green"),
                               na.value = "grey50",
                               limits=c(0, 1), 
                               guide = guide_legend(title = "score"),
                               breaks=c(1,0.5,0),
                               labels=c("best","under\nperforming","fail")) +
          theme(axis.text.x = element_text(angle = 90, 
                                           vjust = 0.5, 
                                           hjust = 1, 
                                           colour=QCM_final.m$axisCol[!duplicated(QCM_final.m$variable2)])) +
          ggtitle("Performance overview") +
          xlab("") +
          ylab("Raw file")
  #print(p)
  return(list(plot = p, table = dcast(QCM_final.m, fc.raw.file ~ variable)))
}  
