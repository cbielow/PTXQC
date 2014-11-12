#'
#' Generate a Heatmap from a list of QC measurements.
#' 
#' Each list entry is a data.frame with two columns. 
#' The first one contains the Raw file name (or the short version).
#' and should be named 'raw.file' (or 'fc.raw.file'). 
#' The second column can have an arbitrary name
#' and contains quality values in the range [0,1]. If values are outside this range, 
#' a warning is issued and values are cut to the nearest allowed value (e.g. '1.2' becomes '1').
#' Columns are ordered by name. Then characters in the column name are modified as follows: "." --> ": " and "_" --> " ".
#' All substrings enclosed by 'X[0-9]*X.' will be removed (can be used for sorting columns).
#' 
#' To judge the quality of each raw file a summary column is added, values being the mean of all other columns per row.
#'        
#' @param QCM List of data.frames, each having a 'raw.file' column and a metric column of arbitrary name
#' @param raw_file_mapping Data.frame with 'from' and 'to' columns for name mapping to unify names from list entries
#' @return A ggplot object for printing
#'
#' @importFrom reshape2 melt
#' @importFrom reshape2 dcast
#'
#' @export
#'
#' @examples
#'   mapping = data.frame(from=c("A.raw","B.raw"), to=c("A","B"))
#'   QC_data = list(somedata  = data.frame(raw.file=c("A.raw", "B.raw"), "X005X.EVD.Some_later" = 1:2),
#'              someother = data.frame(raw.file=c("A.raw", "B.raw"), "X002X.EVD.First_col" = 5:6),
#'              middle = data.frame(raw.file=c("A.raw", "B.raw"), "X003X.EVD.Middle_col" = 3:4))
#'   getQCHeatMap(QC_data, mapping)
#' 
#' 
#' 
getQCHeatMap = function(QCM, raw_file_mapping)
{
  
  QCM_shortNames = lapply(QCM, function(x) {
    if ("raw.file" %in% colnames(x)) {
      x$fc.raw.file = renameFile(x$raw.file, raw_file_mapping)  ## create short name column
      x = x[, !(colnames(x) %in% "raw.file")]  ## remove raw.file column
    }
    if (!("fc.raw.file" %in% colnames(x))) stop("Error in getQCHeatMap(): 'fc.raw.file' column missing from QC measure.")
    return(x)
  })
  
  ## final heat map of QC metrics
  QCM_final = Reduce(function(a,b) merge(a,b,all = TRUE), QCM_shortNames)
  ## add summary column
  QCM_final$X999X.Summary = apply(QCM_final[,!grepl("fc.raw.file", colnames(QCM_final))], 1, function(row) {
    row[is.na(row)] = 0  ## replace NA with 0, since it will bias the mean otherwise
    return(mean(row))
  })
  
  ## reorder file names
  QCM_final = QCM_final[match(raw_file_mapping$to, QCM_final$fc.raw.file), ]
  ## ... fix factor levels
  QCM_final$fc.raw.file = factor(QCM_final$fc.raw.file, levels = QCM_final$fc.raw.file)
  
  ## reorder columns
  QCM_final = QCM_final[, order(colnames(QCM_final))]
    
  QCM_final.m = melt(QCM_final, id.vars="fc.raw.file")
  QCM_final.m$variable = factor(QCM_final.m$variable, ordered = TRUE)
  
  ## some files might not be in the original list (will give NA in table)
  QCM_final.m$value[is.na(QCM_final.m$value)] = 0
  
  if (any(QCM_final.m$value > 1))
  {
    warning("getQCHeatMap() received quality data values larger than one! This can be corrected here, but should be done upstream.")
    QCM_final.m$value[QCM_final.m$value > 1] = 1
  }
  if (any(QCM_final.m$value < 0))
  {
    warning("getQCHeatMap() received quality data values smaller than zero! This can be corrected here, but should be done upstream.")
    QCM_final.m$value[QCM_final.m$value < 0] = 0
  }
  
  ## rename metrics (remove underscores and dots)
  QCM_final.m$variable = gsub("X[0-9]*X\\.", "", as.character(QCM_final.m$variable))  ## prefix sorting
  QCM_final.m$variable = gsub("^\\s+|\\s+$", "", QCM_final.m$variable)             ## trim
  QCM_final.m$variable = gsub("_", " ", QCM_final.m$variable)
  QCM_final.m$variable = gsub("\\.", ": ", QCM_final.m$variable)
  QCM_final.m$variable2 = factor(QCM_final.m$variable, levels = unique(QCM_final.m$variable))
  
  p = ggplot(QCM_final.m, aes_string(y="fc.raw.file", x="variable2")) +
            geom_tile(aes_string(fill = "value"), colour = "white") + 
            scale_fill_continuous(low = "red", high = "green", limits=c(0, 1)) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
            ggtitle("Performance overview") +
            xlab("") +
            ylab("Raw file") +
            scale_y_discrete_reverse(QCM_final.m$fc.raw.file, breaks = ggAxisLabels)
  #print(p)
  return(list(plot = p, table = dcast(QCM_final.m, fc.raw.file ~ variable)))
}  
