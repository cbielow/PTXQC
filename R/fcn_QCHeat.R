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
#' @param lst_qcMetrics List of QCMetric objects
#' @param raw_file_mapping Data.frame with 'from' and 'to' columns for name mapping to unify names from list entries
#' @return A ggplot object for printing
#' 
#' @import ggplot2
#' 
getQCHeatMap = function(lst_qcMetrics, raw_file_mapping)
{
  if (length(lst_qcMetrics) == 0) stop("Heatmap: List of Qc metrics is empty!")
  lst.QCM = lapply(lst_qcMetrics, function(qcm) {
    qcm_sc = qcm$qcScores
    if (plyr::empty(qcm_sc)) return(NULL) ## if metric was not computed, default DF is empty
    if ("raw.file" %in% colnames(qcm_sc)) {
      qcm_sc$fc.raw.file = renameFile(qcm_sc$raw.file, raw_file_mapping)  ## create short name column
      qcm_sc = qcm_sc[, !(colnames(qcm_sc) %in% "raw.file")]  ## remove raw.file column
    }
    if (!("fc.raw.file" %in% colnames(qcm_sc))) {
      cat(paste("columns:", paste0(colnames(qcm_sc), collapse=", ", sep="")))
      stop("Internal error in getQCHeatMap(): 'fc.raw.file' column missing from QC measure.")
    } 
    ## check if fc.raw.filenames are known (e.g. when column was named fc.raw.file but values are from raw.file)
    if (!(all(qcm_sc$fc.raw.file %in% raw_file_mapping$to))) {
      cat("An error occured. Current filename mapping is:\n")
      print(raw_file_mapping)
      stop("Internal error in getQCHeatMap(): 'fc.raw.file' column has invalid entries for metric '", qcm$qcName, "' with fc.raw.file names: [", paste(qcm_sc$fc.raw.file, collapse=","),"]!")
    }
    return(qcm_sc)
  })
  lst.QCM = plyr::compact(lst.QCM) ## remove 'NULL' entries
  ## final heat map of QC metrics
  df.QCM = Reduce(function(a,b) merge(a,b,all = TRUE), lst.QCM)

  ## no metrics...  
  if (is.null(df.QCM))
  {
    p = ggText("HeatMap unavailable", "No metrics were computed (not enough data)");
    return(list(plot = p, table = NULL))
  }
  ## create summary column
  lst_qcMetrics[["qcMetric_AverageQualOverall"]]$setData(df.QCM)
  ## ... add it
  df.AverageQual = lst_qcMetrics[["qcMetric_AverageQualOverall"]]$qcScores
  if (plyr::empty(df.AverageQual)) df.QCMa = df.QCM  else df.QCMa = merge(df.QCM, df.AverageQual)

  ## get order and names for each metric
  df.meta = getMetaData(lst_qcMetrics)
  head(df.meta)
  
  ## reorder rows (file names)
  df.QCMa = df.QCMa[match(raw_file_mapping$to, df.QCMa$fc.raw.file), ]
  ## ... fix factor levels
  df.QCMa$fc.raw.file = factor(df.QCMa$fc.raw.file, levels = df.QCMa$fc.raw.file)
  
  
  ## reorder columns (metrics)
  head(df.QCMa)
  qc_names_all = as.character(df.meta$name) ## contains metrics which have no scores (but the order is ok)
  qc_names_all_scores = qc_names_all[qc_names_all %in% colnames(df.QCMa)]
  df.QCMa = df.QCMa[, c("fc.raw.file", qc_names_all_scores)]
  ## add column numbering (ignore first column, which is 'fc.raw.file')
  df.QCMan = df.QCMa
  idx = 2:(ncol(df.QCMan))
  colnames(df.QCMan)[idx] = paste0(colnames(df.QCMa)[idx], "~\"[", idx-1, "]\"")
  colnames_wNum_map = data.frame(name = colnames(df.QCMa), nameWnum = colnames(df.QCMan))
  
  QCM_final.m = reshape2::melt(df.QCMan, id.vars="fc.raw.file")
  QCM_final.m$variable = factor(QCM_final.m$variable, ordered = TRUE)
  
  ## some files might not be in the original list (will receive 'bad' score in table)
  QCM_final.m$value[is.na(QCM_final.m$value)] = 0
  ## some other files might be missing on purpose
  QCM_final.m$value[QCM_final.m$value == HEATMAP_NA_VALUE] = NA
  
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
  ## color names by category
  ##
  col_axis = rep(c("#606060", "black"), 10) ## grey-black as tic-toc
  head(QCM_final.m)
  QCM_final.m$name = colnames_wNum_map$name[match(QCM_final.m$variable, colnames_wNum_map$nameWnum)]
  QCM_final.m$axisCat = tolower(df.meta$cat[match(QCM_final.m$name, df.meta$name)])
  QCM_final.m$axisCol = col_axis[cumsum(!duplicated(QCM_final.m$axisCat))]
  
  QCM_final.m$variable2 = factor(QCM_final.m$variable, levels = unique(QCM_final.m$variable))
  if (any(is.na(QCM_final.m$value))) QCM_final.m$dummy_col = "NA" ## use color legend for missing values
  
  ## replace the x-axis labels by expressions
  QCM_final.m$variable3 = sapply(as.character(QCM_final.m$variable2), function(x) {
    parse(text=x)
  })
  
  p = ggplot(QCM_final.m, aes(y = .data$fc.raw.file, x = .data$variable2))
  if (any(is.na(QCM_final.m$value)))
  {
    p = p + geom_tile(aes(fill = .data$value, colour = .data$dummy_col)) +
      scale_colour_manual(name="Missing", values=c("NA" = "grey50")) +
      guides(colour = guide_legend(override.aes = list(fill = 'grey50')))
  } else {
    p = p + geom_tile(aes(fill = .data$value))
  }  
  p = p + scale_y_discrete_reverse(QCM_final.m$fc.raw.file, breaks = ggAxisLabels) + 
    scale_x_discrete(labels = QCM_final.m$variable3) +
    scale_fill_gradientn(colours = c("red", "black", "green"),
                         na.value = "grey50",
                         limits=c(0, 1), 
                         guide = guide_legend(title = "score"),
                         breaks=c(1,0.5,0),
                         labels=c("best","under\nperforming","fail")) +
    theme(axis.text.x = suppressWarnings(element_text(angle = 90, 
                                     vjust = 0.5, 
                                     hjust = 1, 
                                     colour=QCM_final.m$axisCol[!duplicated(QCM_final.m$variable2)]))) +
    ggtitle("Performance overview") +
    xlab("") +
    ylab("Raw file")
  #print(p)
  return(list(plot = p, table = reshape2::dcast(QCM_final.m, fc.raw.file ~ variable)))
}  

#'
#' Extract meta information (orderNr, metric name, category)
#' from a list of Qc metric objects
#' 
#' @param lst_qcMetrics List of qcMetrics
#' @return data.frame with columns 'name', 'order' and 'cat' (category)
#'
getMetaData = function(lst_qcMetrics)
{
  df.meta = plyr::ldply(lst_qcMetrics, function(qcm) {
    #qq <<- qcm
    qcm_sc = qcm$qcScores
    if (plyr::empty(qcm_sc)) {
      ## if metric was not computed, default DF is empty
      name = qcm$qcName
    } else {
      # Name of metric (incl. filled placeholders)
      name = grepv("raw.file", colnames(qcm_sc), invert = TRUE)
    }
    data.frame(name = name, order = qcm$orderNr, cat = qcm$qcCat)
  })
  ## order meta
  if (!plyr::empty(df.meta)) df.meta = df.meta[order(df.meta$order), ]
  return(df.meta)
}
