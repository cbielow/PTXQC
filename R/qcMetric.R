
#'
#' Class which can compute plots and generate mzQC output (usually for a single metric).
#'
#' Reference class which is instanciated with a metric description and a
#' worker function (at initialization time, i.e. in the package)
#' and can produce plots and mzQC values (at runtime, when data is provided) using setData().
#' 
#' All derived classes need to implement a 'workerFcn()' function, which returns a list with
#' elements: c("plots", "mzQC", "htmlTable", "qcScores", "title"),
#' where 'plots' is required; all others are optional.
#'
#' @field helpText  Description (lengthy) of the metric and plot elements
#' @field workerFcn Function which generates a result (usually plots). Data is provided using setData().
#' @field plots     List of plots (after setData() was called)
#' @field htmlTable A table for display in the HTML report (preferred over a plot in Html mode)
#' @field qcScores  Data.frame of scores from a qcMetric (computed within workerFcn())
#' @field mzQC      An named list of mzQC MzQCqualityMetric's (named by their fc.raw.file for runQuality or concatenated fc.raw.files for setQualities (e.g. "file 1;file4")) (valid after setData() was called)
#' @field qcCat     QC category (LC, MS, or prep)
#' @field qcName    Name of the qcScore in the heatmap
#' @field orderNr   Column index during heatmap generation and for the general order of plots
#'
#' @exportClass qcMetric
#' @export qcMetric
#' 
#' @examples 
#'
#' require(ggplot2)
#' dd = data.frame(x=1:10, y=11:20)
#' a = qcMetric$new(helpText="small help text", 
#'                  ## arbitrary arguments, matched during setData()
#'                  workerFcn=function(.self, data, gtitle)
#'                  {
#'                    ## usually some code here to produce ggplots
#'                    pl = lapply(1:2, function(xx) {
#'                        ggplot(data) +
#'                          geom_point(aes(x=x*xx,y=y)) +
#'                          ggtitle(gtitle)
#'                      })
#'                      ## add mzQC metric for count of identified clusters
#'                      template_proteinCount = rmzqc::getQualityMetricTemplate("MS:1002406") 
#'                      mzqc = lapply(1:3, function(id){
#'                        out = template_proteinCount$copy();
#'                        out$value = id;
#'                        return(out) })
#'                      names(mzqc) = paste0("file", 1:3);
#'                    return(list(plots = pl, mzQC = mzqc))
#'                  }, 
#'                  qcCat="LC", 
#'                  qcName="MS/MS Peak shape", 
#'                  orderNr = 30)
#' ## test some output
#' a$setData(dd, "my title")
#' a$plots  ## the raw plots
#' a$getPlots(TRUE) ## same as above
#' a$getPlots(FALSE) ## plots without title
#' a$getTitles()  ## get the titles of the all plots
#' a$helpText
#' a$qcName
#' a$mzQC
#' 
#' 
qcMetric = setRefClass("qcMetric",
                       
   fields = list(helpText = "character",         ## identical to helpTextTemplate by default, but can be modified in workerFcn() if desired
                 helpTextTemplate = "character", ## a template (with placeholders for update at runtime) or a static text
                 workerFcn = "function",  ## returns list(plots =, [qcScores=])
                 plots = "list",
                 htmlTable = "character",
                 title = "list",
                 ## the following members are related to the heatmap only
                 qcScores = "data.frame", ## with columns "raw.file", "score"
                 mzQC = "list", ## of MzQCbaseQuality objects
                 qcCat = "character", ## one of "prep", "LC", "MS" or empty (e.g. for PG)
                 qcName = "character", ## expression e.g. "MS^2~ID~Rate"
                 orderNr = "numeric", ## ordering of plots -- also applies to Heatmap; gaps are ignored
                 outData = "list" ## optional auxiliary output data generated in workerFcn
                 ),
   methods = list(
       initialize = function(helpTextTemplate = NA_character_,
                           workerFcn = function(){},
                           qcCat = NA_character_,
                           qcName = NA_character_,
                           orderNr = NaN) {
           .self$helpTextTemplate = helpTextTemplate;
           .self$helpText = helpTextTemplate;
           .self$workerFcn = workerFcn;
           .self$plots = list();  ## obtained from worker
           .self$htmlTable = NA_character_;  ## obtained from worker
           .self$qcScores = data.frame();  ## obtained from worker
           .self$mzQC = list();  ## obtained from worker
           .self$title = list();
           .self$qcCat = qcCat;
           .self$qcName = qcName;
           .self$orderNr = orderNr;
           .self$outData = list();
           return(.self)
       },
       #' @description
       #' Internally calls the workerFcn() , which computes the actual plots metric scores and supporting data (e.g. mzQC metrics) of the derived class; the resulting data is checked and stored in the members of this class
       #' @param df The expected data, usually a data frame. If empty, this function will return immediately without failure.
       #' @param ... Additional arguments passed to the workerFcn()
       #' @return NULL
       setData = function(df, ...) { 
         cat("Starting to work on", gsub("~", " ", .self$qcName), "...\n")
         if (.self$orderNr < 0)
         {
           cat("  Metric disabled. Skipping...\n")
           return (NULL)
         }
         
         if (is.null(df))
         {
           cat("  No data available. Skipping...\n")
           return (NULL)
         }
         
         ## GC stats
         mem_before = gc(verbose = FALSE, reset = TRUE) ## resetting Max to see how much this metric uses
         t_before = proc.time()
         
         r = workerFcn(.self, df, ...)
         
         ## clean memory to get a clean picture of each metrics memory footprint
         ## to enable the user to skip expensive metrics
         mem_after = gc(verbose = FALSE)
         cat(paste0("\nMemory [MB] prior|after|max(diff) : ", 
                    sum(mem_before[, 2]), " | ",
                    sum(mem_after[, 2]), " | ",
                    sum(mem_after[, 6]), " (",
                    round(sum(mem_after[, 6])-sum(mem_before[, 2]),0) ,")\n",
                    "\nDuration: ", round((proc.time() - t_before)[3]), " s\n\n"))
         if (is.null(r)) {
           message(c("Worker of '", .self$qcName, "' returned prematurely due to missing data! Skipping metric!"))
           return(NULL);         
         }
         
         valid_elements = c("plots", "mzQC", "htmlTable", "qcScores", "title")
         required_elements = c("plots")
         
         if (!(all(required_elements %in% names(r)))) stop(c("Worker of '", .self$qcName, "' did not return all required fields (", paste(required_elements, collapse = ","), ")!"))
         if (!inherits(r[["plots"]], "list")) stop(c("Worker of '", .self$qcName, "' did not return plots in list format!"))
         
         if (!all(names(r) %in% valid_elements)) stop(c("Worker of '", .self$qcName, "' return invalid fields (", paste(setdiff(names(r), valid_elements), collapse = ","), ")!"))
         

         lpl = flattenList(r[["plots"]])
         #print(lpl) ## debug
         .self$plots = lpl;
         #.self$plots = r[["plots"]]
         
         if ("mzQC" %in% names(r)) .self$mzQC = r[["mzQC"]];
         if ("htmlTable" %in% names(r)) .self$htmlTable = r[["htmlTable"]];
         if ("qcScores" %in% names(r)) .self$qcScores = r[["qcScores"]];
         if ("title" %in% names(r)) .self$title = r[["title"]]
         
         
         cat("...", gsub("~", " ", .self$qcName), " done\n")
         return(NULL)
       },
       
       getPlots = function(withTitle = TRUE) {
         if (!withTitle) { ## no title
           r = lapply(.self$plots, function(p) {
             ## delete title 
             if (inherits(p, "ggplot")) p = p + ggtitle(NULL)
             return (p)
           })
           return(r)
         };
         return(.self$plots)
       },
       
       getTitles = function(stopOnMissing = TRUE, subtitle_sep = " - ") {
         labels = sapply(1:length(.self$plots), 
                         function(idx) {
                           if (length(.self$title) != 0){
                             return(.self$title[[idx]])
                           }
                           else if ("title" %in% names(.self$plots[[idx]]$labels)){
                             titles = .self$plots[[idx]]$labels$title
                             #title = 'atop("PG: PCA of 'reporter intensity'", scriptstyle("(excludes contaminants)"))'
                             regex = 'atop\\("(.*)", scriptstyle\\("(.*)"\\)\\)'
                             m = regexpr(regex, titles, perl = TRUE)
                             if (m == 1) { ## hit!
                               text = substring(titles, attr(m, "capture.start"), attr(m, "capture.start") + attr(m, "capture.length") - 1)
                               return(paste0(text[1], subtitle_sep, text[2]))
                             }
                             return (titles)
                           } else if (stopOnMissing) {
                             stop(c("getTitles() for ", .self$qcName, ": No title found in ggplot object at index ", idx, "!"))
                           } else return("")
                         })
         return(labels)
       }
         
     )
   ) ## refClass

#########################################################################################


#'
#' Flatten lists of lists with irregular depths to just a list of items,
#' i.e. a list of the leaves (if you consider the input as a tree).
#' 
#' @param x List of 'stuff' (could be lists or items or a mix)
#' @return A flat list
#'
#'
#'
flattenList = function(x) {
  repeat {
    idx_list = sapply(x, function(arg) {return(all(class(arg) == "list"))})
    if (!any(idx_list)) return(x)
    r_list = Reduce(append, x[idx_list])
    items = x[!idx_list]
    x = Reduce(append, list(r_list, items))
  }
}

checkInput = function(required_columns, given_df)
{
  if (!all(required_columns %in% colnames(given_df)))
  {
    warning(paste0("Input check failed: columns '", 
                   paste(setdiff(required_columns, colnames(given_df)), collapse="', '", sep=""),
                   "' are not present in input data!"),
            immediate. = TRUE)
    return (FALSE)
  }
  return (TRUE)
}

qcMetric_AverageQualOverall = 
  setRefClass("qcMetric_AverageQualOverall",
              contains = "qcMetric",
              methods = list(
                initialize=function() {
                  callSuper(
                    helpTextTemplate = "Internal metric to compute the average quality across all other metrics",
                    workerFcn = function(.self, df.QCM)
                    {
                      if (plyr::empty(df.QCM)) stop("AverageQual_qc::workerFcn(): input empty!")
                      lpl = list() ## empty...
                      qcScore = plyr::ddply(df.QCM, "fc.raw.file", function(df.row) {
                        df.row.raw = unlist(df.row[,!grepl("fc.raw.file", colnames(df.row))])
                        df.row.raw[is.infinite(df.row.raw)] = NA  ## mask explicitly missing values, since it will bias the mean otherwise
                        return (data.frame(val = mean(df.row.raw, na.rm = TRUE)))
                      })
                      colnames(qcScore)[colnames(qcScore)=="val"] = .self$qcName
                      
                      return(list(plots = lpl, qcScores = qcScore))
                    },
                    qcCat = "General",
                    qcName = "Average~Overall~Quality",
                    orderNr = 9999
                  )
                  return(.self)
                })
              )

#qcA = qcMetric_AverageQualOverall_$new()
#qcA$setData(data.frame(fc.raw.file = letters[1:5], qual1 = 1:5, qual2 = 5:9))
#qcA$qcScores
#qcA$helpText

