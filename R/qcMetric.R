
#'
#' Class which can compute plots (usually for a single metric).
#'
#' Reference class which is instanciated with a metric description and a
#' worker function (at initialization time, i.e. in the package)
#' and can produce plots (at runtime, when data is provided) using setData().
#'
#' @field helpText Description (lengthy) of the metric and plot elements
#' @field workerFcn Function which generates a result (usually plots). Data is provided using setData().
#' @field plots List of plots (after setData() was called)
#' @field qcScores [placeholder] Data.frame of scores from a qcMetric (computed within workerFcn())
#' @field qcCat [placeholder] QC category (LC, MS, or prep)
#' @field qcName [placeholder] Name of the qcScore in the heatmap
#' @field heatmapOrder [placeholder] column index during heatmap generation
#'
#' @exportClass qcMetric
#' @export qcMetric
#' 
#' @examples 
#'
#' require(ggplot2)
#' dd = data.frame(x=1:10, y=11:20)
#' a = qcMetric$new(helpText="small help text", 
#'                  workerFcn=function(.self, data, gtit) ## arbitrary arguments, matched during setData()
#'                  {
#'                    ## usually some code here to produce ggplots
#'                    pl = lapply(1:2, function(xx) ggplot(data) + geom_point(aes(x=x*xx,y=y)) + ggtitle(gtit))
#'                    return(list(plots = pl))
#'                  }, 
#'                  qcCat="LC", 
#'                  qcName="MS/MS Peak shape", 
#'                  heatmapOrder = 30)
#' ## test some output
#' a$setData(dd, "my title")
#' a$plots  ## the raw plots
#' a$print(TRUE) ## same as above
#' a$print(FALSE) ## plots without title
#' a$getTitles()  ## get the titles of the all plots
#' a$helpText
#' a$qcName
#' 
#' 
#' 
qcMetric = setRefClass("qcMetric",
                       
   fields = list(helpText = "character",
                 workerFcn = "function",  ## returns list(plots =, [qcScores=])
                 plots = "list",
                 ## the following members are related to the heatmap only
                 qcScores = "data.frame", ## with columns "raw.file", "score"
                 qcCat = "character", ## one of "prep", "LC", "MS" or empty (e.g. for PG)
                 qcName = "character", ## expression e.g. "MS^2~ID~Rate"
                 heatmapOrder = "numeric" ## column index in heatmap (gaps are ignored)
                 ),
   methods = list(
       initialize=function(helpText,
                           workerFcn,
                           qcCat = NA_character_,
                           qcName = NA_character_,
                           heatmapOrder = NaN) {
           .self$helpText = helpText;
           .self$workerFcn = workerFcn;
           .self$plots = list();  ## obtained from worker
           .self$qcScores = data.frame();  ## obtained from worker
           .self$qcCat = qcCat;
           .self$qcName = qcName;
           .self$heatmapOrder = heatmapOrder;
           .self
       },
       setData = function(...) { ## fill with MQ data and compute results
         
         r = workerFcn(.self, ...)

         if (!("plots" %in% names(r))) stop(c("Worker of '", .self$qcName, "' did not return valid result format!"))
         if (!class(r[["plots"]]) == "list") stop(c("Worker of '", .self$qcName, "' did not return plots in list format!"))
         .self$plots = r[["plots"]];
         
         if ("qcScores" %in% names(r)) .self$qcScores = r[["qcScores"]];
         
         ##if ("") TODO: extract title?!
         #l = list(...)
         
         return(NULL)
       },
       
       print = function(withTitle = TRUE) {
         if (!withTitle) {
           r = lapply(.self$plots, function(p) p + ggtitle(NULL))
           return(r)
         };
         return(.self$plots)
       },
       
       getTitles = function(stopOnMissing = TRUE) {
         labels = sapply(1:length(.self$plots), 
                         function(idx) {
                           if ("title" %in% names(.self$plots[[idx]]$labels)){
                             return (.self$plots[[idx]]$labels$title)
                           } else if (stopOnMissing) {
                             stop(c("getTitles(): No title found in ggplot object at index ", idx, "!"))
                           } else return("")
                         })
         return(labels)
       }
         
     )
   ) ## refClass






