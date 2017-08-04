
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
#' @field orderNr [placeholder] column index during heatmap generation and for the general order of plots
#'
#' @import methods
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
#'                  workerFcn=function(.self, data, gtit)
#'                  {
#'                    ## usually some code here to produce ggplots
#'                    pl = lapply(1:2, function(xx) {
#'                        ggplot(data) +
#'                          geom_point(aes(x=x*xx,y=y)) +
#'                          ggtitle(gtit)
#'                      })
#'                    return(list(plots = pl))
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
#' 
#' 
#' 
qcMetric = setRefClass("qcMetric",
                       
   fields = list(helpText = "character",         ## identical to helpTextTemplate by default, but can be modified in workerFcn() if desired
                 helpTextTemplate = "character", ## a template (with placeholders for update at runtime) or a static text
                 workerFcn = "function",  ## returns list(plots =, [qcScores=])
                 plots = "list",
                 htmlTable = "character",
                 ## the following members are related to the heatmap only
                 qcScores = "data.frame", ## with columns "raw.file", "score"
                 qcCat = "character", ## one of "prep", "LC", "MS" or empty (e.g. for PG)
                 qcName = "character", ## expression e.g. "MS^2~ID~Rate"
                 orderNr = "numeric", ## ordering of plots -- also applies to Heatmap; gaps are ignored
                 outData = "list" ## optional auxiliary output data generated in workerFcn
                 ),
   methods = list(
       initialize=function(helpTextTemplate = NA_character_,
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
           .self$qcCat = qcCat;
           .self$qcName = qcName;
           .self$orderNr = orderNr;
           .self$outData = list();
           return(.self)
       },
       checkInput = function(required_columns, given_columns)
       {
         if (!all(required_columns %in% given_columns))
         {
           warning(paste0("Input check failed: columns '", 
                          paste(setdiff(required_columns, given_columns), collapse="', '", sep=""),
                          "' are not present in input data!"),
                   immediate. = TRUE)
           return (FALSE)
         }
         return (TRUE)
       },
       setData = function(...) { ## fill with MQ data and compute results
         cat("Starting to work on", gsub("~", " ", .self$qcName), "...\n")
         if (.self$orderNr < 0)
         {
           cat("  Metric disabled. Skipping...\n")
           return(NULL)
         }
         
         ## GC stats
         mem_before = gc(verbose = FALSE, reset = TRUE) ## resetting Max to see how much this metric uses
         t_before = proc.time()
         
         r = workerFcn(.self, ...)
         
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
           warning(c("Worker of '", .self$qcName, "' returned prematurely! Skipping metric!"))
           return(NULL);         
         }
         
         if (!("plots" %in% names(r))) stop(c("Worker of '", .self$qcName, "' did not return valid result format!"))
         if (!class(r[["plots"]]) == "list") stop(c("Worker of '", .self$qcName, "' did not return plots in list format!"))
         

         lpl = flattenList(r[["plots"]])
         #print(lpl) ## debug
         .self$plots = lpl;
         #.self$plots = r[["plots"]]
         
         if ("htmlTable" %in% names(r)) .self$htmlTable = r[["htmlTable"]];
         if ("qcScores" %in% names(r)) .self$qcScores = r[["qcScores"]];
         
         cat("...", gsub("~", " ", .self$qcName), " done\n")
         return(NULL)
       },
       
       getPlots = function(withTitle = TRUE) {
         if (!withTitle) { ## no title
           r = lapply(.self$plots, function(p) {
             ## delete title 
             if ("ggplot" %in% class(p)) p = p + ggtitle(NULL)
             return (p)
           })
           return(r)
         };
         return(.self$plots)
       },
       
       getTitles = function(stopOnMissing = TRUE, subtitle_sep = " - ") {
         labels = sapply(1:length(.self$plots), 
                         function(idx) {
                           if ("title" %in% names(.self$plots[[idx]]$labels)){
                             title = .self$plots[[idx]]$labels$title
                             #title = 'atop("PG: PCA of 'reporter intensity'", scriptstyle("(excludes contaminants)"))'
                             title
                             regex = "atop\\(\"(.*)\", scriptstyle\\(\"(.*)\"\\)\\)"
                             m = regexpr(regex, title, perl = TRUE)
                             if (m == 1) { ## hit!
                               text = substring(title, attr(m, "capture.start"), attr(m, "capture.start") + attr(m, "capture.length") - 1)
                               title = paste0(text[1], subtitle_sep, text[2])
                               title
                             }
                             return (title)
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
    if(!any(idx_list)) return(x)
    r_list = Reduce(append, x[idx_list])
    items = x[!idx_list]
    x = Reduce(append, list(r_list, items))
  }
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
                      if (empty(df.QCM)) stop("AverageQual_qc::workerFcn(): input empty!")
                      lpl = list() ## empty...
                      qcScore = ddply(df.QCM, "fc.raw.file", function(df.row) {
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

