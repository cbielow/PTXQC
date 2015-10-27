
setClass("QCMetric", representation(plotObjs = "list", name = "character"), 
         prototype(name = NA_character_, plotObjs = list()))
asMarkUp <- function(object) 0
setGeneric("asMarkUp")
setMethod("asMarkUp", signature(object = "QCMetric"), function(object) {
  ## meant for ggplot objects
  print(object@plotObjs)
})
m = new("QCMetric", name = "MS/MS")
m@plotObjs



## A proto class to handle a set of objects and provide a simple way of managing them
## (as opposed to using environments or lists directly ...)
##


## CLASS 'ReportData'
ReportData <- proto()

#' Constructor for class 'ReportData'.
#'
#' Simply add items or get them back as an list. 
#'
#' 
#' @name ReportData$new
#' @import proto
#' 
ReportData$new = function(.)
{
  proto(., counter = 0, heatmap = NULL, params = NULL, name_mapping = NULL, metrics = list(),
        allowed_types = c("metric", "heatmap", "params", "name_mapping"))
}


##
## Functions
##

#' Add objects
#' 
#' Also see 'getList'.
#'
#' @param .      A 'this' pointer. Use it to refer/change internal members. It's implicitly added, thus not required too call the function!
#' @param obj   Value to be inserted
#' @param type  Type of object: metric, heatmap, params, name_mapping
#' @param name  String used as object name. Automatically generated if omitted.
#' @return Always TRUE (invisible)
#' 
#' @name ReportData$add
#' 
# (not exported!)
ReportData$add = function(., obj, type = "metric", name = NA)
{
  if (!type %in% .$allowed_types) stop("wrong 'type'")
  
  if (is.na(name))
  {
    .$counter = .$counter + 1
    name = sprintf("v%04d", .$counter) ## pad with 0's, to ensure correct order when ordering by name
  }
  
  if (type == "metric")
  {
    .$metrics[[name]] = obj
  } else if (type == "heatmap")
  {
    .$heatmap = obj
  } else if (type == "params")
  {
    .$params = obj
  } else if (type == "name_mapping")
  {
    .$name_mapping = obj
  }
  #print(obj)
  return (invisible(TRUE))
}


#' Get an object
#' 
#' Also see 'getList'.
#'
#' @param . A 'this' pointer. Use it to refer/change internal members. It's implicitly added, thus not required too call the function!
#' @param name  Name of list entry to return
#' @return The object
#' 
#' @name ReportData$add
#' 
# (not exported!)
ReportData$get = function(., type = "metric", name = NA)
{

  if (!type %in% .$allowed_types) stop("wrong 'type'")
  
  if (type == "metric")
  {
    if (is.na(name))
    {
      return (.$metrics)
    }
    return (.$metrics[[name]])
  } else if (type == "heatmap")
  {
    return (.$heatmap)
  } else if (type == "params")
  {
    return (.$params)
  } else if (type == "name_mapping")
  {
    return (.$name_mapping)
  }
}






