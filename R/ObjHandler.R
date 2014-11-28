## A proto class to handle a set of objects and provide a simple way of managing them
## (as opposed to using environments or lists directly ...)
##


## CLASS 'ObjHandler'
ObjHandler <- proto()

#' Constructor for class 'ObjHandler'.
#'
#' Simply add items or get them back as an list. 
#'
#' 
#' @name ObjHandler$new
#' @import proto
#' 
ObjHandler$new = function(.)
{
  proto(., items = list(), counter = 0)
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
#' @param name  String used as object name. Automatically generated if omitted.
#' @return Always TRUE (invisible)
#' 
#' @name ObjHandler$add
#' 
# (not exported!)
ObjHandler$add = function(., obj, name = NA)
{
  if (is.na(name)) {
    .$counter = .$counter + 1
    name = sprintf("v%04d", .$counter) ## pad with 0's, to ensure correct order when ordering by name
  }
  .$items[[name]] = obj
  print(obj)
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
#' @name ObjHandler$add
#' 
# (not exported!)
ObjHandler$get = function(., name)
{
  return (.$items[[name]])
}


#' 
#' Get the objects as list
#' 
#' Also see 'addObj()'.
#'
#' @param . A 'this' pointer. Use it to refer/change internal members. It's implicitly added, thus not required too call the function!
#' @return A list holding all objects
#' 
#' @name ObjHandler$getList
#' 
# (not exported!)
ObjHandler$getList = function(.)
{
  return (.$items)
}







