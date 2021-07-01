## some generic information on why/how/what?? :)
## 
## + we provide initialize() functions for all RefClasses to enable unnamed construction (shorter syntax)
##
##
##
##
##



##
# Defining this function to enable overload
# e.g. setMethod('asJSON', 'mzQC', function(x, ...) x$toJSON())
#  which allows to use
# jsonlite::toJSON(mzQC$new(content))
asJSON <- jsonlite:::asJSON


#'
#' Tell if a string is undefined (NA or NULL); If yes, and its required by the mzQC standard, we can raise an error
#'
#' You can pass multiple strings, which are all checked. If any of them is undefined, the function returns TRUE
#'
#' @param s A string to be checked for NA/NULL
#' @param ... More strings to be checked
#' @param verbose If TRUE and 's' is NULL/NA, will print the name of the variable which was passed in
#'
#' @examples 
#' isUndefined(NA)       ## TRUE
#' isUndefined(NULL)     ## TRUE
#' isUndefined(NA, NULL) ## TRUE
#' isUndefined("")       ## FALSE
#' isUndefined("", NA)   ## TRUE
#' isUndefined(1)        ## FALSE
#' myVar = NA
#' isUndefined(myVar)    ## TRUE, with warning "variable 'myVar' is NA/NULL!"
#'
#' @export
#'   
isUndefined = function(s, ..., verbose = TRUE)
{
  # anchor
  if (missing(s)) return(FALSE)
  
  r = (is.na(s) || is.null(s))
  name_of_var = deparse(substitute(s))
  # omit the '.self' part of the variable's name
  name_of_var = gsub("^.self\\$", "", name_of_var)
  if (verbose && r) warning(paste0("Variable '", name_of_var, "' is NA/NULL!"), immediate. = TRUE, call. = FALSE)
  ## check remaining args from ... by using '+' (force evaluation)
  return(r + isUndefined(..., verbose = verbose) > 0)
}


#'
#' Checks validity (= completeness) of mzQC objects - or lists (JSON arrays) thereof
#'
#' Note: Returns TRUE for empty lists!
#'
#' You can pass multiple arguments, which are all checked individually.
#' All of them need to be valid, for TRUE to be returned.
#' The reason for combining list support for arguments and ellipsis (...) into this function is that
#' JSON arrays are represented as lists and you can simply pass them as a single argument 
#' (without the need for do.call()) and get the indices of invalid objects (if any).
#' The ellipsis is useful to avoid clutter, 
#' i.e. 
#'      if (!isValidMzQC(a) || !isValidMzQC(b)) doStuff()
#'      is harder to read than
#'      if (!isValidMzQC(a,b)) doStuff()
#'
#' @param x An mzQC refclass (or list of them), each will be subjected to `isValidMzQC()`
#' @param ... Ellipsis, for recursive argument splitting
#'
#' @examples 
#'   isValidMzQC(MzQCcvParameter$new("QC:4000059"))       # FALSE
#'   isValidMzQC(MzQCcvParameter$new("QC:4000059", "Number of MS1 spectra")) # TRUE
#'   isValidMzQC(list(MzQCcvParameter$new("QC:4000059"))) # FALSE
#'   isValidMzQC(list(MzQCcvParameter$new("QC:4000059", "Number of MS1 spectra"))) # TRUE
#'   isValidMzQC(list(MzQCcvParameter$new("QC:4000059", "Number of MS1 spectra")), 
#'               MzQCcvParameter$new()) # FALSE
#'   
#' @export
#'
isValidMzQC = function(x, ...)
{
  # anchor
  if (missing(x)) return(TRUE)
  
  if ("list" %in% class(x)) {
    idx = sapply(x, isValidMzQC)
    if (any(idx == FALSE)) {
      warning(paste0("In list of '", class(x[[1]]), "', the element(s) #[", paste(which(idx == FALSE), collapse = ","), "] is/are invalid."), immediate. = TRUE, call. = FALSE)
    }
    return(all(idx) & isValidMzQC(...))
  }
  r = x$isValid()
  if (r == FALSE)
  {
    warning(paste0("A field in object of type ", class(x), " is invalid."), immediate. = TRUE, call. = FALSE)
  }
  return(r & isValidMzQC(...))
}


#'
#' Allow conversion of plain named lists to mzQC objects
#' 
#' The plain-R representation of your mzQC objects must be wrapped in an outer list,
#' if your mzQC object representation is already a list
#' because upon detecting lists, this function will call 'class$fromData(element)' for every element.
#'
#' @param mzqc_class Prototype of the class to convert 'data' into
#' @param data A datastructure of R lists/arrays as obtained by 'jsonlite::fromJSON()'
#'
#' @examples 
#'  data = MzQCcvParameter$new("acc", "myName", "value")
#'  data_recovered = fromDatatoMzQC(MzQCcvParameter, list(jsonlite::fromJSON(jsonlite::toJSON(data))))
#'  data_recovered
#'
#' @export
#'   
fromDatatoMzQC = function(mzqc_class, data)
{
  if ("list" %in% class(data))
  {
    return(sapply(data, function(x) {
      obj = mzqc_class$new()
      obj$fromData(x)
      obj
    }))
  }
  if (is.na(data) || is.null(data)) return(list())
  obj = mzqc_class$new()
  return(obj$fromData(data))
}

#'
#' converts a NULL to NA_character_; or returns the argument unchanged otherwise
#' 
#' This is useful for missing list elements (which returns NULL), 
#' but when the missing element in refClass should be NA_character_ (and NULL would return an error)
#'
#' @param char_or_NULL A string or NULL
#'
#' @examples 
#'   NULL_to_charNA(NA)   ## NA
#'   NULL_to_charNA(NULL) ## NA_character_
#'   NULL_to_charNA("hi") ## "hi"
#'   
#' @export
#'   
NULL_to_charNA = function(char_or_NULL) {
  if (is.null(char_or_NULL)) return(NA_character_)
  return(char_or_NULL)
}


#'
#' converts a NULL to NA; or returns the argument unchanged otherwise
#' 
#' This is useful for missing list elements (which returns NULL), 
#' but when the missing element in refClass should be NA (and NULL would return an error)
#'
#' @param var_or_NULL A variable of any kind or NULL
#'
#' @examples 
#'   NULL_to_NA(NA)   ## NA
#'   NULL_to_NA(NULL) ## NA
#'   NULL_to_NA("hi") ## "hi"
#'   
#' @export
#'   
NULL_to_NA = function(var_or_NULL) {
  if (is.null(var_or_NULL)) return(NA)
  return(var_or_NULL)
}


#'
#' An mzQC-formatted date+time, as required by the mzQC spec doc
#' 
#' The format is "%Y-%m-%d %H:%M:%S".
#' 
#' @field datetime A correctly formatted date time (use as read-only)
#'
#' @exportClass MzQCDateTime
#' @export MzQCDateTime
#' 
#' @examples
#'    dt1 = MzQCDateTime$new("1900-01-01")
#'    dt2 = MzQCDateTime$new(Sys.time())
#'    ## test roundtrip conversion from/to JSON
#'    dt2$fromData(jsonlite::fromJSON(jsonlite::toJSON(dt1)))
#     dt1$datetime == dt2$datetime    ## TRUE
#'    
#' @export
#'   
MzQCDateTime = setRefClass(
  'MzQCDateTime',
  fields = list(datetime = 'character'),
  methods = list(
    initialize = function(date = as.character(Sys.time()))
    {
      set(date)
    },
    set = function(.self, date)
    {
      .self$datetime = format(as.POSIXct(date), "%Y-%m-%dT%H:%M:%S")  # using ISO8601 format
    },
    isValid = function(.self)
    {
      return(TRUE) ## always valid, because it's designed that way
    },
    toJSON = function(.self, ...)
    {
      if (!isValidMzQC(.self)) stop(paste0("Object of class '", class(.self), "' is not in a valid state for writing to JSON"))
      return(jsonlite::toJSON(.self$datetime))
    },
    fromData = function(.self, data)
    {
      .self$set(data)
    }
  )
)
setMethod('asJSON', 'MzQCDateTime', function(x, ...) x$toJSON(...))


#'
#' An controlled vocabulary document, usually pointing to an .obo file
#'
#' @field name Full name of the controlled vocabulary.
#' @field uri Publicly accessible URI of the controlled vocabulary. 
#' @field version [optional] Version of the controlled vocabulary.
#' 
#' @export MzQCcontrolledVocabulary
#' 
#' @examples 
#'   MzQCcontrolledVocabulary$new(
#'     "Proteomics Standards Initiative Quality Control Ontology",
#'     "https://github.com/HUPO-PSI/mzQC/blob/master/cv/qc-cv.obo",
#'     "1.2.0")
#'    
#' @export
#'   
MzQCcontrolledVocabulary = setRefClass(
  'MzQCcontrolledVocabulary',
  fields = list(name = 'character',
                uri = 'character',
                version = 'character'    # optional
  ),
  methods = list(
    initialize = function(name = NA_character_, uri = NA_character_, version = NA_character_)
    {
      .self$name = name
      .self$uri = uri
      .self$version = version
    },
    isValid = function(.self) {
      if (isUndefined(.self$name, .self$uri)) return(FALSE)
      return(TRUE)
    },
    toJSON = function(.self, ...)
    {
      if (!isValidMzQC(.self)) stop(paste0("Object of class '", class(.self), "' is not in a valid state for writing to JSON"))
      
      r = list("name" = .self$name,
               "uri" = .self$uri,
               "version" = .self$version)
      return (jsonlite::toJSON(r, auto_unbox = TRUE))
    },
    fromData = function(.self, data)
    {
      .self$name = data$name
      .self$uri = data$uri
      .self$version = NULL_to_charNA(data$version)
    }
  )
)
setMethod('asJSON', 'MzQCcontrolledVocabulary', function(x, ...) x$toJSON(...))


#'
#' An controlled vocabulary parameter, as detailed in the OBO file
#'
#' @field accession Accession number identifying the term within its controlled vocabulary (pattern: ^[A-Z]+:[A-Z0-9]+$).
#' @field name Name of the controlled vocabulary term describing the parameter.
#' @field value [optional] Value of the parameter.
#' @field description [optional] Definition of the controlled vocabulary term.
#' 
#' @export MzQCcvParameter
#' 
#' @examples 
#'   MzQCcvParameter$new("QC:4000139",
#'                       "RT acquisition range",
#'                       c(0.2959, 5969.8172))
#'   isValidMzQC(MzQCcvParameter$new("MS:0000000"))
#'    
#' @export
#'   
MzQCcvParameter = setRefClass(
  'MzQCcvParameter',
  fields = list(accession = 'character',
                name = 'character',
                value = 'ANY',              # optional
                description = 'character'   # optional
  ),
  methods = list(
    initialize = function(accession = NA_character_, name = NA_character_, value = NA, description = NA_character_)
    {
      .self$accession = accession
      .self$name = name
      .self$value = value
      .self$description = description
    },
    isValid = function(.self) {
      if (isUndefined(.self$accession, .self$name)) return(FALSE)
      return(TRUE)
    },
    toJSON = function(.self, ...)
    {
      if (!isValidMzQC(.self)) stop(paste0("Object of class '", class(.self), "' is not in a valid state for writing to JSON"))
      
      r = list("accession" = .self$accession,
               "name" = .self$name)
      if (!is.na(.self$description)) r["description"] = .self$description
      if (!is.na(.self$value)) r["value"] = .self$value
      return (jsonlite::toJSON(r, auto_unbox = TRUE))
    },
    fromData = function(.self, data)
    {
      .self$accession = data$accession
      .self$name = data$name
      .self$description = NULL_to_charNA(data$description)
      .self$value = NULL_to_NA(data$value)
    }
  )
)
setMethod('asJSON', 'MzQCcvParameter', function(x, ...) x$toJSON(...))

#'
#' An inputfile within metadata for a run/setQuality
#'
#' @field name The name MUST uniquely match to a location (specified below) listed in the mzQC file.
#' @field location Unique file location, REQUIRED to be specified as a URI. The file URI is RECOMMENDED to be publicly accessible.
#' @field fileFormat A MzQCcvParameter with 'accession' and 'name'.
#' @field fileProperties An array of MzQCcvParameter, usually with 'accession', 'name' and 'value'. Recommended are at least two entries: 
#'        a) Completion time of the input file (MS:1000747) and b) Checksum of the input file (any child of: MS:1000561 ! data file checksum type).
#'
#' @export MzQCinputFile
#' 
#'    
#' @export
#'   
MzQCinputFile = setRefClass(
  'MzQCinputFile',
  fields = list(name = 'character',
                location = 'character',
                fileFormat = 'MzQCcvParameter',
                fileProperties = 'list'         # array of MzQCcvParameter, optional
  ),
  methods = list(
    # defaults are required, otherwise refClasses do not work.
    initialize = function(name = NA_character_, location = NA_character_, fileFormat = MzQCcvParameter$new(), fileProperties = list())
    {
      .self$name = name
      .self$location = location
      .self$fileFormat = fileFormat
      .self$fileProperties = fileProperties
    },
    isValid = function(.self)
    {
      # force evaluation of all fields by '+'
      if (isUndefined(.self$name, .self$location) + !.self$fileFormat$isValid()) return(FALSE)
      return(isValidMzQC(.self$fileProperties)) ## TRUE for empty list, which is ok
    },
    toJSON = function(.self, ...)
    {
      if (!isValidMzQC(.self)) stop(paste0("Object of class '", class(.self), "' is not in a valid state for writing to JSON"))
      # no need to check if optional field fileProperties is present. It will be an (empty) JSON array, which is what we want
      return (jsonlite::toJSON(list(name = .self$name, location = .self$location, fileFormat = .self$fileFormat, fileProperties = .self$fileProperties)))
    },
    fromData = function(.self, data)
    {
      .self$name = data$name
      .self$location = data$location
      .self$fileFormat$fromData(data$fileFormat)
      .self$fileProperties = fromDatatoMzQC(MzQCcvParameter, data$fileProperties) ## for lists, call the free function
    }
  )
)
setMethod('asJSON', 'MzQCinputFile', function(x, ...) x$toJSON(...))

# 
# file_format = MzQCcvParameter$new("MS:1000584", "mzML format")
# nif = MzQCinputFile$new("tmp.mzML", "c:\\", file_format)
# nif
# nif2 = nif
# l2 = list(file_format, file_format)
# nif2$fileProperties = l2
# x = jsonlite::toJSON(nif, pretty = TRUE)
# x
# x2 = jsonlite::toJSON(nif2)
# xdata = jsonlite::fromJSON(x, simplifyDataFrame = FALSE)
# xdata
# class(fromDatatoMzQC(MzQCcvParameter, xdata$fileProperties)) == "list"
# jsonlite::toJSON(xdata, pretty = TRUE, auto_unbox = T)
# isValidMzQC(l2)
# nif$fromData(xdata)




#'
#' Details of the software used to create the QC metrics
#'
#' @field accession Accession number identifying the term within its controlled vocabulary (pattern: ^[A-Z]+:[A-Z0-9]+$).
#' @field name Name of the controlled vocabulary term describing the software tool.
#' @field version Version number of the software tool.
#' @field uri Publicly accessible URI of the software tool or documentation.
#' @field description [optional] Definition of the controlled vocabulary term.
#' @field value [optional] Value of the software tool.
#'
#' @export MzQCanalysisSoftware
#'
#'
#' @export
#'   
MzQCanalysisSoftware = setRefClass(
  'MzQCanalysisSoftware',
  fields = list(accession = 'character',
                name = 'character',
                version = 'character',
                uri = 'character',
                description = 'character',  # optional
                value = 'character'         # optional
  ),
  methods = list(
    # defaults are required, otherwise refClasses do not work.
    initialize = function(accession = NA_character_, 
                          name = NA_character_, 
                          version = NA_character_, 
                          uri = NA_character_, 
                          description = NA_character_, ## optional
                          value = NA_character_        ## optional
                          )
    {
      .self$accession = accession
      .self$name = name
      .self$version = version
      .self$uri = uri
      .self$description = description
      .self$value = value
    },
    isValid = function(.self)
    {
      if (isUndefined(.self$accession, .self$name, .self$version, .self$uri)) return(FALSE)
      return(TRUE)
    },
    toJSON = function(.self, ...)
    {
      if (!isValidMzQC(.self)) stop(paste0("Object of class '", class(.self), "' is not in a valid state for writing to JSON"))
      
      r = list("accession" = .self$accession,
               "name" = .self$name,
               "version" = .self$version,
               "uri" = .self$uri)
      if (!isUndefined(.self$description)) r$description = .self$description
      if (!isUndefined(.self$value)) r$value = .self$value
      return (jsonlite::toJSON(r))
    },
    fromData = function(.self, data)
    {
      .self$accession = data$accession
      .self$name = data$name
      .self$version = data$version
      .self$uri = data$uri
      .self$description = NULL_to_charNA(data$description)
      .self$value = NULL_to_charNA(data$value)
    }
  )
)
setMethod('asJSON', 'MzQCanalysisSoftware', function(x, ...) x$toJSON(...))



#'
#'The metadata for a run/setQuality
#'
#' @field label Unique name for the run (for runQuality) or set (for setQuality).
#' @field inputFiles Array/list of MzQCinputFile objects 
#' @field analysisSoftware Array/list of MzQCanalysisSoftware objects 
#' @field cvParameters [optional] Array of cvParameters objects 
#'
#' @export MzQCmetadata
#' 
#' @export
#'   
MzQCmetadata = setRefClass(
  'MzQCmetadata',
  fields = list(label = 'character',
                inputFiles = 'list',       # array of MzQCinputFile
                analysisSoftware = 'list', # array of MzQCanalysisSoftware
                cvParameters = 'list'      # optional array of MzQCcvParameter
                ),
  methods = list(
    initialize = function(label = NA_character_, inputFiles =list(), analysisSoftware = list(), cvParameters = list())
    {
      .self$label = label
      .self$inputFiles = inputFiles
      .self$analysisSoftware = analysisSoftware
      .self$cvParameters = cvParameters
    },
    isValid = function(.self)
    {
      # force evaluation of all fields by '+'
      if (isUndefined(.self$label) + !isValidMzQC(.self$inputFiles, .self$analysisSoftware, .self$cvParameters)) return(FALSE)
      return(TRUE)
    },
    toJSON = function(.self, ...)
    {
      if (!isValidMzQC(.self)) stop(paste0("Object of class '", class(.self), "' is not in a valid state for writing to JSON"))
      
      r = list("label" = .self$label,
               "inputFiles" = .self$inputFiles,
               "analysisSoftware" = .self$analysisSoftware,
               "cvParameters" = .self$cvParameters)  ## might yield an empty JSON array,but ok
      return (jsonlite::toJSON(r))
    },
    fromData = function(.self, data)
    {
      .self$label = data$label
      .self$inputFiles = fromDatatoMzQC(MzQCinputFile, data$inputFiles)
      .self$analysisSoftware = fromDatatoMzQC(MzQCanalysisSoftware, data$analysisSoftware)
      .self$cvParameters = fromDatatoMzQC(MzQCcvParameter, data$cvParameters)
    }
  )
)
setMethod('asJSON', 'MzQCmetadata', function(x, ...) x$toJSON(...))

################################################################################################################################
#################################################################################################################################'
#' The central class to store QC information
#'
#' @field accession Accession number identifying the term within its controlled vocabulary (pattern: ^[A-Z]+:[A-Z0-9]+$).
#' @field name Name of the controlled vocabulary element describing the metric.
#' @field description [optional] Definition of the controlled vocabulary term.
#' @field value [optional] Value of the metric (single value, n-tuple, table, matrix).
#'        The structure is not checked by our mzQC implementation and must be handled by the caller
#' @field unit [optional] Array of unit(s), stored as MzQcvParameter
#'
#' @export MzQCqualityMetric
#' 
MzQCqualityMetric = setRefClass(
  'MzQCqualityMetric',
  fields = list(accession = 'character',
                name = 'character',      
                description = 'character', # optional
                value = 'ANY',             # optional value of unspecified type
                unit = 'list'              # optional array of MzQCcvParameter
  ),
  methods = list(
    initialize = function(accession = NA_character_, name = NA_character_, description = NA_character_, value = NA, unit = list())
    {
      .self$accession = accession
      .self$name = name
      .self$description = description
      if (!missing(value)) .self$value = value else .self$value = NA  ## need to set as NA explicitly, because the default value 'uninitialized class ANY' cannot be converted to JSON
      .self$unit = unit
    },
    isValid = function(.self)
    {
      if (isUndefined(.self$accession, .self$name)) return(FALSE)
      return(TRUE)
    },
    toJSON = function(.self, ...)
    {
      if (!isValidMzQC(.self)) stop(paste0("Object of class '", class(.self), "' is not in a valid state for writing to JSON"))
      
      r = list("accession" = .self$accession,
               "name" = .self$name,
               "description" = .self$description,
               "value" = .self$value, ## NA is written as "value": [null] and read back as NA
               "unit" = .self$unit)  ## might yield an empty JSON array, but ok
      return (jsonlite::toJSON(r))
    },
    fromData = function(.self, data)
    {
      .self$accession = data$accession
      .self$name = data$name
      .self$description = NULL_to_charNA(data$description)
      if (!is.na(data$value)) .self$value = data$value
      .self$unit = fromDatatoMzQC(MzQCcvParameter, data$unit) ## if data$unit is empty, or NA, the empty list will be returned
    }
  )
)
setMethod('asJSON', 'MzQCqualityMetric', function(x, ...) x$toJSON(...))


#a_qc_metric = MzQCqualityMetric$new("acc", "nnam")
#xq = jsonlite::toJSON(a_qc_metric)
#jsonlite::fromJSON(xq)


#'
#' Base class of runQuality/setQuality
#' 
#' @field metadata The metadata for this run/setQuality
#' @field qualityMetrics Array of MzQCqualityMetric objects 
#'
#' @export MzQCbaseQuality
#' 
#'    
MzQCbaseQuality = setRefClass(
  'MzQCbaseQuality',
  fields = list(metadata = 'MzQCmetadata',
                qualityMetrics = 'list'), # array of MzQCqualityMetric
  methods = list(
    initialize = function(metadata = MzQCmetadata$new(), qualityMetrics = list())
    {
      .self$metadata = metadata
      .self$qualityMetrics = qualityMetrics
    },
    isValid = function(.self)
    {
      if (!isValidMzQC(.self$metadata, .self$qualityMetrics)) return(FALSE)
      return(TRUE)
    },
    toJSON = function(.self, ...)
    {
      if (!isValidMzQC(.self)) stop(paste0("Object of class '", class(.self), "' is not in a valid state for writing to JSON"))
      
      r = list("metadata" = .self$metadata,
               "qualityMetrics" = .self$qualityMetrics)
      return (jsonlite::toJSON(r))
    },
    fromData = function(.self, data)
    {
      .self$metadata = data$metadata
      .self$qualityMetrics = fromDatatoMzQC(MzQCbaseQuality, data$qualityMetrics) ## if data$qualityMetrics is empty, or NA, the empty list will be returned
    }
  )
)
setMethod('asJSON', 'MzQCbaseQuality', function(x, ...) x$toJSON(...))


MzQCrunQuality =  setRefClass(
  "MzQCrunQuality",
  contains = "MzQCbaseQuality"
)
MzQCsetQuality =  setRefClass(
  "MzQCsetQuality",
  contains = "MzQCbaseQuality"
)

###########################################################################

#' Root element of an mzQC document
#' 
#' At least one of runQualities or setQualities MUST be present.
#' 
#' @field version Version of the mzQC format.
#' @field creationDate Creation date of the mzQC file.
#' @field contactName Name of the operator/creator of this mzQC file.
#' @field contactAddress Contact address (mail/e-mail or phone)
#' @field description Description and comments about the mzQC file contents.
#' @field runQualities Array of MzQCrunQuality;
#' @field setQualities Array of MzQCsetQuality
#' @field controlledVocabularies Array of CV domains used (obo files)
#' 
#' @export MzQCmzQC
#' 
#'    
MzQCmzQC = setRefClass(
  'MzQCmzQC',
  fields = list(version = 'character',
                creationDate = 'MzQCDateTime',
                contactName = 'character',            # optional
                contactAddress = 'character',         # optional
                description = 'character',            # optional
                runQualities = 'list',                # either this ... or  (array of MzQCrunQuality)
                setQualities = 'list',                # ... this must be present  (array of MzQCsetQuality)
                controlledVocabularies = 'list'),     # array of MzQCcontrolledVocabulary
  methods = list(
    initialize = function(version = NA_character_, 
                          creationDate = MzQCDateTime$new(), 
                          contactName = NA_character_, 
                          contactAddress = NA_character_, 
                          description = NA_character_,
                          runQualities = list(),
                          setQualities = list(), 
                          controlledVocabularies = list())
    {
      .self$version = version
      .self$creationDate = creationDate
      .self$contactName = contactName
      .self$contactAddress = contactAddress
      .self$description = description
      .self$runQualities = runQualities
      .self$setQualities = setQualities
      .self$controlledVocabularies = controlledVocabularies
    },
    isValid = function(.self)
    {
      # force evaluation using '+'
      if (isUndefined(.self$version) +
          !isValidMzQC(.self$creationDate, .self$runQualities, .self$setQualities, .self$controlledVocabularies)) return(FALSE)
      # at least one must be present
      if (length(.self$runQualities) + length(.self$setQualities) == 0)
      {
        warning("At least one runQuality or setQuality must be present! (currently all empty)", immediate. = TRUE, call. = FALSE)
        return(FALSE)
      }
      return(TRUE)
    },
    toJSON = function(.self, ...)
    {
      if (!isValidMzQC(.self)) stop(paste0("Object of class '", class(.self), "' is not in a valid state for writing to JSON"))
      
      r = list("version" = .self$version,
               "creationDate" = .self$creationDate)
      if (!isUndefined(.self$contactName)) r$contactName = .self$contactName
      if (!isUndefined(.self$contactAddress)) r$contactAddress = .self$contactAddress
      if (!isUndefined(.self$description)) r$description = .self$description
      r$runQualities = .self$runQualities
      r$setQualities = .self$setQualities
      r$controlledVocabularies = .self$controlledVocabularies
      return (jsonlite::toJSON(r))
    },
    fromData = function(.self, data)
    {
      .self$version = data$version
      .self$creationDate = fromDatatoMzQC(MzQCDateTime, data$creationDate)
      .self$contactName = NULL_to_charNA(data$contactName)
      .self$contactAddress = NULL_to_charNA(data$contactAddress)
      .self$description = NULL_to_charNA(data$description)
      .self$runQualities = fromDatatoMzQC(MzQCrunQuality, data$runQualities) ## if data$runQualities is empty, or NA, the empty list will be returned
      .self$setQualities = fromDatatoMzQC(MzQCsetQuality, data$setQualities) ## if data$setQualities is empty, or NA, the empty list will be returned
      .self$controlledVocabularies = fromDatatoMzQC(MzQCcontrolledVocabulary, data$controlledVocabularies) ## if data$controlledVocabularies is empty, or NA, the empty list will be returned
    }
  )
)
setMethod('asJSON', 'MzQCmzQC', function(x, ...) x$toJSON(...))