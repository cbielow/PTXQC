#'
#' Class to read an mzTab file and store the tables internally.
#' 
#' The 'sections' field is initialized after $readMzTab was called.
#' The 'fn_map' fields should be initialized via ...$fn_map$readMappingFile(...) manually if user-defined filename mappings are desired
#' and is automatically updated/queried when $readMzTab is called.
#' 
#' @field sections MzTab sections as list. Valid list entries are: "MTD", "PRT", "PEP", "PSM", "SML", "filename" and "comments"
#' @field fn_map FilenameMapper which can translate raw filenames into something shorter
#'
#'
#'
MzTabReader = setRefClass("MzTabReader",
                       
                       fields = list(sections = "list",
                                     fn_map = "FilenameMapper"
                       ),
                       methods = list(
                         initialize=function() {
                           .self$sections = list();
                           .self$fn_map = FilenameMapper$new();
                           
                           return(.self)
                         },
                         #'
readMzTab = function(.self, file) {
  "Read a mzTab file into a list of 5 data.frames (one df per mzTab section).
   Data.frames in the resulting list are named as follows:
     'MTD', 'PRT', 'PEP', 'PSM', 'SML',.
   Additionally, 'filename' and 'comments' are valid list elements.
  "

  
  ## this implementation is derived from with minor modifications
  ## https://github.com/lgatto/MSnbase/blob/master/R/MzTab.R
  
  lines = readLines(file)
  # remove empty lines
  lines = lines[ nzchar(lines) ]
  
  ## Split on the first two characters (so headers stay in
  ## the same group as table content rows)
  lineType = substring(lines, 1, 2)
  
  ## Could be stricter in the type checking to check that all
  ## three of the first characters match the 10 allowed types
  ## but since it doesn't affect parsing, I don't think it's
  ## worth bothering.
  allowed_types = c("CO", "MT", "PR", "PE", "PS", "SM")
  stopifnot(all(lineType %in% allowed_types))
  linesByType = split(lines, lineType)
  
  ## Comments are easy: just strip the first four characters
  ## from each line.
  comments = substring(linesByType[["CO"]], 5)
  
  ## Parse the other five blocks in a loop, then fix up
  ## metadata afterwards
  res = setNames(
    lapply(
      linesByType[c("MT", "PR", "PE", "PS", "SM")],
      function(x) {
        if (length(x) == 0) return(data.frame())
        ## MTD section has no header...
        if (startsWith(x[1], "MTD")) {
          d = read.delim(text = x,
                         header = FALSE,
                         col.names = c("MTD", "key", "value"),
                         na.strings = c("", "null"),
                         stringsAsFactors = FALSE)
        } else {
          d = read.delim(text = x,
                         header = TRUE,
                         na.strings = c("", "null"),
                         stringsAsFactors = FALSE)
          colnames(d) = make.names(colnames(d), allow_ = FALSE)
        }
        return(d[,-1])
      }),
    c("MTD", "PRT", "PEP", "PSM", "SML"))
  
  ## rewrite MetaData as named vector
  #res[["MTD"]] = setNames(res[["MTD"]][,2], res[["MTD"]][, 1])
  
  ## create Raw filename mapping internally
  mtd = res[["MTD"]]
  idx_run = grep("^ms_run\\[\\d\\]-location", mtd$key, value = FALSE)
  ms_runs = gsub("[.]*-location", "\\1", mtd$key[idx_run])
  raw_filenames = mtd$value[idx_run]
  .self$fn_map$getShortNames(raw_filenames, ms_runs = ms_runs)
  
  
  res[["filename"]] = file
  res[["comments"]] = comments
  .self$sections = res
  return (NULL)
},

getParameters = function()
{
  "Converts internal mzTab metadata section to a two column key-value data.frame similar to MaxQuants parameters.txt."
  
  # just return the whole metadata section for now
  res = .self$sections[["MTD"]]
  setnames(res, old = "key", new = "parameter")
  res = rbind(res, data.frame(parameter= "fasta file", value = paste(basename(unique(.self$sections$PSM$database)), collapse=";")))
  res = res[-(grep("custom",res$parameter)),]
  res[is.na(res)]= "NULL" # temp workaround
  
  return (res)
},

getSummary = function()
{
  "Converts internal mzTab metadata section to a two data.frame with columns 'fc.raw.file', 'ms.ms.identified....'
   similar to MaxQuants summary.txt."
  
  res = .self$fn_map$raw_file_mapping[ , c("from", "to")]
  colnames(res) = c("raw.file", "fc.raw.file")
  
  ## todo: read TIC metadata and attach to current table
  res$ms.ms.identified.... = NA
  
  return (res)
},

## MaxQuant-like representation of PRT table, i.e. augmented this with more columns (or renamed) if a metric requires it
getProteins = function()
{
  "Basically the PRT table ..."
  
  res = .self$sections$PRT

  return ( res )
},

## MaxQuant-like representation of PEP table, i.e. augmented this with more columns (or renamed) if a metric requires it
getEvidence = function()
{
   
  "Basically the PEP table and additionally columns named 'raw.file' and 'fc.raw.file'."
  
  res = .self$sections$PSM
  ## augment PEP with fc.raw.file
  ## The `spectra_ref` looks like ´ms_run[x]:index=y|ms_run´
  ms_runs = sub("[.]*:.*", "\\1", res$spectra.ref)
  res = cbind(res, .self$fn_map$mapRunsToShort(ms_runs))

  #pep=.self$sections$PEP
  #ms_runs = sub("[.]*:.*", "\\1", pep$spectra.ref)
  #pep = cbind(res, .self$fn_map$mapRunsToShort(ms_runs))
  
  if(all(c("opt.global.rt.align", "opt.global.rt.raw") %in% colnames(res))) 
  {
    colnames(res)[colnames(res)=="retention.time"] = "retention.time.pep" #rename existing retention.time column
    colnames(res)[colnames(res)=="opt.global.rt.raw"] = "retention.time"
    colnames(res)[colnames(res)=="opt.global.rt.align"] = "calibrated.retention.time"
    res$retention.time.calibration = res$calibrated.retention.time - res$retention.time 
  }
  else res$retention.time.calibration = NA
  
  res$match.time.difference = NA
  res$type = "MULTI-MSMS"
  
  name = list(opt.global.calibrated.mz.error.ppm = "mass.error..ppm.",
              opt.global.uncalibrated.mz.error.ppm = "uncalibrated.mass.error..ppm.", 
              exp.mass.to.charge = "m.z", 
              opt.global.mass = "mass", 
              opt.global.FWHM = "retention.length",
              opt.global.identified = "identified",
              opt.global.ScanEventNumber = "scan.event.number",
              PSM.ID = "id", 
              opt.global.modified.sequence = "modified.sequence",
              opt.global.is.contaminant = "contaminant",
              opt.global.fragment.mass.error.da = "mass.deviations..da.",
              opt.global.total.ion.count = "total.ion.current",
              opt.global.base.peak.intensity = "base.peak.intensity")
  
  setnames(res, old = names(name), new = unlist(name))
   
  res = aggregate(res[, colnames(res)!="id"], list("id" = res[,"id"]), function(x) {if(length(unique(x)) > 1){ paste0(unique(x), collapse = ".")} else{return (x[1])}})
  
  ## temp workaround
  res = res[!is.na(res$fc.raw.file),]
  
  return ( res )
},

## MaxQuant-like representation of PSM table, i.e. augmented this with more columns (or renamed) if a metric requires it
getMSMSScans = function()
{
  
  "Basically the PSM table and additionally columns named 'raw.file' and 'fc.raw.file'."
  
  res = .self$sections$PSM
  ## augment PSM with fc.raw.file
  ## The `spectra_ref` looks like ´ms_run[x]:index=y|ms_run´
  ms_runs = sub("[.]*:.*", "\\1", res$spectra.ref)
  res = cbind(res, .self$fn_map$mapRunsToShort(ms_runs))

  if(all(c("opt.global.rt.align", "opt.global.rt.raw") %in% colnames(res))) 
  {
    colnames(res)[colnames(res)=="retention.time"] = "retention.time.pep" #rename existing retention.time column
    colnames(res)[colnames(res)=="opt.global.rt.raw"] = "retention.time"
    colnames(res)[colnames(res)=="opt.global.rt.align"] = "calibrated.retention.time"
    res$retention.time.calibration = res$calibrated.retention.time - res$retention.time 
  }
  else res$retention.time.calibration = NA
  
  if("opt.global.ion.injection.time" %in% colnames(res)){setnames(res, old = "opt.global.ion.injection.time", new = "ion.injection.time")}
  
  name = list(opt.global.calibrated.mz.error.ppm = "mass.error..ppm",
              opt.global.uncalibrated.mz.error.ppm = "uncalibrated.mass.error..ppm.", 
              exp.mass.to.charge = "m.z", 
              opt.global.mass = "mass", 
              opt.global.fragment.mass.error.da = "mass.deviations..da.", 
              opt.global.fragment.mass.error.ppm = "mass.deviations..ppm.",
              opt.global.identified = "identified",
              opt.global.ScanEventNumber = "scan.event.number",
              PSM.ID = "id", 
              opt.global.modified.sequence = "modified.sequence",
              opt.global.is.contaminant = "contaminant",
              opt.global.missed.cleavages = "missed.cleavages",
              opt.global.target.decoy = "reverse",
              opt.global.activation.method = "fragmentation")
 
  setnames(res, old = names(name), new = unlist(name))
  
  res$mass.deviations..ppm. = gsub("\\[|\\]", "", res$mass.deviations..ppm.)
  res$mass.deviations..ppm. = gsub(",", ";", res$mass.deviations..ppm.)
  res$mass.deviations..da. = gsub("\\[|\\]", "",  res$mass.deviations..da.)
  res$mass.deviations..da. =  gsub(",", ";", res$mass.deviations..da.)
  res$reverse = (res$reverse=="decoy")

  
  ## temp workaround
  res = res[!is.na(res$contaminant),]
  res = res[order(res$fc.raw.file, res$retention.time), ]
  res = aggregate(res[, colnames(res)!="id"], list("id" = res[,"id"]), function(x) {if(length(unique(x)) > 1){ paste0(unique(x), collapse = ".")} else{return (x[1])}})

  return ( res )
},

getMSMS = function()
{
  res = .self$sections$PSM
  ms_runs = sub("[.]*:.*", "\\1", res$spectra_ref)
  res = cbind(res, mzt$fn_map$mapRunsToShort(ms_runs))
  
  #cbind(.self$sections$PSM$PSM_ID, .self$sections$PSM$PSM_ID$opt.global.missed.cleavages,.self$sections$PSM$PSM_ID$opt.global.target.decoy,mzt$fn_map$mapRunsToShort(ms_runs))
  colnames(res)[colnames(res)=="PSM_ID"] <- "id"
  colnames(res)[colnames(res)=="opt.global.missed.cleavages"]<- "missed.cleavages"
  colnames(res)[colnames(res)=="opt.global.target.decoy"]<- "reverse"
  res$reverse[df=="decoy"]<-TRUE
  res$reverse[df!=TRUE]<-FALSE
  colnames(res)[colnames(res)=="opt.global.target.fragment.mass.error.da"]<- "mass.deviations..da."
  colnames(res)[colnames(res)=="opt.global.target.fragment.mass.error.ppm"]<- "mass.deviations..ppm."
  
  return ( res )
  
}

) # methods
) # class
