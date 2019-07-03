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
  
  f_con = file(mztab_file, open = "r") ## for better error messages
  lines = readLines(f_con)
  close(f_con)
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
                         stringsAsFactors = FALSE,
                         fill = FALSE)
        } else {
          d = read.delim(text = x,
                         header = TRUE,
                         na.strings = c("", "null"),
                         stringsAsFactors = FALSE,
                         fill = FALSE)
          colnames(d) = make.names(colnames(d), allow_ = FALSE)
        }
        return(d[,-1])
      }),
    c("MTD", "PRT", "PEP", "PSM", "SML"))
  
  ## rewrite MetaData as named vector
  #res[["MTD"]] = setNames(res[["MTD"]][,2], res[["MTD"]][, 1])
  
  ## create Raw filename mapping internally
  mtd = res[["MTD"]]
  idx_run = grep("^ms_run\\[\\d*\\]-location", mtd$key, value = FALSE)
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
  
  # copy the whole MTD for now, since R likes shallow copies and we are about to rename columns by reference (see ?setnames)
  res = data.table::copy(.self$sections[["MTD"]])
  res = rbind(res, data.frame(key= "fasta file", value = paste(basename(unique(.self$sections$PSM$database)), collapse=";")))
  res = res[-(grep("custom", res$key)),]
  res[is.na(res)] = "NULL" # temp workaround
  
  data.table::setnames(res, old = "key", new = "parameter") ## todo: remove at some point, since it forces us to use `::copy`

  return (res)
},

getSummary = function()
{
  "Converts internal mzTab metadata section to a two data.frame with columns 'fc.raw.file', 'ms.ms.identified....'
   similar to MaxQuants summary.txt."
  res = .self$fn_map$getRawfm()[ , c("from", "to")]
  colnames(res) = c("raw.file", "fc.raw.file")
  
  ## read all custom entries
  mtd_custom_df = .self$sections$MTD[grep("^custom", .self$sections$MTD$key), ]

  ## ... and subselect the ms2-ID-Rate
  ms2_df = mtd_custom_df[grep("^\\[MS2 identification rate", mtd_custom_df$value), ] 
  res$ms.ms.identified.... = unlist(lapply(gsub(".*, ?(\\d*\\.\\d*)\\]", "\\1", ms2_df$value), as.numeric))

  ## read TIC
  tic_df = mtd_custom_df[grep("total ion current", mtd_custom_df$value),] 
  res$TIC = lapply(strsplit(sub(".* \\[(.*)\\]]", "\\1", tic_df$value), ","), as.numeric)

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
  
  # remove empty PepIDs (with no ConsensusFeature group or unassigned PepIDs (-1); corresponding to unidentfied MS2 scans)
  res = res[!is.na(res$opt.global.cf.id),]
  
  
  ## augment with fc.raw.file
  ## The `spectra_ref` looks like ´ms_run[x]:index=y|ms_run´
  ms_runs = sub("[.]*:.*", "\\1", res$spectra.ref)
  res = cbind(res, .self$fn_map$mapRunsToShort(ms_runs))
  stopifnot(all(!is.na(res$fc.raw.file))) # Spectra-Ref in PSM table not set for all entries
  
  
  if (all(c("opt.global.rt.align", "opt.global.rt.raw") %in% colnames(res))) 
  {
    data.table::setnames(res, old = c("retention.time","opt.global.rt.raw","opt.global.rt.align"), 
                  new = c("retention.time.pep","retention.time","calibrated.retention.time"))
    res$retention.time.calibration = res$calibrated.retention.time - res$retention.time 
  }
  else res$retention.time.calibration = NA
  
  res$match.time.difference = NA
  res$type = "MULTI-MSMS"
  
  name = list(opt.global.calibrated.mz.error.ppm = "mass.error..ppm.",
              opt.global.uncalibrated.mz.error.ppm = "uncalibrated.mass.error..ppm.", 
              exp.mass.to.charge = "m.z", 
              opt.global.mass = "mass", 
              opt.global.identified = "identified",
              opt.global.ScanEventNumber = "scan.event.number",
              PSM.ID = "id", 
              opt.global.modified.sequence = "modified.sequence",
              opt.global.is.contaminant = "contaminant",
              opt.global.fragment.mass.error.da = "mass.deviations..da."
          )
  
  data.table::setnames(res, old = names(name), new = unlist(name))
   
  #res = aggregate(res[, colnames(res)!="id"], list("id" = res[,"id"]), function(x) {if(length(unique(x)) > 1){ paste0(unique(x), collapse = ".")} else{return (x[1])}})

  ## optional in MzTab (depending on which FeatureFinder was used)
  if("opt.global.FWHM" %in% colnames(res)) {
    data.table::setnames(res, old = c("opt.global.FWHM"), new = c("retention.length"))
  }

  
  #ms.ms.count: 
  #1.all different accessions and databases per ID in one row; 2.all IDs only one time; 
  #3.add ms.ms.count (size of groups with same sequence, modified.sequence and charge; 4. set in all groups ms.ms.count all cells but one to NA )
  res_dt = data.table::setDT(res)
  accessions = (res_dt[, .(accession=list(accession)), by=id])$accession
  databases = (res_dt[, .(database=list(database)), by=id])$database
  res_dt = unique(res_dt, by = "id")
  res_dt$proteins = unlist(lapply(accessions, paste, collapse=";"))
  res_dt$protein.group.ids = res_dt$proteins ## todo: these are protein names not IDs, but the code does not care (its not very clean though)
  
  res_dt$database = databases
  res_dt[, ms.ms.count := .N, by = list(raw.file, modified.sequence,charge)]
  toNA = res_dt[, .(toNA = .I[c(1L:.N-1)]), by=list(raw.file, modified.sequence,charge)]$toNA
  res_dt[toNA, ms.ms.count := NA]
  res = as.data.frame(res_dt)
  
  # intensity from PEP to PSM: only labelfree
  pep_df = .self$sections$PEP

  res$pep.id = match(res$spectra.ref, pep_df$spectra.ref, nomatch = NA_integer_)
  res$ms_run_number = as.numeric(sub("\\].*","", sub(".*\\[","", res$spectra.ref)))
  pep_intensity_df = pep_df[, grepl( "peptide.abundance.study.variable." , names(pep_df))]

  res = plyr::ddply(res, "opt.global.cf.id", function(x){
              pep_row = data.table::first(na.omit(x$pep.id))
              x$intensity = as.numeric(pep_intensity_df[pep_row, x$ms_run_number]) 
              return(x)}) 
  
  res$intensity[duplicated(res[,c("opt.global.map.index","opt.global.cf.id")])] = NA

  ## just check if there are no invalid entries
  stopifnot(all(!is.na(res$contaminant)))
  
  return ( res )
},

## MaxQuant-like representation of PSM table, i.e. augmented this with more columns (or renamed) if a metric requires it
getMSMSScans = function(identified_only = FALSE)
{
  "Basically the PSM table and additionally columns named 'raw.file' and 'fc.raw.file'. 
   If identified_only is TRUE, only MS2 scans which were identified (i.e. a PSM) are returned -- this is equivalent to msms.txt in MaxQuant."
  
  res = .self$sections$PSM
  res = as.data.table(res)
  data.table::setkey(res, PSM.ID)
  
  if (identified_only) {
    res = res[opt.global.identified == 1, ]
  }
  
  ## de-duplicate PSM.ID column: take first row for each PSM.ID 
  res_temp = res[!duplicated(res$PSM.ID), ]
  ## ... and append accessions of the complete subset (.SD)
  res_temp$accessions = res[, paste0(.SD$accession, collapse=";"), by = "PSM.ID"]$V1
  res = res_temp
  
  ## Augment fc.raw.file column
  ## ... the `spectra_ref` looks like "ms_run[12]:controllerType=0 controllerNumber=1 scan=25337"
  ms_runs = sub("^(ms_run\\[\\d*\\]):.*", "\\1", res$spectra.ref)
  res = cbind(res, .self$fn_map$mapRunsToShort(ms_runs)) ## gives c("ms_run[1]", "ms_run[2]", ...)

  if(all(c("opt.global.rt.align", "opt.global.rt.raw") %in% colnames(res))) 
  {
    data.table::setnames(res,
             old = c("retention.time","opt.global.rt.raw","opt.global.rt.align"),
             new = c("retention.time.pep","retention.time","calibrated.retention.time"))
    res$retention.time.calibration = res$calibrated.retention.time - res$retention.time 
  }
  else res$retention.time.calibration = NA
  
  if("opt.global.ion.injection.time" %in% colnames(res))
  {
    data.table::setnames(res, old = "opt.global.ion.injection.time", new = "ion.injection.time")
  }
  
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
              opt.global.activation.method = "fragmentation",
              opt.global.total.ion.count = "total.ion.current",
              opt.global.base.peak.intensity = "base.peak.intensity")
 
  data.table::setnames(res, old = names(name), new = unlist(name))
  
  res$mass.deviations..ppm. = substr(res$mass.deviations..ppm., 2, nchar(res$mass.deviations..ppm.) - 2)
  res$mass.deviations..ppm. = gsub(",", ";", res$mass.deviations..ppm., fixed = TRUE)
  res$mass.deviations..da. = substr(res$mass.deviations..da., 2, nchar(res$mass.deviations..da.) - 2)
  res$mass.deviations..da. =  gsub(",", ";", res$mass.deviations..da., fixed = TRUE)
 
  # set reverse to TRUE/FALSE
  res$reverse = (res$reverse == "decoy")
  
  # set identified to MaxQuant values (+/-)
  stopifnot(unique(res$identified) %in% c(0,1)) ## make sure the column has the expected values (0/1)
  res$identified = c("-", "+")[(res$identified==1) + 1]

  ## temp workaround
  #res = res[!is.na(res$contaminant),]
  
  ## order by file and RT
  res = res[order(res$fc.raw.file, res$retention.time), ]
  
  return ( res )
}

) # methods
) # class
