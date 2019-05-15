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
  return (.self$sections[["MTD"]])
},

getSummary = function()
{
  "Converts internal mzTab metadata section to a two data.frame with columns 'fc.raw.file', 'ms.ms.identified....'
   similar to MaxQuants summary.txt."
  res = .self$fn_map$getRawfm()[ , c("from", "to")]
  colnames(res) = c("raw.file", "fc.raw.file")
  
  #read custom entrys
  mtd_custom_df= .self$sections$MTD[grep("custom", .self$sections$MTD$key),]
 
  ##ms2-ID-Rate
  ms2_df=mtd_custom_df[grep("identification", mtd_custom_df$value),] 
  res$ms.ms.identified....=unlist(lapply(lapply(strsplit(gsub("]","",as.character(ms2_df$value)),","), "[[", 4),as.numeric))
  
  ## read TIC
  tic_df=mtd_custom_df[grep("total ion current", mtd_custom_df$value),] 
  res$TIC=lapply(strsplit(sub(".* \\[(.*)\\]", "\\1", tic_df$value), ","), as.numeric)

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
  ## The `spectra_ref` looks like Â´ms_run[x]:index=y|ms_runÂ´
  ms_runs = sub("[.]*:.*", "\\1", res$spectra.ref)
  res = cbind(res, .self$fn_map$mapRunsToShort(ms_runs))

  
  res$match.time.difference = NA
  
  setnames(res, old = c("opt.global.calibrated.mz.error.ppm","opt.global.uncalibrated.mz.error.ppm", "opt.global.activation.method"), new = c("mass.error..ppm.","uncalibrated.mass.error..ppm.","fragmentation"))
  setnames(res, old = c("opt.global.identified","opt.global.ScanEventNumber","PSM.ID", "opt.global.modified.sequence","opt.global.is.contaminant","opt.global.fragment.mass.error.da","opt.global.fragment.mass.error.ppm"), new = c("identified","scan.event.number","id", "modified.sequence","contaminant","mass.deviations..da.","mass.deviations..ppm."))

  if(all(c("opt.global.rt.align", "opt.global.rt.raw") %in% colnames(res))) 
  {
    setnames(res, old = c("retention.time","opt.global.rt.raw","opt.global.rt.align"), new = c("retention.time.pep"," retention.time","calibrated.retention.time"))
    res$retention.time.calibration = res$calibrated.retention.time - res$retention.time 
  }
  else 
  {
    res$retention.time.calibration=NA
  }
  
  if("opt.global.FWHM" %in% colnames(res)) {  setnames(res, old = c("opt.global.FWHM"), new = c("retention.length")) }

  #ms.ms.count: 
  #1.all different accessions and databases per ID in one row; 2.all IDs only one time; 
  #3.add ms.ms.count (size of groups with same sequence, modified.sequence and charge; 4. set in all groups ms.ms.count all cells but one to NA )
  res_dt=setDT(res)
  accessions=(res_dt[, .(accession=list(accession)),by=id])$accession
  databases=(res_dt[, .(database=list(database)),by=id])$database
  res_dt=unique(res_dt, by = "id")
  res_dt$accession=accessions
  res_dt$database=databases
  res_dt[,ms.ms.count:=.N, by=list(raw.file,modified.sequence,charge)]
  toNA=res_dt[, .(toNA = .I[c(1L:.N-1)]), by=list(raw.file,modified.sequence,charge)]$toNA
  res_dt[toNA, ms.ms.count:=NA]
  res=as.data.frame(res_dt)
  
  #intensity from PEP to PSM: only labelfree
  
  pep_df = .self$sections$PEP
  res$pep.id = as.numeric(NA)
  res$intensity = as.numeric(NA)
  res$ms_run_number = as.numeric(NA)
 
   #split data.frame rev in res_df and empty entrys
  empty_entries = res[is.na(res$spectra.ref),]
  res_df = res[!is.na(res$spectra.ref),]

  res_df$pep.id = match(res_df$spectra.ref, pep_df$spectra.ref, nomatch = NA_integer_)
  res_df$ms_run_number = as.numeric(sub("\\].*","", sub(".*\\[","", res_df$spectra.ref)))
  pep_intensity_df = pep_df[ ,grepl( "peptide.abundance.study.variable." , names(pep_df))]

  res_df=ddply(res_df,"opt.global.cf.id",function(x){
             pep_row=first(na.omit(x$pep.id))
             x$intensity=as.numeric(pep_intensity_df[pep_row, x$ms_run_number]) 
             return(x)}) 
  
  res_df$intensity[duplicated(res_df[,c("opt.global.map.index","opt.global.cf.id")])]=NA

  
  #apply empty entrys
  res=rbind(res_df, empty_entries)

  ## temp workaround
  res = res[!is.na(res$fc.raw.file),]
  # res_con=res[which(res$contaminant==1),]
  # print(res_con[,c("spectra.ref","ms_run_number","pep.id","intensity","opt.global.cf.id")])
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

  if("opt.global.ion.injection.time" %in% colnames(res))
  {
    setnames(res, old = c("opt.global.ion.injection.time"), new = c("ion.injection.time"))
  }
  setnames(res, old = c("opt.global.identified","opt.global.ScanEventNumber","PSM.ID", "opt.global.modified.sequence","opt.global.is.contaminant","opt.global.fragment.mass.error.da","opt.global.fragment.mass.error.ppm", "opt.global.missed.cleavages","opt.global.target.decoy"), new = c("identified","scan.event.number","id", "modified.sequence","contaminant","mass.deviations..da.","mass.deviations..ppm.","missed.cleavages","reverse"))
  setnames(res, old = c("opt.global.calibrated.mz.error.ppm","opt.global.uncalibrated.mz.error.ppm","opt.global.activation.method"), new = c("mass.error..ppm","uncalibrated.mass.error..ppm","fragmentation"))

  if(all(c("opt.global.rt.align", "opt.global.rt.raw") %in% colnames(res))) 
  {
    colnames(res)[colnames(res)=="opt.global.rt.raw"] = "retention.time"
    colnames(res)[colnames(res)=="opt.global.rt.align"] = "calibrated.retention.time"
    res$retention.time.calibration = res$calibrated.retention.time - res$retention.time 
  }
  else res$retention.time.calibration = NA
 
 #set reverse to needed values
  res$reverse=(res$reverse=="decoy")


  ## temp workaround
  res = res[!is.na(res$contaminant),]
  return ( res )
}

) # methods
) # class
