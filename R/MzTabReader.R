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

  cat("Reading mzTab '", file, "' ...", sep = "")
  
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
                         na.strings = c("", "null", "not mapped"),
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
  ## in addition to the correct 'ms_run[xx]-location' we also accept 'ms_run[xx]_location' (ProteomeDiscoverer2.2)
  idx_run = grep("^ms_run\\[\\d*\\][_-]location", mtd$key, value = FALSE)
  raw_filenames = mtd$value[idx_run]
  ms_runs = gsub("([.]*)[_-]location", "\\1", mtd$key[idx_run])
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
  if (!is.na(unique(.self$sections$PSM$database))) {
    res = rbind(res, data.frame(key= "fasta file", value = paste(basename(unique(.self$sections$PSM$database)), collapse=";")))
  }
  else {
    res = rbind(res, data.frame(key= "fasta file", value = "NULL"))
  }
  
  res = res[grep("^custom", res$key, invert = TRUE),]
  
  res[is.na(res)] = "NULL" # temp workaround
  
  ## todo: remove at some point, since it forces us to use `::copy`
  renameColumns(res, list(key = "parameter"))

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
  
  if (nrow(mtd_custom_df) == 0) return(NULL)

  ## ... and subselect the ms2-ID-Rate
  ms2_df = mtd_custom_df[grep("^\\[MS2 identification rate", mtd_custom_df$value), ] 
  res$ms.ms.identified.... = unlist(lapply(gsub(".*, ?(\\d*\\.\\d*)\\]", "\\1", ms2_df$value), as.numeric))

  ## read TIC
  tic_df = mtd_custom_df[grep("total ion current", mtd_custom_df$value),]
  if (nrow(tic_df) > 0) res$TIC = lapply(strsplit(sub(".* \\[(.*)\\]]", "\\1", tic_df$value), ","), as.numeric)

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
  "Basically the PSM table and additionally columns named 'raw.file' and 'fc.raw.file'."
  
  res = data.table::as.data.table(.self$sections$PSM)

  ## remove empty PepIDs 
  ## ... unidentfied MS2 scans (NA)      or with no ConsensusFeature group (-1), i.e. unassigned PepIDs
  
  if (all(c("opt.global.cf.id") %in% colnames(res))) {
    res = res[!(is.na(res$opt.global.cf.id) | (res$opt.global.cf.id == -1)),]
    stopifnot(min(res$opt.global.cf.id) >= 0) ## would stop on NA as well
  }
  
  ## augment with fc.raw.file
  ## The 'spectra_ref' looks like 'ms_run[x]:index=y|ms_run'
  res = cbind(res, .self$fn_map$specrefToRawfile(res$spectra.ref))
  stopifnot(all(!is.na(res$fc.raw.file))) # Spectra-Ref in PSM table not set for all entries
  
  
  res$retention.time.calibration = NA
  if (all(c("opt.global.rt.align", "opt.global.rt.raw") %in% colnames(res))) 
  {
    renameColumns(res, list(               retention.time = "retention.time.pep",
                                        opt.global.rt.raw = "retention.time",
                                      opt.global.rt.align = "calibrated.retention.time"
                            ))
    res$retention.time.calibration = res$calibrated.retention.time - res$retention.time 
  }

  renameColumns(res, list(         opt.global.calibrated.mz.error.ppm = "mass.error..ppm.",
                                 opt.global.uncalibrated.mz.error.ppm = "uncalibrated.mass.error..ppm.", 
                                                   exp.mass.to.charge = "m.z", 
                                                      opt.global.mass = "mass", 
                                                opt.global.identified = "identified",
                                           opt.global.ScanEventNumber = "scan.event.number",
                                                               PSM.ID = "id", 
                                         opt.global.modified.sequence = "modified.sequence",
                                            opt.global.is.contaminant = "contaminant",
                                    opt.global.fragment.mass.error.da = "mass.deviations..da.",
                               opt.global.cv.MS.1002217.decoy.peptide = "reverse"
          ))

  if (!"modified.sequence" %in% colnames(res)){
    res$modified.sequence = res$sequence
    warning("modified.sequence is not present in input data, metrics use sequence instead", immediate. = TRUE)
  }
  
  if(!"contaminant" %in% colnames(res)){
    res$contaminant = 0
  }
  # set contaminant to TRUE/FALSE
  res$contaminant = (res$contaminant > 0)

  ## optional in MzTab (depending on which FeatureFinder was used)
  if ("opt.global.FWHM" %in% colnames(res)){
    renameColumns(res, list(opt.global.FWHM = "retention.length"))
  }
  
  ## de-duplicate protein-accession entries
  accessions = res[, .(l_accession = list(accession)), by = id]$l_accession
  res = unique(res, by = "id")
  res$proteins = unlist(lapply(accessions, paste, collapse = ";"))
  
  ## todo: these are protein names not IDs, but the code does not care (its not very clean though)
  res$protein.group.ids = res$proteins 
  
  ## annotate ms.ms.count for identical sequences per rawfile, but only the first member of the group; 
  ## all others get NA to prevent double counting
  res[, ms.ms.count := c(.N, rep(NA, .N-1)), by = list(raw.file, modified.sequence, charge)]

  ## convert values from seconds to minutes for all RT columns
  RTUnitCorrection(res)

  ##
  ## intensity from PEP to PSM: only labelfree ('opt.global.cf.id' links all PSMs belonging to a ConsensusFeature)
  ##
  df_pep = data.frame()
  if ("opt.global.cf.id" %in% colnames(res)){
    df_pep = data.table::as.data.table(.self$sections$PEP)[!is.na(sequence), ]
    renameColumns(df_pep, list(opt.global.modified.sequence = "modified.sequence"))
    ## add raw.file...
    df_pep = cbind(df_pep, .self$fn_map$specrefToRawfile(df_pep$spectra.ref))
    ## .. a unique index
    df_pep$idx = 1:nrow(df_pep)
    ## map from PSM -> PEP row
    ## ... do NOT use spectra.ref since this is ambiguous (IDMapper duplicates MS2 PepIDs to multiple features)
    if (!("opt.global.cf.id" %in% colnames(df_pep)))
    {
      stop("Please re-run the TOPP QualityControl tool to generate an updated version of OpenMS mzTab. Your mzTab is missing some information.")
    } 
    #old: res$pep_idx = match(res$opt.global.feature.id, df_pep$opt.global.feature.id, nomatch = NA_integer_)
    res$pep_idx = match(res$opt.global.cf.id, df_pep$opt.global.cf.id, nomatch = NA_integer_)
    
    res$ms_run_number = as.numeric(gsub("^ms_run\\[(\\d*)\\].*", "\\1", res$ms_run))
    col_abd_df_pep = grepv( "^peptide.abundance.study.variable.", names(df_pep))
    col_RT_df_pep = grepv( "^opt.global.retention.time.study.variable", names(df_pep))
    ## transposed matrix for all abundances (rows = study; cols = pep_idx)
    ## .. this is a significant speedup compared to indexing into .SD[,] in subqueries, since that requires unlist()
    m_pep_abd = t(df_pep[, ..col_abd_df_pep])
    m_pep_rt = t(df_pep[, ..col_RT_df_pep])
    N.studies = length(col_RT_df_pep)
    stopifnot(N.studies == length(col_abd_df_pep))
  }
 

  NA_duplicates = function(vec_abd, idx) {
    ## replaces duplicate indices into the same ms_run intensities with NA
    ## (to avoid counting a feature more than once due to oversampled PSMs from one run assigned to a CF)
    r = vec_abd[idx]
    r[duplicated(idx)] = NA
    return(r)
  }
  
  if (all(c("opt.global.cf.id") %in% colnames(res))) {
    ## assign intensity to genuine PSMs
    res$intensity = NA_real_ ## unassigned PSMs have no MS1 intensity
    res[,
      intensity := NA_duplicates(m_pep_abd[, .SD$pep_idx[1]],
                                 .SD$ms_run_number),
      by = "opt.global.cf.id"]
    summary(res$intensity)
  }
  

  res$is.transferred = FALSE
  res$type = "MULTI-MSMS"

  #set reverse to needed values
  if ("reverse" %in% colnames(res)){
    res$reverse=(res$reverse=="decoy")
  }
  
  ## remove the data.table info, since metrics will break due to different syntax
  class(res) = "data.frame"


  ##
  ## Infer MBR 
  ## --> find all subfeatures in a CF with abundance but missing PSM --> create as MBR-dummy-PSMs
  ##
  ## : df_pep: the PEP table (with some extras)
  ## : res:    the PSM table (with some extras)
  ##     res$pep_idx: row in df_pep (consensusFeature) to which this PSM was assigned
  res_tf = res[NULL,]
  stopifnot(ncol(res_tf) > 0)
  ## skip if MS2 isotope labeling detected
  quant_methods = grep("quantification_method", mzt$sections$MTD$key)
  if (!any(grepl("PRIDE_0000317", .self$sections$MTD$value[quant_methods])) & (!plyr::empty(df_pep)))
  {
    ## iterate though all consensusFeatures .. find subfeatures with intensity but missing PSM
    res_tf_tmp = df_pep[, {#print(idx)
      #idx = 44
      idx_PSM = which(res$pep_idx == idx)
      runs_with_MS2 = unique(res$ms_run_number[idx_PSM]) ## all existing PSMs for this PEP (=consensusFeature)
      runs_wo_MS2 = (1:N.studies)[-runs_with_MS2]
      df = .SD[rep(1, length(runs_wo_MS2)), c("charge", "modified.sequence", "sequence")]
      df$pep_idx = idx
      df$ms_run_number = runs_wo_MS2 ## vector
      df$calibrated.retention.time = m_pep_rt[runs_wo_MS2, idx]
      df$intensity = m_pep_abd[runs_wo_MS2, idx]
      df$protein.group.ids = res$protein.group.ids[idx_PSM[1]] ## use PGI from genuine PSMs
      ## remove all inferred PSMs which do not have an intensity(= MS1 feature)
      df = df[!is.na(df$intensity),]
      df ## return
    }, by = "idx"] ## one PEP row at a time
   
    if(nrow(res_tf_tmp) > 0){
      res_tf = res_tf_tmp
      ## convert from "ms_run_number" to fc.raw.file
      res_tf = cbind(res_tf, .self$fn_map$msrunToRawfile(paste0("ms_run[", res_tf$ms_run_number, "]")))
      res_tf$is.transferred = TRUE
      res_tf$type = "MULTI-MATCH"
      
      ## check: summed intensities should be equal
      if("intensity" %in% colnames(res) && "intensity" %in% colnames(res_tf)){
        stopifnot(sum(df_pep[, ..col_abd_df_pep], na.rm = TRUE) == sum(res$intensity, na.rm = TRUE) + sum(res_tf$intensity, na.rm = TRUE))
      }
    }
  }
  
  ## remove the data.table info, since metrics will break due to different syntax
  class(res_tf) = "data.frame"

  message("Evidence table generated: ", nrow(res), "x", ncol(res), "(genuine); ", nrow(res_tf), "x", ncol(res_tf), "(transferred)")
  
  ## must at least have column names, but can have 0 rows
  stopifnot(ncol(res_tf) > 0)
  
  return (list("genuine" = res, "transferred" = res_tf))
},


getMSMSScans = function(identified_only = FALSE)
{
  "Basically the PSM table (partially renamed columns) and additionally two columns 'raw.file' and 'fc.raw.file'. 
   If identified_only is TRUE, only MS2 scans which were identified (i.e. a PSM) are returned -- this is equivalent to msms.txt in MaxQuant."
  
  res = data.table::as.data.table(.self$sections$PSM)

  stopifnot(all((res$opt.global.identified == 1) == (!is.na(res$sequence))))
  if (identified_only) {
    res = res[!is.na(res$sequence), ] # == NA sequence
  }
  
  ## de-duplicate PSM.ID column: take first row for each PSM.ID 
  res_temp = res[!duplicated(res$PSM.ID), ]
  ## ... and append accessions of the complete subset (.SD)
  res_temp$accessions = res[, paste0(.SD$accession, collapse=";"), by = "PSM.ID"]$V1
  res = res_temp
  
  ## Augment fc.raw.file column
  ## ... the `spectra_ref` looks like "ms_run[12]:controllerType=0 controllerNumber=1 scan=25337"
  res = cbind(res, .self$fn_map$specrefToRawfile(res$spectra.ref))

  ## IDMapper might duplicate PepIDs if two or more features are a good match
  ## ... but we only want each MS2 scan represented once here
  res = unique(res, by = c("fc.raw.file", "spectra.ref"))
  
  if (all(c("opt.global.rt.align", "opt.global.rt.raw") %in% colnames(res))) 
  {
    renameColumns(res, list(     retention.time = "retention.time.pep",
                              opt.global.rt.raw = "retention.time",
                            opt.global.rt.align = "calibrated.retention.time"))
    res$retention.time.calibration = res$calibrated.retention.time - res$retention.time 
  }
  else res$retention.time.calibration = NA
  
  
  name = list(         opt.global.calibrated.mz.error.ppm = "mass.error..ppm",
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
                   opt.global.cv.MS.1002217.decoy.peptide = "reverse",
                             opt.global.activation.method = "fragmentation",
                               opt.global.total.ion.count = "total.ion.current",
                           opt.global.base.peak.intensity = "base.peak.intensity",
                            opt.global.ion.injection.time = "ion.injection.time")
 
  renameColumns(res, name)
 
  if ("mass.deviations..ppm." %in% colnames(res)) {
    res$mass.deviations..ppm. = substr(res$mass.deviations..ppm., 2, nchar(res$mass.deviations..ppm.) - 1)
    res$mass.deviations..ppm. = gsub(",", ";", res$mass.deviations..ppm., fixed = TRUE)
  }
  if ("mass.deviations..da." %in% colnames(res)) {   
    res$mass.deviations..da. = substr(res$mass.deviations..da., 2, nchar(res$mass.deviations..da.) - 1)
    res$mass.deviations..da. =  gsub(",", ";", res$mass.deviations..da., fixed = TRUE)
  }
  
  #set reverse to needed values
  if ("reverse" %in% colnames(res)){
    res$reverse=(res$reverse=="decoy")
  }
  
  if (!"contaminant" %in% colnames(res)){
    res$contaminant = 0
  }
  #set contaminant to TRUE/FALSE
  res$contaminant = (res$contaminant > 0)
  
  #set identified to needed values
  if ("identified" %in% colnames(res)){
    stopifnot(unique(res$identified) %in% c(0,1)) ## make sure the column has the expected values (0/1)
    # set $identified to MaxQuant values (+/-)
    res$identified = c("-", "+")[(res$identified==1) + 1]
  }

  RTUnitCorrection(res)


  ## remove the data.table info, since metrics will break due to different syntax
  class(res) = "data.frame"
  
  ## order by file and specRef as RT proxy (do NOT use RT directly, since it might be NA or non-linearly transformed)
  ## e.g. spectra.ref might be 'ms_run[1]:controllerType=0 controllerNumber=1 scan=13999'
  ##                        or 'ms_run[2]:spectrum=33'
  ##      --> extract scan as numeric, since string compare is insufficient for numbers ("13999" > "140")
  res$scan = as.numeric(gsub(".*index=(\\d*)|.*scan=(\\d*)|.*spectrum=(\\d*)", "\\1\\2\\3", res$spectra.ref))
  stopifnot(all(!is.na(res$scan))) 
  res = res[order(res$fc.raw.file, res$scan), ]
  
  return ( res )
},

RTUnitCorrection = function(dt)
{
  "Convert all RT columns from seconds (OpenMS default) to minutes (MaxQuant default)"
  
  # heuristic to detect presence of unit:seconds; if retention.time has is, we assume that all rt-columns are in seconds
  # retention.time is mandatory for mzTab
  if (max(dt[, "retention.time"], na.rm = TRUE) > 300)
  {
    cn_rt = grepv("retention.time|retention.length", names(dt))
    dt[, c(cn_rt) := lapply(.SD, function(x) x / 60 ), .SDcols = cn_rt]
  }
  
  #dt[, ..cn_rt]
  return(NULL)
},

renameColumns = function(dt, namelist)
{
  "Renames all columns and throws a warning if a column does not exist in the data"
  
  from = names(namelist)
  to = unlist(namelist)
  data.table::setnames(dt, old = from, new = to, skip_absent = TRUE)
  
  existName = to %in% colnames(dt)
  if (!all(existName))
  {
    warning(paste0("Columns\n '", 
                   paste(from[!existName], "' (mzTab name) --> '", to[!existName], collapse="' (internal name),\n '", sep=""),
                   "'\n are not present in input data!"),
            immediate. = TRUE)
  }
}

) # methods
) # class
