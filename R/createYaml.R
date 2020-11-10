#' Creates a yaml file storing the parameters that are used for creating the PTXQC report 
#' and returns these parameters as well as a list of available qc-Metrics objects.
#' 
#' Valid parameters are: 
#'    param_useMQPAR, add_fs_col, id_rate_bad, id_rate_great , pg_ratioLabIncThresh , param_PG_intThresh,
#'    param_EV_protThresh , param_EV_intThresh, param_EV_pepThresh , yaml_contaminants, param_EV_MatchingTolerance,
#'    param_evd_mbr , param_EV_PrecursorTolPPM, param_EV_PrecursorOutOfCalSD , param_EV_PrecursorTolPPMmainSearch, 
#'    param_MSMSScans_ionInjThresh, param_OutputFormats and param_PageNumbers 
#'    
#'    Please provide them as a list() of this format: list$parameter_name
#'
#'
#' @param yc A yaml class object created by YAMLClass$new()
#' @param param list of parameters sorted by names; if empty, will be populated with defaults
#' @param DEBUG_PTXQC print some debugging information; default FALSE
#' @param txt_files list of paths to MaxQuant files; if NULL, it is assumed that the parameters are for mzTab-mode
#' @param metrics list of metric names that should be plotted; if NULL, will be populated with defaults
#' @return list of parameters used for creating the report and list of qc-Metrics objects
#' @export
#'
#'
createYaml <- function(yc, param = list(), DEBUG_PTXQC = FALSE, txt_files = NULL, metrics = NULL){

  ##
  ## YAML default config
  ##
  default_param <- list()
  default_param$param_useMQPAR <- TRUE
  default_param$add_fs_col <- 14 
  default_param$id_rate_bad <- 20
  default_param$id_rate_great <- 35
  default_param$pg_ratioLabIncThresh <- 4
  default_param$param_PG_intThresh <- 25
  default_param$param_EV_protThresh <- 3500
  default_param$param_EV_intThresh <- 23
  default_param$param_EV_pepThresh <- 15000
  default_param$yaml_contaminants <- list("cont_MYCO" = c(name="MYCOPLASMA", threshold=1)) # name (FASTA), threshold for % of unique peptides
  default_param$param_EV_MatchingTolerance <- 1
  default_param$param_evd_mbr <- "auto"
  default_param$param_EV_PrecursorTolPPM <- 20
  default_param$param_EV_PrecursorOutOfCalSD <- 2
  default_param$param_EV_PrecursorTolPPMmainSearch <- NA
  default_param$param_MSMSScans_ionInjThresh <- 10
  default_param$param_OutputFormats <- c("html", "plainPDF")
  default_param$param_PageNumbers <- "on"
  
  
  ##
  ##add missing parameters from default parameter list
  ##
  for(i in c(1:length(default_param))){
    if(!names(default_param)[i] %in% names(param)) param[names(default_param)[i]] <- default_param[i]
  }

  ##
  ## check for invalid parameters
  ##
  for(i in c(length(param):1)){
    if(!names(param)[i] %in% c(NA, names(default_param))) {
      warning(paste0("Invalid parameter detected and removed: ", names(param)[i]))
      param <- param[-i] 
      }
  }
  
  
  ## determines if a local mqpar.xml should be used to grep all YAML parameters whose name starts with "MQpar_" from the
  ## original mqpar.xml instead of the yaml.config. The "MQpar_..." param from the config
  ## will be ignored and the newly written yaml.config will contain the values from mqpar.xml.
  param$param_useMQPAR = yc$getYAML("PTXQC$UseLocalMQPar", param$param_useMQPAR)
  
  param$add_fs_col = yc$getYAML("PTXQC$NameLengthMax_num", param$add_fs_col)
  
  param$id_rate_bad = yc$getYAML("File$Summary$IDRate$Thresh_bad_num", param$id_rate_bad, 0, 100)
  param$id_rate_great = yc$getYAML("File$Summary$IDRate$Thresh_great_num", param$id_rate_great, 0, 100)
  
  param$GL_name_min_length = 8
  
  param$pg_ratioLabIncThresh = yc$getYAML("File$ProteinGroups$RatioPlot$LabelIncThresh_num", param$pg_ratioLabIncThresh)
  ## default median intensity in log2 scale
  param$param_PG_intThresh = yc$getYAML("File$ProteinGroups$IntensityThreshLog2_num", param$param_PG_intThresh, 1, 100)
  
  ## get scoring threshold (upper limit)
  param$param_EV_protThresh = yc$getYAML("File$Evidence$ProteinCountThresh_num", param$param_EV_protThresh, 1, 1e5)
  
  ## default median intensity in log2 scale
  param$param_EV_intThresh = yc$getYAML("File$Evidence$IntensityThreshLog2_num", param$param_EV_intThresh, 1, 100)
  
  ## get scoring threshold (upper limit)
  param$param_EV_pepThresh = yc$getYAML("File$Evidence$PeptideCountThresh_num", param$param_EV_pepThresh, 1, 1e6)
  
  ### warn of special contaminants!
  ## these need to be in FASTA headers (description is not enough)!
  ## syntax:  list( contaminant1 = c(name, threshold), contaminant2 = c(...), ...)
  ##
  ##  if within the YAML file
  ##    SpecialContaminants: no
  ##  is set, then 'yaml_contaminants' will be 'FALSE'
  ##
  ##contaminant_default = FALSE ## to switch it off by default
  
  param$yaml_contaminants = yc$getYAML("File$Evidence$SpecialContaminants", param$yaml_contaminants)
  
  param$param_EV_MatchingTolerance = yc$getYAML("File$Evidence$MQpar_MatchingTimeWindow_num", param$param_EV_MatchingTolerance)
  if (param$param_useMQPAR && !is.null(txt_files)) {
    v = getMQPARValue(txt_files$mqpar, "matchingTimeWindow") ## will also warn() if file is missing
    if (!is.null(v)) {
      param$param_EV_MatchingTolerance = yc$setYAML("File$Evidence$MQpar_MatchingTimeWindow_num", as.numeric(v))
    }
  }
  param$param_evd_mbr = yc$getYAML("File$Evidence$MatchBetweenRuns_wA", param$param_evd_mbr)
  
  param$param_EV_PrecursorTolPPM = yc$getYAML("File$Evidence$MQpar_firstSearchTol_num", param$param_EV_PrecursorTolPPM)
  if (param$param_useMQPAR && !is.null(txt_files)) {
    v = getMQPARValue(txt_files$mqpar, "firstSearchTol") ## will also warn() if file is missing
    if (!is.null(v)) {
      param$param_EV_PrecursorTolPPM = yc$setYAML("File$Evidence$MQpar_firstSearchTol_num", as.numeric(v))
    }
  }
  
  param$param_EV_PrecursorOutOfCalSD = yc$getYAML("File$Evidence$firstSearch_outOfCalWarnSD_num", param$param_EV_PrecursorOutOfCalSD)
  
  ## we do not dare to have a default, since it ranges from 6 - 4.5 ppm across MQ versions
  param$param_EV_PrecursorTolPPMmainSearch = yc$getYAML("File$Evidence$MQpar_mainSearchTol_num", param$param_EV_PrecursorTolPPMmainSearch)
  if (param$param_useMQPAR && !is.null(txt_files)) {
    v = getMQPARValue(txt_files$mqpar, "mainSearchTol") ## will also warn() if file is missing
    if (!is.null(v)) {
      param$param_EV_PrecursorTolPPMmainSearch = yc$setYAML("File$Evidence$MQpar_mainSearchTol_num", as.numeric(v))
    }
  }
  if (is.na(param$param_EV_PrecursorTolPPMmainSearch))
  {
    warning("PTXQC: Cannot draw borders for calibrated mass error, since neither 'File$Evidence$MQpar_mainSearchTol_num' is set nor a mqpar.xml file is present!", immediate. = TRUE)
  }
  
  param$param_MSMSScans_ionInjThresh = yc$getYAML("File$MsMsScans$IonInjectionThresh_num", param$param_MSMSScans_ionInjThresh, 0, 200)
  
  param$param_OutputFormats = yc$getYAML("PTXQC$OutputFormats", param$param_OutputFormats)
  
  param$param_PageNumbers = yc$getYAML("PTXQC$PlainPDF$AddPageNumbers", param$param_PageNumbers)
  
  
  ####
  ####  prepare the metrics
  ####
  lst_qcMetrics = getMetricsObjects(DEBUG_PTXQC)
  df.meta = getMetaData(lst_qcMetrics = lst_qcMetrics)
  df.meta
  ## reorder metrics (required for indexing below!)
  lst_qcMetrics_ord = lst_qcMetrics[df.meta$.id]
  
  ## write/update order from YAML
  for (i in 1:nrow(df.meta))
  {
    #cat(paste("meta id: ", df.meta$.id[i], "\n"))
    pname = paste0("order$", df.meta$.id[i])
    pval = df.meta$order[i]
    
    ##check for shiny metrics input
    if(!is.null(metrics)){
      if(!(df.meta$.id[i] %in% metrics)) pval = -1; ## omit the metric if not listed`
    }
    
    param_v = yc$getYAML(pname, pval) ## read value (if present), or return `pval`
    ## update
    if (is.numeric(param_v)) {
      lst_qcMetrics_ord[[i]]$orderNr = param_v  # for some reason, lst_qcMetrics[[df.meta$.id]] does not work
    } else {
      stop("YAML param '", pname, "' is not numeric (", param_v, "). Please fix the YAML configuration!")
    }
  }
  ## re-read meta (new ordering?)
  df.meta = getMetaData(lst_qcMetrics = lst_qcMetrics)
  ## reorder metrics (again; after param update)
  lst_qcMetrics_ord = lst_qcMetrics[df.meta$.id]
  
  return(list("yc" = yc, "param" = param, "lst_qcMetrics" = lst_qcMetrics_ord))
  
  
}
