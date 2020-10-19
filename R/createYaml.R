#' Creates a yaml file storing the parameters that are used for creating the PTXQC report, 
#' optionally returns these parameters and a list of qc-Metrics
#'
#' @param yc A yaml class object created by YAMLClass$new()
#' @param path path to an (empty) yaml file
#' @param param list of parameters sorted by names (important for shiny application)
#' @param DEBUG_PTXQC default FALSE
#' @param get_parameters Should the parameters and qc-Metrics be returned? Default TRUE
#' @param mztabe_mode default FALSE 
#' @param txt_files list of paths to MaxQuant files
#' @param metrics list of metric names that should be plotted (important for shiny application)
#' @return list of parameters used for creating the report and list of qc-Metrics
#' 
#'
#'
createYaml <- function(yc, path, param = NULL, DEBUG_PTXQC = FALSE, output = TRUE, MZTAB_MODE = FALSE, txt_files = NULL, metrics = NULL){

    ##
    ## YAML default config
    ##
    if(is.null(param)){
      param <- list()
      param$param_useMQPAR <- TRUE
      param$add_fs_col <- 14 
      param$id_rate_bad <- 20
      param$id_rate_great <- 35
      param$pg_ratioLabIncThresh <- 4
      param$param_PG_intThresh <- 25
      param$param_EV_protThresh <- 3500
      param$param_EV_intThresh <- 23
      param$param_EV_pepThresh <- 15000
      param$yaml_contaminants <- list("cont_MYCO" = c(name="MYCOPLASMA", threshold=1)) # name (FASTA), threshold for % of unique peptides
      param$param_EV_MatchingTolerance <- 1
      param$param_evd_mbr <- "auto"
      param$param_EV_PrecursorTolPPM <- 20
      param$param_EV_PrecursorOutOfCalSD <- 2
      param$param_EV_PrecursorTolPPMmainSearch <- NA
      param$param_MSMSScans_ionInjThresh <- 10
      param$param_OutputFormats <- c("html", "plainPDF")
      param$param_PageNumbers <- "on"
      
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
    if (param$param_useMQPAR &! MZTAB_MODE && is.null(param)) {
      v = getMQPARValue(txt_files$mqpar, "matchingTimeWindow") ## will also warn() if file is missing
      if (!is.null(v)) {
        param$param_EV_MatchingTolerance = yc$setYAML("File$Evidence$MQpar_MatchingTimeWindow_num", as.numeric(v))
      }
    }
    param$param_evd_mbr = yc$getYAML("File$Evidence$MatchBetweenRuns_wA", param$param_evd_mbr)
    
    param$param_EV_PrecursorTolPPM = yc$getYAML("File$Evidence$MQpar_firstSearchTol_num", param$param_EV_PrecursorTolPPM)
    if (param$param_useMQPAR & !MZTAB_MODE && is.null(param)) {
      v = getMQPARValue(txt_files$mqpar, "firstSearchTol") ## will also warn() if file is missing
      if (!is.null(v)) {
        param$param_EV_PrecursorTolPPM = yc$setYAML("File$Evidence$MQpar_firstSearchTol_num", as.numeric(v))
      }
    }
    
    param$param_EV_PrecursorOutOfCalSD = yc$getYAML("File$Evidence$firstSearch_outOfCalWarnSD_num", param$param_EV_PrecursorOutOfCalSD)
    
    ## we do not dare to have a default, since it ranges from 6 - 4.5 ppm across MQ versions
    param$param_EV_PrecursorTolPPMmainSearch = yc$getYAML("File$Evidence$MQpar_mainSearchTol_num", param$param_EV_PrecursorTolPPMmainSearch)
    if (param$param_useMQPAR & !MZTAB_MODE && is.null(param)) {
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
    lst_qcMetrics = PTXQC:::getMetricsObjects(DEBUG_PTXQC)
    df.meta = PTXQC:::getMetaData(lst_qcMetrics = lst_qcMetrics)
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
        if(df.meta$.id[i] %in% metrics) yc$setYAML(pname, pval)
        else yc$setYAML(pname, (-1))
      }
      
      param_v = yc$getYAML(pname, pval)
      ## update
      if (is.numeric(param_v)) {
        lst_qcMetrics_ord[[i]]$orderNr = param_v  # for some reason, lst_qcMetrics[[df.meta$.id]] does not work
      } else {
        stop("YAML param '", pname_v, "' is not numeric (", param_v, "). Please fix the YAML configuration!")
      }
    }
    ## re-read meta (new ordering?)
    df.meta = PTXQC:::getMetaData(lst_qcMetrics = lst_qcMetrics)
    ## reorder metrics (again; after param update)
    lst_qcMetrics_ord = lst_qcMetrics[df.meta$.id]
    
    ## write out the final YAML file (so users can disable metrics, if they fail)
    yc$writeYAML(path)
    
    if(output) return(list(param, lst_qcMetrics_ord))
    
    
}
