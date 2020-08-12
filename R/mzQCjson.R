#'
#' Get the information of each cv Term from a obo file.
#' 
#' 
#' 
#' @param file "xxx.obo"
#' @return A list containing cv Term's information
#' 
#' @export
#' 
#' 
#' 
read_cvTerm <- function(file){
  obo <- data.frame(read.csv(file = file, header = FALSE))
  cv_Term <- list(cvRef = list(),accession= list(), name = list(), units = list())
  term_index <- 0
  for (x in obo$V1) {
    if(x == "[Term]"){
      term_index = term_index + 1
      cv_Term$units[[term_index]] = list()
    }
    if(substring(x,1,2)  == "id"){
      cv_Term$cvRef[term_index] = substring(stringr::str_sub(x, stringr::str_locate(x, ":")[1]+2), 1,2)
      cv_Term$accession[term_index] = stringr::str_sub(x, stringr::str_locate(x, ":")[1]+2)
    }
    if(substring(x,1,4)  == "name"){
      cv_Term$name[term_index] = stringr::str_sub(x, stringr::str_locate(x, ":")[1]+2)

    }
    if(substring(x,1,4)  == "is_a"){
      ## Connect multiple unit through commas
      unit <- stringr::str_sub(x, stringr::str_locate(x, ":")[1]+2)
      cv_Term$units[[term_index]] = append(cv_Term$units[[term_index]], rjson::toJSON(list(cvRef = stringr::str_sub(unit, 1, stringr::str_locate(unit, "\\:")[1]-1),
                                                                         accession = stringr::str_sub(unit, 1, stringr::str_locate(unit,"\\ !")[1]-1),
                                                                         name = stringr::str_sub(unit, stringr::str_locate(unit, "\\! ")[1]+2)), indent = 3))
    }
  }
 
  return(cv_Term)
}


#'
#' Get all information needed in metric (cvRef, accession, name, value, unit, type of quality, input_file),
#' type of quality and input_file act on the connection to the other two tables(inputfile_info, analysisSoftware)
#' that belong to metadata
#' 
#' 
#' 
#' @param qcMetric A list containing all metric information
#' @param qc_cv data frame of all cv Parameter
#' @param df_evd data from evidence.txt
#' @param df_msms data from msms.txt
#' @param df_msmsScans data from msmsScans.txt
#' @param d_parAll data from parameters.txt
#' @param df_pg data from proteinGroups.txt
#' @param d_smy data from summary.txt
#' @return data frame of all metric that appears in PTXQC
#' 
#' @export
#' 
#' 
#' 
readvalue <-  function(qcMetric, qc_cv, df_evd, df_msms, df_msmsScans, d_parAll, df_pg, d_smy){
  qc_cv$value <- NA
  qc_cv$quality_type <- NA
  qc_cv$is_raw <- FALSE
  
  new_metric_index = length(qc_cv$name) + 1
  for(klass in qcMetric){
    index <- 1
    for(metric_name in klass$qcCV){
      ##If this metric exists in data, other information about this metric is extracted from data, 
      ##or if not, a new line of data is added to data and a warning is given

      if(metric_name %in% qc_cv$name) {
        locate = which(qc_cv$name == metric_name)
        qc_cv$value[locate] = rjson::toJSON(klass$mzQCdata[index])
        qc_cv$quality_type[locate] = get_raw(klass$raw, df_evd, df_msms, df_msmsScans, d_parAll, df_pg, d_smy)[1]
        qc_cv$is_raw[locate] = get_raw(klass$raw, df_evd, df_msms, df_msmsScans, d_parAll, df_pg, d_smy)[2]
      }else{
        warning("There is no this metric in cv Term of mzQC, please write to the mailing list psidev-qc-dev@lists.sourceforge.net to update mzQC-metrics list.")
        name = metric_name
        value =  rjson::toJSON(klass$mzQCdata[index])
        quality_type = get_raw(klass$raw, df_evd, df_msms, df_msmsScans, d_parAll, df_pg, d_smy)[1]
        is_raw = get_raw(klass$raw, df_evd, df_msms, df_msmsScans, d_parAll, df_pg, d_smy)[2]
        qc_cv[new_metric_index, ] <- c(NA, NA, metric_name, NA, value, quality_type, is_raw)
        new_metric_index = new_metric_index + 1
      }
      index = index + 1
    }
  }
 
  
  ##delete the row where the value of data is NA 
  qc_cv <- qc_cv[! is.na(qc_cv$value), ]
 
  return(qc_cv)
}

#'
#' Reads a string and converts it to data with the same name before getting the data it needs
#' 
#' 
#' 
#' @param raw A String(e.g: "df_evd")
#' @param df_evd data from evidence.txt
#' @param df_msms data from msms.txt
#' @param df_msmsScans data from msmsScans.txt
#' @param d_parAll data from parameters.txt
#' @param df_pg data from proteinGroups.txt
#' @param d_smy data from summary.txt
#' @return a character
#' 
#' 
#' 
#'
get_raw <- function(raw, df_evd, df_msms, df_msmsScans, d_parAll, df_pg, d_smy){
  if(raw == "df_evd") return(get_input(df_evd))
  if(raw == "df_msms") return(get_input(df_msms))
  if(raw == "df_msmsScans") return(get_input(df_msmsScans))
  if(raw == "df_mqpar") return(get_input(d_parAll))
  if(raw == "df_pg") return(get_input(df_pg))
  if(raw == "df_summary") return(get_input(d_smy))
}

#'
#' Check the RAW File information in a data
#' 
#' 
#' 
#' @param input_data A data.frame 
#' @return a character contains type of quality(run-/setQuality) and is_raw(is_raw checks whether the data contains a RAW file)
#' 
#' 
#' 
#'
get_input <- function(input_data){
  is_raw <- FALSE
  if("fc.raw.file" %in% colnames(input_data)) {
    fc <- unique(as.list(input_data$fc.raw.file))
    if(length(fc) > 1){
      quality_type <- "setQuality"
      is_raw <- TRUE
    }else{
      quality_type <- "runQuality"
      is_raw <- TRUE
    }
  }else if("raw.file" %in% colnames(input_data)){
    raw <- unique(as.list(input_data$raw.file))
    if(length(raw) > 1){
      quality_type <- "setQuality"
      is_raw <- TRUE

    }else{
      quality_type <- "runQuality"
      is_raw <- TRUE
    }
  }else{
    quality_type <- "runQuality"
  }
  return(c(quality_type, is_raw))
  
}


#'
#' Extract the data of each row in data and store it in a list
#' 
#' 
#' 
#' @param data A data.frame 
#' @return a list
#' 
#' 
#' 
#' 
transform_data <-function(data){
  trans <- list()
  for (x in (1:length(data$cvRef))) {
    trans[x] = list(data[x, ])
  }
  return(trans)
}



#'
#' Get the Qualities(either runQualities or setQualities) composed of metadata and qualityMetrics
#' 
#' 
#' 
#' @param Qualities the filtered data frame of metric
#' @param data_input File based on which QC metrics were generated
#'                   use analysisSofteware to connect data_input and data_analysis
#' @param data_analysis Software tools used to generate the QC metrics
#' @return a list
#' 
#' @export
#' 
#' 
#' 
getQualities <- function(Qualities, data_input, data_analysis){
  if(length(Qualities[ , 1]) != 0 ){
    if(TRUE %in% Qualities$is_raw){
      metaData <- list(inputFiles = list(inputFile = data_input), analysisSoftware = data_analysis)
      Qualities <- Qualities[! is.na(Qualities$value), ]
      Qualities <- list(metaData = metaData, qualityMetrics = list(qualityMetric = transform_data(data = Qualities[ ,c("cvRef"," accession", "name", "unit", "value")])))
    }else{
      Qualities <- Qualities[! is.na(Qualities$value), ]
      Qualities <- list(metaData = list(), qualityMetrics = list(qualityMetric = transform_data(data = Qualities[ ,c("cvRef"," accession", "name", "unit", "value")])))
    }
  }else{
    Qualities <- list()
  }
  return(Qualities)
}


#'
#' Write all mzQC related information to the JSON file
#' 
#' 
#' 
#' @param version version of mzQC
#' @param creationDate JSON file generation time
#' @param runQualities data frame of runQuality elements
#' @param setQualities data frame of setQuality elements
#' @param controlledVocabularies use to refer to the source of the used cv terms in the qualityMetric objects
#' @param file output file path
#' @param data_input File based on which QC metrics were generated
#'                   use analysisSofteware to connect data_input and data_analysis
#' @param data_analysis Software tools used to generate the QC metrics
#' @return a list
#' 
#' @export
#' 
#' 
#' 
createmzQC <- function(version, creationDate, runQualities, setQualities, controlledVocabularies, file, data_input, 
                       data_analysis){

  mzQC <- list(version = version, creationDate =creationDate, 
               runQualities = getQualities(runQualities,data_input,data_analysis),
               setQualities = getQualities(setQualities,data_input,data_analysis), 
               controlledVocabularies = controlledVocabularies)
    
  cat(rjson::toJSON((list(mzQC = mzQC)), indent = 3), file = file)
  
  return(file)
}
