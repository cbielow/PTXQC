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
readValue <-  function(qcMetric, qc_cv){
  qc_cv$value <- NA
  qc_cv$quality_type <- NA
  qc_cv$rawfile <- NA
  
  new_metric_index = length(qc_cv$name) + 1
  for(klass in qcMetric){
    index <- 1
    for(metric_name in klass$qcCV){
      ##If this metric exists in data, other information about this metric is extracted from data, 
      ##or if not, a new line of data is added to data and a warning is given

      if(metric_name %in% qc_cv$name) {
        locate = which(qc_cv$name == metric_name)
        qc_cv$value[locate] = rjson::toJSON(klass$mzQCdata[index])
        qc_cv$quality_type[locate] = klass$quality_type[index]
        qc_cv$rawfile[locate] = klass$raw[index]
      }else{
        warning("There is no this metric in cv Term of mzQC, please write to the mailing list psidev-qc-dev@lists.sourceforge.net to update mzQC-metrics list.")
        name = metric_name
        value =  rjson::toJSON(klass$mzQCdata[index])
        quality_type = klass$quality_type[index]
        rawfile = klass$raw[index]
        qc_cv[new_metric_index, ] <- c(NA, NA, metric_name, NA, value, quality_type, rawfile)
        new_metric_index = new_metric_index + 1
      }
      index = index + 1
    }
  }
 
  
  ##delete the row where the value of data is NA 
  qc_cv <- qc_cv[! is.na(qc_cv$value), ]
  qc_cv$quality_type <- unlist(qc_cv$quality_type)
  qc_cv$rawfile <- unlist(qc_cv$rawfile)
 
  return(qc_cv)
}

#'
#' Extract raw file information from the source data
#' 
#' 
#' 
#' @param input_rawData A list containing all the source data
#' @param rawData_name A list containing the names of all the source data
#' @return a data.frame
#' 
#' 
#' 
#' 
get_inputInfo <- function(input_rawData, rawData_name){
  inputfile_info <- list()
  for (index in 1:length(input_rawData)) {
    inputfile_info[[index]] <- list(rawData = "", fc.raw.file = list(), raw.file = list())
    inputfile_info[[index]]$rawData = rawData_name[index]
    x <- input_rawData[index][[1]]
    if("fc.raw.file" %in% colnames(x)) {
      inputfile_info[[index]]$fc.raw.file = list(name = stringr::str_sub(x$fc.raw.file,0, stringr::str_locate(x$fc.raw.file, "_0")[1]-1)[1],
                                             location = stringr::str_sub(x$fc.raw.file,0, stringr::str_locate(x$fc.raw.file, "_0")[1]-1)[1],
                                             fileformat = list(accession = "", name = "raw file"))
    }
    if("raw.file" %in% colnames(x)){
      inputfile_info[[index]]$raw.file = list(name = stringr::str_sub(x$raw.file,0, stringr::str_locate(x$raw.file, "_0")[1]-1)[1],
                                          location = stringr::str_sub(x$raw.file,0, stringr::str_locate(x$raw.file, "_0")[1]-1)[1],
                                          fileformat = list(accession = "", name = "raw file"))
    }
    if(rawData_name[index] == "d_parAll") inputfile_info[[index]]$raw.file = list(name = "parameter.txt",
                                                                                  location = "parameter.txt",
                                                                                  fileformat = list(accession = "", name = "raw file"))
    if(rawData_name[index] == "df_pg") inputfile_info[[index]][c("fc.raw.file", "raw.file")] = inputfile_info[[1]][c("fc.raw.file", "raw.file")]
  }
  return(inputfile_info)
}


#'
#' By searching the name of raw file, it is to extract the information we needed from the table that generated by function get_inputInfo()
#' 
#' 
#' 
#' @param inputfiles_info  return value of function named get-infileInfo()
#' @param inputfile_name name of a single source data
#' @param rawData_name A list containing the names of all the source data
#' @return a list
#' 
#' 
#' 
#' 
get_inputfiles <- function(inputfiles_info, inputfile_name, rawData_name){
  inputfiles <- inputfiles_info[[which(inputfile_name == rawData_name)]][c("fc.raw.file","raw.file")]
  inputfiles <- lapply(inputfiles, function(x){
    if(length(x) == 0) return(NA)
    else return(x)
  })
  inputfiles <- inputfiles[ ! is.na(inputfiles)]
  return(inputfiles)
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
#' @param is_run Determine whether it is runQualities
#' @return a list
#' 
#' @export
#' 
#' 
#' 
getQualities <- function(Qualities, data_input, data_analysis,rawData_name){
  if(length(Qualities[ , 1]) != 0 ){
    
    Qualities <- Qualities[! is.na(Qualities$value), ]
    Qualities <- lapply(transform_data(data = Qualities), function(x){
      Quality <- list(metadata = list(inputFiles = get_inputfiles(data_input, x$raw, rawData_name), analysisSoftware = data_analysis),
                      qualityMetrics = list(qualityMetric = x[ ,c("cvRef"," accession", "name", "unit", "value")]))
      return(Quality)
    })
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
                       data_analysis,rawData_name){

  mzQC <- list(version = version, creationDate =creationDate, 
               runQualities = getQualities(runQualities,data_input,data_analysis,rawData_name),
               setQualities = getQualities(setQualities,data_input,data_analysis,rawData_name), 
               controlledVocabularies = controlledVocabularies)
    
  cat(rjson::toJSON((list(mzQC = mzQC)), indent = 3), file = file)
  
  return(file)
}
