library(stringr)
library(readr)

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
      cv_Term$units[term_index] = ""
    }
    if(substring(x,1,2)  == "id"){
      cv_Term$cvRef[term_index] = substring(str_sub(x, str_locate(x, ":")[1]+2), 1,2)
      cv_Term$accession[term_index] = str_sub(x, str_locate(x, ":")[1]+2)
    }
    if(substring(x,1,4)  == "name"){
      cv_Term$name[term_index] = str_sub(x, str_locate(x, ":")[1]+2)

    }
    if(substring(x,1,4)  == "is_a"){
      ## Connect multiple unit through commas
      cv_Term$units[term_index] = paste(cv_Term$units[term_index], str_sub(x, str_locate(x, ":")[1]+2), seq=",")

    }
  }
  ##Divide the units joined by commas into a list
  cv_Term$units <- lapply(cv_Term$units, function(x){
    if(x != ""){
      
      units <- unlist(strsplit(x , split = ","))
      data <- (lapply(units, function(y){
        cvRef <- str_sub(y, 1, str_locate(y, "\\:")[1]-1)
        accession <- str_sub(y, 1, str_locate(y,"\\ !")[1]-1)
        name <- str_sub(y, str_locate(y, "\\! ")[1]+2)
        return(list(cvRef, accession, name))
      }))
      return(rjson::toJSON(data))
      
    }else{
      data <- NA
      return(data)
    }
  })
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
#' @param data data frame of all cv Parameter
#' @return data frame of all metric that appears in PTXQC
#' 
#' @export
#' 
#' 
#' 
readvalue <-  function(qcMetric, data){
  runQuality = FALSE
  setQuality = FALSE
  new_metric = length(data$cvRef) + 1
  for(klass in qcMetric){
    index <- 1
    for(metric_name in klass$qcCV){
      ##If this metric exists in data, other information about this metric is extracted from data, 
      ##or if not, a new line of data is added to data and a warning is given
      if(metric_name %in% data$name) {
        if(! is.na(data$value[which(metric_name == data$name)])) warning("Repeat metric! Each metric has only one. Please compare the data entered before with the existing data and Retain original data.")
        else {
          locate = which(data$name == metric_name)
          data$value[locate] = rjson::toJSON(klass$mzQCdata[index])
          data$quality_type[locate] = klass$quality_type[index]
          data$input_file[locate] = klass$input_file[index]
          
        }
      }else{
        warning("There is no this metric in cv Term of mzQC, please write to the mailing list psidev-qc-dev@lists.sourceforge.net to update mzQC-metrics list.")
        name = metric_name
        value =  rjson::toJSON(klass$mzQCdata[index])
        quality_type = klass$quality_type[index]
        input_file = klass$input_file[index]
        data[new_metric, ] <- c(NA, NA, metric_name, NA, value, quality_type, input_file)
        new_metric = new_metric + 1
      }
      index = index + 1
    }
  }
  ##Check whether "..Ecoli_01" or "..Ecoli_02" is included in metric's value
  for(klass in qcMetric){                                                                                                                 index <- 1
    for(mzQC in klass$mzQCdata){
      if("fc.raw.file" %in% colnames(mzQC)){
        if("..Ecoli_01" %in% mzQC$fc.raw.file || "..Ecoli_02" %in% mzQC$fc.raw.file){  
          if(unlist(klass$quality_type[index] == "runQuality")) runQuality = TRUE
          else setQuality = TRUE                                                                                                          }
      }
    index = index + 1
  }                                                                                                                               }
  
  ##delete the row where the value of data is NA 
  data <- data[! is.na(data$value), ]
  if(runQuality == TRUE) data[length(data$quality_type) + 1, ] <- c(NA, NA, NA, NA, NA, "runQuality", "Raw file: Ecoli")
  if(setQuality == TRUE) data[length(data$quality_type) + 1, ] <- c(NA, NA, NA, NA, NA, "setQuality", "Raw file: Ecoli")
  
  return(data)
}

#'
#' Extract the data of each row in data and store it in a list
#' 
#' 
#' 
#' @param data A list 
#' @return a list
#' 
#' @export
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
#' Get the Metadata composed of inputFiles and analysisSoftware
#' 
#' 
#' 
#' @param input Raw file of metric data
#' @param data_input File based on which QC metrics were generated
#'                   use analysisSofteware to connect data_input and data_analysis
#' @param data_analysis Software tools used to generate the QC metrics
#' @return a list
#' 
#' @export
#' 
#' 
#' 
getmetaData <- function(input, data_input, data_analysis){
  inputFile <- list()
  analysis <- list()
  input <- unique(input)
  for (x in (1:length(input))) {
    locate <- which(input[x] == data_input$name)
    inputFile[x] <- list(data_input[locate, ][1:3])
    analysis[x] <- list(data_analysis[which(data_input[locate,]$analysisSoftware == data_analysis$name), ])
  }
  analysis <- unique(analysis)
  
  metaData = list(inputFiles = list(inputFile = inputFile), analysisSoftware = analysis)
  return(metaData)
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
#' @param set judge whether the quality is setQuality
#' @return a list
#' 
#' @export
#' 
#' 
#' 
getQualities <- function(Qualities, data_input, data_analysis, set){
  if(length(Qualities[ , 1]) != 0 ){
    if(set == TRUE){
      input <- as.list(Qualities$input_file)
      input_file <- unlist(lapply(input, function(x){
        if (grepl(",", x)){
          input[which(x == as.list(Qualities$input_file))] <- NULL
          input <- unlist(c(input, unlist(strsplit(x , split = ","))))
        }
      }))
    }else{
      input_file <- as.list(Qualities$input_file)
    }
    metaData <- getmetaData(input = input_file, data_input, data_analysis)
    Qualities <- Qualities[! is.na(Qualities$value), ]
    Qualities <- list(metaData = metaData, qualityMetrics = list(qualityMetric = transform_data(data = Qualities[1:5])))
    return(Qualities)
  }
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
createmzQC <- function(version, creationDate, runQualities, setQualities, controlledVocabularies, file, data_input, data_analysis){
  ##if data frame of runQualities is empty, set runQualities to list()
  if(length(runQualities[ ,1]) == 0){
    mzQC <- list(version = version, creationDate =creationDate, runQualities =list(),
                 setQualities = getQualities(setQualities, data_input,data_analysis, set = TRUE), controlledVocabularies = controlledVocabularies)
    
  }else if(length(setQualities[ ,1])== 0){
    ##if data frame of setQualities is empty, set setQualities to list()
    mzQC <- list(version = version, creationDate =creationDate, runQualities = getQualities(runQualities,data_input,data_analysis, set = FALSE),
                 setQualities = list(), controlledVocabularies = controlledVocabularies)
    
  }else if(length(runQualities[ ,1]) == 0 && setQualities[ ,1] == 0) {
    mzQC <- list(version = version, creationDate =creationDate, runQualities = list(),
                 setQualities = list(), controlledVocabularies = controlledVocabularies)
    
  }else{
    mzQC <- list(version = version, creationDate =creationDate, runQualities = getQualities(runQualities,data_input,data_analysis, set = FALSE),
                 setQualities = getQualities(setQualities,data_input,data_analysis, set = TRUE), controlledVocabularies = controlledVocabularies)
    
  }
  cat(rjson::toJSON((list(mzQC = mzQC)), indent = 3), file = file)
  
  return(file)
}
