########################################################################################################################
#'
#' Class which holds the data from a .mzTab.
#' 
#' Memberfunction read.mzTab() takes path to .mzTab file
#' Data from .mzTab sections is stored in seperate data.frames
#'
#' Usage: foo <- mzTab$new()
#'        foo <- foo$read.mzTab("file.mztab")
#'        file.mztab <- foo$ALL
#'        metadata_from_file.mzTab <- foo$MTD 
#'        

mzTab = setRefClass("mzTab",
                    fields = list(ALL="data.frame",
                                  MTD="data.frame", ##
                                  PRT="data.frame", ##
                                  PEP="data.frame", ##
                                  PSM="data.frame", ##
                                  SML="data.frame", ##
                                  rawFileMapping="data.frame" ##
                    ),
                    methods = list(
                      initialize = function(){
                        .self$ALL = data.frame();
                        .self$MTD = data.frame();
                        .self$PRT = data.frame();
                        .self$PEP = data.frame();
                        .self$PSM = data.frame();
                        .self$SML = data.frame();
                        .self$rawFileMapping <- data.frame();
                      },
                      #' function to read mzTab
                      read.mzTab = function(file){
                        
                        #get max number of columns
                        nr <- max(count.fields(file, sep = "\t", quote = "'"))
                        cn <- paste0("c", 1:nr)
                        data <- read.table(file=file, col.names = cn , sep="\t", fill=TRUE, na.strings=c("null","NA", ""), stringsAsFactors = FALSE)
                        
                        #.self$ALL <- data
                        
                        .self$MTD <- data[which(data[,1]=="MTD"),]
                        
                        .self$PRT <- data[which(data[,1]=="PRT"),]
                        if (length(.self$PRT) != 0) {
                          colnames(.self$PRT) <- as.character(data[which(data[,1]=="PRH"),])
                        }
                        
                        .self$PEP <- data[which(data[,1]=="PEP"),]
                        if (length(.self$PEP) != 0) {
                          colnames(.self$PEP) <- as.character(data[which(data[,1]=="PEH"),])
                        }
                        
                        .self$PSM <- data[which(data[,1]=="PSM"),]
                        if (length(.self$PSM) != 0) {
                          colnames(.self$PSM) <- as.character(data[which(data[,1]=="PSH"),])
                        }
                        
                        .self$SML <- data[which(data[,1]=="SML"),]
                        if (length(.self$SML) != 0) {
                          colnames(.self$SML) <- as.character(data[which(data[,1]=="SMH"),])
                        }
                        
                      },
                      #### function to generate evidence dataframe (df_evd) from an mzTab object ---------------------
                      get_df_evd = function(.self,add_fs_col){
                        
                     
                        
                        # the base for the df_evd is the PSM section of the mzTab file
                        # helper function 
                        pasteUnique <- function(vector){
                          return(paste(unique(vector),collapse=";"))
                        }
                        
                        if("PSM_ID" %in% colnames(.self$PSM)){
                          df_evd <- .self$PSM
                          df_evd <- aggregate(.~PSM_ID, pasteUnique, na.action = na.pass, data = df_evd[which(colnames(df_evd)!="")])
                        }else{
                          warning("Dataframe for EVD Metrics could not be created. Column named PSM_ID was not found.")
                        }
                        
                        #### rename the columns (column names determined by metric functions) -------------------------
                        
                        df_evd <- rename(df_evd, c("sequence" = "modified.sequence","accession" = "proteins","PSM_ID" = "id","spectra_ref" = "raw.file",
                                                   "exp_mass_to_charge" = "m.z", "retention_time" = "retention.time", "search_engine_score[1]" = "score") )
                        
                        #### set datatypes for cloumns ------------------------
                        
                        try(df_evd$retention.time <- as.numeric(df_evd$retention.time))
                        try(df_evd$id <- as.numeric(df_evd$id))
                        try(df_evd$m.z <- as.numeric(df_evd$m.z))
                        try(df_evd$charge <- as.numeric(df_evd$charge))
                        try(df_evd$score <- as.numeric(df_evd$score))
                        
                        #### add "extra" columns needed for the evd metrics but not immediate present in the mzTabfile -----------
                        if("proteins" %in% colnames(df_evd))
                        {
                          "Problem: mzTab unterscheidet nicht zwischen leading proteins und proteins, leading proteins ist im moment einfach das 
                           erste Protein aus Proteins"
                          df_evd$leading.proteins =  sapply(X = strsplit(x = df_evd$proteins ,split = ";"),FUN = function(list){return(list[[1]][1])} )
                          # protein name info is extracted from the PRT section of the mzTab file
                          protein.names = sapply(X = .self$PSM$accession,
                                                      FUN = function(accession){return(.self$PRT$description[which(.self$PRT$accession == accession)[1]])},
                                                      USE.NAMES = TRUE)
                          # map protein names onto df_evd
                          df_evd$protein.names = sapply(X = df_evd$proteins,
                                                        FUN = function(protein){paste0(protein.names[unlist(strsplit(protein,";"))],collapse = ";")},
                                                        USE.NAMES = FALSE)
                          rm(protein.names)
                        }
                        
                        if(all(c("charge","m.z") %in% colnames(df_evd)))
                        {
                          df_evd$mass = df_evd$charge * df_evd$m.z
                        }
                        
                        if("raw.file" %in% colnames(df_evd))
                        {
                          df_evd$raw.file = gsub(pattern = ":.*$",replacement = "",x = df_evd$raw.file)
                          # test cases had file names with characters that cause proplems with regular expressions (e.g. [,(...)
                          df_evd$raw.file =  gsub(pattern = "\\[",replacement = " ",x = df_evd$raw.file)
                          df_evd$raw.file =  gsub(pattern = "\\]",replacement = " ",x = df_evd$raw.file)
                          
                          df_evd$fc.raw.file = .self$get_file_mapping(add_fs_col,df_evd,.self)
                          
                        }
                        
                        return(df_evd)
                      },
                      #### function to generate summary dataframe (df_smy) from an mzTab object ---------------------
                      get_df_smy = function(.self,add_fs_col){},
                      #### function to generate parAll dataframe (df_parAll) from an mzTab object ---------------------
                      get_df_parAll = function(.self){},
                      #### function to generate protein groups dataframe (df_pg) from an mzTab object ---------------------
                      df_pg = function(.self){},
                      #### function to generate msms dataframe (df_msms_s) from an mzTab object ---------------------
                      df_msms_s = function(.self){},
                      #### function to generate msms dataframe (df_msmsScan_h) from an mzTab object ---------------------
                      df_msmsScan_h = function(.self){},
                      get_file_mapping = function(add_fs_col,df,.self){
                          # code für Mapping aus dem MQReader (Leicht abgeändert)
                          if (add_fs_col & "raw.file" %in% colnames(df))
                          {
                            ## check if we already have a mapping
                            if (nrow(.self$rawFileMapping) == 0)
                            {
                              .self$rawFileMapping = getShortNames(unique(df$raw.file), add_fs_col)
                              ## indicate to outside that a new table is ready
                              ##.self$mapping.creation = .$getMappingCreation()['auto']
                            }
                            cat(paste0("Adding fc.raw.file column ..."))
                            
                            ## do the mapping
                            df$fc.raw.file = as.factor(.self$rawFileMapping$to[match(df$raw.file, .self$rawFileMapping$from)])
                            ## check for NA's
                            if (any(is.na(df$fc.raw.file)))
                            { ## if mapping is incomplete
                              #missing = unique(df_evd$raw.file[is.na(df_evd$fc.raw.file)])
                              # if (.$mapping.creation == .$getMappingCreation()['user'])
                              # {
                              #   ## the user has re-run MaxQuant with more Raw files,
                              #   ## but the old _filename_sort.txt file was used to read the (now incomplete mapping)
                              #   warning("Incomplete mapping file '", .$external.mapping.file, "'.\nAugmenting shortened Raw files:\n  " %+%
                              #             paste(missing, collapse="\n  ", sep="") %+% ".\nEdit the table if necessary and re-run PTXQC.")
                              #   ## augment
                              #   addon = .$getShortNames(missing, add_fs_col, nrow(.$raw_file_mapping) + 1)
                              #   .$raw_file_mapping = rbind(.$raw_file_mapping,
                              #                              addon)
                              # } else {
                              #   stop("Hithero unknown Raw files: " %+% paste(missing, collapse=", ", sep="") %+% " occurred in file '" %+% file %+% "' which were not present in previous txt files.")
                              # }
                            }
                            cat(paste0(" done\n"))
                          }
                        return(df$fc.raw.file)
                        }
                    )
)