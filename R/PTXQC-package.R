#' PTXQC: A package for computing Quality Control (QC) metrics for Proteomics (PTX)
#'
#' The following sections describe the main components:
#' 
#' @section Input:
#' Valid input data are either the files from MaxQuant's .txt folder (all versions from MaxQuant >= 1.0 upwards are supported)
#' or a single mzTab file. All mzTab files will work, but most metrics can be obtained from OpenMS' mzTab as produced
#' by the QualityControl TOPP tool (from OpenMS 2.5 onwards).
#' 
#' @section Important functions:
#' The central function of this package is called \code{\link{createReport}} and it accepts either MaxQuant or mzTab data, along with 
#' a configuration (optional).
#' There is a parser for mzTab \code{\link{MzTabReader}} and MaxQuant txt files \code{\link{MQDataReader}}, as well as a plethora of QC metrics
#' derived from a common \code{\link{qcMetric}} class and scoring functions \code{qual...}, e.g. \code{\link{qualGaussDev}}.
#' 
#' @section Configuration:
#' The user can modify the behaviour of PTXQC, e.g. to enable/disable certain metrics or change scoring thresholds, via a YAML object/file.
#' By default a Yaml file is written automatically side-by-side to the input files upon running PTXQC for the first time on a particular input.
#' A custom Yaml object can be passed to the main \code{\link{createReport}} function for customization. 
#' Use \code{yaml::yaml.load_file(input = 'myYAML.yaml')} to load an existing file and pass the Yaml object along.
#' 
#' @section Output:
#' Either a PDF and/or Html report which contains QC plots and a description of the metrics.
#'
"_PACKAGE"
#' @name PTXQC
#' 
#' @import data.table
#' @import ggplot2
#' @import ggdendro
#' @import grid
#' @import gridExtra
#' @import grDevices
#' @import gtable
#' @import knitr
#' @import methods
#' @import plyr
#' @import R6
#' @import R6P
#' @import RColorBrewer
#' @rawNamespace import(reshape2, except = c(dcast, melt))
#' @import rmarkdown
#' @importFrom seqinr circle
#' @import stats
#' @import utils
#' @import UpSetR
#' @import xml2
#' @import yaml
#' @importFrom rlang .data 
#' 
NULL