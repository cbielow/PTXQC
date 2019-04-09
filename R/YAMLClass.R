
#yamlConfig = "foo:\n bar: 456 \n"
#yamlObj = yaml.load(yamlConfig)
#eval(parse(text="yamlObj$bla$foo = 2"))
#yamlObj
#useP = getYAML(yamlObj, "SEC_Parameters$use", TRUE)

#' 
#' Query a YAML object for a certain parameter.
#' 
#' If the object has the param, then return it.
#' If the param is unknown, create it with the given default value and return the default.
#' 
#' @field yamlObj A Yaml object as created by \code{\link[yaml]{yaml.load}}
#'
#' @importFrom yaml as.yaml
#' 
#' @exportClass YAMLClass
#' @export YAMLClass
#' 
#' @examples 
#'     yc = YAMLClass$new(list())
#'     val = yc$getYAML("cat$subCat", "someDefault")
#'     val  ## someDefault
#'     val = yc$setYAML("cat$subCat", "someValue")
#'     val  ## someValue
#'     yc$getYAML("cat$subCat", "someDefault") ## still 'someValue' (since its set already)
#' 
YAMLClass = setRefClass(
  "YAMLClass",
  
  fields = list(yamlObj = "list" # A Yaml object as created by \code{\link[yaml]{yaml.load}}
  ),
  
  methods = list(
    
    ##
    ##  ctor
    ##
    initialize = function(yamlObj = list()) {
      .self$yamlObj = yamlObj;
      return(.self)
    },

    getYAML = function(param_name, default)
    {
      "Query this YAML object for a certain parameter and return its value. If it does not exist it is created with a default value."
      cat(paste0("YAML: ", param_name, " def: ", paste(default, sep="", collapse=",")))
      pval = eval(parse(text=paste0(".self$yamlObj$", param_name))) 
      if (is.null(pval))
      { ## param not known yet --> add
        expr = paste0(".self$yamlObj$", param_name, " = ", quote(default))
        eval(parse(text=expr)) 
        cat("\n")
        return (default)
      } else {
        cat(paste0(" || new val: ", paste(pval, sep="", collapse=","), "\n"))
        return (pval)
      }
    },
    
    setYAML = function(param_name, value)
    {
      "Set a YAML parameter to a certain value. Overwrites the old value or creates a new entry if hithero unknown."
      expr = paste0(".self$yamlObj$", param_name, " = ", quote(value))
      eval(parse(text=expr)) 
      return (value)
    },
    
    
    writeYAML = function(filename)
    {
      "Write YAML config (including some documentation) to a YAML file. Returns TRUE on success (always), unless writing the file generates an error."
      yaml.user.warning = 
        "# This is a configuration file for PTXQC reporting.
# One such file is generated automatically every time a report PDF is created.
# You can make a copy of this file, then modify its values and use the copy as an input to another round of report generation,
# e.g., to exclude/include certain plots or change some global settings.
#
# It is recommended to work on a copy of this file, such that this original configuration reflects the content of the report.
# 
# Note that upon report generation, this file will be overwritten again (with potentially new values, if and depending on the YAML config you provided).
#
# This file has a certain structure, which should be *retained* when editing.
# Note that each parameter level has two more spaces for indentation
# A value is assigned using a colon followed by a space, i.e. ': '
#
# The values that can be taken by a parameter are indicated by the SUFFIX of the name.
# '<name>_num' :  a numeric value is expected
# '<name>_wA'  :  'no', 'yes', 'auto' (see below)
# '<name>      :  'no', 'yes'
#
# By default (no special ending) parameters are BINARY (yes or no).
# Parameters ending in 'wA' (== 'withAuto') offer a heuristic to decide if something should be plotted.
# Finally, numerical parameters... well... are numerical (usually integer).
#
#
# The subsection 'order' enables reordering of metrics by assigning other numbers.
# The order affects both the order of columns in the heatmap, and the order of plots in the report.
# If you set a metrics order number to a negative value, the metric will not be computed!
# This is useful to switch of uninteresting or broken metrics (until a fix is rolled out).
#
# Parameters whose name starts with 'MQpar_' should be matched to the parameters set in MaxQuant.
# E.g. 'File$Evidence$MQpar_MatchingTimeWindow_num' is the matching (not alignment!) tolerance
#      for Match-between-Runs, and has a default of 1 (minute).
#      Older MQ version prior to MaxQuant 1.4 allowed 2 minutes, very recent ones use 0.7 minutes.
# The parameter 'PTXQC$UseLocalMQPar' controls if the mqpar.xml file will be looked up or not. 
# It needs to be present within this txt folder (i.e. you need to copy it here!).
# If the mqpar.xml file is not found, a warning will be issued and the internal PTXQC default will
# be used (bearing the risk of being inaccurate).
#
# Furthermore, there is the SpecialContaminants section, where you can trigger the generation of a plot
# which just shows the fraction of proteins containing a certain string.
# By default we search for 'MYCOPLASMA' in the protein identifier, i.e. your txt should be derived from a 
# run using a database containing these identifiers (the FASTA description is not enough).
# Each of the special contaminants requires has its own section (e.g. 'cont_MYCO:') with two parameters: 
# - a string          = name within the protein identifier
# - an integer number = a threshold in % which will be plotted to visually ease interpretation
# If you do not want any contaminants search, just set
#     SpecialContaminants: no
#
# With the exception(!) of extra 'SpecialContaminants':
#   Do NOT add extra values (since they are ignored anyway and might even destroy the integrity of this configuration file).
#   Only modify existing values, but not their names. I.e. only change 123 but not 'test'
#   test: 123
#
#
"
      cat(paste0(yaml.user.warning, as.yaml(.self$yamlObj)), file=filename)
      
      return (TRUE);
    }
    
  )
)