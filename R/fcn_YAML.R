
#yamlConfig = "foo:\n bar: 456 \n"
#yamlObj = yaml.load(yamlConfig)
#eval(parse(text="yamlObj$bla$foo = 2"))
#yamlObj
#useP = getYAML(yamlObj, "SEC_Parameters$use", TRUE)

#' Query a YAML object for a certain parameter.
#' 
#' If the object has the param, then return it.
#' If the param is unknown, create it with the given default value and return the default.
#' 
#' @note This function has a side-effect: it updates the parameter passed as 'config' in the calling environment!!!
#' 
#' @param config     A Yaml object as created by \code{\link[yaml]{yaml.load}}
#' @param param_name A string which holds the param name, i.e. referring to a variable within 'config', e.g. "proteinGroup$doPlot"
#' @param default    If seeked value is unknown, this will be used as new value
#' 
#' @return The stored value (might be the default, if no value was present)
#'
#' @export
#'
getYAML = function(config, param_name, default)
{
  orig_var = deparse(substitute(config)) ## do this BEFORE accessing 'config'
  #print(orig_var)
  pval = eval(parse(text=paste0("config$", param_name))) 
  if (is.null(pval))
  { ## param not known yet --> add
    expr = paste0("config$", param_name, " = ", quote(default))
    #cat("parsing expr: ", expr)
    eval(parse(text=expr)) 
    #print("new obj:");  #print(config)
    assign(x=orig_var, value=config, pos=parent.frame(1))
    return (default)
  } else
  {
    return (pval)
  }
}


#' Set a YAML parameter to a certain value.
#' 
#' @note This function has a side-effect: it updates the parameter passed as 'config' in the calling environment!!!
#' 
#' @param config     A Yaml object as created by \code{\link[yaml]{yaml.load}}
#' @param param_name A string which holds the param name, i.e. referring to a variable within 'config', e.g. "proteinGroup$doPlot"
#' @param value      New value to set
#' 
#' @return The value
#'
#' @export
#'
setYAML = function(config, param_name, value)
{
  orig_var = deparse(substitute(config)) ## do this BEFORE accessing 'config'
  #print(orig_var)
  expr = paste0("config$", param_name, " = ", quote(value))
  #cat("parsing expr: ", expr)
  eval(parse(text=expr)) 
  #print("new obj:");  #print(config)
  assign(x=orig_var, value=config, pos=parent.frame(1))
  cat("Setting parameter '", param_name, "' using external data to '", value, "'.\n", sep="")
  return (value)
}


#'
#' Write YAML config (including some documentation) to a YAML file
#' 
#' @param filename File to write
#' @param yaml_obj A nested list object with configuration parameters for the report.
#'                 Useful to switch off certain plots or skip entire sections.
#'                 Will be converted via 'as.yaml()'
#' @return TRUE on success
#' 
#' @importFrom yaml as.yaml
#' 
writeYAML = function(filename, yaml_obj)
{
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
# 
#
# With the exception(!) of extra 'SpecialContaminants':
#   Do NOT add extra values (since they are ignored anyway and might even destroy the integrity of this configuration file).
#   Only modify existing values, but not their names. I.e. only change 123 but not 'test'
#   test: 123
#
#
"
  cat(paste0(yaml.user.warning, as.yaml(yaml_obj)), file=filename)
  
  return (TRUE);
}