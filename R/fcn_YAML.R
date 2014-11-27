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
# Finally, numerical parameters... well.. are numerical (usually integer).
#
# Parameters whose name contains 'MQ_' should be matched to the parameters set in MaxQuant.
# E.g. 'File$Evidence$MQ_MatchingTolerance_num' is the matching (not alignment!) tolerance
#      for Match-between-Runs, and has a default of 1 (minute).
#      Older MQ version prior to MQ1.4 allowed 2 minutes.
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