##
## Minimal usage example after checkout of PTXQC/mzTab_support git repo
##
## install Pandoc, see https://github.com/cbielow/PTXQC/tree/mzTab_support#installation
##

if (0) {
  ## just do this once
  install.packages(c("rmarkdown", "knitr", "ggplot2", "reshape", "plyr", "yaml", "RColorBrewer", 
                     "proto", "ggdendro", "gtable", "grid", "gridExtra", "kableExtra", "data.table"))
}
library(rmarkdown)
library(knitr)
library(ggplot2)
library(reshape2)
library(plyr)
library(yaml)
library(RColorBrewer)
library(proto)
library(ggdendro)
library(gtable)
library(grid)
library(gridExtra)
library(kableExtra)
library(data.table)

## your path
basedir = "C:/dev/PTXQC/R"  ## this is the path to the git repo. subdirectories are, e.g. ./R, or ./inst
setwd("C:/dev")

##
## source all R files in ./R
##
all.sources = list.files(paste0(basedir, "/R"), pattern = "*.R", full.names = TRUE)
source(paste0(basedir, "/qcMetric.R"))
for (s in all.sources) source(s, keep.source = TRUE)

yamlfile = "C:/Users/Jule/Documents/Bioinformatik/Softwarepraktikum-OpenMS/R/report_v0.92.6.yaml"
yaml_obj = list()
yaml_obj = yaml.load_file(input = yamlfile)
#yaml_obj$PTXQC$OutputFormats = "plainPDF" ## use this if you do not have Pandoc (omits Html report; just PDF)
yaml_obj

mztab_file = NULL
report_filenames = NULL
txt_folder = NULL
DEBUG_PTXQC = TRUE
mztab_file = "C:/Users/Jule/Documents/Bioinformatik/Softwarepraktikum-OpenMS/R/crap.mzTab"

createReport(txt_folder = NULL, mztab_file = mztab_file,  yaml_obj)

