## load packages
library(PTXQC)
library(yaml)
## the next library() is needed to prevent a spurious error in certain R versions (might be a bug in R or a package)
## error message is:
##    Error in Scales$new : could not find function "loadMethod"
library(methods)

argv = commandArgs(TRUE)
#argv = c('Z:\\projects\\QC\\PTXQC\\data\\_txt_withMBR_withFractions_MQ15')
cat("Command line args are:\n")
cat(paste(argv, collapse="\n", sep=""))
cat("\n")

if(!(length(argv) %in% 1:2))
{
  stop("Wrong number of parameters!\n",
    "Received: \n - ", paste(argv, collapse="\n - ", sep=""),
    "\n\nUsage: <thisScript.R> <PATH_TO_TXT> [<PATH_TO_YAML_CONFIG>]\n");
}

PATH_TO_TXT = argv[1]
fi = file.info(PATH_TO_TXT)
if (is.na(fi$isdir) || !fi$isdir)
{
  stop(paste0("Argument '", PATH_TO_TXT, "' is not a valid directory\n"));
}

YAML_CONFIG = list()
if (length(argv)==2 && nchar(argv[2])>0)
{ ## YAML was passed via command line
  cat("\nUsing YAML config provided via command line ...\n")
  YAML_CONFIG = yaml.load_file(input = argv[2])
} else {
  ## use a YAML config inside the target directory if present
  rprt_fns = getReportFilenames(PATH_TO_TXT)
  if (file.exists(rprt_fns$yaml_file))
  {
    cat("\nUsing YAML config already present in target directory ...\n")
    YAML_CONFIG = yaml.load_file(input = rprt_fns$yaml_file)
  }
}
## use YAML_CONFIG to get output-filenames (which contains the log-file name)
yc = YAMLClass$new(YAML_CONFIG)
use_extended_reportname = yc$getYAML("PTXQC$ReportFilename$extended", TRUE)
rprt_fns = getReportFilenames(PATH_TO_TXT, use_extended_reportname)

sink(rprt_fns$log_file, split = TRUE) ## log output to file
output_files = try(createReport(PATH_TO_TXT, NULL, YAML_CONFIG, rprt_fns))
sink() ## undo sink()



