## load packages
require(PTXQC)
require(yaml)
## the next require() is needed to prevent a spurious error in certain R versions (might be a bug in R or a package)
## error message is:
##    Error in Scales$new : could not find function "loadMethod"
require(methods)

argv = commandArgs(TRUE)
#argv = c('C:\\projects\\QC\\data\\txt_SILAC')
#cat("Command line args are:\n")
#cat(paste(argv, collapse="\n", sep=""))
#cat("\n")

if(!(length(argv) %in% 1:2))
{
  stop("Wrong number of parameters!\n",
    "Received: \n - ", paste(argv, collapse="\n - ", sep=""),
    "\n\nUsage: <thisScript.R> <PATH_TO_TXT> [<PATH_TO_YAML_CONFIG>]\n");
}

PATH_TO_TXT = argv[1]
if (!file.info(PATH_TO_TXT)$isdir)
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
  fh_out = getReportFilenames(PATH_TO_TXT)
  if (file.exists(fh_out$yaml_file))
  {
    cat("\nUsing YAML config already present in target directory ...\n")
    YAML_CONFIG = yaml.load_file(input = fh_out$yaml_file)
  }
}

r = try(createReport(PATH_TO_TXT, YAML_CONFIG))




