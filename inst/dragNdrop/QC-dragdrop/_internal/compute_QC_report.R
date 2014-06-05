## load packages
require(PTXQC)
#install.packages("yaml")
require(yaml)
## the next require() is needed to prevent a spurious error in certain R versions (might be a bug in R or a package)
## error message is:
##    Error in Scales$new : could not find function "loadMethod"
require(methods)

## the next require() should not be needed, since PTXQC imports it, but on some systems it seems that a subfunction 
## dispatch within 'directlabels' is not working properly. If 'directlabels' is attached, all is well. So ...
require(directlabels)


argv = commandArgs(TRUE)
#argv = c('C:\\projects\\QC\\data\\txt_SILAC', 'OFF')
#cat("Command line args are:\n")
#cat(paste(argv, collapse="\n", sep=""))
#cat("\n")

if(!(length(argv) %in% 1:2))
{
  stop("Wrong number of parameters!\nUsage: <thisScript.R> <PATH_TO_TXT> [<PATH_TO_YAML_CONFIG>]\n");
}

PATH_TO_TXT = argv[1]
if (!file.info(PATH_TO_TXT)$isdir)
{
  stop(paste0("Argument '", PATH_TO_TXT, "' is not a valid directory\n"));
}

YAML_CONFIG = list()
if (length(argv)==2 && nchar(argv[2])>0)
{
  YAML_CONFIG = yaml.load_file(input = argv[2])
}

r = createReport(PATH_TO_TXT, YAML_CONFIG)


