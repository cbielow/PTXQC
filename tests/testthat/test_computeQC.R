library(PTXQC)

## the next require() is needed to prevent a spurious error in certain R versions (might be a bug in R or a package)
## error message is:
##    Error in Scales$new : could not find function "loadMethod"
require(methods)

context("fcn_computeQC.R")

test_that("createReport", {
  ## this is a rather lengthy function, and its hard to test in all its granlarity (hence we test
  ## the functions within in other tests), but at least we want to see that it produces
  ## results -- even though its hard to test if the resulting PDF report is exactly what we expect
  ## without going into mindboggling image-comparison tests...
  
  
  ## get some data
  local_zip = tempfile(fileext=".zip")
  ## getting the correct URL is tricky, since download methods like 'curl' will not follow
  ## redirects, but download the server message instead. Try 'curl' etc on the command line, to see
  ## what's going on and find the right URL (for GitHub: we need https and raw.github....)
  target_url = "https://raw.githubusercontent.com/cbielow/PTXQC_data/master/txt_Ecoli.zip"
  tryCatch({
    dl = download.file(target_url, destfile = local_zip)
  }, error = function(err) {
    ## in case of error, try with Curl
    dl = download.file(target_url, destfile = local_zip, method='curl') ## for Linux/MacOSX
  })

  unzip(local_zip, exdir = tempdir()) ## extracts content
  txt_folder = file.path(tempdir(), "txt")
  yaml_obj = list() ## so special config...
  
  r = createReport(txt_folder, yaml_obj)
  expect_equal(c("yaml_file", "heatmap_values_file", "R_plots_file", "filename_sorting", "stats_file",         
                 "report_file_simple", "report_file_extended", "report_file", "report_file_extension"), names(r))
  rep_file = paste0(r[["report_file"]], r[["report_file_extension"]])
  
  expect_equal(file.exists(rep_file), TRUE)
  expect_gt(file.info(rep_file)$size, 100*1024) ## ~119kb PDF
  
  expect_equal(file.exists(r[["heatmap_values_file"]]), TRUE)
  d_heatmap = read.delim(r[["heatmap_values_file"]])
  expect_equal(dim(d_heatmap), c(2, 22)) ## two files, 22 metrics
  expect_equal(as.character(d_heatmap$fc.raw.file), c("..Ecoli_01", "..Ecoli_02"))
  
  expect_equal(file.exists(r[["filename_sorting"]]), TRUE)
  d_filenamesort = read.delim(r[["filename_sorting"]], comment.char="#")
  expect_equal(dim(d_filenamesort), c(2, 3)) ## two files, three columns
  expect_equal(as.character(d_filenamesort$new.Name), c("..Ecoli_01", "..Ecoli_02"))
  
  
  unlink(local_zip) ## delete zip
  unlink(txt_folder, recursive = TRUE) ## delete txt-folder
  
})