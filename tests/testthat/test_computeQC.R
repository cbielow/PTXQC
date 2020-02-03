library(PTXQC)

## the next require() is needed to prevent a spurious error in certain R versions (might be a bug in R or a package)
## error message is:
##    Error in Scales$new : could not find function "loadMethod"
require(methods)

context("createReport.R")

test_that("createReport", {
  ## this is a rather lengthy function, and its hard to test in all its granularity (hence we test
  ## the functions within in other tests), but at least we want to see that it produces
  ## results -- even though its hard to test if the resulting PDF report is exactly what we expect
  ## without going into mindboggling image-comparison tests...
  
  
  ## get some data
  local_zip = tempfile(fileext=".zip")
  ## getting the correct URL is tricky, since download methods like 'curl' will not follow
  ## redirects, but download the server message instead. Try 'curl' etc on the command line, to see
  ## what's going on and find the right URL (for GitHub: we need https and raw.github....)
  target_url = "https://raw.githubusercontent.com/cbielow/PTXQC_data/master/txt_Ecoli.zip"
  dl = NULL
  tryCatch({
    dl = download.file(target_url, destfile = local_zip, quiet = TRUE)
  }, silent = TRUE, warning = function(w) { }, error = function(err) {
    ## in case of error, try with Curl
    tryCatch({
      dl = download.file(target_url, destfile = local_zip, method='curl', quiet = TRUE) ## for Linux/MacOSX
    }, silent = TRUE, warning = function(w) { }, error = function(err2) {
      print("Internet down. Aborting test gracefully")
      return(); ## fail gracefully
    })
  })
  if (is.null(dl)) return()
  
  unzip(local_zip, exdir = tempdir()) ## extracts content
  txt_folder = file.path(tempdir(), "txt")
  yaml_obj = list() ## no special config...
  
  r = createReport(txt_folder, NULL, yaml_obj)
  expect_equal(c("yaml_file", "heatmap_values_file", "R_plots_file", "filename_sorting", "stats_file",         
                 "log_file", "report_file_prefix", "report_file_PDF", "report_file_HTML"), names(r))
  rep_files = c(r[["report_file_PDF"]], r[["report_file_HTML"]])
  
  print(list.files(path = txt_folder))
  
  for (f in rep_files)
  {
    cat("Checking file ", f, "\n")
    # HTML file might not exist if PANDOC is not installed
    if (file.exists(f)) expect_equal(file.info(f)$size > 100*1024, TRUE) ## ~119kb PDF & HTML
  }
  
  expect_equal(file.exists(r[["heatmap_values_file"]]), TRUE)
  d_heatmap = read.delim(r[["heatmap_values_file"]])
  expect_equal(dim(d_heatmap)[1], 2) ## two files
  expect_equal(dim(d_heatmap)[2] >= 22, TRUE) ## 22 (or more) metrics
  expect_equal(as.character(d_heatmap$fc.raw.file), c("..Ecoli_01", "..Ecoli_02"))
  
  expect_equal(file.exists(r[["filename_sorting"]]), TRUE)
  d_filenamesort = read.delim(r[["filename_sorting"]], comment.char="#")
  expect_equal(dim(d_filenamesort), c(2, 3)) ## two files, three columns
  expect_equal(as.character(d_filenamesort$new.Name), c("..Ecoli_01", "..Ecoli_02"))
  
  
  unlink(local_zip) ## delete zip
  unlink(txt_folder, recursive = TRUE) ## delete txt-folder
  
})
