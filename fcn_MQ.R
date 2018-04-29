###
### Author: Chris Bielow
###
###

#' Boxplots - one for each condition (=column) in a data frame.
#' 
#' Given a data.frame with two/three columns in long format (name, value, [contaminant]; in that order), each group (given from 1st column)
#' is plotted as a bar.
#' Contaminants (if given) are separated and plotted as yellow bars.
#' 
#' Boxes are shaded: many NA or Inf lead to more transparency. Allows to easily spot sparse groups
#' 
#' 
#' @param data    Data frame in long format with numerical expression data
#' @param log2    Apply log2 to the data (yes/no)
#' @param ylab    Label on Y-axis
#' @param mainlab Main title
#' @param sublab  Sub title
#' @param boxes_per_page  Maximum number of boxplots per plot. Yields multiple plots if more groups are given.
#' @param abline          Draw a horizontal green line at the specified y-position (e.g. to indicate target median values)
#' @param coord_flip      Exchange Y and X-axis for better readability
#' @param names           An optional data.frame(long=.., short=..), giving a renaming scheme (long->short) for the 'name' column
#' @return List of ggplot objects
#' 
#' @import ggplot2
#' @importFrom plyr ddply
#' @importFrom grDevices boxplot.stats
#' 
#' @export
#' 
boxplotCompare = function(data, 
                          log2 = TRUE,
                          ylab = "intensity",
                          mainlab = ylab,
                          sublab = "",
                          boxes_per_page = 30,
                          abline = NA,
                          coord_flip = TRUE,
                          names = NA)
{
  
  if (ncol(data) == 2) {
    data$contaminant = FALSE ## add a third column, if missing
  }
  colnames(data) = c("group", "value", "contaminant")
  
  if (log2) {
    data$value = log2(data$value)
  }
  
  ## shorten group names
  if (class(names)=='data.frame')
  {
    stopifnot(sort(colnames(names)) == c("long", "short"))
    data$group = names$short[match(data$group, names$long)]
    if (any(is.na(data$group))) 
    {
      print(names)
      stop("Group renaming is incomplete! Aborting...")
    }
  }
  ## make it a factor
  if (!("factor" %in% class(data$group))) data$group = factor(data$group)
  
  ## actual number of entries in each column (e.g. LFQ often has 0)
  ncol.stat = ddply(data, "group", function(x){
    notNA = sum(!is.infinite(x$value) & !is.na(x$value));
    data.frame(n = nrow(x), notNA = notNA, newname = paste0(x$group[1], " (n=", notNA, ")"))})
  head(ncol.stat)
  ## rename (augment with '#n')
  data$group2 = ncol.stat$newname[match(data$group, ncol.stat$group)]

  ## check (ddply makes  ncol.stat$newname a factor with matching levels)
  stopifnot(class(data$group2) == "factor")
  #stopifnot(all(as.numeric(data$group) == as.numeric(data$group2))) ## as.numeric() is meaningless on factors with different levels
  data$group = data$group2
  
  ## remove -inf and NA's
  data = data[!is.infinite(data$value) & !is.na(data$value), ]
  
  groups = unique(data$group);
  ## add color for H vs L (if SILAC)
  cols = c("sample" = "black", 
           "sample (light)" = "black", 
           "sample (medium)" = "blue", 
           "sample (heavy)" = "green",
           "contaminant" = "yellow")
  cat_names = names(cols)
  cat = factor(cat_names, levels=cat_names)
  data$cat = cat[1]
  if (sum(grepl("^[^HLM]", groups )) == 0 || sum(grepl("^intensity\\.[hlm]\\.", groups )) > 0)
  { ## all start with either L, M or H
    data$cat = cat[2]
    data$cat[grepl("^M", data$group) | grepl("^intensity\\.m\\.", data$group)] = cat[3]
    data$cat[grepl("^H", data$group) | grepl("^intensity\\.h\\.", data$group)] = cat[4]
  }
  data$cat[data$contaminant] = cat[5]
  
  ## compute global y-limits (so we can fix it across plots)
  ylims = boxplot.stats(data$value)$stats[c(1, 5)]
  ## make sure to inlude abline (if existing)
  if (!is.na(abline))
  {
    ylims = c(ylims[1], max(ylims[2], abline))
  }
  
  fcn_boxplot_internal = function(data, abline = NA) 
  {
    #require(ggplot2)
    pl = ggplot(data=data, aes_string(x = "group", y = "value", fill = "cat")) + ## do not use col="cat", since this will dodge bars and loose scaling
      geom_boxplot(varwidth = TRUE) +
      xlab("") + 
      ylab(ylab) +
      ylim(ylims) +
      scale_alpha(guide = FALSE) +
      scale_fill_manual(values=cols, name = "Category") + 
      scale_color_manual(values=cols, name = "Category") + 
      theme(axis.text.x = element_text(angle=90, vjust = 0.5)) +
      theme(legend.position=ifelse(length(cols)==1, "none", "right")) +
      addGGtitle(mainlab, sublab) + 
      scale_x_discrete_reverse(unique(data$group))
    
    if (!is.na(abline))
    {
      pl = pl + geom_abline(alpha = 0.5, intercept = abline, slope = 0, colour = "green")
    }
    if (coord_flip == TRUE)
    {
      pl = pl + coord_flip()
    }
    #print(pl)
    return(pl)
  }
  lpl = byXflex(data = data, indices = data$group, subset_size = boxes_per_page, sort_indices = TRUE, FUN = fcn_boxplot_internal, abline)
  return (lpl)
}


#'
#' Extract fragment mass deviation errors from a data.frame from msms.txt
#' 
#' Given a data.frame as obtainable from a msms.txt with 
#'  - a 'mass.analyzer' column which contains only a single value for the whole column
#'  - a 'mass.deviations..da.' and (if available) 'mass.deviations..ppm.'
#'  - a 'masses' column (only required if 'mass.deviations..ppm.' is unavailable and the mass.analyzer
#'    indicates hig-res data)
#' 
#' Mass deviations are extracted from the columns, e.g. each cell containing values separated by
#' semicolons is split into single values. The appropriate unit is chosen (Da or ppm, depending on
#' ITMS or FTMS data). Also the fragmentation type can be used: CID indicates ITMS, HCD to FTMS.
#' This is not 100% safe, but older MQ versions do not report the mass analyzer properly.
#' 
#' Sometimes, peptides are identified purely based on MS1, i.e. have no fragments. These will be ignored.
#' 
#' If ppm mass deviations are not available, errors in Da will be converted to ppm using the corresponding mass values.
#' 
#' @param x Data frame in long format with numerical expression data
#' @return  Data frame with mass errors ('msErr') and their 'unit' (Da or ppm) or NULL (if no fragments were given)
#' 
#' @export
#' 
getFragmentErrors = function(x)
{
  ## require only one mass analyzer type:
  stopifnot(length(unique(x$mass.analyzer))==1)
  stopifnot(all(c("mass.analyzer", "mass.deviations..da.") %in% colnames(x)))
  
  convert_Da2PPM = FALSE
  if (grepl("ITMS|TOF|CID", x$mass.analyzer[1]) & ("mass.deviations..da." %in% colnames(x)))
  {
    ms2_unit = "[Da]"; ms2_col = "mass.deviations..da."
  } else if (grepl("FTMS|HCD", x$mass.analyzer[1]) & ("mass.deviations..ppm." %in% colnames(x))) {
    ms2_unit = "[ppm]"; ms2_col = "mass.deviations..ppm."
  } else if (grepl("FTMS|HCD", x$mass.analyzer[1]) & ("mass.deviations..da." %in% colnames(x))) {
    ## we know its high resolution, but this MQ version only gave us Dalton mass deviations
    ## --> convert back to ppm
    ms2_unit = "[ppm]"; ms2_col = "mass.deviations..da."
    convert_Da2PPM = TRUE
  } else {
    ## fallback if mass.analyzer is 'Unknown' (e.g. for mzXML input)
    if ("mass.deviations..ppm." %in% colnames(x)) {
      ms2_unit = "[ppm]"; ms2_col = "mass.deviations..ppm."
    } else {
      ms2_unit = "[Da]"; ms2_col = "mass.deviations..da."
    }
  }

  ## sometimes, peptides are identified purely based on MS1, i.e. have no fragments
  ## if only those peptides are present here in 'x', then 'err' below will be empty
  ## and no data.frame can be constructed, so...
  if (all(nchar(x[, ms2_col]) == 0)) return(NULL)
    
  err = as.numeric(unlist(strsplit(x[, ms2_col], split=";", fixed = TRUE)))
  if (convert_Da2PPM) {
    stopifnot("masses" %in% colnames(x))
    mass = unlist(strsplit(x$masses, split=";", fixed = TRUE))
    err = err / as.numeric(mass) * 1e6
  }
  
  if ((ms2_unit == "[ppm]") & (median(abs(err)) > 10)) {
    cat(paste0("MS/MS fragment error seems rather large ", median(abs(err)), ". Reporting in [Da]...\n"))
    # heuristic: ppm errors seem to be way to big. Use 'Da' instead.
    x$mass.analyzer = "ITMS"
    return (getFragmentErrors(x))
  }
  
  return(data.frame(msErr = err, unit = ms2_unit))
}

#'
#' Detect (and fix) MaxQuant mass recalibration columns, since they
#' sometimes report wrong values.
#' 
#' Returns a list of items for both diagnostics and possibly a fixed evidence data.frame.
#' Also two strings with messages are returned, which can serve as user message for
#' pre and post calibration status. 
#' 
#' @param df_evd Evidence data.frame with columns ()
#' @param df_idrate Data.frame from summary.txt, giving ID rates for each raw file (cols: "ms.ms.identified....", "fc.raw.file"). Can also be NULL.
#' @param tolerance_sd_PCoutOfCal Maximal standard deviation allowed before considered 'failed'
#' @param low_id_rate Minimum ID rate in Percent before a Raw file is considered 'failed'
#' @return list of data (stats, affected_raw_files, df_evd, recal_message, recal_message_post)
#'  
#'
fixCalibration = function(df_evd, df_idrate = NULL, tolerance_sd_PCoutOfCal = 2, low_id_rate = 1)
{
  
  stopifnot(c("fc.raw.file", "mass", "charge", "m.z", "mass.error..ppm.", "uncalibrated.mass.error..ppm.") %in% colnames(df_evd))
  
  ## heuristic to determine if the instrument is completely out of calibration, 
  ## i.e. all ID's are false positives, since the Precursor mass is wrong
  ## -- we use the SD; if larger than 2ppm, 
  ## then ID's are supposedly random
  ## -- alt: we use the 1%-to-99% quantile range: if > 10ppm
  ## -- uninformative for detection is the distribution (it's still Gaussian for a strange reason)
  MS1_decal_smr = ddply(df_evd, "fc.raw.file", function(x) 
    data.frame(n = nrow(x), 
               sd = round(sd(x$mass.error..ppm., na.rm = TRUE), 1), 
               range = diff(quantile(x$mass.error..ppm., c(0.01, 0.99), na.rm = TRUE)),
               decal = (median(abs(x$uncalibrated.mass.error..ppm.), na.rm = TRUE) > 1e3),
               hasMassErrorBug = FALSE,
               hasMassErrorBug_unfixable = FALSE)
  )
  ## additionally use MS2-ID rate (should be below 1%)
  if (is.null(df_idrate)) {
    MS1_decal_smr$ms.ms.identified.... = 0 ## upon no info: assume low IDrate
  } else {
    MS1_decal_smr = merge(MS1_decal_smr, df_idrate)
  }
  MS1_decal_smr$lowIDRate = MS1_decal_smr$ms.ms.identified.... < low_id_rate
  
  recal_message = ""
  recal_message_post = ""
  ## check each raw file individually (usually its just a few who are affected)
  if (any(MS1_decal_smr$decal, na.rm = TRUE))
  {
    recal_message = "MQ bug: data rescued"
    recal_message_post = 'MQ bug: data cannot be rescued'
    
    MS1_decal_smr$hasMassErrorBug[ MS1_decal_smr$fc.raw.file %in% MS1_decal_smr$fc.raw.file[MS1_decal_smr$decal > 0] ] = TRUE
    
    ## re-compute 'uncalibrated.mass.error..ppm.' and 'mass.error..ppm.'
    df_evd$theomz = df_evd$mass / df_evd$charge + 1.00726
    df_evd$mass.error..ppm.2 = (df_evd$theomz - df_evd$m.z) / df_evd$theomz * 1e6
    df_evd$uncalibrated.mass.error..ppm.2 = df_evd$mass.error..ppm.2 + df_evd$uncalibrated...calibrated.m.z..ppm.
    
    ## check if fix worked
    de_cal2 = ddply(df_evd, "fc.raw.file", .fun = function(x)  data.frame(q = (median(abs(x$uncalibrated.mass.error..ppm.2), na.rm = TRUE) > 1e3)))
    if (any(de_cal2$q, na.rm = TRUE))
    { ## fix did not work
      MS1_decal_smr$hasMassErrorBug_unfixable[ MS1_decal_smr$fc.raw.file %in% de_cal2$fc.raw.file[de_cal2$q] ] = TRUE
      recal_message = "m/z recalibration bugfix applied but failed\n(there are still large numbers)"
    }
    
    idx_overwrite = (df_evd$fc.raw.file %in% MS1_decal_smr$fc.raw.file[MS1_decal_smr$decal > 0])
    ## overwrite original values
    df_evd$mass.error..ppm.[idx_overwrite] = df_evd$mass.error..ppm.2[idx_overwrite]
    df_evd$uncalibrated.mass.error..ppm.[idx_overwrite] = df_evd$uncalibrated.mass.error..ppm.2[idx_overwrite]
  }
  
  MS1_decal_smr$outOfCal = (MS1_decal_smr$sd > tolerance_sd_PCoutOfCal) & 
    (MS1_decal_smr$lowIDRate) & 
    (MS1_decal_smr$sd < 100)  ## upper bound, to distinguish from MQ bug (which has much larger SD's)
  ## NA's are possible if only a single peptide was seen (yielding MS1_decal_smr$sd == NA)
  MS1_decal_smr$outOfCal[is.na(MS1_decal_smr$outOfCal)] = FALSE
  
  ## report too small search tolerance
  if (any(MS1_decal_smr$outOfCal)) recal_message = "search tolerance too small"
  
  
  ## Raw files where the mass error is obviously wrong (PSM's not substracted etc...)
  affected_raw_files = MS1_decal_smr$fc.raw.file[MS1_decal_smr$outOfCal | MS1_decal_smr$hasMassErrorBug]
  
  
  return(list(stats = MS1_decal_smr,
              affected_raw_files = affected_raw_files,
              df_evd = df_evd[, c("mass.error..ppm.", "uncalibrated.mass.error..ppm.", "fc.raw.file")],
              recal_message = recal_message,
              recal_message_post = recal_message_post
              ))
}
