#install.packages('Biobase')
#source("http://bioconductor.org/biocLite.R")
#biocLite("Biobase")


#' Removes the longest common prefix (LCP) from a vector of strings.
#' 
#' If only a single string is given, no shortening will be done.
#' 
#' @param x Vector of strings with common prefix
#' @return Shortened vector of strings
#' 
#' @importFrom Biobase lcPrefix
#' 
#' @export
delLCP <- function(x)
{
  x = as.character(x)
  #require(Biobase)
  if (length(x) <= 1) return(x)
  
  lcp = nchar(lcPrefix( x ))   # shorten string (remove common prefix)
  x = sapply(x, substring, lcp+1)
  return (x)
}

#' Removes the longest common suffix (LCS) from a vector of strings.
#' 
#' @param x Vector of strings with common suffix
#' @return Shortened vector of strings
#' 
#' @importFrom Biobase lcSuffix
#' 
#' @export
delLCS <- function(x)
{
  #require(Biobase)
  lcs = nchar(lcSuffix( x ))   # shorten string (remove common suffix)
  x = sapply(x, function(v) substr(v, 1, nchar(v)-lcs))
  return (x)
}

#' Count the number of chars of the longest common prefix
#' 
#' @param x Vector of strings with common prefix
#' @return Length of LCP
#' 
#' @importFrom Biobase lcPrefix
#' 
#' @export
lcpCount <- function(x)
{
  #require(Biobase)
  lcp = nchar(lcPrefix( x ))   # number of common prefix chars
  return (lcp)
}

#' Count the number of chars of the longest common suffix
#' 
#' @param x Vector of strings with common suffix
#' @return Length of LCS
#' 
#' @importFrom Biobase lcSuffix
#' 
#' @export
lcsCount <- function(x)
{
  #require(Biobase)
  lcs = nchar(lcSuffix( x ))   # number of common suffix chars
  return(lcs)
}

#' Compute longest common substring of two strings.
#' 
#' Implementation is very inefficient (dynamic programming in R)
#' --> use only on small instances
#' 
#' @param s1 String one
#' @param s2 String two
#' @return String containing the longest common substring
#' 
#' @export
LCS <- function(s1, s2) 
{
  v1 <- unlist(strsplit(s1,split=""))
  v2 <- unlist(strsplit(s2,split=""))
  
  num <- matrix(0,nchar(s1), nchar(s2))    
  maxlen <- 0
  pstart = 0
  
  for (i in 1:nchar(s1)) {
    for (j in 1:nchar(s2)) {
      if (v1[i] == v2[j]) {
        if ((i==1) || (j==1)) { 
          num[i,j] <- 1
        } 
        else {
          num[i,j] <- 1+num[i-1,j-1]
        }
        if (num[i,j] > maxlen) {
          maxlen <- num[i,j]
          pstart = i-maxlen+1
        }
      }
    }
  }
  
  ## return the substring found
  return (substr(s1, pstart, pstart+maxlen-1))
}

#' Find longest common substring from 'n' strings.
#' 
#' Warning: heuristic! This is not guaranteed to find the best solution, since its done pairwise with the shortest input string as reference.
#' 
#' @param strings A vector of strings in which to search for LCS
#' @param min_LCS_length Minimum length expected. Search is aborted if result is shorter
#' @return longest common substring (or "" if shorter than \code{min_LCS_length})
#' 
#' @examples
#' LCSn(c("1_abcde...", "2_abcd...", "x_abc..."))  ## result: "_abc"
#' 
#' @export
#' 
LCSn = function(strings, min_LCS_length=7)
{
  ## abort if there is no chance of finding a suitably long substring
  if (min(nchar(strings)) < min_LCS_length) return("");
  
  if (length(strings) <= 1) return (strings)
  
  ## apply LCS to all strings, using the shortest as reference
  idx_ref = which(nchar(strings)==min(nchar(strings)))[1]
  strings_other = strings[-idx_ref]  
  r = unique(sapply(strings_other, LCS, strings[idx_ref]))
  r
  ## if only one string remains, we're done
  if (length(r) == 1) return (r[1])
  ## if its more, call recursively until a solution is found
  return (LCSn(r, min_LCS_length))
}
# LCSn(c("16_IMU008_CISPLA_E5_R11", "48_IMU008_CISPLA_P4_E7_R31", "60_IMU008_CISPLA_E7_R11"), 3)



#' Removes common substrings in a set of strings.
#' 
#' Usually handy for plots, where condition names should be as concise as possible.
#' E.g. you do not want names like 
#' 'TK20130501_H2M1_010_IMU008_CISPLA_E3_R1.raw' and 
#' 'TK20130501_H2M1_026_IMU008_CISPLA_E7_R2.raw' 
#' but rather 'TK.._010_I.._E3_R1.raw' and
#'            'TK.._026_I.._E7_R2.raw'
#'            
#' If multiple such substrings exist, the algorithm will remove the longest first and iterate
#' a number of times (two by default) to find the second/third etc longest common substring.
#' Each substring must fulfill a minimum length requirement - if its shorter, its not considered worth removing
#' and the iteration is aborted.
#' 
#' @param strings          A vector of strings which are to be shortened
#' @param infix_iterations Number of successive rounds of substring removal
#' @param min_LCS_length   Minimum length of the longest common substring (default:7, minimum: 6)
#' @return A list of shortened strings, with the same length as the input                       
#'
#' @examples
#' require(PTXQC)
#' simplifyNames(c('TK20130501_H2M1_010_IMU008_CISPLA_E3_R1.raw', 'TK20130501_H2M1_026_IMU008_CISPLA_E7_R2.raw'), infix_iterations = 2)
#' # --> "TK.._010_I.._E3_R1.raw","TK.._026_I.._E7_R2.raw"
#' 
#' try(simplifyNames(c("bla", "foo"), min_LCS_length=5))
#' # --> error, since min_LCS_length must be >=6
#' 
#' @export
#' 
simplifyNames = function(strings, infix_iterations=2, min_LCS_length=7)
{
  if (min_LCS_length<6) stop( "simplifyNames(): param 'min_LCS_length' must be 6 at least.")
  
  for (it in 1:infix_iterations)
  {
    lcs = LCSn(strings, min_LCS_length=min_LCS_length)
    if (nchar(lcs)==0) return (strings)
    ## replace infix with '..'
    strings = sub(paste0("(.*)", lcs, "(.*)"), paste0("\\1", substring(lcs,1,2), "..", substring(lcs, nchar(lcs)-1), "\\2"), strings)
  }
  return (strings)
}

#' Shorten a string to a maximum length and indicate shorting by appending '..'
#' 
#' Some axis labels are sometimes just too long and printing them will either
#' squeeze the actual plot (ggplot) or make the labels disappear beyond the margins (graphics::plot)
#' One ad-hoc way of avoiding this is to shorten the names, hoping they are still meaningful to the viewer.
#' 
#' @param x                Vector of input strings
#' @param max_len          Maximum length allowed
#' @param verbose          Print which strings were shortened
#' @param allow_duplicates If shortened strings are not discernible any longer, consider the short version valid (not the default), otherwise (default) return the full string (--> no-op)
#' @return  A vector of shortened strings
#' 
#' @examples
#' r = shortenStrings(c("gamg_101", "gamg_101230100451", "jurkat_06_100731121305", "jurkat_06_1"))
#' all(r == c("gamg_101", "gamg_101230100..", "jurkat_06_1007..", "jurkat_06_1"))
#' 
#' @export
#' 
shortenStrings = function(x, max_len = 20, verbose = TRUE, allow_duplicates = FALSE)
{
  idx = nchar(x) > max_len
  if (any(duplicated(x))) stop("Duplicated input given to 'shortenStrings'!\n  ", paste(x, collapse="\n  "))
  
  xr = x
  xr[idx] = paste0(substr(x[idx], 1, max_len-2), "..")
  
  if (!allow_duplicates & any(duplicated(xr)))
  {
    ## try cutting the prefix instead of the suffix...
    xr[idx] = paste0("..", sapply(x[idx], function(s) substr(s, nchar(s)-max_len+2, nchar(s))))
    if (any(duplicated(xr)))
    {
      stop("Duplicated output produced by 'shortenStrings'!\n  ", paste(xr, collapse="\n  "), "\nPlease rename your items (make them shorter and less similar)")
    }
  }
  
  if (verbose & any(idx))
  {
    cat("The following labels will be shortened to ease plotting:\n")
    cat(paste0("  ", paste(x[idx], collapse="\n  ")))
    cat("\n")
  }

  return (xr)
}

#' Compute shortest prefix length which makes all strings in a vector uniquely identifyable.
#'
#' If there is no unique prefix (e.g. if a string is contained twice), then the length
#' of the longest string is returned, i.e. if the return value is used in a call to substr, nothing happens
#' e.g.  substr(x, 1, supCount(x)) == x
#' 
#' @param x        Vector of strings
#' @param prefix_l Starting prefix length, which is incremented in steps of 1 until all prefixes are unique (or maximum string length is reached)
#' @return Integer with minimal prefix length required
#' 
#' @examples
#' supCount(c("abcde...", "abcd...", "abc..."))  ## 5
#'
#' x= c("doubled", "doubled", "aLongDummyString")
#' all( substr(x, 1, supCount(x)) == x )   ## TRUE (since due to duplicated entries there is no prefix which makes them unique)
#' 
#' @export
#' 
supCount <- function(x, prefix_l=1)
{
  max_length = max(nchar(x))
  while (prefix_l < max_length)
  {
    dups = duplicated(sapply(x, substr, 1, prefix_l), NA)
    if (sum(dups)==0)
    {
      break;
    }
    prefix_l = prefix_l + 1
  }
  return (prefix_l)
}

#' Reverse a string.
#' 
#' @param x String to be reversed
#' @return Reversed string
#' 
#' @export
#' 
strRev <- function(x) 
{
  return (sapply(lapply(strsplit(x, NULL), rev), paste, collapse=""))
}


#' Re-estimate a new set size to split a number of items into equally sized sets.
#'
#' This is useful for plotting large datasets where multiple pages are needed.
#' E.g. you know that you need 101 barplots, but you only want to fit about 25 per page.
#' Naively one would now do five plots, with the last one only containing a single barplot.
#' Using this function with correctSetSize(101, 25) would tell you to use 26 barplots per page,
#' so you end up with four plots, all roughly equally filled.
#' It also works the other extreme case, where your initial size is chosen slightly too high, e.g.
#' Sets of size 5 for just 8 items is too much, because we can reduce the set size to 4 and still
#' need two sets but now they are much more equally filled (correctSetSize(8, 5) == 4).
#' 
#' We allow for up to set sizes of 150\% from default, to avoid the last set being sparse (we remove it and distribute to the other bins)
#' 150%\ oversize is the extreme case, which only happens with sets of size two. With more sets the overhead is much smaller (1/X).
#' Once the number of sets is fixed, we distribute all items equally.
#' 
#' E.g. 6 items & initial_set_size=5, would result in 2 bins (5 items, 1 item), but we'd rather have one bin of 6 items
#' or 8 items & initial_set_size=5, would result in 2 bins (5+3 items), since the last set is more than half full, but we'd rather have 4+4
#'
#' @param item_count Known number of items which need to assigned to sets
#' @param initial_set_size Desired number of items a single set should hold
#' @return re-estimated set size which a set should hold in order to avoid underfilled sets
#' 
#' @examples
#'  stopifnot(
#'    correctSetSize(8, 5) == 4
#'  )
#'  stopifnot(
#'    correctSetSize(101, 25) == 26
#'  )
#'  
#' @export
#'  
correctSetSize = function(item_count, initial_set_size)
{
  blocks = seq(from = 1, to = item_count, by = initial_set_size)
  blockcount = length(blocks)
  lastblocksize = item_count - tail(blocks, n = 1) + 1
  #cat(paste("Naively, last set has size", lastblocksize, "\n"))
  if (lastblocksize < initial_set_size * 0.5 & blockcount > 1)
  { ## last block is not full, reduce block count if more than one block
    blockcount = blockcount - 1
    #cat(paste("reducing number of sets to", blockcount, "\n"))
  }
  ## distribute equally among fixed number of sets
  set_size = ceiling(item_count / blockcount)
  return (set_size)
}


#' Calls FUN on a subset of data in blocks of size 'subset_size' of unique indices.
#' 
#' One subset consists of 'subset_size' unique groups and thus of all rows which
#' in 'data' which have any of these groups.
#' The last subset might have less groups, if the number of unique groups is not dividable by subset_size.
#' 
#' FUN is applied on each subset.
#' 
#' @param data         Data.frame whose subsets to use on FUN
#' @param indices      Vector of group assignments, same length as nrow(data)
#' @param subset_size  Number of groups to use in one subset
#' @param FUN          Function applied to subsets of data
#' @param sort_indices Sort groups (by their sorted character(!) names) before building subsets
#' @param ...          More arguments to FUN
#'
#' @return list of function result (one entry for each subset)
#' 
#' @examples
#'  byX(data.frame(d=1:10), 1:10, 2, sum)
#'  
#' @export  
#' 
byX <- function(data, indices, subset_size = 5, FUN, sort_indices = TRUE, ...)
{
  stopifnot(subset_size > 0)
  
  groups = unique(indices)
  if (sort_indices)
  {
    #cat(paste0("Sorting indices (", length(groups), ")...\n"))
    groups = factor(sort(as.character(groups)))
    #cat(paste0(groups))
    #cat()
  }
  blocks = seq(from = 1, to = length(groups), by = subset_size)
  result = lapply(blocks, function(x) {
    #cat(paste("block", x, " ... "))
    range = x:(min(x+subset_size-1, length(groups)))
    lns = indices %in% groups[range];
    subset = droplevels(data[lns, , drop=F])
    rownames(subset) = rownames(data)[lns]
    #cat(paste("call\n"))
    r = FUN(subset, ...)    
    flush.console()
    return (r)
  })
  return (result)
}


#' Same as 'byX', but with more flexible group size, to avoid that the last group has only a few entries (<50\% of desired size).
#' 
#' The 'subset_size' param is internally optimized using \code{\link{correctSetSize}} and
#' then \code{\link{byX}} is called.
#' 
#' @param data         Data.frame whose subset to use on FUN
#' @param indices      Vector of group assignments, same length as nrow(data)
#' @param subset_size  Ideal number of groups to use in one subset -- this can be changed internally, from 75\%-150\%
#' @param FUN function Applied to subsets of data
#' @param sort_indices Groups are formed by their sorted character(!) names
#' @param ...          More arguments to FUN
#'
#' @return list of function result (one entry for each subset)
#' 
#' @examples 
#'  stopifnot(
#'    byXflex(data.frame(d=1:10), 1:10, 2, sum, sort_indices = FALSE) ==
#'    c(3, 7, 11, 15, 19)
#'  )
#'  
#' @export
#'  
byXflex <- function(data, indices, subset_size = 5, FUN, sort_indices = TRUE, ...)
{
  stopifnot(subset_size>0)
  #indices=1:24
  #subset_size=5
  groups = unique(indices)
  subset_size = correctSetSize(length(groups), subset_size)
  return (byX(data, indices, subset_size, FUN, sort_indices, ...))
}

#' Assign set numbers to a vector of values.
#'
#' Each set has size set_size (internally optimized using \code{\link{correctSetSize}}), holding values from 'values'.
#' This gives n such sets and the return value is just the set index for each value.
#' 
#' @param values       Vector of values
#' @param set_size     Number of distinct values allowed in a set
#' @param sort_values  Before assigning values to sets, sort the values?
#'
#' @return Vector (same length as input) with set numbers
#' 
#' @examples
#'  library(PTXQC)
#'  assignBlocks(c(1:11, 1), set_size = 3, sort_values = FALSE)
#'  ## --> 1 1 1 2 2 2 3 3 3 4 4 1
#'  
#' @export
#'  
assignBlocks = function(values, set_size = 5, sort_values = TRUE)
{
  groups = unique(values)
  ## correct
  set_size = correctSetSize(length(groups), set_size)
  ## sort
  if (sort_values)
  {
    #cat(paste0("Sorting values (", length(groups), ")..."))
    groups = factor(sort(as.character(groups)))
  }
  
  blocks = seq(from=1, to=length(groups), by=set_size)
  result = rep(NA, length(values))
  for (x in 1:length(blocks))
  {
    range = blocks[x]:(min(blocks[x]+set_size-1, length(groups)))
    lns = values %in% groups[range];
    result[lns] = x
  }
  #cat("done\n")
  return (result)
}




#' Grep with values returned instead of indices
#' 
#' @param reg  regex param
#' @param data container 
#' @param ... other params forwarded to grep()
#'
#' @return values of data which matched the regex
#' 
#' @examples 
#'   grepv("x", c("abc", "xyz"))
#'   ## --> "xyz"
#'   
#' @export
#'    
grepv = function(reg, data, ...)
{
  idx = grep(reg, data, ...)
  return (data[idx])  
}

#' paste with tab as separator
#' 
#' @param ... Arguments forwarded to paste()
#' @return return value of paste()
#' 
#' @examples
#'   pastet("tab","separated")
#'   ## --> "tab\tseparated"
#'   
#' @export
#'    
pastet = function(...)
{
  paste(..., sep="\t")
}

#' paste with newline as separator
#' 
#' @param ... Arguments forwarded to paste()
#' @return return value of paste()
#' 
#' @examples
#'   pasten("newline","separated")
#'   ## --> "newline\nseparated"
#'   
#' @export
#'    
pasten = function(...)
{
  paste(..., sep="\n")
}

#' Coefficient of variation (CV)
#' 
#' Computes sd(x) / mean(x)
#' 
#' @param x Vector of numeric values
#' @return CV
CV = function(x){sd(x) / mean(x)} 

#' Relative standard deviation (RSD)
#' 
#' Simply \code{\link{CV}}*100
#' 
#' @param x Vector of numeric values
#' @return RSD
#' 
RSD = function(x){CV(x)*100}

#' Replace 0 with NA in a vector
#' 
#' @param x A numeric vector
#' @return Vector of same size as 'x', with 0's replaced by NA
#' 
del0 = function(x)
{
  x[x==0] = NA
  return(x)
}


#' Find the local maxima in a vector of numbers.
#' 
#' A vector of booleans is returned with the same length as input (omitting NA's)
#' which contains TRUE when there is a maximum.
#' Simply sum up the vector to get the number of maxima.
#' 
#' @param x           Vector of numbers
#' @param thresh_rel  Minimum relative intensity to maximum intensity of 'x' required
#'                    to be a maximum (i.e., a noise threshold). Default is 20\%.
#' @return Vector of bool's, where TRUE indicates a local maximum.
#' 
#' @examples
#'     r = getMaxima(c(1,0,3,4,5,0))                                
#'     all(r == c(1,0,0,0,1,0))
#'     
#' @export     
#'     
getMaxima = function(x, thresh_rel = 0.2)
{
  pos = rep(FALSE, length(x))
  x = na.omit(x)
  thresh_abs = max(x) * thresh_rel
  last = x[1]
  up = TRUE ## we start going up (if the next point is lower, the first point is a maximum)
  for (i in 2:length(x))
  {
    if (last < x[i] && !up) 
    { # ascending in down mode
      up = !up
    } else if (last > x[i] && up)
    { # going down in up mode
      up = !up
      if (x[i-1]>=thresh_abs) pos[i-1] = TRUE
    }
    last = x[i]
  }
  if (up) pos[length(x)] = TRUE ## if we ended going up, the last point is a maximum
  return (pos)
}

#' Mosaic plot of two columns in long format.
#' 
#' Found at http://stackoverflow.com/questions/19233365/how-to-create-a-marimekko-mosaic-plot-in-ggplot2
#' Modified (e.g. to pass R check)
#' 
#' Returns a ggplot object without printing it.
#' This allows to further modify the plot (e.g. rename axis, color schemes etc)
#' 
#' @param var1 Vector of factors
#' @param var2 Vector of factors
#' @return ggplot object
#' 
#' 
#' @export
#' 
mosaicPlot = function(var1, var2)
{
  lev_var1 = length(levels(var1))
  
  pl.data = as.data.frame(prop.table(table(var1, var2)))
  pl.data$margin_var1 = prop.table(table(var1))
  pl.data$var2_height = pl.data$Freq / pl.data$margin_var1
  pl.data$var1_center = c(0, cumsum(pl.data$margin_var1)[(1:lev_var1) -1]) + pl.data$margin_var1 / 2
  
  pl = 
    ggplot(pl.data, aes_string(x = "var1_center", y = "var2_height")) +
      geom_bar(stat = "identity", aes_string(width = "margin_var1", fill = "var2"), color = "black")  +
      geom_text(aes_string(label = "as.character(var1)", x = "var1_center", y = 1.05)) 
  return (pl)
}

#' A string concatenation function, more readable than 'paste()'.
#' 
#' @param a Char vector
#' @param b Char vector
#' @return Concatenated string (no separator)
#'
#' @export
#'
`%+%` <- function(a, b)
{
  return (paste(a, b, sep=""))
} 

#'
#' Given a vector of (short/long) filenames, translate to the (long/short) version
#' 
#' @param f_names Vector of filenames
#' @param mapping A data.frame with from,to columns
#' @return A vector of translated file names as factor (ordered by mapping!)
#'
#' @export
#'
renameFile = function(f_names, mapping)
{
  if (all(f_names %in% mapping$from)) {      ## from -> to
    f_new = mapping$to[match(f_names, mapping$from)]
    f_new = factor(f_new, levels=mapping$to)
  } else if (all(f_names %in% mapping$to)) { ## to -> from
    f_new = mapping$from[match(f_names, mapping$to)]
    f_new = factor(f_new, levels=mapping$from)
  } else {
    stop("Error in renameFile: filenames cannot be found in mapping!")
  }
  
  return (f_new)
} 

#'
#' Add the value of a variable to an environment (fast append)
#'
#' The environment must exist, and its name must be given as string literal in 'env_name'!
#' The value of the variable 'v' will be stored under the name given in 'v_name'.
#' If 'v_name' is not given, a variable name will be created by increasing an internal counter
#' and using the its value padded with zeros as name (i.e., "0001", "0002" etc).
#' 
#' @param env_name String of the environment variable
#' @param v Value to be inserted
#' @param v_name String used as variable name. Automatically generated if omitted.
#' @return Always TRUE
#'
#'
appendEnv = function(env_name, v, v_name = NULL)
{
  e = eval.parent(as.name(env_name), n=3)
  e$.counter <- e$.counter + 1
  
  e[[sprintf("%04d",e$.counter)]] <- v ## pad with 0's, to ensure correct order when calling ls() on env_name
  return(TRUE)
}

#'
#' Thin out a data.frame by removing rows with similar values in a certain column.
#' 
#' All values in the column 'filterColname' are assigned to bins if width 'binsize'.
#' Only one value per bin is retained. All other rows are removed and the reduced
#' data frame will all its columns is returned.
#' 
#' @param data The data.frame to be filtered
#' @param filterColname Name of the filter column as string
#' @param binsize Width of a bin
#'
#'
thinOut = function(data, filterColname, binsize)
{
  nstart = nrow(data)
  ## first remove duplicates (they cannot possibly pass the filter)
  data = data[!duplicated(data[, filterColname]), ]
  ##
  data = data[order(data[, filterColname]), ]
  local_deltas = c(0, diff(data[, filterColname]))
  ld_cs = cumsum(local_deltas)
  data = data[!duplicated(round(ld_cs / binsize)),]
  #print("Saved " %+% round(100 - nrow(x)/nstart*100) %+% "% data")
  return(data)
}

#'
#' Apply 'thinOut' on all subsets of a data.frame, split by a batch column
#' 
#' The binsize is computed from the global data range of the filter column by dividing the range
#' into binCount bins.
#' 
#' @param data The data.frame to be split and filtered(thinned)
#' @param filterColname Name of the filter column as string
#' @param batchColname Name of the split column as string
#' @param binCount Number of bins in the 'filterColname' dimension.
#'
#' @importFrom plyr ddply
thinOutBatch = function(data, filterColname, batchColname, binCount = 1000)
{
  binsize = (max(data[, filterColname], na.rm=T) - min(data[, filterColname], na.rm=T)) / binCount
  r = ddply(data, batchColname, thinOut, filterColname, binsize)
  return (r)  
}


#'
#' Assign a relative abundance class to a set of (log10) abundance values
#' 
#' Abundances (should be logged already) are grouped into different levels,
#' starting from the smallest values ("low") to the highest values ("high").
#' Intermediate abundances are either assigned as "mid", or "low-mid".
#' If the range is too large, only "low" and "high" are assigned, the intermediate values
#' are just numbers.
#' 
#' Example:
#'  getAbundanceClass(c(12.4, 17.1, 14.9, 12.3)) ## --> factor(c("low", "high", "mid", "low"))
#' 
#' 
#' @param x Vector of numeric values (in log10)
#' @return Vector of factors corresponding to input with abundance class names (e.g. low, high)
#'
getAbundanceClass = function(x) {
  r = data.frame(x = x)
  r$x_diff = round(r$x - min(r$x))
  r$x_diff_fac = as.numeric(factor(r$x_diff)) ## make a factor, to get ordering (lowAbd=1, ...)
  ## assign names to abundance classes
  cc = length(unique(r$x_diff))
  if (cc==1) lvls = "mid"
  if (cc==2) lvls = c("low", "high")
  if (cc==3) lvls = c("low", "mid", "high")
  if (cc==4) lvls = c("low", "low-mid", "mid-high", "high")
  if (cc>4) {
    lvls = 1:max(r$x_diff_fac)
    lvls[1] = "low"
    lvls[length(lvls)] = "high"
  }
  r$x_cl = lvls[r$x_diff_fac]
  r$x_cl = factor(r$x_cl, levels=lvls, ordered=T)
  r
  return(r$x_cl)
}
