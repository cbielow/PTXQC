
#'
#' Quality metric with linear response to input, reaching the maximum score at the given threshold.
#' 
#' Ranges between 0 (worst score) and 1 (best score).
#' Useful for performance measures where reaching a certain reference threshold 't'
#' will be enough to reach 100\%.
#' The input range from [0, t] is scored from 0-100\%.
#' 
#' @param x Numeric value(s) between [0, inf]
#' @param t Threshold value, which indicates 100\%
#' @return Value between [0, 1]
#' 
#' 
qualLinThresh = function(x, t = 1)
{
  x[is.na(x)] = 0; ## replace NA with 0
  if (any(x < 0)) stop("qualLinThresh(): negative input 'x' not allowed!")
  if (t < 0) stop("qualLinThresh(): negative threshold 't' not allowed!")
  
  q = pmin(1, x / t)
  return (q)
}


#'
#' Quality metric for 'centeredness' of a distribution around zero.
#' 
#' Ranges between 0 (worst score) and 1 (best score).
#' A median of zero gives the best score of 1.
#' The closer the median is to the most extreme value of the distribution, the smaller the score (until reaching 0).
#' Can be used for calibrated mass errors, as a measure of how well they are centered around 0.
#' E.g. if the median is 0.1, while the range is [-0.5,0.5], the score will be 0.8 (punishing the 20% deviation).
#' If the range of data is asymmetric, e.g. [-1.5,-0.5] and does not include zero, the score cannot reach 1, 
#' since the median can never be zero.
#' 
#' @param x Numeric values (e.g. ppm errors)
#' @return Value between [0, 1]
#' 
#' 
qualCentered = function(x)
{
  q = 1 - (abs(median(x)) / max(abs(x)))
  return (q)
}

#'
#' Quality metric for 'centeredness' of a distribution around zero with a user-supplied range threshold.
#' 
#' Ranges between 0 (worst score) and 1 (best score).
#' The best score is achieved when the median of 'x' is close to the center of the interval (given by 'tol').
#' If median of 'x' is close to the border (on either side), the score decreases linearly to zero.
#' Can be used for uncalibrated mass errors, as a measure of how well they are centered around 0.
#' 
#' @param x  Vector of values (hopefully in interval [-tol, tol])
#' @param tol Border of interval (must be positive)
#' @return Value between [0, 1]
#' 
#' 
qualCenteredRef = function(x, tol)
{
  x = na.omit(x)
  if (length(x) == 0) stop("qualCenteredRef(): input is empty or consists of NA only!")
  if (tol <= 0) stop("qualCenteredRef(): negative or zero interval border not allowed!")
  m = median(x, na.rm=T)
  if (abs(m) > tol) warning("qualCenteredRef(): Median of x is outside of interval. Score will be set to 0.")
  q = 1 - (abs(m) / tol)
  q = max(0, q) ## avoid negative scores if abs(x)>tol
  return (q)
}


#'
#' Quality metric which measures the absolute distance from median.
#' 
#' Ranges between 0 (worst score) and 1 (best score).
#' Input must be between [0,1].
#' Deviations from the median of the sample represent the score for each sample point.
#' 
#' @param x A vector numeric values between [0,1]
#' @return A vector of the same size as x, with quality values between [0, 1]
#'
qualMedianDist = function(x)
{
  if (any(x < 0 | x > 1)) stop("qualMedianDist(): x values out of range [0,1]!")
  q = 1 - abs(x - median(x, na.rm=T))
  return (q)
}


#'
#' Compute deviation from uniform distribution
#' 
#' Ranges between 0 (worst score) and 1 (best score).
#' Input 'x' is a vector of counts (or probabilities) for equally spaced bins in a histogram.
#' A uniform distribution (e.g. c(3,3,3) will get a score of 1. The worst possible case (e.g. c(4,0,0)), will get a score of 0,
#' and a linear increasing function (e.g. c(1,2,3)) will get something in between (0.585 here)
#' 
#' In addition, bin values can be weighted (e.g. by their confidence). The total sum of weights is normalized to 1 internally.
#' 
#' The distance function used is the square root of the absolute difference between a uniform distribution and the input 'x'
#' (summed for each element of 'x').
#' This distance is normalized to the worst possible input (e.g. one bin with 100% counts, all other bins being empty).
#' 
#' @param x Vector of numeric intensity/count values (e.g. ID's per RT bin); bins are assumed to have equal widths
#' @param weight Vector of weights for values in 'x' (same length as 'x').
#' @return Value between [0, 1]
#' 
#' @export
#' 
#' @examples 
#'  stopifnot(qualUniform(c(3,3,3))==1)
#'  stopifnot(qualUniform(c(4,0,0))==0)         
#'
#'  ## how 'uniform' is a vector where only a single index has weight?-- answer: very
#'  stopifnot(qualUniform(c(4,0,0), c(1,0,0))==1)   
#'  stopifnot(qualUniform(c(4,0,0), c(0,1,0))==1)     
#'  stopifnot(qualUniform(c(0,4,0))==0)              
#'  stopifnot(abs(qualUniform(c(3,2,1))-0.58578) < 0.0001)
#'  stopifnot(abs(qualUniform(c(1,2,3))-0.58578) < 0.0001)
#'  stopifnot(qualUniform(c(1,2,3), c(0,1,0))==1)   
#'  stopifnot(abs(qualUniform(c(1,2,3))-0.58578) < 0.0001)
#'  stopifnot(abs(qualUniform(c(1,2,3), c(0,1,1))- 0.590316) < 0.0001)
#'  stopifnot(abs(qualUniform(c(2,3), c(1,1))-0.552786) < 0.0001)
#'  stopifnot(abs(qualUniform(1:120)-0.38661) < 0.0001)
#'  
qualUniform = function(x, weight=vector())
{
  if (any(is.na(x) | (x<0))) stop("Error: qualUniform() received negative or missing data!")
  
  if(length(weight)==0) {
    weight=rep(1, length(x)) 
  } else if (length(weight)!=length(x)) {
    stop("Error: qualUniform() received weights which do not match length of input!")
  }
  weight = weight/sum(weight)  ## normalize to 1
  x = x/sum(x)

  n = length(x)
  #y = 1/n ## expected intensity (no weights)
  y = sum(x*weight) ## expected weighted intensity
  y
  p = 1/2
  ## for worst possible score, assume bin with highest weight has all the counts
  ## i.e. error is 1-expected == 1-y for 100% bin, and abs(0-y)=y for all other bins of zero intensity
  idx_maxWeight = which(weight == max(weight))[1]
  worst = ((1-y)^p)*weight[idx_maxWeight] + (y^p)*sum(weight[-idx_maxWeight])
  ## worst ==0, iff weight-vector has only one non-zero entry, e.g. c(0,1,0)
  #worst = ((1-y)^p) + (n-1)*((1/n)^p)
  worst  
  ## score of current distribution
  sc = sum( (abs(x-y)^p) * weight)
  sc
  q = ifelse(worst==0, 1, (worst - sc) / worst)
  return (q)    
}


#'
#' Test for uniform distribution using Kolmogorov-Smirnoff
#'
#' If only 'x' is given, a one sample KS test is done (min and max are taken from 'x').
#' If both 'x' and 'y' are given, a two-sample KS test is conducted.
#' 
#' Thoughts: this test has multiple problems:
#'      - is much too stringent, leading to p-values  ~ 0 very quickly
#'      - looking at the 'D' statistic instead is also not good:
#'            - using a Gaussian, centered at the middle of the data range, gives
#'              D~0.17, i.e. q=1-D=0.83, which seems to high
#'            - if only one bin dominates, its position strongly influences D:
#'              e.g. 
#'                (all data at the right (or left)):
#'                ks.test(c(100,100,rep(0,100)), y="punif", min=min(x), max=max(x)) ==> d=0.98 (bad fit, thus good metric)
#'                (all data at the center):
#'                ks.test(c(0,100,rep(50,100)), y="punif", min=min(x), max=max(x))  ==> d=0.57 (similarly bad fit, but metric just does not reflect it)
#'
#' @param x Vector of values from distribution 1
#' @param y Vector of values from distribution 2 (or assumed uniform if omitted)
#' @return p.value of KS test
#' 
qualUnifKS = function(x, y = NULL)
{
  if (is.null(y)) return (ks.test(x, y="punif", min=min(x), max=max(x))$p.value)
  else return (ks.test(x, y)$p.value)
}

#'
#' Discete ChiSquare-test on raw (unbinned) data
#' 
#' Data are first assigned to bins of equal width using 'hist()'.
#' 
#' If only 'x' is given, its distribution of compared to a uniform distribution.
#' If 'y' is given as well, 'x' vs. 'y' is compared (after binning).
#' If counts of 'x' and 'y' differ, the binned counts of the vector with higher counts 
#' are rescaled by N/M to match the smaller one (this might lead to bias).
#' 
#' 
#'
#'
# qualUniform_C2 = function(x, y = NULL, bin_count = 30)
# {
#   
#   if (!is.null(y)) {
#     xy = c(x,y)
#     h = hist(c(x,y), breaks = bin_count, plot = FALSE)
#     hx = hist(x, breaks = h$breaks, plot = FALSE)
#     hy = hist(y, breaks = h$breaks, plot = FALSE)
#     if (length(x) != length(y)) { ## rescale
#       fac = length(x) / length(y)
#     }
#   }
#   else xy = x;
#   
#   
#   
#   chisq.test
#   
# }

#' 
#' Compute position of each element in a vector, assuming the rest of the data is a Gaussian
#' 
#' For each index i, compute a Gaussian G from x \\ x_i and assess the prob of x_i in G.
#' The best result is 1, i.e. x_i is in the center of all other values.
#' The worst results is 0 (for very far outliers), i.e. the value is very far from the median.
#' 
#' Problem: magic constant in here...
#' 
#' @param x Numeric vector of values (e.g. counts for charge state of 2)
#' @return Numeric vector of qualities, same length of 'x', ranging from 0 (bad agreement) to 1 (perfect)
#' 
qualGauss = function(x)
{
  r = rep(NA, length(x))
  for (i in 1:length(x))
  {
    med = median(x[-i])
    gsd = 0.05 * med; # constant 5%, otherwise some datapoint is always bad;  old: "mad(x[-i], center = med)"
    p_best = dnorm(med, mean = med, sd = gsd)
    p_obs =  dnorm(x[i], mean = med, sd = gsd)
    r[i]  = p_obs/p_best  ## ratio of current position divided by best possible position (=center)
  }
  return(r)
}



#' 
#' Compute probability of Gaussian (mu=m, sd=s) at a position 0, with reference 
#' to the max obtainable probability of that Gaussian at its center.
#' 
#' Measure for centeredness around 0.
#' Highest score is 1, worst score is 0.
#' 
#' @param mu Center of Gaussian
#' @param sd SD of Gaussian
#' @return quality, ranging from 0 (bad agreement) to 1 (perfect, i.e. centered at 0)
#' 
qualGaussDev = function(mu, sd)
{
  q = dnorm(0, mean = mu, sd = sd) / dnorm(mu, mean = mu, sd = sd)
  return (q)
}

#'
#' Score an empirical density distribution of values, where the best possible distribution is right-skewed.
#' 
#' The score is computed according to
#' 
#' q = ((N-1) - sum_i( ((N-i-1)*x_i) ) / (N-1)
#' 
#' Scores range from 0 (worst), to 1 (best).
#' E.g. c(0,0,0,16) would yield a score of 1.
#' c(16,0,0,0,0) gives a score of 0.
#' 
#' @param x Vector of numeric values (e.g. height of histogram bins)
#' @param N Length of x (just a precaution currently)
#' @return Quality score in the range of [0,1]
#' 
#' @export
#' 
#' @examples
#'  qualHighest(c(0,0,0,16), 4)   ## 1
#'  qualHighest(c(16,0,0,0), 4)   ## 0 
#'  qualHighest(c(1,1,1,1), 4)    ## 0.5
#'  qualHighest(c(0,16,0,0), 4)   ## 1/3
#'  
qualHighest = function(x, N)
{
  if (length(x)!=N) stop("Error in qualHighest(): length of x not equal to N")
  
  # normalize
  x = x / sum(x)
  p = 1 ## exponent needs to be >=1, using larger values will only mildly penalize intermediate values
  worst = (N-1)^p
  penalty = sum( (x*(N:1 - 1))^p )
  q = (worst - penalty) / worst
  return (q)
}



#'
#' Given a vector of values sampled from a mixture of two Gaussians, compute the area contribution
#' of the higher/more abundant (=target) Gaussian.
#' 
#' Returns a score between 1 (Target only) and 0 (Decoy only). Higher scores are better.
#' We assume that the target distribution is much more peaked and narrow, whereas the decoy (false matches) are 
#' uniformly distributed across all bins (the Gaussian has a large SD).
#' If there is not enough data, NA is returned.
#' 
#' @note This is pre-tuned for Match-time-differences
#'
#' @param x Vector of numeric values
#' @param debug Set to TRUE to generate a plot, showing the fitted Gaussians to the data
#' @return Ratio of area of the target distribution vs. whole area under the two Gaussians
#'
#' @importFrom mixtools normalmixEM2comp
#' 
#' 
qualGauss2Mix = function(x, debug = FALSE)
{
  d = na.omit(x)
  f = file()
  sink(file=f) ## silence output of 'normalmixEM2comp'
  fit2 = normalmixEM2comp(d, lambda = c(0.3,0.6), mu = c(mean(d),mean(d)), c(1.5,0.3), maxit = 10000)
  sink() ## undo silencing
  close(f)
  if (inherits(fit2, "try-error"))
  { ## not enough data?!
    return (NA);
  }
  if (debug)
  {
    fit2$mu
    fit2$sigma
    fit2$lambda
    h=hist(d, 100, plot=F)
    x = seq(-4,4,by=0.1)
    plot(h$mids, h$density)
    lines(x, dnorm(x, fit2$mu[1], fit2$sigma[1]) *fit2$lambda[1], col="red")
    lines(x, dnorm(x, fit2$mu[2], fit2$sigma[2]) *fit2$lambda[2], col="green")
    lines(x, dnorm(x, fit2$mu[1], fit2$sigma[1]) *fit2$lambda[1] + dnorm(x, fit2$mu[2], fit2$sigma[2]) *fit2$lambda[2])
  }
  ## look at scaling of the one with better SD (order defined above)
  q = fit2$lambda[2] / sum(fit2$lambda)
  #q
  return (q)
}



#'
#' From a list of vectors, compute all vs. all Kolmogorov-Smirnoff distance statistics (D)
#' 
#' ... and report the row of the matrix which has maximum sum (i.e the best "reference" distribution).
#' The returned data.frame has as many rows as distributions given and two columns.
#' The first column 'name' gives the name of the list element, the second column 'ks_best' gives '1-statistic' of the
#' Kolmogorov-Smirnoff test to the "reference" distribution (which was picked by maximising the sum of 'ks_best').
#' Thus, the row with a 'ks_best' of 1 is the reference distribution.
#' 
#' @param x List of vectors, where each vector holds a distribution
#' @return A data.frame with ks-test values of the "reference" to all other distributions (see Details)
#'
#'
qualBestKS = function(x) {
  
  if (class(x) != "list") stop("Function bestKS() expects a list!")
  
  ## result matrix
  rr = matrix(0, nrow=length(x), ncol=length(x))
  ## compute upper diagonal
  for (i in 1:length(x))
  {
    for (j in (i+1):length(x))
    {
      if (j>length(x)) next;
      rr[i,j] = 1 - ks.test(x[[i]], x[[j]])$statistic
    }
  }
  ## fill diagonal with 1's
  diag(rr) = 1
  ## add up all values for each reference distribution (equivalent to just rowSums on full matrix)
  rr_sums = rowSums(rr, na.rm=T) + colSums(rr, na.rm=T)
  ## pick best
  r_max = which.max(rr_sums)
  ## get values
  v = pmax(rr[r_max,], rr[, r_max])
  
  return(data.frame(name = names(x), ks_best = v))
         
}



qualMBR = function(d_evd)
{
  d_evd$hasMTD = !is.na(d_evd$Match.time.difference)
  wronglyInferred = ddply(d_evd[d_evd$Type!="MSMS",], c("Raw.file", "Modified.sequence", "Charge"), function(x)
  {
    ratio =  NA
    if (nrow(x)==2 & sum(x$hasMTD)==1) ratio = x$Intensity[!x$hasMTD] / x$Intensity[x$hasMTD]
    #r = rep(NA, nrow(x))
    ## mixed state
    #if (sum(x$hasMTD) >0 & sum(x$hasMTD)<nrow(x))  {
    #  r[!x$hasMTD]
    #} 
    return(data.frame(nNative = sum(!x$hasMTD), nMatched = sum(x$hasMTD), ratio = ratio))
  })
  head(wronglyInferred)
  mbr_score = ddply(wronglyInferred, "Raw.file", function(wronglyInferred)
  {
    ddt = table(wronglyInferred[,c("nMatched", "nNative")])
    ## add nMatched=1 row, if not present
    ddt = rbind(ddt, 0)
    ## native split frequency:
    n.perf = sum(ddt[,colnames(ddt)==1])
    n.all = sum(ddt[,!colnames(ddt)==0])
    corr.nat = n.perf/n.all ## 0.959
    ## inferred splits:
    i.perf = ddt[1,colnames(ddt)==1] + max(0,ddt[2,colnames(ddt)==0])
    i.all = sum(ddt)
    i.perf/i.all  ## 0.923
    
    ## correctly inferred %
    corr.inf = max(0,ddt[2,colnames(ddt)==0]) / (sum(ddt[,colnames(ddt)==0]) + sum(ddt[-1,colnames(ddt)==1]))
    ## 83%
    wrong.inf = 1-corr.inf
    ##
    wrong.inf.ref = wrong.inf * corr.nat
    data.frame(corr.nat = corr.nat, wrong.inf.ref = wrong.inf.ref)
  })
  return (list(mbr_score = mbr_score, ratio = wronglyInferred$ratio))
}

