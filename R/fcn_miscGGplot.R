
#'
#' Plot a text as graphic using ggplot2.
#' 
#' @param title The title of the plot
#' @param text Centered text, can contain linebreaks
#' @param col Colour of text (excluding the title)
#' @return ggplot object
#' 
#' @import ggplot2
#' @importFrom grid unit
#'
ggText = function(title, text, col = "black") {
  pl = ggplot(data.frame(text = text, ypos=1, xpos=1), 
              aes_string(x = "xpos", y = "ypos"))  +
    geom_text(aes_string(label = "text"), colour = col, family="mono") +
    theme_bw() +
    theme(plot.margin = unit(c(1,1,1,1), "cm"), line = element_blank(), axis.title = element_blank(), panel.border = element_blank(),
          axis.text = element_blank(), strip.text = element_blank(), legend.position="none") +
    ggtitle(title)
  return(pl)
}


#'
#' Inverse the order of items on the x-axis (for discrete scales)
#'
#' @param values The vector of values as given to the x aestetic
#' @param ... Other arguments forwarded to 'scale_y_discrete()'
#' @return ggplot object, concatenatable with '+'
#'
scale_x_discrete_reverse = function(values, ...)
{
  if (!("factor" %in% class(values))) stop("Cannot use scale_x_discrete_reverse() on non-factor()")
  values = droplevels(values)
  return (scale_x_discrete(limits = rev(levels(values)), ... ))
}

#'
#' Inverse the order of items on the y-axis (for discrete scales)
#'
#' @param values The vector of values as given to the y aestetic
#' @param ... Other arguments forwarded to 'scale_y_discrete()'
#' @return ggplot object, concatenatable with '+'
#'
scale_y_discrete_reverse = function(values, ...)
{
  if (!("factor" %in% class(values))) stop("Cannot use scale_y_discrete_reverse() on non-factor()")
  values = droplevels(values)
  return (scale_y_discrete(limits = rev(levels(values)), ... ))
}

#'
#' Function to thin out the number of labels shown on an axis in GGplot
#' 
#' By default, 20 labels (or up to 40 see below) are shown.
#' If the number of items is less than twice the number of desired labels,
#' all labels will be shown (to avoid irregular holes for some labels).
#' I.e. if n=20, and x has 22 entries, there would be only two labels removed, giving a very irregular picture.
#' It only becomes somewhat regular if after any label there is at least one blank, i.e. at most
#' half the entries are labeled.
#' #' 
#' Example: 
#'   ## p is any ggplot object
#'   p + scale_y_discrete(breaks = ggAxisLabels) 
#'   ## customize 'n'
#'   my.ggAxisLabels = function(x) ggAxisLabels(x, n = 4)
#'   p + scale_y_discrete(breaks = my.ggAxisLabels) 
#' 
#' @param x Vector of labels (passed by GGplot)
#' @param n Number of labels to show
#' @return Shortened version of 'x'
#'
#'
ggAxisLabels = function(x, n = 20)
{
  ## avoid irregular holes in labels if |x| < 2n (some neighbouring labels would be present, while others not)  
  ## giving patterns like 1101011010110. So we might as well use all labels, ie. 11111111111111
  if (length(x) <= 2*n) n = n*2
  
  breaks = x[seq(1, length(x), length.out = n)]
  names(breaks) = x[breaks]
  return (breaks)
}


#' Add title and subtitle to a ggplot
#' 
#' Found in http://www.antoni.fr/blog/?p=39 .. whewww... modified a little
#' 
#' @param pl GGplot object
#' @param main String for main title
#' @param sub String for sub title
#' @return A ggplot object
#' 
#' @import ggplot2
#' 
#' @export
#' 
addGGtitle <- function(pl, main, sub=""){
  #require(ggplot2)
  if(sub==""){
    pl <- pl + ggtitle(eval(parse(text=paste("expression(atop(\"",main,"\",","))", sep=""))))
  }else{
    pl <- pl + ggtitle(eval(parse(text=paste("expression(atop(\"",main, "\",", " atop(\"", sub , "\",\"\")))", sep=""))))
  }
  return (pl)
}

