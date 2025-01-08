
#'
#' Plot a text as graphic using ggplot2.
#' 
#' @param title The title of the plot
#' @param text Centered text, can contain linebreaks
#' @param col Colour of text (excluding the title)
#' @return ggplot object
#' 
#'
ggText = function(title, text, col = "black") {
  pl = ggplot(data.frame(text = text, ypos=1, xpos=1), 
                       aes(x = .data$xpos, y = .data$ypos))  +
    geom_text(aes(label = .data$text), colour = col, family="mono") +
    theme_bw() +
    theme(plot.margin = grid::unit(c(1,1,1,1), "cm"), line = element_blank(), axis.title = element_blank(), panel.border = element_blank(),
                   axis.text = element_blank(), strip.text = element_blank(), legend.position="none") +
    ggtitle(title)
  return(pl)
}

#'
#' Augment a ggplot with footer text
#' 
#' @param gg_obj ggplot2 object to be printed
#' @param bottom_left Footer text for bottom left side
#' @param bottom_right Footer text for bottom right side
#' @return -
#'
printWithFooter = function(gg_obj, bottom_left = NULL, bottom_right = NULL) 
{
  if ("gtable" %in% class(gg_obj)) plot(gg_obj) else print(gg_obj)
  
  if (!is.null(bottom_left))
  {
    label = grid::textGrob(bottom_left,
                           x = 0.02,  # left side
                           y = 0.0,   # bottom
                           just = "left", 
                           hjust = NULL,
                           vjust = -.5,
                           gp = grid::gpar(fontsize=7, col="#333333"))  
    grid::grid.draw(label)
  }
  if (!is.null(bottom_right))
  {
    label = grid::textGrob(bottom_right,
                           x = 0.98,  # right side
                           y = 0.0,   # bottom
                           just = "right", 
                           hjust = NULL,
                           vjust = -.5,
                           gp = grid::gpar(fontsize=7, col="#333333"))
    grid::grid.draw(label)
  }
}

#'
#' Inverse the order of items on the x-axis (for discrete scales)
#'
#' @param values The vector of values as given to the x aestetic
#' @param ... Other arguments forwarded to 'scale_y_discrete()'
#' @return ggplot object, concatenatable with '+'
#'
#' @import ggplot2
#' 
scale_x_discrete_reverse = function(values, ...)
{
  if (inherits(values, "factor")) values = droplevels(values)
  return (scale_x_discrete(limits = rev(levels(values)), ... ))
}

#'
#' Inverse the order of items on the y-axis (for discrete scales)
#'
#' @param values The vector of values as given to the y aestetic
#' @param ... Other arguments forwarded to 'scale_y_discrete()'
#' @return ggplot object, concatenatable with '+'
#' 
#' @import ggplot2
#'
scale_y_discrete_reverse = function(values, ...)
{
  if (inherits(values, "factor")) values = droplevels(values)
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


#' 
#' Distribute a set of points with fixed y-values on a stretch of the x-axis.
#' 
#' #' 
#' Usage: 
#'   ggplot(...) + geom_X(...) + pointsPutX(...)
#'   
#' @param x_range [min,max] valid range of x-values
#' @param x_section [min,max] fraction in which to distribute the values (in [0,1] for min,max, e.g. c(0.03,0.08) for 3-8\%)
#' @param y Y-values
#' @param col Colour of the points (used as argument to aes(colour=))
#' @return ggplot object with new geom_point 
#' 
#' @import ggplot2
#' @export
#' 
pointsPutX = function(x_range, x_section, y, col = NA){
  #require(ggplot2)
  x_dist = dist(x_range)
  x_pos = x_range[1] + x_section*x_dist
  d = data.frame(x_ = seq(x_pos[1], x_pos[2], length.out=length(y)), y_ = y, col_= col)
  
  if (is.na(d$col_[1])) pl = geom_point(data = d, aes(x = .data$x_, y = .data$y_))
  else pl = geom_point(data = d, aes(x = .data$x_, y = .data$y_, col = .data$col_))
  
  return (pl)
}


#'
#' A blank theme (similar to the deprecated theme_blank())
#' 
#' @return A ggplot2 object, representing an empty theme
#' 
#' @import ggplot2
#' @export
#' 
theme_blank = function()
{
  theme_blank = theme_bw()
  theme_blank$line = element_blank()
  theme_blank$rect = element_blank()
  theme_blank$strip.text = element_blank()
  theme_blank$axis.title = element_blank()
  return (theme_blank)
}

#'
#' Return color brew palettes, but fail hard if number of requested colors
#' is larger than the palette is holding.
#' 
#' Internally calls 'brewer.pal(n, palette)', checking 'n' beforehand.
#' 
#' @param n Number of colours
#' @param palette Name of palette (e.g. "set1")
#' @return character vector of colors
#' 
brewer.pal.Safe = function(n = 3, palette = "Set1")
{
  idx = which(rownames(RColorBrewer::brewer.pal.info) == palette)
  if (length(idx) != 1) stop("Palette ", palette," unknown!")
  if (RColorBrewer::brewer.pal.info$maxcolors[idx] < n) stop("Palette ", palette, " provides ", RColorBrewer::brewer.pal.info$maxcolors[idx], " colors, but not ", n, " as requested!")
  ## avoid warning about less than 3 colors, if n<3
  return (RColorBrewer::brewer.pal(max(3,n), palette)[1:n])
}
