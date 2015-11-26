#' Create a principal component analysis (PCA) plot for the first two dimensions.
#' 
#' @param data               Matrix(!) where each row is one high-dimensional point, with ncol dimensions, 
#'                           e.g. a mouse as an array of proteinexpressions
#'                           rownames(data) give classes for colouring (can be duplicates in matrices, as opposed to data.frames)
#' @param do_plot            Show PCA plot? if ==2, then shows correlations plot as well
#' @param connect_line_order Connect points by lines, the order is given by this vector.
#'                           Default: NA (no lines)
#' @param gg_layer           More parameters added to a ggplot object (ggplot(x) + gg_layer)                 

#' @return [invisible] Named list with "PCA": The PCA object as returned by \code{\link[stats]{prcomp}}, access $x for PC values
#'                                 and "plots": list of plot objects (one or two)
#' 
#' @import ggplot2
#' @import stats
#' @importFrom seqinr circle
#' 
#' @export
#' 
getPCA = function(data, do_plot = TRUE, connect_line_order = NA, gg_layer)
{
  #require(ggplot2)
  #require(stats)
  #require(seqinr)
  
  pc = prcomp(data, scale. = TRUE)
  useOrd = !is.na(connect_line_order[1])
  # create data frame with scores
  scores = as.data.frame(pc$x)
  scores$class = rownames(data)
  rownames(scores) = paste0(rownames(scores),1:nrow(scores))
  show.Line = !is.na(connect_line_order[1])
  if (show.Line)
  {
    scores$pathX = scores$PC1[connect_line_order]
    scores$pathY = scores$PC2[connect_line_order]
  }
  if (do_plot & show.Line)
  {
    scores$ord = connect_line_order / max(connect_line_order) * 10
  } else
  {
    scores$ord = as.numeric(factor(rownames(data)))
  }
  
  lpl = list();
  
  # plot of observations
  if (do_plot)
  {
    
    var_PC1 = round(pc$sdev[1]^2 / sum(pc$sdev^2) * 100)
    var_PC2 = round(pc$sdev[2]^2 / sum(pc$sdev^2) * 100)
    
    pl = ggplot(data=scores, aes_string(x = "PC1", y = "PC2"));
    if (show.Line) pl = pl + geom_path(aes_string(x = "pathX", y = "pathY", color = "ord", alpha = "0.5"))
    pl = pl + scale_colour_gradient(low="red", high="darkgreen")
    pl = pl + 
              gg_layer + 
              geom_point(aes_string(colour = "ord"), size = 1) +
              geom_text(size = 3, angle=0, aes_string(label = "class", colour = "ord", vjust = "1")) +
              geom_hline(yintercept=0, colour="gray65") +
              #theme(panel.background=element_rect("black")) +
              geom_vline(xintercept=0, colour="gray65") +
              theme_bw() + 
              xlab(paste0("PC #1 (", var_PC1, "% var)")) +
              ylab(paste0("PC #2 (", var_PC2, "% var)"))
    lpl[[1]] = pl
  }
  
  correlations = NA
  
  if (do_plot==2)
  {
    circle <- function(center=c(0,0), npoints=100)
    {
      r = 1
      tt = seq(0, 2*pi, length=npoints)
      xx = center[1] + r * cos(tt)
      yy = center[1] + r * sin(tt)
      return(data.frame(x = xx, y = yy))
    }
    corcir = circle(c(0,0), npoints = 100)
    # create data frame with correlations between variables and PCs
    correlations = as.data.frame(cor(data, pc$x))
    # data frame with arrows coordinates
    arrows = data.frame(x1=rep(0, ncol(data)), y1=rep(0, ncol(data)),
                        x2=correlations$PC1, y2=correlations$PC2)
    
    
    lpl[[2]] = 
      ggplot() +
      geom_path(data=corcir, aes_string(x = "x", y = "y"), colour="gray65") +  ## open circles
      geom_segment(data=arrows, aes_string(x = "x1", y = "y1", xend = "x2", yend = "y2"), colour="gray65") +
      geom_text(data=correlations, aes_string(x  = "PC1", y = "PC2", label = "rownames(correlations)")) +
      geom_hline(yintercept = 0, colour = "gray65") +
      geom_vline(xintercept = 0, colour = "gray65") +
      xlim(-1.1,1.1) + ylim(-1.1,1.1) +
      labs(x="pc1 axis", y="pc2 axis") +
      ggtitle("Circle of correlations")

  }
  return (list("PCA" = invisible(pc), "plots" = lpl, "correlations" = correlations))
}




