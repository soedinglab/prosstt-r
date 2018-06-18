#' Translates numbers to colors for plotting.
#'
#' Bins the input number values in as many bins as there are available colors
#' and returns a vector of the appropriate color for each input value.
#'
#' @param x the input vector
#' @param pal the color palette
#' @param limits the limits within x moves. If `NULL`, then defaults to 
#' `range(x)`
#'
#' @return a vector of the appropriate palette color for every number in x
#'
#' @export
map2color <- function(x, pal, limits=NULL) {
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

#' Creates ggplot2 colors
#'
#' Creates colors from the default ggplot2 palette.
#'
#' @param n how many colors to generate
#'
#' @return the n first colors from the ggplot2 palette
#'
#' @examples
#' gg_color_hue(3)
#'
#' @export
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#' Plots a simulation in 2D
#'
#' Arranges simulation cells according to branch and pseudotime and 
#' colors them by their predicted branch or pseudotime.
#'
#' @param cell_params data frame that contains branch and pseudotime labels
#' @param branch_orientation the coordinates of each branch
#' @param prediction predicted pseudotime or branch labels
#' @param pcex controls cell size at plotting
#' @param plot_title the title of the plot
#' @param col_pal a color palette; will be generated if absent
#' @param time_step the step size for the x axis label (pseudotime)
#'
#' @return None
#'
#' @export
plot_flat_tree <- function(cell_params, branch_orientation, prediction, pcex=1,
                           plot_title = "", col_pal = NA, time_step = 50) {
  times <- cell_params$pseudotime
  branches <- cell_params$branches + 1
  branch_names <- sort(unique(branches))

  if (is.na(col_pal)) {
    n <- length(unique(prediction))
    if (n < 25) { # we are in branch territory
      col_pal <- LSD::distinctcolors(n, show = FALSE, bw = TRUE)
      cols <- col_pal[as.factor(prediction)]
    } else { # we are in pseudotime territory
      col_pal <- viridis::viridis(100)
      cols <- map2color(prediction, col_pal)
    }
  }

  par(mar = c(2., 1., 2., 1.))
  plot(times, branches, ylim=c(min(branch_orientation$to), max(branch_orientation$to)),
       type="n", axes = FALSE, main=plot_title, xlab = "pseudotime", ylab = "")
  time_labs <- seq(min(min(times), 0), max(times)+1, by = time_step)
  axis(1, at=time_labs, labels = time_labs)
  for (i in seq_along(branch_names)) {
    b <- branch_names[i]
    x <- times[branches == b]
    poss_y <- seq(branch_orientation$from[i], branch_orientation$to[i], length.out=length(x))
    y <- poss_y[x - min(x) + 1]
    points(x, y, pch="|", col=cols[branches == b], cex=pcex)
  }
}
