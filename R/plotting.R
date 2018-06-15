library(ggplot2)

map2color <- function(x, pal, limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

plot_flat_tree <- function(job, branch_orientation, prediction, pcex=1,
                           plot_title = "", col_pal = NA, time_step = 50) {
  params <- read.table(file = paste(job, "cellparams.txt", sep = "_"), sep = "\t", header = T, row.names = 1)
  times <- params$pseudotime
  branches <- params$branches + 1
  branch_names <- sort(unique(branches))

  if (is.na(col_pal)) {
    n <- length(unique(prediction))
    if (n < 25) { # we are in branch territory
      col_pal <- distinctcolors(n, show = FALSE, bw = TRUE)
      cols <- col_pal[as.factor(prediction)]
    } else { # we are in pseudotime territory
      col_pal <- viridis(100)
      cols <- map2color(prediction, col_pal)
    }
  }

  par(mar = c(2., 1., 2., 1.))
  plot(times, branches, ylim=c(min(branch_orientation$to), max(branch_orientation$to)),
       type="n", axes = FALSE, main=plot_title, xlab = "pseudotime", ylab = "")
  time_labs <- seq(min(min(times), 0), max(times), by = time_step)
  axis(1, at=time_labs, labels = time_labs)
  for (i in seq_along(branch_names)) {
    b <- branch_names[i]
    x <- times[branches == b]
    poss_y <- seq(branch_orientation$from[i], branch_orientation$to[i], length.out=length(x))
    y <- poss_y[x - min(x) + 1]
    points(x, y, pch="|", col=cols[branches == b], cex=pcex)
  }
}
