library(readr)
library(igraph)
library(viridis)
library(LSD)

read_prosstt_topology <- function(job) {
  params <- read_delim(paste(job, "_params.txt", sep=""), ":", col_names = FALSE, trim_ws = TRUE)
  top_index <- which(params$X1 == "topology")
  group <- strsplit(params$X2[top_index], "],")[[1]]
  processed <- gsub("\\[|\\]|\\s", "", group, perl = TRUE)
  aslist <- strsplit(processed, ",")
  topology <- matrix(0, ncol = 2, nrow = length(aslist))
  for (i in seq_along(aslist)) {
    topology[i, ] <- as.numeric(aslist[[i]])
  }
  topology
}

time_ordered_branch <- function(b, branches, times) {
  x <- which(branches == b)
  y <- times[x]
  x[order(y)]
}

fan_out <- function(n) {
  candidates <- seq(-(n-1) %/% 2, n %/% 2)
  if (n %% 2 == 1) {
    candidates
  } else { # leave out the zero
    candidates[- (n+1)%/%2]
  }
}

to_flip <- function(duplics) {
  tmp <- which(duplics, arr.ind = TRUE)
  keep <- (tmp[,"row"] < tmp[,"col"]) #& (abs(tmp[,"row"] - tmp[,"col"]) == 1)
  tmp[keep,]
}

branch_order <- function(x, zone, parents, branch_orientation) {
  index <- seq_along(x)
  duplics <- outer(x, x, "==")
  diag(duplics) <- FALSE
  flip_pairs <- to_flip(duplics)
  res <- order(x)
  if (is.null(dim(flip_pairs))) {
    where <- match(flip_pairs, res)
    # are the parents of these in the correct order?
    end_one <- branch_orientation[parents[zone[where[1]]], ]$to
    end_two <- branch_orientation[parents[zone[where[2]]], ]$to
    res[where] <- res[where][order(c(end_one, end_two))]
  } else {
    for(i in seq_along(flip_pairs[,1])) {
      where <- match(flip_pairs[i, ], res)
      # are the parents of these in the correct order?
      end_one <- branch_orientation[parents[zone[where[1]]], ]$to
      end_two <- branch_orientation[parents[zone[where[2]]], ]$to
      res[where] <- res[where][order(c(end_one, end_two))]
    }
  }
  res
}

flat_simulation <- function(job, mode = "prosstt") {
  params <- read.table(file = paste(job, "cellparams.txt", sep = "_"), sep = "\t", header = T, row.names = 1)
  if (mode == "prosstt") {
    branches <- params$branches + 1
    topology <- read_prosstt_topology(job) + 1
  } else {
    branches <- params$branches
    topology <- read_prosstt_topology(job)
  }

  N <- length(branches)

  branch_names <- sort(unique(branches))
  num_branches <- length(unique(branches))
  branch_orientation <- data.frame(from = rep(0, num_branches),
                                  to = rep(0, num_branches),
                                  zone = rep(0, num_branches))
  g <- graph_from_edgelist(topology)
  branch_orientation$zone <- distances(g)[1,]
  bla <- table(topology[,1])
  children <- rep(0, num_branches)
  indices <- match(as.numeric(names(bla)), branch_names)
  children[indices] <- bla
  timezones <- sort(unique(branch_orientation$zone))
  parents <- rep(0, max(topology))
  for (i in seq(1, max(topology))) {
    where <- which(topology[,2] == i)
    if (length(where) > 0) {
      parents[i] <- topology[where, 1]
    }
  }

  root <- which(sapply(sapply(V(g), function(x) neighbors(g,x, mode="in")), length) == 0)
  depth_first_search <- bfs(g, root)

  for (n_from in depth_first_search$order) {
    b_from <- which(branch_names == n_from)
    offsets <- fan_out(children[b_from])
    destinations <- topology[,2][topology[,1] == n_from]
    for (i in seq_along(destinations)) {
      n_to <- destinations[i]
      b_to <- which(branch_names == n_to)
      branch_orientation$from[b_to] <- branch_orientation$to[b_from]
      branch_orientation$to[b_to] <- branch_orientation$from[b_to] + offsets[i]
      # plot_flat_tree(job, branch_orientation, params$branches, plot_title = paste(b_from, "to", b_to))
    }
  }

  for (z in timezones[-1]) {
  # for (z in 1:3) {
    zone <- which(branch_orientation$zone == z)
    ord <- order(branch_orientation$to[zone])
    zone <- zone[ord]
    overlaps <- branch_orientation$to[zone]
    if (any(duplicated(overlaps))) {
      replacement <- seq(min(overlaps), min(overlaps) + length(overlaps), length.out = length(overlaps))
      ord <- branch_order(overlaps, zone, parents, branch_orientation)
      branch_orientation$to[zone][ord] <- replacement
    }
    # plot_flat_tree(job, branch_orientation, params$branches, plot_title = paste("zone", z))
  }


  for (n_from in depth_first_search$order) {
    b_from <- which(branch_names == n_from)
    destinations <- topology[,2][topology[,1] == n_from]
    for (i in seq_along(destinations)) {
      n_to <- destinations[i]
      b_to <- which(branch_names == n_to)
      branch_orientation$from[b_to] <- branch_orientation$to[b_from]
    }
  }
  branch_orientation
}

