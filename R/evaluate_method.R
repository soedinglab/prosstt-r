#' Evaluate a lineage tree prediction by all available methods
#'
#' Run all implemented methods for branch assignment and pseudotime prediction
#' evaluation.
#'
#' @param method the name of the assessed trajectory inference algorithm
#' @param mbranches the predicted branch assignments
#' @param mtimes the predicted pseudotime assignments
#' @param cell_params a data frame that includes branch and pseudotime labels
#' @param par_loc location of the PROSSTT parameter file
#' @param returnNA if set to TRUE will return an empty result matrix
#'
#' @return an one-row matrix with the evaluation results for all methods
#'
#' @export
evaluate_method <- function(method, mbranches, mtimes, cell_params, par_loc, returnNA=FALSE) {
  funcnames <- c("rand index", "matthews corr. coef.",
                 "F1 measure", "jaccard index", "fowkles-mallows index", "adjusted MI")

  timenames <- c("goodman-kruskall (unweighted)", "goodman-kruskall (weighted)",
                 "kendall index (unweighted)", "kendall index (weighted)")

  res <- data.frame(matrix(NA, nrow = length(funcnames)+length(timenames)+1, ncol = 1))
  rownames(res) <- c(funcnames, timenames, "branches")
  colnames(res) <- method
  # print(res)

  if (returnNA) { return(t(res)) }

  mstatus <- assign_status(mbranches, cell_params$branches+1)

  res[1, 1] <- randInd_manual(mstatus)
  res[2, 1] <- matthews_cor(mstatus)
  res[3, 1] <- f_measure(mstatus)
  res[4, 1] <- jaccard(mstatus)
  res[5, 1] <- fowkles_mallows(mstatus)
  res[6, 1] <- adjusted_mi(prediction = mbranches, truth = cell_params$branches+1)

  # longest path goodman kruskal:
  lp_id <- get_lpgk_indices(par_loc, cell_params)
  ltimes <- cell_params$pseudotime[lp_id]
  mtimes <- mtimes[lp_id]

  res[7, 1] <- goodman_kruskal_index(ltimes, mtimes, weighted=FALSE)
  res[8, 1] <- goodman_kruskal_index(ltimes, mtimes, weighted=TRUE)
  res[9, 1] <- kendall_index(ltimes, mtimes, weighted=FALSE)
  res[10, 1] <- kendall_index(ltimes, mtimes, weighted=TRUE)

  # add number of predicted branches
  num_branches <- length(levels(as.factor(mbranches)))
  res[length(funcnames) + length(timenames) + 1,] <- num_branches

  return(t(res))
}

#' Identifies the longest path in the simulation
#'
#' Finds the indices of the cells on the longest path in the simulated lineage
#' tree.
#'
#' @param par_loc location of the PROSSTT parameter file
#' @param cell_params a data frame that includes branch and pseudotime labels
#'
#' @return the indices of the cells on the longest path
#'
#' @export
#'
#' @importFrom igraph graph_from_adjacency_matrix shortest_paths
get_lpgk_indices <- function(par_loc, cell_params) {
  cell_params$branches <- cell_params$branches + 1
  time <- cell_params$pseudotime - min(cell_params$pseudotime) + 1

  StartCell <- rownames(cell_params)[which(time == min(time))][1]
  EndCell <- rownames(cell_params)[which(time == max(time))][1]

  StartBranch <- cell_params[StartCell,]$branches
  EndBranch <- cell_params[EndCell,]$branches

  results <- system(sprintf("awk '{if($0~/topology/) print}' %s", par_loc), intern = T)
  digits <- gregexpr("[0-9]+", results)
  Branches <- as.numeric(regmatches(results, digits)[[1]])
  # if the simulation contains a linear topology then the line will read
  # topology: []
  # of course this means that digits will not contain any numbers,
  # so the regex will not match and we have to account for that
  if (length(Branches) == 0) {
    Branches <- 1
  } else {
    Branches = Branches + 1
  }

  # Create adjacency matrix and calculate shortest path in between start and end branches
  BranchConnections <- matrix(0, max(Branches), max(Branches))
  for (i in seq(from=1, to=length(Branches), by=2)) {
    BranchConnections[Branches[i], Branches[i+1]]=1
    BranchConnections[Branches[i+1], Branches[i]]=1
  }

  Graph_Branches <- igraph::graph_from_adjacency_matrix(BranchConnections)
  LongestPath <- unlist(igraph::shortest_paths(Graph_Branches, from = StartBranch, to = EndBranch)$vpath)

  LongestPathCellsID <- which(cell_params$branches %in% LongestPath)
  return(LongestPathCellsID)
}

#' Appends a result to a file
#'
#' Appends an one-row matrix to the end of a file with a matrix of the same
#' column number (with same column names). Useful for benchmarking multiple
#' methods and keeping all results in one file.
#'
#' @param res the result matrix (produced by `evaluate_method()`)
#' @param out the location of the file to append to
#'
#' @return None
#'
#' @export
#'
#' @importFrom utils read.table write.table
read_write_output <- function(res, out) {
  if (!file.exists(out)) {
    write.table(res, out)
  } else {
    old <- read.table(out, check.names=FALSE)
    res <- rbind(old, res)

    write.table(res, out, append=F)
  }
}
