#' Normalizes expression matrix per cell.
#'
#' Normalizes the expression matrix per row (cell) by dividing with the average
#' library size (row sum) and multiplying with the original row sum.
#'
#' @param X the expression matrix to be normalized
#'
#' @return normalized expression matrix
#'
#' @export
norm_lib_size <- function(X) {
  sum1 <- apply(X, 1, sum)
  scalings <- sum1 / mean(sum1)
  X <- (1 / scalings) * X
  return (X)
}

#' Calculates the expected mutual information
#'
#' Calculates the mutual information background, assuming a hypergeometric 
#' distribution. Code comes from Stack Overflow 
#' (\url{http://stackoverflow.com/a/38010916/7195775}).
#'
#' @param s1 counts per cluster for the prediction
#' @param s2 counts per cluster for the truth
#' @param l1 length of s1
#' @param l2 length of s2
#' @param n number of samples
#'
#' @return expected Mutual Information
#'
#' @export
expected_mi <- function(s1, s2, l1, l2, n){    # expected mutual information
  s_emi <- 0
  for(i in 1:l1){
    for (j in 1:l2){
      min_nij <- max(1,s1[i]+s2[j]-n)
      max_nij <- min(s1[i],s2[j])
      n.ij <- seq(min_nij, max_nij)   # sequence of consecutive numbers
      t1<- (n.ij / n) * log((n.ij * n) / (s1[i]*s2[j]))
      t2 <- exp(lfactorial(s1[i]) + lfactorial(s2[j]) + lfactorial(n - s1[i]) + lfactorial(n - s2[j]) - lfactorial(n) - lfactorial(n.ij) - lfactorial(s1[i] - n.ij) - lfactorial(s2[j] - n.ij) - lfactorial(n - s1[i] - s2[j] + n.ij))
      emi <- sum(t1*t2)
      s_emi <- s_emi + emi
    }
  }
  return(s_emi)
}

#' Calculates adjusted (normalized) Mutual Information
#'
#' Calculates the mutual information corrected for a hypergeometric-distributed
#' background. Code comes from \url{http://stackoverflow.com/a/38010916/7195775}.
#'
#' @param prediction the predicted cluster assignments for each cell
#' @param truth the true cluster assignments for each cell
#'
#' @return the adjusted Mutual Information for the prediction/truth pair
#'
#' @export
#'
#' @importFrom infotheo mutinformation
adjusted_mi <- function(prediction, truth){
  status <- assign_status(prediction, truth)
  if (all(status == 0)) {
    return (NaN)
  }

  prediction <- as.numeric(prediction)
  truth <- as.numeric(truth)

  # apparently NA values in prediction are just ignored D:
  # so we need to compensate
  prediction[is.na(prediction)] <- 0
  prediction[prediction==0] <- max(prediction) + 1

  prediction = as.factor(prediction)
  truth = as.factor(truth)

  s1 <- tabulate(prediction)
  s2 <- tabulate(truth)
  l1 <- length(s1)
  l2 <- length(s2)
  N <- length(prediction)
  tij <- table(prediction, truth)
  mi <- infotheo::mutinformation(prediction, truth)

  # avoid log(0)
  ls1 <- s1/N
  ls1[ls1==0] <- 1
  ls1 <- log(ls1)

  ls2 <- s2/N
  ls2[ls2==0] <- 1
  ls2 <- log(ls2)

  h1 <- -sum(s1*ls1)/N
  h2 <- -sum(s2*ls2)/N

  emi <- expected_mi(s1, s2, l1, l2, N) # EMI Expected MI
  ami <- (mi-emi)/(max(h1,h2) - emi)  #AMI Adjusted MI
  return(ami)
}

#' Calculates the truth table for a branch assignment prediction.
#'
#' Goes over all pairs of points in the prediction and truth and classifies
#' them:
#' - TP: pair of points in same cluster in prediction, truth
#' - TN: pair of points in diff. cluster in prediction, truth
#' - FN: pair of points in diff. cluster in prediction, same in truth
#' - FP: pair of points in same cluster in prediction, diff. in truth
#' NA values in the predictions are considered FNs.
#'
#' @param prediction predicted cluster assignments
#' @param truth true cluster assignments
#'
#' @return the truth table for this prediction given the truth: (TP, TN, FN, FP)
#'
#' @export
assign_status <- function(prediction, truth) {
  if (!all(is.na(prediction))) { # this will happen when a method is unable to run
    bu <- outer(prediction, prediction, "==")
    bv <- outer(truth, truth, "==")
    TP <- length(intersect(which(bu), which(bv))) - length(prediction)
    TN <- length(intersect(which(!bu), which(!bv)))
    FP <- length(intersect(which(bu), which(!bv)))
    FN <- length(intersect(which(!bu), which(bv)))
    FN <- FN + length(which(is.na(bu)))
    return(c(TP, TN, FN, FP)/2)
  }
  return(c(0,0,0,0))
}

#' Calculates the (unweighted) Rand index
#'
#' Calculates the (unweighted) Rand index:
#' \deqn{\frac{\mathrm{TP} + \mathrm{TN}}{\mathrm{TP} + \mathrm{FP} + \mathrm{TN} + \mathrm{FN}}}
#'
#' @param status a truth table
#'
#' @return the Rand index of the input truth table
#'
#' @export
randInd_manual <- function(status) {
  TP <- status[1]
  TN <- status[2]
  FN <- status[3]
  FP <- status[4]

  res <- (TP+TN) / (TP+FP+TN+FN)
  return(res)
}

#' Calculates the Matthews Correlation Coefficient.
#'
#' Calculates the Matthews Correlation Coefficient:
#' \deqn{\frac{\mathrm{TP}\mathrm{TN} - \mathrm{FP}\mathrm{FN}}{\sqrt{(\mathrm{TP}+\mathrm{FP}) (\mathrm{TP}+\mathrm{FN}) (\mathrm{TN}+\mathrm{FP}) (\mathrm{TN}+\mathrm{FN})}}}
#'
#' @param status a truth table
#'
#' @return the Matthews Correlation Coefficient of the input truth table
#'
#' @export
matthews_cor <- function(status) {
  TP <- status[1]
  TN <- status[2]
  FN <- status[3]
  FP <- status[4]

  res <- (TP*TN - FP*FN) / sqrt((TP+FP) * (TP+FN) * (TN+FP) * (TN+FN))
  return(res)
}

#' Calculates the F measure
#'
#' Calculates the F measure, which measures the effectiveness of retrieval with
#' respect to a user who attaches \code{b} times as much importance to recall
#' \code{R} as precision \code{P}:
#' 
#' \deqn{\frac{(b^2 + 1)\mathrm{P} \mathrm{R}}{b^2 \mathrm{P} + \mathrm{R}}}
#'
#' @param status a truth table
#' @param b determines the balance between precision and recall
#'
#' @return the F measure of the input truth table
#'
#' @export
f_measure <- function(status, b=1) {
  TP <- status[1]
  TN <- status[2]
  FN <- status[3]
  FP <- status[4]

  P <- TP / (TP+FP)
  R <- TP / (TP+FN)
  res <- (b^2 + 1) * P * R / (b^2 * P + R)
  return(res)
}

#' Calculates the Jaccard index.
#'
#' Calculates the Jaccard index, also known as Intersection over Union:
#'
#' \deqn{\frac{\mathrm{TP}}{\mathrm{TP} + \mathrm{FP} + \mathrm{FN}}}
#'
#' @param status a truth table
#'
#' @return the Jaccard index of the input truth table
#'
#' @export
jaccard <- function(status) {
  TP <- status[1]
  TN <- status[2]
  FN <- status[3]
  FP <- status[4]

  res <- TP / (TP+FP+FN)
  return(res)
}

#' Calculates the Fowkles-Mallows index.
#'
#' Calculates the Fowkles-Mallows index:
#'
#' \deqn{\sqrt{\frac{\mathrm{TP}}{\mathrm{TP} + \mathrm{FP}} \frac{\mathrm{TP}}{\mathrm{TP} + \mathrm{FN}} }}
#'
#' @param status a truth table
#'
#' @return the Fowkles-Mallows index of the input truth table
#'
#' @export
fowkles_mallows <- function(status) {
  TP <- status[1]
  TN <- status[2]
  FN <- status[3]
  FP <- status[4]

  P <- TP / (TP+FP)
  R <- TP / (TP+FN)
  res <- sqrt(P*R)
  return(res)
}

# Calculates the indices for the Goodman-Kruskal and Kendall indices.
#
# Finds positive (S+) and negative (S-) pairs of cells, where S+ are cells
# with the same order in both sets and S- are cells with opposite order.
# Neutral pairs are ignored. Instead of only counting, the distance of the points
# is taken into account. The indices are symmetric, so the input order of
# prediction and truth don't matter.
#
# @param a the first set
# @param b the second set
#
# @return the score for each pair
#
# @export
goodman_kruskal_weights <- function(a, b) {
  l <- length(a)

  a_width <- max(a) - min(a)
  b_width <- max(b) - min(b)

  wA <- matrix(0, ncol = l, nrow = l)
  wB <- matrix(0, ncol = l, nrow = l)

  if (a_width>0) {
    wA <- outer(a, a, '-') / a_width
  }
  if (b_width>0) {
    wB <- outer(b, b, '-') / b_width
  }

  W <- array(0, dim=c(l, l, 2))
  sA <- sign(wA)
  sB <- sign(wB)

  same <- (sA*sB == 1)
  opp <- (sA*sB == -1)
  null <- (sA == 0 & sB == 0)

  W[,,1] = wA / wB
  W[,,2] = wB / wA

  res <- array(0, dim=c(l,l))

  res[same] = pmin(W[,,1], W[,,2])[same]
  res[opp] = pmax(W[,,1], W[,,2])[opp]
  res[null] = 1

  return(res)
}

# Calculates the signs for the Goodman-Kruskal and Kendall indices.
#
# Finds positive (S+) and negative (S-) pairs of cells, where S+ are cells
# with the same order in both sets and S- are cells with opposite order.
# Neutral pairs are ignored. The indices are symmetric, so the input order of
# prediction and truth don't matter.
#
# @param a the first set
# @param b the second set
#
# @return the score for each pair
#
# @export
goodman_kruskal_signs <- function(a, b) {
  l <- length(a)

  wA <- matrix(0, ncol = l, nrow = l)

  wA <- sign(outer(a, a, '-'))
  wB <- sign(outer(b, b, '-'))

  uneq <- (wB != 0)
  null <- (wA == 0 & wB == 0)

  res <- array(0, dim=c(l,l))

  res[uneq] = (wA/wB)[uneq]
  res[null] = 1

  return(res)
}

#' Calculates the Goodman-Kruskal index
#'
#' Calculates the Goodman-Kruskal index, a comparison of two ordered sequences,
#' essentially the fraction of concordant pairs in the two sets if neutral pairs
#' are ignored.
#'
#' @param a ordering of the first sequence
#' @param b ordering of the second sequence
#' @param weighted calculate the weighted version, which takes distances between
#' points into account.
#'
#' @return the Goodman-Kruskal of the input truth table
#'
#' @export
goodman_kruskal_index <- function(a, b, weighted=TRUE) {
  if (length(a) == 1) {
    return(NaN)
  }

  if (weighted) {
    W <- goodman_kruskal_weights(a, b)
  } else {
    W <- goodman_kruskal_signs(a, b)
  }

  tri <- upper.tri(W, diag=FALSE)

  up <- sum(W[tri])
  down <- sum(abs(W)[tri])
  return(up/down)
}

#' Calculates the Kendall index
#'
#' Calculates the Kendall index, a comparison of two ordered sequences,
#' essentially the fraction of concordant pairs in the two sets. Very similar to
#' the Goodman-Kruskal index, but it also considers neutral pairs.
#'
#' @param a ordering of the first sequence
#' @param b ordering of the second sequence
#' @param weighted calculate the weighted version, which takes distances between
#' points into account.
#'
#' @return the Kendall index of the input truth table
#'
#' @export
kendall_index <- function(a, b, weighted=TRUE) {
  if (length(a) == 1) {
    return(NaN)
  }

  if (weighted) {
    W <- goodman_kruskal_weights(a, b)
  } else {
    W <- goodman_kruskal_signs(a, b)
  }

  tri <- upper.tri(W, diag=FALSE)

  up <- sum(W[tri])
  n <- dim(W)[1]
  down <- n*(n-1)/2

  return(up/down)
}
