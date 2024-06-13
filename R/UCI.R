#' UCI: U-statistic permutation test of conditional independence
#'
#' Given a three-way contingency table, the function \code{UCI} computes two test statistics for conditional independence
#' depending on \code{weight} option. When \code{weight} is \code{TRUE}, it calculates the test statistic \eqn{T_W^\dagger}
#' based on weighted U-statistics proposed by \insertCite{kim2022;textual}{UCI}.
#' On the other hand, when \code{weight} is \code{FALSE}, it calculates the test statistic \eqn{T} based on
#' unweighted U-statistics proposed by \insertCite{canonne2018testing;textual}{UCI}. See \insertCite{kim2022;textual}{UCI} for more details.
#' Both test statistics are calibrated by the Monte Carlo permutation method, yielding valid p-values.
#' Note that the computation of \eqn{T_W^\dagger} requires the knowledge of \code{l1} (the dimension of \eqn{X}) and
#' \code{l2} (the dimension of \eqn{Y}). When these are unknown, then the function estimates \code{l1}
#' by the length of the unique elements of \eqn{X} in the dataset (and \code{l2} is similarly estimated).
#'
#' @param dat a three-way frequency table
#' @param l1 the dimension of \eqn{X}
#' @param l2 the dimension of \eqn{Y}
#' @param B the number of random permutations
#' @param weight a logical indicating whether to consider weights; If \code{weight} is \code{TRUE}, the function uses
#' \eqn{T_W^\dagger} based on weighted U-statistics as a test statistic. Otherwise, the function uses \eqn{T} based on unweighted U-statistics
#' as a test statistic. See \insertCite{kim2022;textual}{UCI} for details.
#'
#' @return Returns the p-value for \code{UCI}-test and the value of the test statistic
#' @export
#' @references{\insertAllCited{}}
#' @encoding UTF-8
#' @importFrom stats r2dtable
#' @importFrom Rdpack reprompt
#' @seealso \code{\link{chi2test}} and \code{\link{gtest}}
#'
#' @examples
#' ## sample size and dimensions
#' n = 100; l1 = 5; l2 = 5; d = 5
#'
#' ## Generate categorical data
#' x = sample(1:l1, n, replace = TRUE)
#' y = sample(1:l2, n, replace = TRUE)
#' z = sample(1:d, n, replace = TRUE)
#'
#' ## Generate the corresponding three-way frequency table
#' freq = table(x, y, z)
#'
#' ## Implement the test
#' # weighted statistic
#' UCI(freq, l1, l2, 199, TRUE)
#' # unweighted statistic
#' UCI(freq, 199, FALSE)
UCI <- function(dat = dat, l1 = 0, l2 = 0, B = 199, weight = TRUE) {
  # sample size
  sig <-apply(dat, 3, sum)

  # Permuted data
  Perm <- apply(dat, 3, function(x) {

    rs <- rowSums(x)
    cs <- colSums(x)

    if(min(length(rs), length(cs)) > 1) {
      return(r2dtable(B, rs, cs))
    } else if(length(rs) == 1) {
      M <- list(as.matrix(cs,ncol = 1))
      return(rep(M, B))
    } else {
      M <- list(as.matrix(rs,nrow = 1))
      return(rep(M, B))
    }
  })

  if(weight == TRUE) {
    # weighted U-statistic approach
    # when l1 and l2 are not specified
    # take the union of X over bins and set l1 to be the cardinality of X
    # similarly, take the union of Y over bins and set l2 to be the cardinality of Y
    if(l1 == 0 | l2 == 0) {
      l1 <- length(unique(c(apply(dat, 3, rownames))))
      l2 <- length(unique(c(apply(dat, 3, colnames))))
    }

    # compute weights
    w1 <- apply(cbind(sig, l1), 1, min)
    w2 <- apply(cbind(sig, l2), 1, min)
    w <- sqrt(w1 * w2)

    # original statistic
    Ustats <- apply(dat, 3, Ustat, l1 = l1, l2 = l2, weight = TRUE)

    Perm_M_Ustats <- sapply(1:dim(dat)[3], function(k) unlist(lapply(Perm[[k]], Ustat, l1 = l1, l2 = l2, weight = TRUE)))
    Perm_M_Ustats_c <- rbind(Ustats, Perm_M_Ustats)
    CI_Ustats <- apply(Perm_M_Ustats_c, 1, function(z) {
      sum(z * sig * w)
    })

    pvalU <- 1 / length(CI_Ustats) * sum(CI_Ustats >= CI_Ustats[1])

  } else {
    # unweighted U-statistic approach
    # this approach does not use l1 and l2
    Ustats <- sapply(dat, Ustat, weight = FALSE)

    Perm_M_Ustats <- sapply(1:dim(dat)[3], function(k) unlist(lapply(Perm[[k]], Ustat, weight = FALSE)))
    Perm_M_Ustats_c <- rbind(Ustats, Perm_M_Ustats)
    CI_Ustats <- apply(Perm_M_Ustats_c, 1, function(z) {
      sum(z * sig)
    })

    pvalU <- 1 / length(CI_Ustats) * sum(CI_Ustats >= CI_Ustats[1])
  }

  return(list(pvalue = pvalU, statistic = as.numeric(CI_Ustats[1])))
}

# test statistic
Ustat <- function(freq, l1 = 0, l2 = 0, weight = TRUE) {
  # sample size
  n <- sum(freq)

  # intrinsic (empirical) dimensions of X and Y
  l <- dim(freq)
  l1_int <- l[1]
  l2_int <- l[2]

  if(n >= 4) {
    # row sum / column sum
    freq_x <- rowSums(freq)
    freq_y <- colSums(freq)

    if(weight == TRUE) {

      # when l1 and l2 are not specified
      if(l1 == 0 | l2 == 0) {
        l1 <- l1_int
        l2 <- l2_int
      }

      # weights (1/eta_q * 1/upsilon_r) based on the true dimensions
      eta <- 1 + l1 * freq_x / n
      upsilon <- 1 + l2 * freq_y / n
      W <- 1 / (eta %*% t(upsilon))

      A1 <- sum((freq^2 - freq) * W)
      A2 <- sum((freq_x^2 -freq_x) / eta) * sum((freq_y^2 - freq_y) / upsilon)
      A3 <- sum(freq * (freq_x %*% t(freq_y) - freq_x %*% t(rep(1,l2_int)) - rep(1,l1_int) %*% t(freq_y) + matrix(1,l1_int,l2_int)) * W)

      return(1/n/(n-3) * (A1 + 1/(n-1)/(n-2) * A2 - 2/(n-2) * A3))

    } else {
      # identity (uniform) weights
      A1 <- sum((freq^2 - freq))
      A2 <- sum(freq_x^2 -freq_x) * sum(freq_y^2 - freq_y)
      A3 <- sum(freq * (freq_x %*% t(freq_y) - freq_x %*% t(rep(1,l2_int)) - rep(1,l1_int) %*% t(freq_y) + matrix(1,l1_int,l2_int)))

      return(1/n/(n-3) * (A1 + 1/(n-1)/(n-2) * A2 - 2/(n-2) * A3))
    }
  } else {
    # When the sample size is less than 4, return zero
    return(0)
  }
}
