#' Permutation-based \eqn{G}-test for conditional independence
#'
#' Given a three-way contingency table, the function \code{gtest} computes
#' the likelihood test statistic (or \code{G}-test statistic) for conditional independence and its
#' corresponding permutation p-value via Monte Carlo simulation. Unlike the asymptotic test calibrated by a limiting
#' chi-squared distribution, the resulting test based on the permutation p-value is exact in any finite sample scenarios.
#' See Section 4.3 of \insertCite{kim2022;textual}{UCI} for an explicit formula of the test statistic.
#'
#' @param dat a three-way frequency table
#' @param B the number of random permutations
#'
#' @return Returns the p-value for permutation-based \eqn{G}-test and the value of the test statistic
#' @export
#' @references{\insertAllCited{}}
#' @encoding UTF-8
#' @importFrom stats r2dtable
#' @importFrom Rdpack reprompt
#' @seealso \code{\link{UCI}} and \code{\link{gtest}}
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
#' gtest(freq, 199)
gtest <- function(dat = dat, B = 199) {
  # Calculate the permutation p-value
  # sample size
  sig <- apply(dat, 3, sum)

  # unconditional statistics
  gstats <- apply(dat, 3, Gstat)

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

  Perm_M_gstats <- sapply(1:dim(dat)[3], function(k) unlist(lapply(Perm[[k]], Gstat)))
  Perm_M_gstats_c <- rbind(gstats, Perm_M_gstats)
  CI_gstats <- rowSums(Perm_M_gstats_c)

  pvalg <- 1 / length(CI_gstats) * sum(CI_gstats >= CI_gstats[1])

  return(list(pvalue = pvalg, statistic = as.numeric(CI_gstats[1])))
}

# test statistic
Gstat <- function(freq) {

  # the dimension of the frequency table
  dim_freq <- dim(freq)
  l1 <- dim_freq[1]
  l2 <- dim_freq[2]

  # sample size
  n <- sum(freq)

  if(n >= 2) {
    # row sum / column sum
    freq_x <- rowSums(freq)
    freq_y <- colSums(freq)

    E <- freq_x %*% t(freq_y) / n
    G <- 2 * sum(log((freq / E)^freq))

    return(G)

    # when the sample size is less than two, return zero
  } else {
    return(0)
  }
}
