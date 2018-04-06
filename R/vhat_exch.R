#' Variance computation for linear regression of network response
#'
#' Stand-alone estimation of exchangeable variance matrix based on residuals and design matrix
#'
#' @export
#' @keywords external
#'
#' @param e Vector of residuals, of length \eqn{d}. Column-wise unfolding of adjacency matrix without diagonal entries (self-loops). 
#' @param X Matrix of covariates from regression, must have \eqn{d} rows. 
#' @param directed Optional logical indicator of whether input data is for a directed network, default is \code{TRUE}. Undirected data format is lower triangle of adjacencey matrix. 
#' @param nodes Optional \eqn{d \times 2} matrix indicating the (directed) relation pairs to which each entry in \eqn{e} and each row in \eqn{X} corresponds. If not input, complete network observation is assumed and the size \eqn{d} and \code{directed} must correspond to an appropriate network of size \eqn{n}. 
#'
#'
#' @details This function takes \eqn{X} and \eqn{e} values computes the variance-covariante matrix of \eqn{\hat{\beta}} that resulted in the residuals \eqn{e = Y - X \hat{\beta}} assuming that the errors are exchangeable, as based on Marrs (2017).
#'
#' @return 
#' A an object of class \code{vhat} containing summary information:
#' \item{vhat}{Estimated variance-covariance matrix of cofficient estimates \eqn{\hat{\beta}}.}
#' \item{phi}{Vector of variance-covariance parameter estimates.}
#' \item{corrected}{Logical of whether variance-covariance matrix was corrected from negative definite to positive semi-definite.}
#'
#' @seealso \code{\link{lmnet}}, \code{\link{inputs_lmnet}} 
#' 
#' @references Marrs, F. W., McCormick, T. H., & Fosdick, B. K. (2017). Standard errors for regression on relational data with exchangeable errors. arXiv preprint arXiv:1701.05530.
#'
#' @examples
#' set.seed(1)
#' n <- 10
#' d <- n*(n-1)
#' X <- cbind(1, rnorm(d), sample(c(0,1), d, replace=TRUE))
#' e <- rnorm(d)
#' vhat_exch(e,X)
#' 
vhat_exch <- function(e, X, directed=TRUE, nodes=NULL)
{
  subtract <- NULL
  
  #### Data check and formats
  temp <- node_preprocess(e, X, directed, nodes)
  e <- temp$Y  ;  X <- temp$X  ;  row_list <- temp$row_list  ;  n <- temp$n  
  rm(temp)
  ####
  
  ### Calculate variance estimate
  XX <- solve(crossprod(X))
  meat <- meat.E.row(row_list, X, e)
  phi <- meat$phi
  v0 <- make.positive.var( XX %*% meat$M %*% XX )
  
  vhat <- v0$V
  colnames(vhat) <- rownames(vhat) <- colnames(X)
  flag <- as.logical(v0$flag == 1)
  
  listout <- list(vhat=vhat, phi=phi, corrected=flag)
  class(listout) <- "vexch"
  return(listout)
}


#' Print generic for vexch object
#' @keywords internal
print.vexch <- function(x)
{
  print(signif(x$vhat, 4))
}

#' Summary generic for vexch object
#' @keywords internal
summary.vexch <- function(x)
{
  out <- x
  class(out) <- "summary.vexch"
  out
}

#' Print generic for summary.vexch object
#' @keywords internal
print.summary.vexch <- function(x)
{
  cat("\nEstimated variance-covariance matrix of coefficients:\n")
  print(signif(x$vhat, 4))
  cat("\nWas above matrix corrected to be nonnegative definite?\n")
  print(x$corrected)
  cat("\nEstimated parameters in variance-covariance matrix of the errors:\n")
  print(signif(x$phi, 4))
}

