#' Variance computation for linear regression of network response
#'
#' Stand-alone estimation of exchangeable variance matrix based on residuals and design matrix.
#'
#' @export
#' @keywords external
#'
#' @param e Optional vector of residuals, of length \eqn{d}. Column-wise unfolding of adjacency matrix without diagonal entries (self-loops). 
#' @param X Optional matrix of covariates from regression, must have \eqn{d} rows. 
#' @param directed Optional logical indicator of whether input data is for a directed network, default is \code{TRUE}. Undirected data format is lower triangle of adjacencey matrix. 
#' @param nodes Optional \eqn{d \times 2} matrix indicating the (directed) relation pairs to which each entry in \eqn{e} and each row in \eqn{X} corresponds. If not input, complete network observation is assumed and the size \eqn{d} and \code{directed} must correspond to an appropriate network of size \eqn{n}. 
#' @param type Optional string indicating whether the `meat' in the sandwich variance estimator is estimated using exchangeable theory (see Marrs et. al. (2017)) or using dyadic clustering (Fafchamps and Gubert (2007)). 
#' @param tmax Optional numeric of third dimension of relational data array, default is \code{1}, i.e. a relational matrix. Currently only accepts \code{tmax = 1}.
#' @param fit Optional fitted model object. One of either \code{fit} or the pair \code{(e, X)} must be specified. Defaults to \code{fit} if both entered. Designed around `lmnet' class but may work for others, such as `lm'
#'
#' @details This function takes \eqn{X} and \eqn{e} values computes the variance-covariance matrix of \eqn{\hat{\beta}} that resulted in the residuals \eqn{e = Y - X \hat{\beta}} assuming that the errors are exchangeable, as based on Marrs et. al. (2017) when \code{type = "exchangeable"}. When \code{type = "dyadic clustering"}, the theory from Fafchamps and Gubert (2007) is implemented. 
#'
#' @return 
#' A an object of class \code{vhat} containing summary information:
#' \item{vhat}{Estimated variance-covariance matrix of cofficient estimates \eqn{\hat{\beta}}.}
#' \item{phi}{Vector of variance-covariance parameter estimates.}
#' \item{corrected}{Logical of whether variance-covariance matrix was corrected from negative definite to positive semi-definite.}
#' \item{type}{See inputs.}
#' \item{tmax}{See inputs.}
#' 
#' @seealso \code{\link{lmnet}}, \code{\link{inputs_lmnet}} 
#' 
#' @aliases vhat_exch
#' 
#' @references Marrs, F. W., McCormick, T. H., & Fosdick, B. K. (2017). Standard errors for regression on relational data with exchangeable errors. arXiv preprint arXiv:1701.05530.
#' @references Fafchamps, M., & Gubert, F. (2007). Risk sharing and network formation. American Economic Review, 97(2), 75-79.
#'
#' @examples
#' set.seed(1)
#' n <- 10
#' d <- n*(n-1)
#' X <- cbind(1, rnorm(d), sample(c(0,1), d, replace=TRUE))
#' e <- rnorm(d)
#' vnet(e=e,X=X)
#' 
vnet <- function(e=NULL, X=NULL, directed=TRUE, nodes=NULL, type="exchangeable", tmax=1, fit=NULL)
{
  #### Data check and formats
  subtract <- NULL
  directed <- as.logical(directed)
  type <- as.character(type)
  tmax <- as.numeric(tmax)
  if(tmax > 1){stop("tmax > 1 not currently supported")}
  
  if(!is.null(nodes)){
    nodes <- as.matrix(nodes)
  }
  
  if(!is.null(fit)){
    X <- as.matrix(model.matrix(fit))
    e <- as.numeric(resid(fit))
    if(class(fit) == "lmnet"){
      bread <- fit$bread
      W <- fit$W
    } else {
      warning("fit is not of class lmnet")
      bread <- vcov(fit)
      W <- diag(nrow(X))
    }
  } else {
    W <- diag(nrow(X))
  }
  
  if(is.null(fit) & is.null(e)){stop("One of fit or e must be input")}
  if(is.null(fit) & is.null(X)){stop("One of fit or X must be input")}
  
  temp <- node_preprocess(as.numeric(e), as.matrix(X), as.logical(directed), nodes)
  e <- temp$Y  ;  X <- temp$X  ;  row_list <- temp$row_list  ;  n <- temp$n  
  rm(temp)
  gc()
  
  if(is.null(fit)){
    bread <- solve(crossprod(X))
  }
  ####
  
  ### Calculate variance estimate
  if(strtrim(type, 1) == "E" | strtrim(type, 1) == "e"){
    type <- "exchangeable"
    meat <- meat.E.row(row_list, as.matrix( W %*% X ), as.numeric(e) )
    phi <- meat$phi
  } else if (strtrim(type, 1) == "D" | strtrim(type, 1) == "d") {
    type <- "dyadic clustering"
    meat <- meat.DC.row(row_list, as.matrix( W %*% X ), as.numeric(e) )
    phi <- NA
  } else {
    stop("type must be one of exchangeable or dyadic clustering")
  }

  v0 <- make.positive.var( bread %*% meat$M %*% bread )
  
  vhat <- v0$V
  colnames(vhat) <- rownames(vhat) <- colnames(X)
  flag <- as.logical(v0$flag == 1)
  
  listout <- list(vhat=vhat, phi=phi, corrected=flag, type=type)
  class(listout) <- "vnet"
  return(listout)
}


#' Print S3 generic for vnet object
#' @export
#' @param x vnet object
#' @param ... ignored
print.vnet <- function(x, ...)
{
  print(signif(x$vhat, 4))
}

#' Summary S3 generic for vnet object
#' @export
#' @param object vnet object
#' @param ... ignored
summary.vnet <- function(object, ...)
{
  out <- object
  class(out) <- "summary.vnet"
  out
}

#' Print S3 generic for summary.vnet object
#' @export
#' @param x summary.vnet object
#' @param ... ignored
print.summary.vnet <- function(x, ...)
{
  if(nchar(as.character(x$type))>1){
    cat("Type of `meat' in sandwich variance estimate:\n")
    print(x$type)
  }
  cat("\nEstimated variance-covariance matrix of coefficients:\n")
  print(signif(x$vhat, 4))
  cat("\nWas above matrix corrected to be nonnegative definite?\n")
  print(x$corrected)
  cat("\nEstimated parameters in variance-covariance matrix of the errors:\n")
  print(signif(x$phi, 4))
}

