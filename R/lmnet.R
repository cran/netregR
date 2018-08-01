#' Linear regression for network response
#'
#' This function takes \eqn{X} and \eqn{Y} values and fits the multiple linear regression \eqn{Y = X \beta + \epsilon} and returns standard errors.
#'
#' @export
#' @keywords external
#'
#' @param Y Vector of relations to be regress, of length \eqn{d}. Column-wise vectorization of adjacency matrix without diagonal entries (self-loops).
#' @param X Matrix of covariates to be regressed upon, including intercept if intercept is desired, must have \eqn{d} rows. Ordering of rows should match \code{Y} and optional input \code{nodes}.
#' @param directed Optional logical indicator of whether input data is for a directed network, default is \code{TRUE}. Undirected data format is lower triangle of adjacencey matrix. 
#' @param tmax Optional numeric of third dimension of relational data array, default is \code{1}, i.e. a relational matrix. 
#' @param nodes Optional \eqn{d \times 2} matrix indicating the (directed) relation pairs to which each entry in \eqn{Y} and each row in \eqn{X} corresponds. If not input, complete network observation with column-wise vectorization of adjacency matrix without diagonal entries (self-loops) is assumed. The size \eqn{d} and \code{directed} must correspond to an appropriate network of size \eqn{n}. 
#' @param reweight Optional logical indicator of whether iteratively reweighted least squares should be used to compute estimate of \eqn{\beta}. Default is \code{FALSE}.
#' @param type Optional character specifying degree of exchangeability of third dimension of array (when present, i.e. in temporal relational arrays). Default is \code{exchangeable}, and the remaining option is \code{independent}. Truncated inputs are accepted. See details below.
#' @param tol Optional numeric, tolerance of stopping criteria of iteratively reweighted least squares estimate of \eqn{\beta}. Default is \code{tol=1e-6}.
#' @param maxit Optional numeric, maximum number of iterations for iteratively reweighted least squares estimate of \eqn{\beta}. Default is \code{maxit=1e4}.
#' @param ndstop Optional logical indicator of whether negative definite weighting matrix in iteratively reweighted least squares should stop the descent. Default is \code{TRUE}.
#' @param verbose Optional logical indicator of whether information from iteratively reweighted least squares estimate of \eqn{\beta} should be printed. Default is \code{FALSE}.
#' 
#' 
#' @details 
#' This function takes \eqn{X} and \eqn{Y} values and fits the multiple linear regression \eqn{Y = X \beta + \epsilon} by ordinary least squares or iteratively reweighted least squares as indicated by the input. The covariance structure is exchangeable from that of Marrs et. al. (2017). The standard errors and test statistics are based on the same paper.    
#' 
#' The three dimensional relational array case, i.e. temporal relational data, requires a specification of the type of exchangeability in this third dimension. We may assume that different time periods are independent. On the other hand, we might assume each repeated observation is exchangeable (for example decomposing trade networks into sectors of trade: goods vs. services). See Figure 6a of Marrs et. al. (2017) for the exchangeable case and the surrounding discussion for the independent case.
#' 
#' 
#' @return 
#' \item{fit}{An \code{lmnet} object containing summary information.}
#'
#' @seealso \code{\link{vhat_exch}}, \code{\link{inputs_lmnet}} 
#' 
#' @references Marrs, F. W., Fosdick, B. K., & McCormick, T. H., (2017). Standard errors for regression on relational data with exchangeable errors. arXiv preprint arXiv:1701.05530.
#'
#' @examples
#' set.seed(1)
#' n <- 10
#' d <- n*(n-1)
#' X <- cbind(1, rnorm(d), sample(c(0,1), d, replace=TRUE))
#' betatrue <- rep(1,3)
#' Y <- X %*% betatrue + rnorm(d)
#' fit <- lmnet(Y,X)
#' fit
#' fit2 <- lmnet(Y,X,reweight=TRUE)
#' fit2
#' 
lmnet <- function(Y, X, directed=TRUE, tmax=1, nodes=NULL, reweight=FALSE, type="exchangeable", tol=1e-6, maxit=1e4, ndstop=TRUE, verbose=FALSE)
{
  
  #### Data check and formats
  tmax <- as.numeric(tmax)
  directed <- as.logical(directed)
  
  if(tmax == 1){
    temp <- node_preprocess(Y,X,directed,nodes)
  } else {
    temp <- node_preprocess_time(Y,X,directed,nodes,tmax,type,subtract=NULL)
  }
  Y <- temp$Y  ;  X <- temp$X  ;  missing <- temp$missing  ;  row_list <- temp$row_list  ;  dyads <- temp$dyads  ;  n <- temp$n  ;  type <- temp$type
  rm(temp)
  reweight <- as.logical(reweight)
  tol <- as.numeric(tol)
  maxit <- as.numeric(maxit)
  verbose <- as.logical(verbose)
  if(sum(is.na(X))!=0){warning("NAs in X; no action taken.")}
  ####
  
  
  if(missing & tmax > 1){
    stop("Missing data not yet implemented for temporal data")
  }
  
  
  #### Initial fits
  fit <- lm(Y ~ X - 1)
  beta_ols <- coef(fit)
  
  X <- model.matrix(fit)
  p <- ncol(X)
  XX <- solve(crossprod(X))
  
  e <- Y - X %*% beta_ols
  meat <- meat.E.row(row_list, X, e)
  phi_ols <- meat$phi
  v0 <- make.positive.var( XX %*% meat$M %*% XX )
  
  Vhat_ols <- v0$V
  Vflag_ols <- v0$flag
  ####
  
  
  
  #### GEE fit if desired
  if(reweight){  # compute GEE fit
    if(tmax ==1){
      fit_weighted <- GEE.est(row_list, Y, X, n, directed, tol.in=tol, beta_start=beta_ols, missing=missing, dyads=dyads, ndstop=ndstop, verbose=verbose)
    } else if (tmax > 1){
      fit_weighted <- GEE_est_time(Y, X, n, tmax,  directed, type, write_dir=NULL, missing=missing, tol.in=tol, maxit=maxit, verbose=verbose) 
    } else {
      stop("tmax must be a positive integer")
    }
    
    beta_weighted <- fit_weighted$beta
    v0 <- make.positive.var( solve( fit_weighted$bread ) )
    e <- fit_weighted$residuals
      
    Vhat_weighted <- v0$V
    Vflag_weighted <- v0$flag
    phiout <- as.numeric(fit_weighted$phi)
    nit <- as.numeric(fit_weighted$nit)
    conv <- as.logical(fit_weighted$convergence)
    
    if(!conv){warning("Iteratively reweighted least squares procedure stopped based on maximum number of iterations (did not converge)\n")}
    
    betaout <- beta_weighted
    Vout <- Vhat_weighted
    flagout <- as.logical(Vflag_weighted == 1)
    
    bread = Vout
    W <- fit_weighted$W
  
  } else {
    beta_weighted <- Vhat_weighted <- Vflag_weighted  <- nit <- tol <- conv <- NA
    
    betaout <- beta_ols
    Vout <- Vhat_ols
    flagout <- as.logical(Vflag_ols == 1)
    phiout=phi_ols
    W <- diag(nrow(X))
    
    bread <- XX 

  }
  ####
  
  df <- nrow(X) - length(betaout) - 1
  if(length(betaout) == ncol(X)){names(betaout) <- colnames(X)}
  
  fitout <- list(call=match.call(), coefficients=betaout, residuals=e, vcov=Vout, fitted.values=X %*% betaout,
                 df=df, sigma=sqrt(phiout[1]), 
                 reweight=reweight,
                 corrected=flagout, phi_hat=phiout, nit=nit, converged=conv, 
                 X=X, nodes=nodes,
                 bread=bread, W=W,
                 tmax=tmax, type=type, ndstop=ndstop)
  class(fitout) <- "lmnet"
  return(fitout)
}


#' Print S3 generic for class lmnet
#' @export
#' @param x lmnet object
#' @param ... ignored
print.lmnet <- function(x, ...)
{
  cat("\nCall: \n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
  cat("\n")
}


#' Coef S3 generic for class lmnet
#' @export
#' @param object lmnet object
#' @param ... ignored
coef.lmnet <- function(object, ...)
{
  object$coefficients
}


#' vcov S3 generic for class lmnet
#' @export
#' @param object lmnet object
#' @param ... ignored
vcov.lmnet <- function(object, ...)
{
  object$vcov
}


#' Summary S3 generic for class lmnet
#' @export
#' @param object lmnet object
#' @param ... ignored
summary.lmnet <- function(object, ...)
{
  x <- object
  out <- matrix(coef(x), ncol=1)
  out <- cbind(out, sqrt(diag(vcov(x))))
  out <- cbind(out, out[,1] / out[,2])
  out <- cbind(out, 1-pt(abs(out[,3]), df=x$df))
  rownames(out) <- names(coef(x))
  colnames(out) <- c("Estimate", "Std. Error", "t value", "Pr(|t| > 0)")             
  
  listout <- list(coefficients=out, call=x$call)
  class(listout) <- "summary.lmnet"
  return(listout)
}


#' Print S3 generic for class summary.lmnet
#' @export
#' @param x summary.lmnet object
#' @param ... ignored
print.summary.lmnet <- function(x, ...)
{
  cat("\nCall:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  printCoefmat(x$coefficients)
}


#' Plot S3 generic for class lmnet
#' @export
#' @param x lmnet object
#' @param ... ignored
plot.lmnet <- function(x, ...)
{
  hist(scale(resid(x)), freq=F, xlab="standardized residuals", main="")
  
  plot(fitted.values(x), scale(resid(x)), xlab="fitted values", ylab="standardized residuals", main="")
  
  qqnorm(scale(resid(x)), main="Normal Q-Q Plot for residuals")
  abline(0,1, col="red")
}


#' model.matrix S3 generic for class lmnet
#' @export
#' @param object lmnet object
#' @param ... ignored
model.matrix.lmnet <- function(object, ...)
{
  object$X
}
