#' Input preprocessing 
#'
#' Prepare covariates and optional response in adjacency matrix form. If undirected, the values are drawn from the lower triangle of the adjacency matrices.
#'
#' @export
#' @keywords external
#'
#' @param Xlist List of \eqn{n \times n \times tmax} matrices, possibly containing response matrix labeled `Y'. Diagonals (self-loops) are ignored.
#' @param Y Optional \eqn{n \times n \times tmax} response matrix. NAs in this matrix will be automatically removed. Diagonals (self-loops) are ignored.
#' @param directed Optional logical indicator of whether input data is for a directed network, default is \code{TRUE}. Undirected data format is lower triangle of adjacencey matrix. 
#' @param add_intercept Optional logical indicator of whether intercept should be added to X, default is \code{TRUE}.
#' @param time_intercept Optional logical indicator of whether separate intercept should be added to X for each observation of the relational matrix, default is \code{FALSE}.
#'
#' @details This function takes a list of network covariates (in adjacency matrix form) and prepares them for the regression code \code{lmnet}. Accomodates 3-dimensional relational arrays with \code{tmax} repeated observations of the network (over time or context). Typical network data with a single observation may be input as matrices, i.e. \code{tmax = 1}.
#'
#' @return 
#' A list of:
#' \item{Y}{Vector of responses (column-wise vectorization order) of appropriate length.}
#' \item{X}{Matrix of covariates (column-wise vectorization order) of appropriate size.}
#' \item{nodes}{2-column matrix (or 3-column for repeated observations) indicating directed relation pairs to which each entry in \eqn{Y} and each row in \eqn{X} corresponds.}
#'
#' @seealso \code{\link{lmnet}}, \code{\link{vhat_exch}} 
#'
#' @examples
#' # tmax = 1
#' set.seed(1)
#' n <- 10
#' Xlist <- list(matrix(rnorm(n^2),n,n), matrix(sample(c(0,1), n^2, replace=TRUE),n,n))
#' Xlist$Y <- matrix(rnorm(n^2), n, n)
#' Xlist$Y[1:5] <- NA
#' r <- inputs_lmnet(Xlist)
#' r
#' lmnet(r$Y,r$X,nodes=r$nodes)
#' 
#' # tmax = 4
#' set.seed(1)
#' n <- 10
#' tmax <- 4
#' X1 <- array(rnorm(n^2*tmax),c(n,n,tmax))
#' X2 <- array(sample(c(0,1), n^2*tmax, replace=TRUE), c(n,n,tmax))
#' Xlist <- list(X1, X2)
#' Xlist$Y <- array(rnorm(n^2)*tmax, c(n, n, tmax))
#' Xlist$Y[1:5] <- NA
#' r <- inputs_lmnet(Xlist)
#' head(r$nodes)
#' 
inputs_lmnet <- function(Xlist, Y=NULL, directed=TRUE, add_intercept=TRUE, time_intercept=FALSE)
{
  add_intercept <- as.logical(add_intercept)
  directed <- as.logical(directed)
  
  if(class(Xlist) != "list"){stop("Xlist must be a list")}
  yin <- !is.null(Y)
  ylist <- ("Y" %in% names(Xlist)) | ("y" %in% names(Xlist))
  if(yin & ylist){warning("Y found in Xlist and Y input, defaulting to the first one in Xlist")}
  if(yin &! ylist ){
    # n <- nrow(yin)
    n <- dim(Y)[1]
  } 
  if(ylist){
    iy <- which(names(Xlist) %in% c("Y", "y"))
    if(length(iy) > 1){warning("Multiple entries in Xlist are named 'Y', taking the first one by default")}
    Y <- Xlist[[ iy[1] ]]
    Xlist <- Xlist[-iy]
    n <- dim(Y)[1]
  }

  
  # if( !all(c(unlist(sapply(Xlist, function(z) dim(as.matrix(z))))) == n) ){stop("All entries in Xlist must be square of same size")}
  if( !all(c(unlist(sapply(Xlist, function(z) dim(z)[1:2]) )) == n) ){stop("All entries in Xlist must be of same size")}
  
  n <- dim(Xlist[[1]])[1]
  d <- (1 + directed)/2*n*(n-1)
  
  tcheck <- length(dim(Y)) > 2
  if(tcheck){
    if(length(dim(Y)) > 3){stop("Y must be a 2- or 3-mode array")}
    tmax <- dim(Y)[3]
    Yout <- rep(NA, d*tmax)
    for(t in 1:tmax){
      Yout[1:d + d*(t-1)] <- vec.net(Y[,,t], directed)
    }
  } else {
    if(length(dim(Y)) < 2){stop("Y must be a 2- or 3-mode array")}
    tmax <- 1
    Yout <- vec.net(Y, directed)
  }
  
  
  # Build X matrix
  p <- length(Xlist)  # number of columns in X
  if(tcheck){
    Xmat <- matrix(NA, d*tmax, p)
    for(t in 1:tmax){
      for(j in 1:p){
        Xmat[1:d + d*(t-1),j] <- vec.net(Xlist[[j]][,,t], directed)
      }
    }
    if(add_intercept){
      if(time_intercept){   # separate intercept for each time
        Xmat <- cbind(matrix(0, nrow(Xmat), tmax), Xmat)
        for(t in 1:tmax){
          Xmat[1:d + d*(t-1), t] <- 1
        }
      } else {
        Xmat <- cbind(1, Xmat)
      }
    }
  } else {
    Xmat <- matrix(NA, d, p)
    for(j in 1:p){
      Xmat[,j] <- vec.net(Xlist[[j]], directed)
    }
    if(add_intercept){Xmat <- cbind(1, Xmat)}
  }

  
  # Make node set
  if(tcheck){
    nodes <- node_gen_time(n, directed, tmax)
  } else{
    nodes <- node.gen(n, directed)
  }
  
  # remove any NAs in Y
  remove <- which(is.na(Yout))
  if(length(remove) > 0){
    Yout <- Yout[-remove]
    Xmat <- Xmat[-remove,]
    nodes <- nodes[-remove,]
  }
  
  return(list(Y=as.numeric(Yout), X=as.matrix(Xmat), nodes=as.matrix(nodes)))
}

