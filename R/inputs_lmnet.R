#' Input preprocessing 
#'
#' Prepare covariates and optional response in adjacency matrix form. If undirected, the values are drawn from the lower triangle of the adjacency matrices.
#'
#' @export
#' @keywords external
#'
#' @param Xlist List of \eqn{n \times n} matrices, possibly containing response matrix labeled `Y'. Diagonals (self-loops) are ignored.
#' @param Y Optional \eqn{n \times n} response matrix. NAs in this matrix will be automatically removed. Diagonals (self-loops) are ignored.
#' @param directed Optional logical indicator of whether input data is for a directed network, default is \code{TRUE}. Undirected data format is lower triangle of adjacencey matrix. 
#' @param add_intercept Optional logical indicator of whether intercept should be added to X, default is \code{TRUE}.
#'
#' @details This function takes a list of network covariates (in adjacency matrix form) and prepares them for the regression code \code{lmnet}.
#'
#' @return 
#' A list of:
#' \item{Y}{Vector of responses (column-wise vectorization order) of appropriate length.}
#' \item{X}{Matrix of covariates (column-wise vectorization order) of appropriate size.}
#' \item{nodes}{2-column matrix indicating directed relation pairs to which each entry in \eqn{Y} and each row in \eqn{X} corresponds.}
#'
#' @seealso \code{\link{lmnet}}, \code{\link{vhat_exch}} 
#'
#' @examples
#' set.seed(1)
#' n <- 10
#' Xlist <- list(matrix(rnorm(n^2),n,n), matrix(sample(c(0,1), n^2, replace=TRUE),n,n))
#' Xlist$Y <- matrix(rnorm(n^2), n, n)
#' Xlist$Y[1:5] <- NA
#' r <- inputs_lmnet(Xlist)
#' r
#' lmnet(r$Y,r$X,nodes=r$nodes)
#' 
inputs_lmnet <- function(Xlist, Y=NULL, directed=TRUE,  add_intercept=TRUE )
{
  if(class(Xlist) != "list"){stop("Xlist must be a list")}
  yin <- !is.null(Y)
  ylist <- ("Y" %in% names(Xlist)) | ("y" %in% names(Xlist))
  if(yin & ylist){warning("Y found in Xlist and Y input, defaulting to the Y in Xlist")}
  if(yin &!ylist ){
    n <- nrow(yin)
  } 
  if(ylist){
    iy <- which(names(Xlist) %in% c("Y", "y"))[1]
    n <- nrow(as.matrix(Xlist[[iy]]))
  }
  if( !all(c(unlist(sapply(Xlist, function(z) dim(as.matrix(z))))) == n) ){stop("All entries in Xlist must be square of same size")}
  
  
  n <- nrow(Xlist[[1]])
  d <- (1 + directed)/2*n*(n-1)
  
  # Build Y vector
  if(ylist){
    iy <- which(names(Xlist) %in% c("Y", "y"))
    if(length(iy) > 1){warning("Multiple entries in Xlist are named 'Y', taking the first one by default")}
    Y <- as.matrix(Xlist[[ iy[1] ]])
    Xlist <- Xlist[-iy]
  }
  Yout <- vec.net(Y, directed)
  
  # Build X matrix
  p <- length(Xlist)  # number of columns in X
  Xmat <- matrix(NA, d, p)
  for(j in 1:p){
    Xmat[,j] <- vec.net(Xlist[[j]], directed)
  }
  if(add_intercept){Xmat <- cbind(1, Xmat)}
  
  # Make node set
  nodes <- node.gen(n, directed)
  
  # remove any NAs in Y
  remove <- which(is.na(Yout))
  if(length(remove) > 0){
    Yout <- Yout[-remove]
    Xmat <- Xmat[-remove,]
    nodes <- nodes[-remove,]
  }
  
  return(list(Y=as.numeric(Yout), X=as.matrix(Xmat), nodes=as.matrix(nodes)))
}

