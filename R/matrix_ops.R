# Matrix operation functions


#' Generate positive definite phi set
#'
#' @export
#' @keywords external
#' 
#' @param n Number of actors in the network, scalar numeric.
#' @param seed Optional numeric seed to set, default is \code{NULL}.  
#' @param phi6 Optional logical indicator of whether sixth parameter \eqn{\phi_6} should be considered nonzero. Default is \code{FALSE}.
#' @param directed Optional logical indicator of whether input data is for a directed network, default is \code{TRUE}. Undirected data format is lower triangle of adjacencey matrix. 
#'
#'
#' @details This function generates a set of 5 (or 6, as appropriate) parameters that corresponds to positive definite exchangeable covariance matrix for a network of size \code{n}. See Marrs et. al. (2017). 
#'
#' @return 
#' \item{phi}{Vector of parameters.}
#' 
#' @seealso \code{\link{build_exchangeable_matrix}}, \code{\link{invert_exchangeable_matrix}} 
#' 
#' @references Marrs, F. W., Fosdick, B. K., & McCormick, T. H.,  (2017). Standard errors for regression on relational data with exchangeable errors. arXiv preprint arXiv:1701.05530.
#'
#' @examples
#' rphi(10, seed=1)
#' 
rphi <- function(n, seed=NULL, phi6=FALSE, directed=TRUE)
{
  if(is.numeric(seed)){ set.seed(seed)}
  n <- as.numeric(n)[1]
  phi6 <- as.logical(phi6)
  directed <- as.logical(directed)
  
  if(directed){
    mineig <- -1
    count <- 0
    while(mineig <= 0){
      
      count <- count + 1
      phi <- phi0 <- c(1, runif(5,-1,1))
      if(!phi6){
        phi[6] <- 0
      }
      # C <- build_phi_matrix(n,phi)
      # C0 <- build_phi_matrix(n,phi0)
      
      # mineig <- min(eigen(C)$values)
      #min(min(eigen(C)$values), min(eigen(C0)$values))
      mineig <- min( eigen_exch(n, phi, TRUE, FALSE)$uniquevals )
      
    }
    # S <- Sigma.ind(n,T)
    # if(phi6){
    #   Omega <- Reduce("+", lapply(1:6, function(z) phi[z]*S[[z]]))
    # } else {
    #   Omega <- Reduce("+", lapply(1:5, function(z) phi[z]*S[[z]]))
    # }
    # # Omega0 <- Reduce("+", lapply(1:6, function(z) phi0[z]*S[[z]]))
    # if( min(eigen(Omega)$values) <= 0){ 
    #   warning("C eigens did not carry over to Omega")
    #   if(is.numeric(seed)){seed <-  seed + 1} else {seed <- sample(1:1000, 1)}
    #   phi <- rphi(n, seed, directed)
    # }
    
  } else {   # undirected
    phi <- c(1,0,0)
    phi[2] <- runif(1, -.5/(n-2), .5)
  }
  
  return(phi)
}



#' Build an exchangeable matrix of sparseMatrix class
#'
#' @export
#' @keywords external
#' 
#' @param n Number of actors in the network, scalar numeric.
#' @param phi Appropriate-length vector of parameters, must be length 5 or 6 for directed=\code{TRUE} or length 2 or 3 for directed=\code{FALSE}.
#' @param directed Optional logical indicator of whether input data is for a directed network, default is \code{TRUE}. Undirected data format is lower triangle of adjacencey matrix. 
#' @param dyads Optional numeric vector of dyads to subset the matrix to. 
#'
#' @details This function builds a covariance matrix in the exchangeable class from the vector of parameters input. See Marrs et.al. (2017).
#'
#' @return 
#' \item{out}{Exchangeable matrix.}
#' 
#' @seealso \code{\link{rphi}}, \code{\link{invert_exchangeable_matrix}} 
#' 
#' @references Marrs, F. W., Fosdick, B. K., & McCormick, T. H., (2017). Standard errors for regression on relational data with exchangeable errors. arXiv preprint arXiv:1701.05530.
#'
#' @examples
#' n <- 5
#' build_exchangeable_matrix(n, rphi(n, seed=1))
#' 
build_exchangeable_matrix <- function(n, phi, directed=TRUE, dyads=NULL)
{
  n <- as.numeric(n)[1]
  phi <- as.numeric(phi)
  directed <- as.logical(directed)
  
  if( !(length(phi) %in% c(5,6)) & directed){stop("phi must be length 5 or 6 for directed=TRUE")}
  if( !(length(phi) %in% c(2,3)) & !directed){stop("phi must be length 2 or 3 for directed=FALSE")}
  
  S <- Sigma.ind(n, directed)
  if(directed){
    L <- max(5, length(phi))
    if(L == 6){
      L <- 5 + 1*(phi[6]!=0)
    }
  } else {
    L <- max(2, length(phi))
    if(L == 3){
      L <- 2 + 1*(phi[3]!=0)
    }
  }
  d <- nrow(S[[1]])
  out <- Reduce("+", lapply(1:L, function(z) phi[z]*S[[z]]))
  if(!is.null(dyads)){   # subset to dyads of interest (and ordering) as appropriate
    dyads <- as.numeric(dyads)
    out <- out[dyads, dyads]
  }
  
  return(out)
}


#' Invert an exchangeable matrix 
#'
#' @export
#' @keywords external
#' 
#' @param n Number of actors in the network, scalar numeric.
#' @param phi Appropriate-length vector of parameters, must be length 5 or 6 for directed=\code{TRUE} or length 2 or 3 for directed=\code{FALSE}.
#' @param directed Optional logical indicator of whether input data is for a directed network, default is \code{TRUE}. Undirected data format is lower triangle of adjacencey matrix. 
#'
#' @details This function inverts a covariance matrix of the exchangeable class in a manner much faster than the direct inverse, and the computational cost does not scale with n. See Marrs et. al. (2017). This approach will only work for complete networks. 
#'
#' @return 
#' \item{out}{Parameters of inverted matrix of exchangeable class.}
#' 
#' @seealso \code{\link{rphi}}, \code{\link{build_exchangeable_matrix}}
#' 
#' @references Marrs, F. W., Fosdick, B. K., & McCormick, T. H., (2017). Standard errors for regression on relational data with exchangeable errors. arXiv preprint arXiv:1701.05530.
#'
#' @examples
#' n <- 10
#' phi <- rphi(n, seed=1)
#' p <- invert_exchangeable_matrix(n, phi)
#' I1 <- build_exchangeable_matrix(n, phi) %*% build_exchangeable_matrix(n, p)
#' range(I1 -  diag(n*(n-1)))   # it works
#'
invert_exchangeable_matrix <- function(n, phi, directed=TRUE)
{
  n <- as.numeric(n)[1]
  phi <- as.numeric(phi)
  directed <- as.logical(directed)
  
  if( !(length(phi) %in% c(5,6)) & directed){stop("phi must be length 5 or 6 for directed=TRUE")}
  if( !(length(phi) %in% c(2,3)) & !directed){stop("phi must be length 2 or 3 for directed=FALSE")}
  
  C <- build_phi_matrix(n, phi, directed)
  
  return(solve(C, c(1, rep(0, nrow(C) - 1))))
}
