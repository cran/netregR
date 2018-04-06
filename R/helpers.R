#' Perform GEE estimate / IRWLS of coefficients
#'
#' @keywords internal
GEE.est <- function(row.list, Y, X, n, directed=T, beta_start=NULL, missing=F, dyads=NULL, tol.in=1e-6, maxit=1e4, verbose=F) #, node.mat=1) 
{
  # Calculate GEE estimate of regressors beta for continuous data
  # Return beta, residuals, weight matrix, # iterations, objective function Q, and actual tolerance
  # 
  
  
  t1 <- proc.time()[3]
  
  if(is.null(beta_start)){
    beta_start <- coef(lm(Y ~ X - 1))   # initialize at OLS
  } 
  beta_weighted <- beta_start
  
  count.in <- 0
  d <- n*(n-1)/2*(1 + directed)
  
  XX <- chol2inv(chol(crossprod(X)))
  e <- Y - X %*% beta_start
  Q.old <- Q.new <- sum(e^2)
  MSE <- mean(e^2)
  delta.loop <- tol.in*100
  
  while ((delta.loop) > tol.in & count.in < maxit) {
    count.in <- count.in + 1
    
    phi <- param_est(row.list, e)   
    # calculate_matrix_params(ilist, e.test, t.max, type)  # unique parameters in Var(Y)
    if(!missing){
      inv_phi <- invert_exchangeable_matrix(n, c(phi,0), directed)
      # W <- build_exchangeable_matrix(n, inv_phi, directed, dyads)
      XWX <- meatABC(row.list, inv_phi, X, X, directed)
      XWY <- meatABC(row.list, inv_phi, X, matrix(Y, ncol=1), directed)
      
    } else {  # must invert directly with missingness
      Winv <- build_exchangeable_matrix(n, phi, directed, dyads)
      W <- solve(Winv)  
      XWX <- crossprod(X,W) %*% X
      XWY <- crossprod(X,W) %*% Y
      
    }
    
    
    beta_new <- as.vector( solve(XWX, XWY) )
    e <- Y - X %*% beta_weighted
    # Q.new <- as.numeric(-determinant(as.matrix(W))$modulus + crossprod(e, W) %*% e )  # want this to be small!
    Q.new <- mean( (beta_new - beta_weighted)^2 )   # want this to be small!
    
    if(verbose & count.in%%100 == 0){
      cat('Weighted iteration: \t\t', count.in, '\n')
      cat('Change in criterion: \t', abs(delta.loop), '\n')
      cat('Elapsed time: \t\t', proc.time()[3]- t1, 'sec \n')
      print(warnings())
      cat('************************************ \n')
    }
    delta.loop <- as.numeric( -(Q.new - Q.old) )      # decrease in criterion
    if(delta.loop > 0){   # if continuing to decrease
      delta.loop <- delta.loop / Q.old * 100
    }
    
    # Update for next loop
    Q.old <- Q.new # new baseline
    if((delta.loop) > 0 & count.in < maxit){  # R finishes the loop, so only update beta if it's a good choice, ow use the last one
      # also keeps e, W, etc. consistent
      beta_weighted <- beta_new
    }
    
  }
  
  if(verbose){
    cat("\n Weighted estimation complete \n")
  }
  output <- list(beta_weighted, e, XWX, phi, count.in, as.logical(count.in < maxit))
  names(output) <- c('beta','residuals','bread', 'phi', 'nit', 'convergence')
  
  return(output)
}



#' Pre-processes data for ordering etc.
#' #'
#' @keywords internal
node_preprocess <- function(Y, X, directed, nodes, subtract=NULL)
{
  #### Preprep
  directed <- as.logical(directed)
  
  Y <- as.vector(as.numeric(Y))
  if(is.null(dim(X))){  # if X is a vector
    X <- matrix(X, ncol=1)
  } else {
    X <- as.matrix(X)  # matrix
  }
  
  
  d <- length(Y)
  if(nrow(X) != d){stop("X and Y must have the same number of observations")}
  
  remove <- which(is.na(Y))
  if(length(remove) > 0){
    Y <- Y[-remove]
    X <- X[-remove,]
    nodes <- nodes[-remove,]
  }
  ####
  
  
  if(is.null(nodes)){
    missing <- F
    cc <- 4 + 4*(directed==F)
    n <- (1+sqrt(1+cc*d))/2
    if(n != round(n)){stop("Size of Y and X must be valid for a network of some size; assuming complete network at this point")}
    
    node_list <- node.set(n, directed) 
    row_list <- lapply(node_list, function(z) cbind( dyad(z[,1], z[,2], n, directed), dyad(z[,3], z[,4], n, directed)) )
    dyads <- 1:d
    
  } else {
    nodes <- as.matrix(nodes)
    n <- max(nodes[,1:2])   # number of nodes
    if(!directed){   # if undirected, sort dyads into lower triangle
      nodes[,1:2] <-  t(apply(nodes[,1:2], 1, sort, decreasing=T))  # lower triangle dyads
      nodes <- unique(nodes)   # unique nodes
      if(nrow(nodes) != nrow(X)){stop("Number of relation pairs in nodes does no match number of rows of X. Check for repeated relation pairs in nodes input.")}
    }
    
    dyads <- dyad(nodes[,1], nodes[,2], n, directed)
    dmax <- n*(n-1)*(1 + directed)/2
    
    if(length(dyads) == dmax){  # just reorder
      missing <- F
      reorder <- order(nodes[,2], nodes[,1])   # columnwise unfolding
      
      # proceed as usual
      Y <- Y[reorder]
      X <- X[reorder,]
      node_list <- node.set(n, directed) 
      row_list <- lapply(node_list, function(z) cbind( dyad(z[,1], z[,2], n, directed), dyad(z[,3], z[,4], n, directed)) )
      
    } else {  #build row_list based on nodes input
      missing <- T
      if(is.null(subtract)){
        subtract <- length(dyads)/max(dyads) > .5   # are over half the dyads present?
      }
      row_list <- row_list_missing(nodes, dyads, directed, subtract)   # row list based on overlaps
    }
  }
  
  return(list(X=X, Y=Y, dyads=dyads, missing=missing, row_list=row_list, n=n))
}



#' Generate node sets of various overlapping dyad pairs
#'
#' @keywords internal
node.set <- function(n.tot, directed=T)
{  
  # Generate node sets of various overlapping dyad pairs
  # Return list of each set of nodes, with null for first set (diagonal)
  
  # if(!directed){stop("node.set() not yet coded for undirected")}
  
  if(directed){
    nodes.1 <- cbind(rep(1:n.tot, each=n.tot)[-seq(1,n.tot^2,n.tot+1)], rep(1:n.tot, n.tot)[-seq(1,n.tot^2,n.tot+1)],
                     rep(1:n.tot, each=n.tot)[-seq(1,n.tot^2,n.tot+1)], rep(1:n.tot, n.tot)[-seq(1,n.tot^2,n.tot+1)])
    
    nodes.2 <- nodes.3 <- nodes.4 <- nodes.5 <-  c()
    
    for (i in 1:n.tot){ 
      
      # ij,ji
      if (i<n.tot){
        c1 <- rep(i,(n.tot-i))
        #       c2 <- (1:n.tot)[-i]
        c2 <- ((i+1):n.tot)
        
        nodes.2 <- rbind(nodes.2,cbind(c1,c2,c2,c1))
      }
      
      # ij,il  ;  ij,kj
      c1 <- rep(i,(n.tot - 1)*(n.tot-2)/2)
      c2 <- rep( (1:(n.tot-1))[-i], times=(n.tot-2):(1 - 1*(i==n.tot) ) )
      c3 <- c1
      c4.mat <- outer(1:n.tot,rep(1,n.tot))[-i,-i]
      c4 <- c4.mat[lower.tri(diag(n.tot-1))]
      
      nodes.3 <- rbind(nodes.3, cbind(c1,c2,c3,c4))
      nodes.4 <- rbind(nodes.4, cbind(c2,c1,c4,c3))
      
      # ij,jl and ij,ki
      nodes.5 <- rbind(nodes.5, cbind(c1,c2,c4,c3), cbind(c2,c1,c3,c4))   
    }
    return(list(n1=nodes.1,n2=nodes.2,n3=nodes.3,n4=nodes.4,n5=nodes.5))
    
  } else { # undirected
    
    node.list <- node.set(n.tot, directed=T)
    node.list <- lapply(node.list, function(z) z[z[,2] < z[,1] & z[,4] < z[,3], ])   # keep lower tri only
    node.list[[2]] <- unique( rbind(node.list[[2]], node.list[[3]], node.list[[4]], node.list[[5]]) )
    
    node.list <- node.list[1:2]
    return(node.list)
  }
  
}


#' Generate row list based on nodes input with missingness
#'
#' @keywords internal
row_list_missing <- function(nodes, dyads, directed, subtract)
{
  n <- max(nodes[,1:2])  # number of nodes
  dmax <- n*(n-1)*(1 + directed)/2   # total number of dyads for complete data
  
  if(subtract){   # if build row list by removing dyads
    node_list <- node.set(n, directed)
    dyad_list <- lapply(node_list, function(z) cbind( dyad(z[,1], z[,2], n, directed), dyad(z[,3], z[,4], n, directed)) )
    
    dremove <- setdiff(1:dmax, unique(c(dyads)))
    dout <- lapply(dyad_list, function(z) z[-unique(which( apply(z, 1, function(x) any(x %in% dremove)), arr.ind=T)), ])
    dout1 <- lapply(dout, function(z) cbind(match(matrix(z, ncol=2)[,1], dyads),  match(matrix(z, ncol=2)[,2], dyads) ))
    
  } else {
    
    if(directed){
      node_list <- dyad_list <- vector("list", 5)
      node_list[[1]] <- cbind(nodes, nodes)
      node_list[[2]] <- node_list[[3]] <- node_list[[4]] <- node_list[[5]] <- matrix(0, 0, 4)
      
      dyads_flipped <- dyad(nodes[,2], nodes[,1], n, directed)
      i2b <- which(sapply(dyads_flipped, function(z) z %in% dyads))
      if(length(i2b) > 0){
        d2 <- matrix( cbind(dyads[i2b], dyads_flipped[i2b]), ncol=2)
        for(i in 1:(length(i2b)/2)){
          d2 <- matrix( d2[ -which(d2[,2] == d2[i,1]), ], ncol=2)
        }
        r2 <- match(d2[,1], dyads)
        node_list[[2]] <- matrix( cbind(nodes[r2,], nodes[r2, c(2,1)]), ncol=4)
      }
      
      
      possible_nodes <- sort(unique(c(nodes[,1:2])))
      for(i in possible_nodes){
        i3 <- matrix( combine(which(nodes[,1] == i)), ncol=2)
        if(length(i3) > 0){
          node_list[[3]] <- rbind(node_list[[3]], cbind(matrix(nodes[i3[,1],], ncol=2), matrix(nodes[i3[,2],], ncol=2)) )
        }
        
        i4 <- matrix( combine(which(nodes[,2] == i)), ncol=2)
        if(length(i4) > 0){
          node_list[[4]] <- rbind(node_list[[4]], cbind( matrix(nodes[i4[,1],], ncol=2), matrix(nodes[i4[,2],], ncol=2) ) )
        }
        
        i5 <- matrix( combine(which(nodes[,1] == i), which(nodes[,2] == i)), ncol=2)
        if(length(i5) > 0){
          temp <- matrix( cbind(nodes[i5[,1],], nodes[i5[,2],]), ncol=4)
          rem <- which(temp[,1] == temp[,4] & temp[,2] == temp[,3])
          if(length(rem) > 0){
            temp <- temp[-rem,]
          }
          node_list[[5]] <- rbind(node_list[[5]], temp)
        }
      }
      
    } else { # undirected
      node_list <- dyad_list <- vector("list", 2)
      node_list[[1]] <- cbind(nodes, nodes)
      
      node_list[[2]] <- matrix(0,0,4)
      possible_nodes <- sort(unique(c(nodes[,1:2])))
      for(i in possible_nodes){
        i2 <- matrix( combine( c(which(nodes[,1] == i), which(nodes[,2] == i)) ), ncol=2)
        if(length(i2) > 0){
          temp <- matrix( cbind(matrix(nodes[i2[,1],], ncol=2), matrix(nodes[i2[,2],], ncol=2)), ncol=4) 
          node_list[[2]] <- rbind(node_list[[2]], temp)
        }
      }
    }
    
    node_list <- lapply(node_list, unique)
    
    dyad_list <- lapply(node_list, function(z) cbind( dyad(z[,1], z[,2], n, directed), dyad(z[,3], z[,4], n, directed)) )
    dout1 <- lapply(dyad_list, function(z) cbind(match(z[,1], dyads),  match(z[,2], dyads) ))
  }
  
  return(dout1)
}



#' Calculate E meat using rows of X, e
#'
#' @keywords internal
meat.E.row <- function(row.list, X, e)
{
  # Build meat for exchangeable variances estimator
  # Takes in design matrix X.1, residuals e.1, and list of overlapping nodes 
  
  
  # sizes
  p <- ncol(X)
  
  # Estimates
  phi <- rep(0, length(row.list))
  phi[1] <- mean(e^2, na.rm=T)
  meat.1 <- phi[1]*crossprod(X)
  
  meat.2 <- matrix(0,p,p)
  for (i in 2:length(row.list)){
    r <- matrix( row.list[[i]], ncol=2)
    if(nrow(r) > 0){
      phi[i] <- mean(e[ r[,1] ]*e[ r[,2] ], na.rm=T)
      
      meat.2 <- meat.2 + crossprod(X[r[,1], ,drop=F], X[r[,2], ,drop=F])*phi[i]  + crossprod(X[r[,2], ,drop=F], X[r[,1], ,drop=F])*phi[i]
    }
  }
  
  meat.out <- meat.2 + meat.1
  
  return(list(M=meat.out, phi=phi))
}


#' Calculate parameter estimates using rows of e
#'
#' @keywords internal
param_est <- function(row.list, e)
{
  
  # Estimates
  phi <- rep(0, length(row.list))
  phi[1] <- mean(e^2, na.rm=T)
  
  for (i in 2:length(row.list)){
    r <- row.list[[i]]
    phi[i] <- mean(e[ r[,1] ]*e[ r[,2] ], na.rm=T)
  }
  
  return(phi)
}


#' Replace negative eigenvalues with zeros in variance matrix
#'
#' @keywords internal
make.positive.var <- function(V.test)
{
  # Make variance matrix positive definite by zeroing out negative eigenvalues
  eig.V <- eigen(as.matrix(V.test), symmetric=T)
  eig.flag <- sum(eig.V$values < 0)                # creat flag if correction made
  
  if (eig.flag>0) { 
    eig.V$values[which(eig.V$values<0)] <- 0         # zero out negative eigenvalues 
    V.corrected <- eig.V$vectors %*% diag(eig.V$values) %*% t(eig.V$vectors)
  } else { 
    V.corrected <- V.test
  }
  
  return(list(V=V.corrected, flag=1*(eig.flag>0) ))
}



#' Vectorize a network matrix (without diagonal)
#'
#' @keywords internal
vec.net <- function(A, directed=T)
{
  # convert matrix A to vector
  # unfold along columns, omitting diagonal
  # return n(n-1) vector or n(n-1)/2 vector depending on directed option input
  # Uses lower triangular unfolding for undirected case
  
  n <- nrow(A)
  if (directed==T){ 
    vec.out <- as.vector(A)[-seq(1,n^2,n+1)] 
  } else if (directed==F) {
    vec.out <- A[lower.tri(A)]
  } else {stop('Logical input only for directed input')}
  return(vec.out)
}


#' Matricize a network vector (without diagonal)
#'
#' @keywords internal
mat.net <- function(V, directed=T)
{
  # convert vector V to a matrix
  # build along columns, omitting diagonal
  # return n(n-1) vector or n(n-1)/2 vector depending on directed option input
  # Uses lower triangular unfolding for undirected case
  
  if (directed==T){   
    d <- length(as.vector(V))    
    n <- floor((1+sqrt(1+4*d))/2)
    
    A <- matrix(1:n^2,n,n)
    ind <- vec.net(A, directed)            # indices in matrix
    Mat.out <- matrix(0,n,n)
    Mat.out[ind] <- V
  } else if (directed==F){
    d <- length(as.vector(V))    
    n <- floor((1+sqrt(1+8*d))/2)
    
    Mat.out <- matrix(0,n,n)
    Mat.out[lower.tri(Mat.out)] <- V
    Mat.out <- Mat.out + t(Mat.out)
  } else {stop('Logical only for directed input')}
  return(Mat.out)
}


#' Generate node pairs for complete network
#'
#' @keywords internal
node.gen <- function(n, directed=T)
{
  # Generate pairs of nodes for n.code given network size and directed vs. undirected
  # for columnwise unfolding of matrix
  
  if (directed==T){
    mat.out <- c()
    for (i in 1:n){
      mat.out <- rbind(mat.out, cbind((1:n)[-i], rep(i,n-1)))
    }
  } else if (directed==F){
    mat.out <- c()
    for (i in 1:(n-1)){
      mat.out <- rbind(mat.out, cbind( (i+1):n, rep(i,n-i)))
    }
  } else (stop('Choose T or F logical for directed'))
  
  return(mat.out)
}


#' Find all possible combinations of elements in two vectors, or all combinations of all elements in one without repeats
#'
#' @keywords internal
combine <- function(v1, v2=NA)
{
  # Find all possible combinations of elements in two vectors
  # No pairs (a,a)
  # return Nx2 matrix of all pairs
  # If v2 = NULL, do all possible combinations of elements in v1
  #
  # Make sure don't have repeats in v1 or v2, not coded to remove repeats
  #
  v1 <- sort(v1)
  l1 <- length(v1)
  
  if(is.na(v2[1])){  # if null v2 or not input
    c1 <- rep(v1, times=(l1-1):0)    # first column
    v.mat <- outer(v1,rep(1,l1))
    c2 <- v.mat[lower.tri(v.mat)]    # second column
    m.out <- cbind(c1,c2) 
  } else {    # all combinations of v1 and v2
    l2 <- length(v2)
    v2 <- sort(v2)
    c1 <- rep(v1, each=l2)    # first column
    c2 <- rep(v2, times=l1)   # second column
    m.out <- cbind(c1,c2)
    remove <- which(m.out[,1]==m.out[,2], arr.ind=T)
    if(length(remove) > 0){
      m.out <- m.out[-remove,]   # remove any duplicate rows
    }
  }
  
  return(m.out)
}


#' Generate list indicator matrix of overlapping dyads
#'
#' @keywords internal
Sigma.ind <- function(n.tot, directed, sta_flag=F)
{  
  # Generate dyad sets of various overlapping dyad pairs
  # Return list of each set of dyads
  # Map dyads from matricization numbering to actual rows through d.map
  #    may include zeros
  
  
  if (directed==T){  
    nodes.1 <- cbind(rep(1:n.tot, each=n.tot)[-seq(1,n.tot^2,n.tot+1)], rep(1:n.tot, n.tot)[-seq(1,n.tot^2,n.tot+1)],
                     rep(1:n.tot, each=n.tot)[-seq(1,n.tot^2,n.tot+1)], rep(1:n.tot, n.tot)[-seq(1,n.tot^2,n.tot+1)])
    
    nodes.2 <- nodes.3 <- nodes.4 <- nodes.5 <- nodes.6 <- c()
    
    for (i in 1:n.tot){ 
      
      # ij,ji
      if (i<n.tot){
        c1 <- rep(i,(n.tot-i))
        #       c2 <- (1:n.tot)[-i]
        c2 <- ((i+1):n.tot)
        
        nodes.2 <- rbind(nodes.2,cbind(c1,c2,c2,c1))
      }
      
      # ij,il  ;  ij,kj
      c1 <- rep(i,(n.tot - 1)*(n.tot-2)/2)
      c2 <- rep( (1:(n.tot-1))[-i], times=(n.tot-2):(1 - 1*(i==n.tot) ) )
      c3 <- c1
      c4.mat <- outer(1:n.tot,rep(1,n.tot))[-i,-i]
      c4 <- c4.mat[lower.tri(diag(n.tot-1))]
      
      nodes.3 <- rbind(nodes.3, cbind(c1,c2,c3,c4))
      nodes.4 <- rbind(nodes.4, cbind(c2,c1,c4,c3))
      
      # ij,jl and ij,ki
      nodes.5 <- rbind(nodes.5, cbind(c1,c2,c4,c3), cbind(c2,c1,c3,c4))   
      n.list.tmp <- list(n1=nodes.1, n2=nodes.2, n3=nodes.3, n4=nodes.4, n5=nodes.5)
    }
    
  } else if (directed==F) {
    nodes.1 <- nodes.3 <- nodes.4 <- nodes.5 <- c()
    
    rows <- outer(1:n.tot, rep(1,n.tot))
    cols <- outer(rep(1,n.tot), 1:n.tot)
    r.1 <- rows[lower.tri(rows)]
    c.1 <- cols[lower.tri(cols)]
    
    nodes.1 <- cbind(r.1,c.1,r.1,c.1)    # diagonal/variance
    for (i in 1:n.tot){
      if(i>2){ 
        c.1 <- rep(1:(i-2), (i-2):1)
        c.mat <- outer(1:(i-1), rep(1,(i-1)))
        c.2 <- c.mat[lower.tri(c.mat)]
        nodes.3 <- rbind(nodes.3, cbind(rep(i, length(c.1)),c.1,rep(i, length(c.1)),c.2))
        #         nodes.5 <- rbind(nodes.5, cbind(rep(i, length(c.1)),c.1,rep(i, length(c.1)),c.2))
      }
      if (i < (n.tot - 1)){
        r.1 <- rep( (i + 1):(n.tot-1), (n.tot-i-1):1 ) 
        r.mat <- outer((1+i):n.tot, rep(1,(n.tot - i)))
        r.2 <- r.mat[lower.tri(r.mat)]
        nodes.4 <- rbind(nodes.4, cbind(r.1, rep(i, length(r.1)), r.2, rep(i, length(r.1))))
        nodes.5 <- rbind(nodes.5, cbind(r.1, rep(i, length(r.1)), r.2, r.1))
      }
    }
    n.list.tmp <- list(n1=nodes.1, n2=rbind(nodes.3,nodes.4,nodes.5))
  } else {stop('Logical only for directed input')}
  
  
  # Form dyad list 
  d.tot <- n.tot*(n.tot-1)/2*(1 + directed)
  S.list <- vector('list', length(n.list.tmp))
  d.list <- vector('list', length(n.list.tmp))
  names(S.list) <- paste('S',1:length(n.list.tmp),sep='')
  for (i in 1:length(n.list.tmp)){
    nodes.tmp <- n.list.tmp[[i]]
    d1 <- dyad(nodes.tmp[,1], nodes.tmp[,2], n.tot, directed)
    d2 <- dyad(nodes.tmp[,3], nodes.tmp[,4], n.tot, directed)
    #     mat <- cbind(d1, d2)
    #     mat <- apply(mat,1,sort)
    L1 <- (d1<=d2)
    L2 <- d2<d1
    s1 <- d1*L1 + d2*L2
    s2 <- d1*L2 + d2*L1
    #     S.list[[i]] <- cbind(d1,d2)
    #     S.list[[i]] <- sparseMatrix(i=mat[,1], j=mat[,2], x=1, dims=c(d.tot, d.tot))  #, symmetric=T)
    S.list[[i]] <- sparseMatrix(i=s1, j=s2, x=1, dims=c(d.tot, d.tot), symmetric=T)
    d.list[[i]] <- cbind(s1,s2)
  }
  
  if(sta_flag == T & directed==T){   # change listing of 5th overlap 
    tmp <- matrix(1, n.tot, n.tot)
    diag(tmp) <- 0
    nodemat <- which(tmp == 1, arr.ind=T)   # columnwise unfolding dyad list
    i1 <- which(S.list[[5]]==1, arr.ind=T)
    nodes.5 <- cbind(nodemat[i1[,1], ], nodemat[i1[,2], ])
    i5 <- which(nodes.5[,1] == nodes.5[,4])
    S.list[[5]] <- S.list[[6]] <- sparseMatrix(i=1, j=1, x=0, dims=c(d.tot, d.tot), symmetric=T)
    S.list[[5]][i1[i5, ]] <- 1
    S.list[[6]][i1[-i5, ]] <- 1
    
  } else if (sta_flag == T & directed==F){
    stop("stationary listing not impememented for undirected")
  }
  
  S.list[[length(S.list) + 1]] <- 1*(Reduce("+", S.list) == 0)   # all the rest (zeros)
  
  return(S.list)
}


#' Build intermediate C(phi,n) matrix in inversion of Exchangeable variance matrix
#'
#' @keywords internal
build_phi_matrix <- function(n, phi, directed=T, sta_flag=F )
{
  # Build intermediate phi/p matrix in inversion of E/E matrix
  # A(n, phi) %*% p = c(1,rep(0, 5))
  #
  if(length(phi) != 6 & directed){stop("wrong length of phi")}
  
  if(sta_flag == F & directed==T){
    
    A <- matrix(NA, 6, 6)
    A[1,] <- c(phi[1], phi[2], (n-2)*phi[3], (n-2)*phi[4], 2*(n-2)*phi[5], (n-3)*(n-2)*phi[6])
    A[2,] <- c(phi[2], phi[1], (n-2)*phi[5], (n-2)*phi[5], (n-2)*phi[3] + (n-2)*phi[4], (n-3)*(n-2)*phi[6])
    A[3,] <- c(phi[3], phi[5], phi[1] + (n-3)*phi[3], phi[5] + (n-3)*phi[6], phi[2] + phi[4] + (n-3)*phi[5] + (n-3)*phi[6],
               (n-3)*(phi[4] + phi[5] + (n-4)*phi[6]))
    A[4,] <- c(phi[4], phi[5], phi[5] + (n-3)*phi[6], phi[1] + (n-3)*phi[4], phi[2] + phi[3] + (n-3)*phi[5] + (n-3)*phi[6],
               (n-3)*(phi[3] + phi[5] + (n-4)*phi[6]))
    A[5,] <- c(phi[5], phi[4], phi[2] + (n-3)*phi[5], phi[3] + (n-3)*phi[6], phi[1] + phi[5] + (n-3)*phi[4] + (n-3)*phi[6],
               (n-3)*(phi[3] + phi[5] + (n-4)*phi[6]))
    # A[5,] <- c(phi[5], phi[3], phi[4] + (n-3)*phi[6], phi[2] + (n-3)*phi[5], phi[1] + phi[5] + (n-3)*phi[3] + (n-3)*phi[6],
    # (n-3)*(phi[4] + phi[5] + (n-4)*phi[6]))
    A[6,] <- c(phi[6], phi[6], phi[4] + phi[5] + (n-4)*phi[6], phi[3] + phi[5] + (n-4)*phi[6], phi[3] + phi[4] + 2*phi[5] + 2*(n-4)*phi[6],
               phi[1] + phi[2] + (n-4)*(phi[3] + phi[4] + 2*phi[5] + (n-5)*phi[6]))
    
  } else if (sta_flag == F & directed==T) {
    A <- matrix(NA, 7, 7)
    A[1,] <- c(phi[1], phi[2], (n-2)*phi[3], (n-2)*phi[4], (n-2)*phi[6], (n-2)*phi[5], (n-3)*(n-2)*phi[7])
    
    A[2,] <- c(phi[2], phi[1], (n-2)*phi[6], (n-2)*phi[5], (n-2)*phi[3], (n-2)*phi[4], (n-3)*(n-2)*phi[7])
    
    A[3,] <- c(phi[3],   phi[5],   phi[1] + (n-3)*phi[3], 
               phi[6] + (n-3)*phi[7],   phi[4] + (n-3)*phi[7],   phi[2] + (n-3)*phi[5],
               (n-3)*(phi[4] + phi[6] + (n-4)*phi[7]))
    
    A[4,] <- c(phi[4],   phi[6],   phi[5] + (n-3)*phi[7], 
               phi[1] + (n-3)*phi[4],   phi[2] + (n-3)*phi[6],   phi[3] + (n-3)*phi[7],
               (n-3)*(phi[3] + phi[5] + (n-4)*phi[7]))
    
    A[5,] <- c(phi[5],   phi[3],  phi[4] + (n-3)*phi[7], 
               phi[2] + (n-3)*phi[5],   phi[1] + (n-3)*phi[3],   phi[6] + (n-3)*phi[7],
               (n-3)*(phi[4] + phi[6] + (n-4)*phi[7]))
    
    A[6,] <- c(phi[6],   phi[4],   phi[2] + (n-3)*phi[6], 
               phi[3] + (n-3)*phi[7],   phi[5] + (n-3)*phi[7],   phi[1] +  (n-3)*phi[4],
               (n-3)*(phi[3] + phi[5] + (n-4)*phi[7]))
    
    A[7,] <- c(phi[7],   phi[7],   phi[4] + phi[5] + (n-4)*phi[7], 
               phi[3] + phi[6] + (n-4)*phi[7],    phi[4] + phi[5] + (n-4)*phi[7],    phi[3] + phi[6] + (n-4)*phi[7],
               phi[1] + phi[2] + (n-4)*(phi[3] + phi[4] + phi[5] + phi[6] + (n-5)*phi[7]))
    
  } else if (directed == F){
    
    if(length(phi) != 3){stop("wrong length of phi")}
    C <- matrix(0, 3, 3)
    C[1,] <- c(phi[1],  2*(n-2)*phi[2],  .5*(n-2)*(n-3)*phi[3])
    C[2,] <- c(phi[2],  phi[1] + (n-2)*phi[2] + (n-3)*phi[3], (n-3)*phi[2] + (.5*(n-2)*(n-3) - (n-3))*phi[3] )
    C[3,] <- c(phi[3],  4*phi[2] + (2*n - 8)*phi[3],  phi[1] + (2*n - 8)*phi[2] + (.5*(n-2)*(n-3) - ( 2*(n-2) - 4) - 1)*phi[3])
    A <- C
  }
  
  return(A)
}




#' Dyad map from nodes i,j --> dyad d
#'
#' @keywords internal
dyad <- function(i.in, j.in, n.tot, directed=T) 
{
  # transform n indices to dyad indices
  if (directed==T){ 
    dyad.result <- ((i.in-1) + (j.in-1)*(n.tot-1) + (j.in > i.in)) * (i.in != j.in)
  } else if (directed==F) {
    #     if (j.in>i.in){stop('Only lower triangular matrix for undirected network')}
    ij <- cbind(i.in, j.in)
    ij <- t(apply(ij, 1, sort, decreasing=T))
    j.1 <- ij[,2]    # want invariance to ordering
    i.1 <- ij[,1]    # want invariance to ordering
    dyad.result <- ( i.1 - 1 -.5*j.1^2 + (n.tot - 1/2)*j.1 - n.tot + 1  ) * (i.1 != j.1)
  } else {stop('Logical T/F only for directed input')}
  return(dyad.result)
}





#' Compute symmetric square root of A, assuming it is real, symmetric, positive definite
#'
#' @keywords internal
symm_square_root <- function(A, prec=15)
{
  e1 <- eigen(A)
  vals <- round(e1$val, prec)   # get rid of barely negative zeros
  r1 <- e1$vec %*% diag(sqrt(vals)) %*% t(e1$vec)   # symmetric square root 
  r1
}


#' Matrix product of A^TBC where B is a short list of parameters
#'  A and C must be matrices
#'  B is parameterized by phi, row.list, and assumed symmetric without repeats, with phi[1] along diagonal
#'
#' @keywords internal
meatABC <- function(row.list, phi, A, C, directed)
{
  L <- 2 + 3*directed
  phi6 <- length(phi) > L
  outX <- crossprod(A,C)
  out <- phi[1]*outX
  
  for(i in 2:L){
    r <- row.list[[i]]
    X1 <- crossprod(A[r[,1],], C[r[,2],]) + crossprod(A[r[,2],], C[r[,1],])
    out <- out + phi[i]*( X1 )
    if(phi6){
      outX <- outX + X1
    }
  }
  
  if(phi6){
    XX <- tcrossprod(apply(A, 2, sum), apply(C, 2, sum))
    out <- out + (XX - outX)*phi[L+1]
  }
  
  return(out)
}


