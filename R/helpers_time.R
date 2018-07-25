# Helper functions to accomodate temporal networks
#
#




#' Pre-processes data for ordering, FOR TEMPORAL DATA, etc.
#'
#' TO-DO:
#'    missing data
#' 
#' @keywords internal
node_preprocess_time <- function(Y, X, directed, nodes, tmax, type, subtract=NULL)
{
  if(!is.null(nodes)){
    nodes <- as.matrix(nodes)
    if(ncol(nodes) != 3){stop("nodes must have 3 columns when tmax > 1")}
    tmax1 <- max(nodes[,3])
    if(tmax1 != tmax){stop("Third dimension of relational array in nodes must match input tmax")}
  }
  
  
  type <- as.character(type)
  if( !(strtrim(type, 1) %in% c("e", "E", "I", "i")) ){stop("Invalid input for type.")}
  if( tolower(strtrim(type, 1)) ==  "e"){ errortype <- "EE" }
  if( tolower(strtrim(type, 1)) ==  "i"){ errortype <- "EI" }
  
  
  #### Preprep
  directed <- as.logical(directed)
  
  Y <- as.vector(as.numeric(Y))
  if(is.null(dim(X))){  # if X is a vector
    X <- matrix(X, ncol=1)
  } else {
    X <- as.matrix(X)  # matrix
  }
  
  
  d <- length(Y) / tmax
  if(nrow(X) != length(Y)){stop("X and Y must have the same number of observations")}
  
  remove <- which(is.na(Y))
  if(length(remove) > 0){
    Y <- Y[-remove]
    X <- X[-remove,]
    missing <- TRUE
    if(!is.null(nodes)){
      nodes <- nodes[-remove,]
    } else {
      nodes <- node_gen_time(n, directed, tmax)[-remove,]
    }
  }
  ####
  
  
  if(is.null(nodes)){
    missing <- FALSE
    
    cc <- 4 + 4*(directed==FALSE)
    n <- (1+sqrt(1+cc*d))/2
    if(n != round(n)){stop("Size of Y and X must be valid for a network of some size; assuming complete network at this point")}
    row_list <- row_list_time(n, directed, tmax, errortype)
    dyads <- 1:(d*tmax)
    
    
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
    
    if(length(dyads) == dmax*tmax){  # just reorder
      missing <- FALSE
      reorder <- order(nodes[,3], nodes[,2], nodes[,1])   # columnwise unfolding
      
      # proceed as usual
      Y <- Y[reorder]
      X <- X[reorder,]
      row_list <- row_list_time(n, directed, tmax, errortype)
      
    } else {  #build row_list based on nodes input
      
      stop("Missing data not yet implemented for temporal data")
      # missing <- TRUE
      # if(is.null(subtract)){
      #   subtract <- length(dyads)/max(dyads) > .5   # are over half the dyads present?
      # }
      # row_list <- row_list_missing(nodes, dyads, directed, subtract)   # row list based on overlaps
    }
  }
  
  return(list(X=X, Y=Y, dyads=dyads, missing=missing, row_list=row_list, n=n, type=errortype))
}





#' Perform GEE estimate / IRWLS of coefficients for temporal data
#' 
#' TO DO:
#'  eigenvalues / ndstop
#'  fix for missingness
#'
#' @keywords internal
GEE_est_time <- function(Y, X, n, tmax,  directed=TRUE, type="EE", write_dir=NULL, missing=FALSE, tol.in=1e-6, maxit=1e4, ndstop=TRUE, verbose=FALSE) 
{
  
  if(missing){stop("GEE.est.time(.) not yet coded for missingness")}
  # Calculate GEE estimate of regressors beta for continuous data
  # Return beta, residuals, weight matrix, # iterations, objective function Q, and actual tolerance
  # 
  Y.in <- Y
  X.in <- X
  wo <- !is.null(write_dir)
  if(wo){
    dir.create(write_dir)
  }
  
  t1 <- proc.time()[3]
  
  count.in <- 0
  d <- d.tot <- n*(n-1)/2*(1 + directed)
  
  XX <- chol2inv(chol(crossprod(X)))
  
  # Overlapping pair indicators
  S.list <- Sigma.ind(n, directed, sta_flag=FALSE)
  # W.ind <- Sigma.ind(n, directed, sta_flag=TRUE)   # 7 parameter version
  # ilist <- lapply(S.list, find_ones)
  # d.tot <- nrow(S.list[[1]])
  node_list <- node.set(n, directed) 
  ilist <- lapply(node_list, function(z) cbind( dyad(z[,1], z[,2], n, directed), dyad(z[,3], z[,4], n, directed)) )
  ilist[length(ilist) + 1] <- NA
  
  # Initial fit
  beta.0 <- XX %*% t(X.in) %*% Y.in
  se.ols <- sqrt(diag(XX * sum((Y.in - X.in %*% beta.0)^2)/(nrow(X.in) - ncol(X.in))))
  if(wo){
    write.table(beta.0, file.path(write_dir, "beta_ols.txt"), row.names=F, col.names=F)
  }
  # write.table(se.ols, file.path(write_dir, "beta_olsSE.txt"), row.names=F, col.names=F)
  if(verbose){
    cat("OLS fit complete \n")
  }
  
  eols <- Y.in - X.in %*% beta.0
  # meat_DC <- meatDC2(node.mat, X.in, eols)
  # Vdc <- XX %*% meat_DC %*% XX
  # Vdc <- make.positive.var(Vdc)
  # se_dc <- sqrt(diag(Vdc$V))
  # write.table(se_dc, file.path(write_dir, paste0("se_dc_t", tmax,".txt")), row.names=F, col.names=F)
  # cat("DC SE complete \n")
  
  # # Start looping
  # if(nchar(beta_start) == 0){
  #   # if(type != "STA") {
  #   beta.gls <- beta.0
  # } else {
  #   
  #   beta.gls <- as.vector(read.table(beta_start, header=F)[,1])
  #   # # For stationary fit, first do EE fit!!
  #   # cat('************************************ \n')
  #   # cat("Fitting GEE with EE/6a before doing STA/6b \n")
  #   # temp <- GEE.est(Y.in, X.in, n.tot, t.max, tol.in, directed, type="EE", write_dir)
  #   # beta.gls <- temp$beta
  #   # rm(temp)
  #   # gc()
  #   # cat("Finished preliminary EE/6a fit \n")
  #   # cat('************************************ \n\n')
  # }
  
  beta.gls <- beta.old <- beta.0
  e.test <- Y - X %*% beta.gls
  Q.old <- Q.new <- sum(beta.gls^2)   # sum(e.test^2)
  delta.loop <- tol.in*100
  params <- NA
  
  if(wo){
    out <- data.frame(matrix(c(count.in, Q.new), ncol=2))
    names(out) <-  c("iteration", "criterion")
    write.table(out, file.path(write_dir, "GEE_criterion.txt"), row.names = F, sep="\t", quote=F)
  }
  # if(verbose){
    # cat("Iteration:  ", count.in, " Criterion:  ", Q.new, "\n")
  # }
  
  
  while (abs(delta.loop) > tol.in & count.in < maxit) {
    count.in <- count.in + 1
    
    params <- calculate_matrix_params(ilist, e.test, tmax, type)  # unique parameters in Var(Y)
    invParams <- calculate_parameter_inverse(n, tmax, params, type)  # calculate inverse parameters of Var(Y)^{-1}
    
    mine <- min( eigen_exch_time(n, tmax, params, directed, type) )
    if(mine < -.Machine$double.eps^.5){
      warning("Negative definite weight matrix encountered")
      cat("\nValue: ", mine, "\n")
      if(ndstop){
        break()
      }
    }
    
    # Following if-statement calculates GEE criterion, beta, etc. based on type of matrix
    # Avoid building large inverse matrix, W = Var(Y)^{-1} 
    # find bread = (X^T W X)^{-1} and X^T W Y WITHOUT building W
    #
    if(type == "EI"){
      
      Wmini <- Reduce("+", lapply(1:6, function(x) invParams$p[x]*S.list[[x]]))
      Wmini <- as.matrix(Wmini)

      XWX <- Reduce("+", lapply(1:tmax, function(x) crossprod(X.in[1:d.tot + (x-1)*d.tot,], Wmini %*% X.in[1:d.tot + (x-1)*d.tot,]) ))
      
      bread <- solve(XWX)
      
      XWY <- Reduce("+", lapply(1:tmax, function(x) crossprod(X.in[1:d.tot + (x-1)*d.tot,], Wmini %*% Y.in[1:d.tot + (x-1)*d.tot]) ))
      
      beta.gls <- as.vector( solve(XWX, XWY) )  # new beta fit
      e.test <- as.vector( Y.in - X.in %*% beta.gls )   # residuals
      # eWe <- lapply(1:tmax, function(x) crossprod(e.test[1:d.tot + (x-1)*d.tot],Wmini %*% e.test[1:d.tot + (x-1)*d.tot]) )
      # Q.new <- Reduce("+", eWe)
      
      
    } else if (type == "EE"){
      A <- as.matrix( Reduce("+", lapply(1:6, function(x) invParams$p[x]*S.list[[x]])) )
      B <- as.matrix( Reduce("+", lapply(1:6, function(x) invParams$p[6 + x]*S.list[[x]])) )
      
      ABX <- lapply(1:tmax, function(x) A %*% X.in[1:d.tot + (x-1)*d.tot,]
                    + B %*% Reduce("+", lapply((1:tmax)[-x], function(t) X.in[1:d.tot + (t-1)*d.tot,])) )
      XWX <- Reduce("+", lapply(1:tmax, function(x) crossprod(X.in[1:d.tot + (x-1)*d.tot,], ABX[[x]])))
      bread <- solve(XWX)
      rm(ABX)
      
      ABY <- lapply(1:tmax, function(x) A %*% Y.in[1:d.tot + (x-1)*d.tot]
                    + B %*% Reduce("+", lapply((1:tmax)[-x], function(t) Y.in[1:d.tot + (t-1)*d.tot])) )
      XWY <- Reduce("+", lapply(1:tmax, function(x) crossprod(X.in[1:d.tot + (x-1)*d.tot,], ABY[[x]])))
      rm(ABY)
      
      beta.gls <- as.vector( solve(XWX, XWY) )  # new beta fit
      e.test <- as.vector( Y.in - X.in %*% beta.gls )   # residuals
      
      # ABe <- lapply(1:tmax, function(x) A %*% e.test[1:d.tot + (x-1)*d.tot]
                    # + B %*% Reduce("+", lapply((1:tmax)[-x], function(t) e.test[1:d.tot + (t-1)*d.tot])) )
      # Q.new <- Reduce("+", lapply(1:tmax, function(x) crossprod(e.test[1:d.tot + (x-1)*d.tot], ABe[[x]])))
      # rm(ABe)
      # gc()
      
    } else if (type == "STA"){
      stop("Not yet implemented")
      # 
      # blocklist <- build_blocklist(tmax, type=invParams$type)
      # Psi_blocks <- matrix(NA, tmax, tmax)
      # for(i in 1:length(blocklist)){
      #   Psi_blocks[matrix(blocklist[[i]], ncol=2)] <- i
      # }
      # 
      # Alist <- lapply(1:length(blocklist), function(y) Reduce("+", lapply(1:7, function(x) invParams$p[x + 7*(y-1)]*W.ind[[x]])))   # list of all blocks
      # 
      # WX <- lapply(1:tmax, function(x) Reduce("+", 
      #                                         lapply(1:tmax, function(t) Alist[[ Psi_blocks[x,t] ]] %*% X.in[1:d.tot + (t-1)*d.tot, ])) )   # list entries are block columns in WX
      # XWX <- Reduce("+", lapply(1:tmax, function(x) crossprod(X.in[1:d.tot + (x-1)*d.tot,], WX[[x]]))  )
      # bread <- solve(XWX)
      # rm(WX)
      # 
      # 
      # WY <- lapply(1:tmax, function(x) Reduce("+", 
      #                                         lapply(1:tmax, function(t) Alist[[ Psi_blocks[x,t] ]] %*% Y.in[1:d.tot + (t-1)*d.tot])) )   # list entries are block columns in WX
      # XWY <- Reduce("+", lapply(1:tmax, function(x) crossprod(X.in[1:d.tot + (x-1)*d.tot,], WY[[x]])) )
      # rm(WY)
      # 
      # beta.gls <- as.vector( solve(XWX, XWY) )  # new beta fit
      # e.test <- as.vector( Y.in - X.in %*% beta.gls )   # residuals
      # 
      # We <- lapply(1:tmax, function(x) Reduce("+", 
      #                                         lapply(1:tmax, function(t) Alist[[ Psi_blocks[x,t] ]] %*% e.test[1:d.tot + (t-1)*d.tot])) )
      # Q.new <- Reduce("+", lapply(1:tmax, function(x) crossprod(e.test[1:d.tot + (x-1)*d.tot], We[[x]])))
      # rm(We)
      # gc()
      
    } else {
      stop("Not yet implemented")
    }
    
    Q.new <- delta.loop <- mean( (beta.gls - beta.old)^2 )   # want this to be small!
    
    
    if(wo){
      write.table(cbind(count.in, as.numeric(Q.new)), file.path(write_dir, "GEE_criterion.txt"), row.names = F, sep="\t", append=T, col.names=F)  # write out criterion
      write.table(beta.gls, file=file.path(write_dir, paste('beta_loop',count.in,'_t', tmax, '.txt',sep='')), row.names=F, col.names=F, sep='\t')
    }
    
    if(verbose){
      cat('\n************************************ \n')
      cat('GEE iteration: \t\t', count.in, '\n')
      cat('Change in criterion: \t', abs(delta.loop), '\n')
      cat('Elapsed time: \t\t', proc.time()[3]- t1, 'sec \n')
      cat("Warnings: \n")
      print(warnings())
      cat('\n************************************ \n')
    }
    # delta.loop <- as.numeric(Q.old - Q.new)
    
    # Update for next loop
    Q.old <- Q.new # new baseline
    beta.old <- beta.gls
  }
  
  if(wo){
    write.table(sqrt(diag(bread)), file=file.path(write_dir, paste('se_t', tmax, '.txt',sep='')), row.names=F, col.names=F, sep='\t')
  }
  if(verbose){
    cat("\nGEE estimation complete \n")
  }
  e.test <- as.vector( Y.in - X.in %*% beta.gls )   # residuals
  
  
  if(!missing){
    W <- list(Omega1=NA, Omega2=NA)
    W$Omega1 <- Reduce("+", lapply(1:6, function(x) invParams$p[x]*S.list[[x]]))
    if(type == "EE"){
      W$Omega2 <- Reduce("+", lapply(1:6, function(x) invParams$p[x+6]*S.list[[x]]))
    }
  }
  
  output <- list(beta.gls, e.test, XWX, params, count.in, as.logical(count.in < maxit), W)
  names(output) <- c('beta','residuals', 'bread', 'phi', 'nit', 'convergence', 'W')
  
  
  # output <- list(beta.gls, e.test, bread, invParams=invParams$p, count.in, delta.loop)
  # names(output) <- c('beta','resid','bread','W', 'n', 'tol')
  return(output)
}



#' Make complete node indices for temporal relational data
#' 
#' @keywords internal
node_gen_time <- function(n, directed, tmax)
{
  ntemp <- node.gen(n, directed)
  d <- nrow(ntemp)
  nodesout <- matrix(0, d*tmax, 3)
  for(t in 1:tmax){
    nodesout[1:d + d*(t-1),] <- cbind(ntemp, t)
  }
  return(nodesout)
}


#' Make row list for complete temporal relational data
#' 
#' @keywords internal
row_list_time <- function(n, directed, tmax, type)
{
  n <- as.numeric(n)
  directed <- as.logical(directed)
  d <- n*(n-1)/2*(1 + directed)
  
  node_list <- node.set(n, directed) 
  row_list_temp <- lapply(node_list, function(z) cbind( dyad(z[,1], z[,2], n, directed), dyad(z[,3], z[,4], n, directed)) )
  dyads <- 1:(d*tmax)
  
  if(type == "EE"){
    row_list <- vector("list", length(row_list_temp)*2)
    row_list[[1]] <- cbind(1:(d*tmax), 1:(d*tmax))
    
    for(k in 2:length(row_list_temp)){
      for(t in 1:tmax){
        row_list[[k]] <- rbind(row_list[[k]], cbind( row_list_temp[[k]][,1] + (t-1)*d, row_list_temp[[k]][,2] + (t-1)*d))
      }
    }
    
    
    for(s in 1:(tmax-1)){
      for(t in (s+1):tmax){
        for(k in 1:length(row_list_temp)){
          j <- k + length(row_list_temp)
          row_list[[j]] <- rbind(row_list[[j]], cbind( row_list_temp[[k]][,1] + (s-1)*d, row_list_temp[[k]][,2] + (t-1)*d))
          if(directed){
            row_list[[j]] <- rbind(row_list[[j]], cbind( row_list_temp[[k]][,2] + (s-1)*d, row_list_temp[[k]][,1] + (t-1)*d)) 
          }
        }
      }
    }
  
  } else if (type == "EI"){
    row_list <- vector("list", length(row_list_temp))
    row_list[[1]] <- cbind(1:(d*tmax), 1:(d*tmax))
    
    for(k in 2:length(row_list_temp)){
      for(t in 1:tmax){
        row_list[[k]] <- rbind(row_list[[k]], cbind( row_list_temp[[k]][,1] + (t-1)*d, row_list_temp[[k]][,2] + (t-1)*d))
      }
    }
  }

  return(row_list)
}

#' calculate parameter estimates for different types of matrices, i.e. 6a, 6b, etc.
#' 
#' @keywords internal
calculate_matrix_params <- function(ilist, e, tmax, type="EI")
{ # calculate parameter estimates for different types of matrices
  # EI is 6a with \Omega_2 = 0
  # EE is 6a
  # STA is 6b
  # HET is 6c
  # 
  d.tot <- length(e)/tmax
  
  ilist <- ilist[1:(length(ilist) - 1)]   # remove zero-entry list
  blocklist <- build_blocklist(tmax, type)
  
  # Calculate parameters for first 5 entries or first 2 (directed vs. undirected)
  params <- c(sapply(blocklist, function(x) sapply(ilist, param_est_single_ilist, x, e, d.tot)))
  
  # Add in zeros for 6th parameters
  param_mat <- matrix(params, nrow=length(ilist))
  param_mat <- rbind(param_mat, 0)
  params <- c(param_mat)
  
  return(params)
}

#' build list of time blocks that are correlated based on the maximum time and type of temporal model
#' #' 
#' @keywords internal
build_blocklist <- function(tmax, type="EI")
{  # build list of time blocks that are correlated based on the maximum time and type of temporal model
  # EI is 6a with \Omega_2 = 0
  # EE is 6a
  # STA is 6b
  # HET is 6c
  # 
  if(type == "EI"){
    blocklist <- vector("list", 1)    # EE matrix
    blocklist[[1]] <- cbind(1:tmax, 1:tmax)   # on-diagonal
  } else if (type == "EE"){
    blocklist <- vector("list", 2)    # EE matrix
    blocklist[[1]] <- cbind(1:tmax, 1:tmax)   # on-diagonal
    blocklist[[2]] <- combine(1:tmax)  # off-diagonal
  } else if (type == "STA"){
    blocklist <- lapply(1:tmax, function(x) cbind(1:(tmax - x + 1), x:tmax))   # stationary block mat
  } else if (type == "STAINV") {
    Psi_blocks <- diag(c(1:(tmax/2),ceiling(tmax/2):1))  # block matrix of Omega^{-1} for STA
    Psi_blocks[lower.tri(Psi_blocks)] <- max(Psi_blocks) +1:(tmax*(tmax-1)/2)
    Psi_blocks <- Psi_blocks + t(Psi_blocks*lower.tri(Psi_blocks, diag=F))
    blocklist <- lapply(1:max(Psi_blocks), function(x) which(Psi_blocks == x, arr.ind=T))
  } else if (type == "HET") {
    blockmat <- rbind(cbind(1:tmax, 1:tmax), combine(1:tmax))
    blocklist <- lapply(1:nrow(blockmat), function(x) matrix(blockmat[x, ], ncol=2))   # same as stationary
  } else if (type == "HET2") {
    blockmat <- rbind(cbind(1:tmax, 1:tmax), combine(1:tmax), combine(1:tmax)[,c(2,1)])
    blocklist <- lapply(1:nrow(blockmat), function(x) matrix(blockmat[x, ], ncol=2))   # same as stationary
  } else {
    stop("Invalid type")
  }
  
  return(blocklist)
}

#' Invert matrix parameters based on inputs.
#' 
#' @keywords internal
calculate_parameter_inverse <- function(n, tmax, params, type="EI")
{
  # Invert matrix parameters based on inputs. 
  #
  # Returns: 
  #   p, vector double, is parameters of inverse matrix, length 6*# unique blocks
  #   type, the type of matrix that should be built using build_matrix function      
  #
  # Args:
  #   n, integer, is network size (nodes)
  #   t, integer, is number of ``time'' periods
  #   phi, vector double, is parameters STA matrix to invert.
  
  if(tmax==1 | type=="EI"){ 
    A <- build_phi_matrix(n, params[1:6])
    return(list(p = solve(A, c(1, rep(0,5))), type="EI"))
  }    # single time period solution
  
  else if (type == "EE" & tmax!=1) {
    A <- build_phi_matrix(n, params[1:6])
    B <- build_phi_matrix(n, params[7:12])
    C <- matrix(NA, 12, 12)
    C[1:6, 1:6] <- A
    C[1:6, 7:12] <- (tmax-1)*B
    C[7:12, 1:6] <- B
    C[7:12, 7:12] <- A + (tmax-2)*B
    
    return(list(p=solve(C, c(1, rep(0,11))), type="EE"))
    
  } else if (type == "STA" & tmax!=1) {  
    # Stationary, needs 7-parameter solution
    nt <- tmax*(tmax - 1)/2 + ceiling(tmax/2)   # number of unique blocks in inverse
    
    numcol <- nt*7 - tmax   # number of unique parameters in inverse
    
    # Stationary, needs 7-parameter solution
    pmat <- matrix(params, 6, tmax)   # blow up parameters to 7 of them
    pmat <- rbind(pmat, pmat[6,])
    pmat[6,] <- pmat[5,]
    
    blocklist <- build_blocklist(tmax, "HET2")
    nblocks <- length(blocklist)
    blockmat <- matrix(unlist(blocklist), ncol=2, byrow=2)
    
    Cmat <- matrix(0, 7*nblocks, 7*nblocks)
    Dvec <- rep(0, 7*nblocks)
    
    Alist7 <- lapply(1:tmax, function(x) build_phi_matrix(n, as.vector(pmat[,x]), sta_flag=T))   # find all phi matrices
    
    Oblocks <- diag(tmax)
    Oblocks <- abs(row(Oblocks) - col(Oblocks)) + 1  # block matrix of Omega's for STA
    
    Psi_blocks <- diag(1:tmax)  # block matrix of Omega^{-1} for STA
    Psi_blocks[blockmat] <- 1:nblocks
    
    rc <- blockmat  # each row is row in O and column in Psi to multiply and get an equation
    
    # pair_mat <- NULL
    for(i in 1:nrow(rc)){
      x <- rc[i,1]
      y <- rc[i,2]
      Dvec[1 + (i-1)*7] <- x == y   # if along diagonal, should be Identity, o.w. zero
      
      for(j in 1:tmax){
        
        Cmat[(i - 1)*7 + 1:7, (Psi_blocks[j, y]-1)*7 + 1:7] <- Alist7[[ Oblocks[x, j] ]]
        
      }
    }
    Dvec <- Dvec*1
    
    
    return(list(p=solve(Cmat, Dvec), type="HET2"))
    
  } else if (type == "HET" & tmax!=1) {
    stop("HET inverse not yet implemented")
    
  }
  
  
}

#' Given matrix of time blocks and a particular exchangeable parameter set (within each block),
#' calculate a single parameter/phi. ASSUMES NO MISSING DATA
#' 
#' @keywords internal
param_est_single_ilist <- function(ilist, blockmat, e, d.tot)
{

  #
  # e, vector double, residuals
  # ilist, matrix of 2 columns, list of i^th parameter entries (within block)
  # blockmat, matrix of 2 columns, time blocks to consider
  # 
  
  partial_sum <- sapply(1:nrow(blockmat), function(x) 
    sum(e[(blockmat[x,1]-1)*d.tot + ilist[,1]]*e[(blockmat[x,2]-1)*d.tot + ilist[,2]]) )
  
  return(sum(partial_sum)/(nrow(ilist)*nrow(blockmat)))
}


#' Compute eigenvalues of covariance matrices of jointly exchangeable errors with repeated observations
#' 
#' @keywords internal
eigen_exch_time <- function(n, tmax, params, directed, type)
{
  if(type == "EI"){
    e <- eigen_exch(n, params[1:6], directed, calcall=FALSE)$uniquevals
  } else if (type == "EE"){
    e1 <- eigen_exch(n, params[1:6] - params[7:12], directed, calcall=FALSE)$uniquevals
    e2 <- eigen_exch(n, params[1:6] + (tmax-1)*params[7:12], directed, calcall=FALSE)$uniquevals
    e <- c(e1,e2)
  } else {
    stop("eigen_exch_time(.) Only coded for types EI and EE")
  }
  
}