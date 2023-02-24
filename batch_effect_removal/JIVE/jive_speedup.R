library(tictoc)
library(r.jive)
library(RSpectra)
library(Rcpp)

sourceCpp("matrix_multiplication.cpp", showOutput = FALSE)
CORES <- parallel::detectCores()

####################################################################################################

jive_v2 <- function (data, rankJ = 1, rankA = rep(1, length(data)), method = "perm", 
                     dnames = names(data), conv = "default", maxiter = 1000, 
                     scale = TRUE, center = TRUE, orthIndiv = TRUE, est = TRUE, 
                     showProgress = TRUE) 
{
  # START PART 1
  l <- length(data)
  n <- c()
  d <- c()
  for (i in 1:length(data)) {
    n[i] <- nrow(data[[i]]) * ncol(data[[i]])
    d[i] <- nrow(data[[i]])
  }
  for (i in 1:l) {
    temp <- SVDmiss(data[[i]], ncomp = min(ncol(data[[i]]), 
                                           nrow(data[[i]])))[[1]]
    data[[i]] <- temp$u %*% diag(x = temp$d) %*% t(temp$v)
  }
  centerValues <- list()
  scaleValues <- c()
  for (i in 1:l) {
    if (center) {
      centerValues[[i]] <- apply(data[[i]], 1, mean, na.rm = T)
      data[[i]] <- data[[i]] - matrix(rep(centerValues[[i]], 
                                          ncol(data[[i]])), nrow = nrow(data[[i]]))
    }
    if (!center) {
      centerValues[[i]] <- rep(0, d[i])
    }
    if (scale) {
      scaleValues[i] <- norm(data[[i]], type = "f") * 
        sqrt(sum(n))
      data[[i]] <- data[[i]]/scaleValues[i]
    }
    if (!scale) {
      scaleValues[i] <- 1
    }
  }
  if (conv == "default") {
    conv = 10^(-6) * norm(do.call(rbind, data), type = "f")
  }
  orig <- data
  # END PART 1
  
  if (method == "given") {
    # START PART 2
    if (est) {
      u <- list()
      for (i in 1:l) {
        if (nrow(data[[i]]) > ncol(data[[i]])) {
          temp <- svdwrapper(data[[i]], nu = ncol(data[[i]]), 
                             nv = ncol(data[[i]]))
          data[[i]] <- diag(x = temp$d[1:ncol(data[[1]])], 
                            nrow = ncol(data[[1]])) %*% t(temp$v[, 1:ncol(data[[1]])])
          u[[i]] <- temp$u
        }
        else {
          u[[i]] <- diag(1, nrow(data[[i]]))
        }
      }
    }
    # END PART 2
    if (showProgress) {
      cat("Running JIVE algorithm for ranks:\njoint rank:", 
          rankJ, ", individual ranks:", rankA, "\n")
      cat("----------------------------------------------\n")
    }
    temp <- jive_iter(data, rankJ, rankA, conv = conv, maxiter = maxiter, 
                      orthIndiv = orthIndiv, showProgress = showProgress)
    joint <- temp$joint
    individual <- temp$individual
    if (est) {
      for (i in 1:l) {
        joint[[i]] <- u[[i]] %*% joint[[i]]
        individual[[i]] <- u[[i]] %*% individual[[i]]
      }
    }
  }

  if (is.null(dnames)) {
    for (i in 1:l) {
      names(orig)[[i]] <- paste("Source", i, sep = "_")
    }
  }
  else {
    names(orig) <- dnames
  }
  result <- list(data = orig, joint = joint, individual = individual, 
                 rankJ = rankJ, rankA = rankA, method = method)
  if (method == "bic") {
    result$bic.table <- bic.table
  }
  if (method == "perm") {
    result$converged <- converged
  }
  result$scale <- list(center, scale, centerValues, scaleValues)
  names(result$scale) <- c("Center", "Scale", "Center Values", 
                           "Scale Values")
  class(result) <- "jive"
  return(result)
}


jive_iter <- function (data, rankJ = 1, rankA = rep(1, length(data)), conv = 1e-06, 
                           maxiter = 1000, orthIndiv = TRUE, showProgress = TRUE) 
{
  # START PART 3
  l <- length(data) # Number of data sets being integrated
  A <- list() 
  for (i in 1:l) {
    A[[i]] <- matrix(0, nrow(data[[i]]), ncol(data[[i]]))
  }
  Xtot <- do.call(rbind, data)
  Jtot <- matrix(-1, nrow(Xtot), ncol(Xtot))
  Atot <- matrix(-1, nrow(Xtot), ncol(Xtot))
  nrun = 0
  converged = F
  if (orthIndiv) {
    Vind <- list()
  }
  while (nrun < maxiter & !converged) {
    cat(paste0("Iteration ", nrun, " begin\n"))
    tic(paste0("Iteration ", nrun))
    Jlast <- Jtot
    Alast <- Atot
    
    tic("  Joint")
    if (rankJ > 0) {
      temp <- Xtot - Atot
      
      # s <- svdwrapper(temp, nu = rankJ, nv = rankJ)
      s <- svds(temp, rankJ) #########################################################################
      
      # Jtot <- s$u[, 1:rankJ] %*% diag(x = s$d[1:rankJ], nrow = rankJ) %*% t(s$v[, 1:rankJ])
      # Jtot <- Tcrossprod(Tcrossprod(s$u, diag(s$d)), (s$v)) ##########################################
      Jtot <- eigenMapMatMult2(eigenMapMatMult2(s$u, diag(s$d), n_cores = CORES), t(s$v), n_cores = CORES)
      
      V <- s$v
    } else {
      Jtot <- matrix(0, nrow(Xtot), ncol(Xtot))
      V <- matrix(0, ncol(Xtot), rankJ)
    }
    temp <- Jtot
    J <- list()
    for (i in 1:l) {
      J[[i]] <- temp[1:nrow(data[[i]]), ]
      temp <- temp[-(1:nrow(data[[i]])), ]
    }
    toc()
    
    tic("  Indiv")
    A <- list()
    for (i in 1:l) {
      if (rankA[i] > 0) {
        
        # temp <- (data[[i]] - J[[i]]) %*% (diag(ncol(Xtot)) - V %*% t(V))
        # temp <- Tcrossprod(data[[i]] - J[[i]], diag(ncol(Xtot)) - Tcrossprod(V, V)) ##################
        temp <- eigenMapMatMult2(data[[i]] - J[[i]], diag(ncol(Xtot)) - eigenMapMatMult2(V, t(V), n_cores = CORES), n_cores = CORES)
        
        if (orthIndiv & nrun > 0) {
          for (j in (1:l)[-i]) {
            
            # temp <- temp %*% (diag(ncol(Xtot)) - Vind[[j]] %*% t(Vind[[j]]))
            # temp <- Tcrossprod(temp, diag(ncol(Xtot)) - Tcrossprod(Vind[[j]], Vind[[j]])) ############
            temp <- eigenMapMatMult2(temp, diag(ncol(Xtot)) - eigenMapMatMult2(Vind[[j]], t(Vind[[j]]), n_cores = CORES), n_cores = CORES)
            
          }
        }
        
        # s <- svdwrapper(temp, nu = rankA[i], nv = rankA[i])
        s <- svds(temp, rankA[i]) ####################################################################
        
        if (orthIndiv) {
          Vind[[i]] <- s$v[, 1:rankA[i]]
        }
        
        # A[[i]] <- s$u[, 1:rankA[i]] %*% diag(x = s$d[1:rankA[i]], nrow = rankA[i]) %*% t(s$v[, 1:rankA[i]])
        # A[[i]] <- Tcrossprod(Tcrossprod(s$u[, 1:rankA[i]], diag(x = s$d[1:rankA[i]], nrow = rankA[i])), s$v[, 1:rankA[i]]) ##################################
        A[[i]] <- eigenMapMatMult2(eigenMapMatMult2(s$u[, 1:rankA[i]], diag(x = s$d[1:rankA[i]], nrow = rankA[i]), n_cores = CORES), t(s$v[, 1:rankA[i]]), n_cores = CORES)
        
      }
      else {
        A[[i]] <- matrix(0, nrow(data[[i]]), ncol(data[[i]]))
        if (orthIndiv) {
          Vind[[i]] <- matrix(0, ncol(Xtot), rankA[[i]])
        }
      }
    }
    toc()
    
    if (orthIndiv & nrun == 0) {
      for (i in 1:l) {
        for (j in (1:l)[-i]) {
          
          # A[[i]] <- A[[i]] %*% (diag(ncol(Xtot)) - Vind[[j]] %*% t(Vind[[j]]))
          # A[[i]] <- Tcrossprod(A[[i]], diag(ncol(Xtot)) - Tcrossprod(Vind[[j]], Vind[[j]])) ##########
          A[[i]] <- eigenMapMatMult2(A[[i]], diag(ncol(Xtot)) - eigenMapMatMult2(Vind[[j]], t(Vind[[j]]), n_cores = CORES), n_cores = CORES)
          
        }
      }
      for (i in 1:l) {
        if (rankA[i] > 0) {
          
          # s <- svdwrapper(A[[i]], nu = rankA[i], nv = rankA[i])
          s <- svds(A[[i]], rankA[i]) ################################################################
          
          Vind[[i]] <- s$v[, 1:rankA[i]]
        }
      }
    }
    Atot <- do.call(rbind, A)
    
    normJ <- norm(Jtot - Jlast, type = "f")
    normA <- norm(Atot - Alast, type = "f")
    if (normJ <= conv & normA <= conv) {
      converged <- T
    }
    nrun = nrun + 1
    cat("    Convergence Target: ", conv, "\n")
    cat("    normJ: ", normJ, ", normA: ", normA, "\n")
    toc()
    cat("----------------------------------------------\n")
  }
  if (showProgress) {
    if (converged) {
      cat(paste("JIVE algorithm converged after ", nrun, 
                " iterations.\n"))
    }
    else {
      cat(paste("JIVE algorithm did not converge after ", 
                nrun, " iterations.\n"))
    }
  }
  return(list(data = data, joint = J, individual = A, rankJ, 
              rankA, method = "given"))
}

# tic("JIVE update")
# JIVE_results <- jive_v2(sim1_data, rankJ = 5, rankA = c(7, 7), method = "given", maxiter = 2000)
# toc()
# 
# saveRDS(JIVE_results, file = "data/JIVE_dataset1_speedup2.rds")
