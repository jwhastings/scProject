jive_iter <- function (data, rankJ = 1, rankA = rep(1, length(data)), conv = 1e-06, 
                       maxiter = 1000, orthIndiv = TRUE, showProgress = TRUE, CORES = 1) 
{
    l <- length(data)
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
      s <- svds(temp, rankJ)
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
        temp <- eigenMapMatMult2(data[[i]] - J[[i]], diag(ncol(Xtot)) - eigenMapMatMult2(V, t(V), n_cores = CORES), n_cores = CORES)
        
        if (orthIndiv & nrun > 0) {
          for (j in (1:l)[-i]) {
            temp <- eigenMapMatMult2(temp, diag(ncol(Xtot)) - eigenMapMatMult2(Vind[[j]], t(Vind[[j]]), n_cores = CORES), n_cores = CORES)
          }
        }
        s <- svds(temp, rankA[i])
        if (orthIndiv) {
          Vind[[i]] <- s$v[, 1:rankA[i]]
        }
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
          A[[i]] <- eigenMapMatMult2(A[[i]], diag(ncol(Xtot)) - eigenMapMatMult2(Vind[[j]], t(Vind[[j]]), n_cores = CORES), n_cores = CORES)
        }
      }
      for (i in 1:l) {
        if (rankA[i] > 0) {
          s <- svds(A[[i]], rankA[i])
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
