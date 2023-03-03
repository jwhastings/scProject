jive.iter <- function (data, rankJ = 1, rankA = rep(1, length(data)), conv = 1e-06, 
                       maxiter = 1000, orthIndiv = TRUE, showProgress = TRUE) 
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
    Jlast <- Jtot
    Alast <- Atot
    timeJ <- system.time({
      if (rankJ > 0) {
        temp <- Xtot - Atot
        s <- svdwrapper(temp, nu = rankJ, nv = rankJ)
        Jtot <- s$u[, 1:rankJ] %*% diag(x = s$d[1:rankJ], 
                                        nrow = rankJ) %*% t(s$v[, 1:rankJ])
        V <- s$v[, 1:rankJ]
      }
      else {
        Jtot <- matrix(0, nrow(Xtot), ncol(Xtot))
        V <- matrix(0, ncol(Xtot), rankJ)
      }
      temp <- Jtot
      J <- list()
      for (i in 1:l) {
        J[[i]] <- temp[1:nrow(data[[i]]), ]
        temp <- temp[-(1:nrow(data[[i]])), ]
      }
    })
    timeA <- system.time({
      A <- list()
      for (i in 1:l) {
        if (rankA[i] > 0) {
          temp <- (data[[i]] - J[[i]]) %*% (diag(ncol(Xtot)) - 
                                              V %*% t(V))
          if (orthIndiv & nrun > 0) {
            for (j in (1:l)[-i]) {
              temp <- temp %*% (diag(ncol(Xtot)) - Vind[[j]] %*% 
                                  t(Vind[[j]]))
            }
          }
          s <- svdwrapper(temp, nu = rankA[i], nv = rankA[i])
          if (orthIndiv) {
            Vind[[i]] <- s$v[, 1:rankA[i]]
          }
          A[[i]] <- s$u[, 1:rankA[i]] %*% diag(x = s$d[1:rankA[i]], 
                                               nrow = rankA[i]) %*% t(s$v[, 1:rankA[i]])
        }
        else {
          A[[i]] <- matrix(0, nrow(data[[i]]), ncol(data[[i]]))
          if (orthIndiv) {
            Vind[[i]] <- matrix(0, ncol(Xtot), rankA[[i]])
          }
        }
      }
    })
    if (orthIndiv & nrun == 0) {
      for (i in 1:l) {
        for (j in (1:l)[-i]) {
          A[[i]] <- A[[i]] %*% (diag(ncol(Xtot)) - Vind[[j]] %*% 
                                  t(Vind[[j]]))
        }
      }
      for (i in 1:l) {
        if (rankA[i] > 0) {
          s <- svdwrapper(A[[i]], nu = rankA[i], nv = rankA[i])
          Vind[[i]] <- s$v[, 1:rankA[i]]
        }
      }
    }
    Atot <- do.call(rbind, A)
    if (norm(Jtot - Jlast, type = "f") <= conv & norm(Atot - 
                                                      Alast, type = "f") <= conv) {
      converged <- T
    }
    nrun = nrun + 1
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
