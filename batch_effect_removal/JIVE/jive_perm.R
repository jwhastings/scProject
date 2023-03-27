jive_perm <- function (data, nperms = 100, alpha = 0.05, est = TRUE, conv = 1e-06, 
          maxiter = 1000, orthIndiv = TRUE, showProgress = TRUE, CORES = CORES) 
{
  nrun <- 0
  Jperp <- list()
  Aperp <- list()
  for (i in 1:length(data)) {
    Jperp[[i]] <- matrix(0, nrow = nrow(data[[i]]), ncol = ncol(data[[i]]))
    Aperp[[i]] <- matrix(0, nrow = nrow(data[[i]]), ncol = ncol(data[[i]]))
  }
  last <- rep(-2, length(data) + 1)
  current <- rep(-1, length(data) + 1)
  while (!isTRUE(all.equal(last, current)) & nrun < 10) {
    last <- current
    if (showProgress) {
      if (nrun == 0) {
        cat("Estimating  joint and individual ranks via permutation...\n")
      }
      else {
        cat("Re-estimating  joint and individual ranks via permutation...\n")
      }
    }
    tic("Estimating joint ranks")
    full <- list()
    for (i in 1:length(data)) {
      full[[i]] <- data[[i]] - Aperp[[i]]
    }
    n <- ncol(full[[1]])
    
    # actual <- svdwrapper(do.call(rbind, full), nu = 0, nv = 0)$d
    tic("  SV for joint (actual)")
    actual <- eigenBDCSVD(do.call(rbind, full), n_cores = 8)$d
    toc(quiet = !showProgress)
    
    perms <- matrix(NA, nperms, min(n, sum(unlist(lapply(data, nrow)))))
    for (i in 1:nperms) {
      temp <- list()
      for (j in 1:length(data)) {
        temp[[j]] <- full[[j]][, sample(1:n, n, replace = F)]
      }

      # perms[i, ] <- svdwrapper(do.call(rbind, temp), nu = 0, nv = 0)$d
      tic(paste0("  SV for joint permutation ", i, " of ", nperms))
      perms[i, ] <- eigenBDCSVD(do.call(rbind, temp), n_cores = CORES)$d
      toc(quiet = !showProgress)
    }
    rankJ <- 0
    for (i in 1:n) {
      if (actual[i] > quantile(perms[, i], 1 - alpha)) {
        rankJ <- rankJ + 1
      }
      else {
        break
      }
    }
    rankJ <- max(rankJ, last[1])
    toc(quiet = !showProgress)
    
    tic("Estimating Individual Ranks")
    rankA <- c()
    for (i in 1:length(data)) {
      ind <- data[[i]] - Jperp[[i]]
      ev_len <- min(dim(ind))
      
      tic(paste0("  SV for individual ", i, " (actual)"))
      # actual <- svdwrapper(ind, nu = 0, nv = 0)$d
      actual <- eigenBDCSVD(ind, n_cores = CORES)$d
      toc(quiet = !showProgress)
      
      perms <- matrix(NA, nperms, min(n, nrow(data[[i]])))
      for (k in 1:nperms) {
        perm <- t(ind)
        pind <- order(c(col(perm)), runif(length(perm)))
        perm <- matrix(perm[pind], nrow = nrow(ind), ncol = n, byrow = TRUE)
        ev_len <- min(dim(perm))
        
        tic(paste0("  SV for individual ", i, " permutation ", k, " of ", nperms))
        # perms[k, ] <- svdwrapper(perm, nu = 0, nv = 0)$d
        perms[k, ] <- eigenBDCSVD(perm, n_cores = CORES)$d
        toc(quiet = !showProgress)
      }
      rankA[i] <- 0
      for (j in 1:n) {
        if (actual[j] > quantile(perms[, j], 1 - alpha)) {
          rankA[i] <- rankA[i] + 1
        }
        else {
          break
        }
      }
    }
    toc()
    current <- c(rankJ, rankA)
    if (!isTRUE(all.equal(last, current))) {
      dataR <- list()
      if (est) {
        u <- list()
        for (i in 1:length(data)) {
          if (nrow(data[[i]]) > ncol(data[[i]])) {
            temp <- svds(data[[i]], k = ncol(data[[i]]))
            dataR[[i]] <- diag(x = temp$d[1:ncol(data[[1]])], nrow = ncol(data[[1]])) %*% t(temp$v[, 1:ncol(data[[1]])])
            u[[i]] <- temp$u
          }
          else {
            u[[i]] <- diag(1, nrow(data[[i]]))
            dataR[[i]] <- data[[i]]
          }
        }
      }
      else {
        dataR <- data
      }
      if (showProgress) {
        cat("Running JIVE algorithm for ranks:\njoint rank:", 
            rankJ, ", individual ranks:", rankA, "\n")
      }
      tempjive <- jive_iter(dataR, rankJ, rankA, conv = conv, 
                            maxiter = maxiter, orthIndiv = orthIndiv, showProgress = showProgress, CORES = CORES)
      Jperp = tempjive$joint
      Aperp = tempjive$individual
      if (est) {
        for (i in 1:length(data)) {
          Jperp[[i]] <- u[[i]] %*% Jperp[[i]]
          Aperp[[i]] <- u[[i]] %*% Aperp[[i]]
        }
      }
    }
    nrun <- nrun + 1
  }
  converged <- ifelse(nrun == 10, F, T)
  if (showProgress) {
    cat("Final joint rank:", rankJ, ", final individual ranks:", 
        rankA, "\n")
  }
  return(list(data = data, joint = Jperp, individual = Aperp, 
              rankJ = rankJ, rankA = rankA, converged = converged))
}
