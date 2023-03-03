bic.jive <- function (data, n = unlist(lapply(data, ncol)) * unlist(lapply(data, nrow)),
                      d = unlist(lapply(data, nrow)), conv = 1e-06, maxiter = 1000, 
                      orthIndiv = TRUE, showProgress = TRUE) 
{
  l <- length(data)
  nc <- ncol(data[[1]])
  lambda <- log(sum(n))
  sse <- c()
  Xtot <- do.call(rbind, data)
  bic.improve <- T
  rankJ <- 0
  rankA <- rep(0, l)
  if (showProgress) {
    cat("Running JIVE algorithm for ranks:\njoint rank:", 
        rankJ, ", individual ranks:", rankA, "\n")
  }
  current <- jive.iter(data, rankJ, rankA, conv, maxiter, 
                       orthIndiv, showProgress = showProgress)
  for (k in 1:length(data)) {
    sse[k] <- norm(data[[k]] - current$joint[[k]] - current$individual[[k]], 
                   type = "f")^2
  }
  p.jive <- 0
  current.bic <- sum(n * log(sse/n)) + p.jive * lambda
  bic.table <- c(rankJ, rankA, current.bic)
  while (bic.improve) {
    bic.improve <- F
    temp <- list()
    if (showProgress) {
      cat("Running JIVE algorithm for ranks:\njoint rank:", 
          rankJ + 1, ", individual ranks:", rankA, "\n")
    }
    temp[[1]] <- jive.iter(data, rankJ + 1, rankA, conv, 
                           maxiter, orthIndiv, showProgress = showProgress)
    for (k in 1:length(data)) {
      sse[k] <- norm(data[[k]] - temp[[1]]$joint[[k]] - 
                       temp[[1]]$individual[[k]], type = "f")^2
    }
    p.jive <- sum(sum(d):(sum(d) - (rankJ + 1) + 1)) + sum(nc:(nc - 
                                                                 (rankJ + 1) + 1)) + pjsum(d, rankA) + pjsum(rep(nc, 
                                                                                                                 length(data)) - (rankJ + 1), rankA)
    bic <- sum(n * log(sse/n)) + p.jive * lambda
    bicvec <- bic
    bic.table <- rbind(bic.table, c(rankJ + 1, rankA, bic))
    for (i in 1:l) {
      tempR <- rankA
      tempR[i] <- tempR[i] + 1
      if (tempR[i] < min(n, nrow(data[[i]]))) {
        if (showProgress) {
          cat("Running JIVE algorithm for ranks:\njoint rank:", 
              rankJ, ", individual ranks:", tempR, "\n")
        }
        temp[[i + 1]] <- jive.iter(data, rankJ, tempR, 
                                   conv, maxiter, orthIndiv, showProgress = showProgress)
        for (k in 1:length(data)) {
          sse[k] <- norm(data[[k]] - temp[[i + 1]]$joint[[k]] - 
                           temp[[i + 1]]$individual[[k]], type = "f")^2
        }
        p.jive <- ifelse(rankJ == 0, 0, sum(sum(d):(sum(d) - 
                                                      rankJ + 1)) + sum(nc:(nc - rankJ + 1))) + 
          pjsum(d, tempR) + pjsum(rep(nc, length(data)) - 
                                    rankJ, tempR)
        bic <- sum(n * log(sse/n)) + p.jive * lambda
      }
      else {
        bic <- NA
      }
      bicvec <- c(bicvec, bic)
      bic.table <- rbind(bic.table, c(rankJ, tempR, bic))
    }
    lowest.bic <- temp[[which.min(bicvec)]]
    if (min(bicvec, na.rm = T) < current.bic) {
      bic.improve <- T
      current <- lowest.bic
      current.bic <- min(bicvec, na.rm = T)
      if (which.min(bicvec) == 1) {
        rankJ <- rankJ + 1
      }
      else {
        rankA[which.min(bicvec) - 1] <- rankA[which.min(bicvec) - 
                                                1] + 1
      }
    }
  }
  return(list(data = data, joint = current$joint, individual = current$individual, 
              rankJ = rankJ, rankA = rankA, bic.table))
}
