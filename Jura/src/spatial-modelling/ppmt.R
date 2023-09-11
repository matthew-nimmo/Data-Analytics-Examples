library(fastICA)
library(gstat)
library(sp)

nscore <- function(x) {
  y <- qqnorm(x, plot.it=FALSE)$x
  trn.table <- data.frame(x=sort(x), nscore=sort(y))
  y <- list(nscore=y, trn.table=trn.table)
  return (y)
}

nscore2 <- function(x) {
  y <- qqnorm(x, plot.it=FALSE)$x
  return (y)
}

nsbacktr <- function(scores, nscore, tails="none", draw=FALSE) {
  if(tails == "separate") { 
    mean.x <- mean(nscore$trn.table$x)
    small.x <- nscore$trn.table$x < mean.x
    large.x <- nscore$trn.table$x > mean.x
    small.sd <- sqrt(sum((nscore$trn.table$x[small.x]-mean.x)^2) / (length(nscore$trn.table$x[small.x])-1))
    large.sd <- sqrt(sum((nscore$trn.table$x[large.x]-mean.x)^2) / (length(nscore$trn.table$x[large.x])-1))
    min.x <- mean(nscore$trn.table$x) + (min(scores) * small.sd)
    max.x <- mean(nscore$trn.table$x) + (max(scores) * large.sd)
    # check to see if these values are LESS extreme than the
    # initial data - if so, use the initial data.
    #print(paste('lg.sd is:',large.sd,'max.x is:',max.x,'max nsc.x is:',max(nscore$trn.table$x)))
    if(min.x > min(nscore$trn.table$x)) {
      min.x <- min(nscore$trn.table$x)
    }
    if(max.x < max(nscore$trn.table$x)) {
      max.x <- max(nscore$trn.table$x)
    }
  }

  if(tails == "equal") { # assumes symmetric distribution around the mean
    mean.x <- mean(nscore$trn.table$x)
    sd.x <- sd(nscore$trn.table$x)
    min.x <- mean(nscore$trn.table$x) + (min(scores) * sd.x)
    max.x <- mean(nscore$trn.table$x) + (max(scores) * sd.x)
    # check to see if these values are LESS extreme than the
    # initial data - if so, use the initial data.
    if(min.x > min(nscore$trn.table$x)) {
      min.x <- min(nscore$trn.table$x)
    }
    if(max.x < max(nscore$trn.table$x)) {
      max.x <- max(nscore$trn.table$x)
    }
  }

  if(tails == "none") {   # No extrapolation
    min.x <- min(nscore$trn.table$x)
    max.x <- max(nscore$trn.table$x)
  }

  min.sc <- min(scores)
  max.sc <- max(scores)
  x <- c(min.x, nscore$trn.table$x, max.x)
  nsc <- c(min.sc, nscore$trn.table$nscore, max.sc)

  if(draw) {
    plot(nsc, x, main='Transform Function')
  }

  back.xf <- approxfun(nsc, x) # Develop the back transform function
  val <- back.xf(scores)

  return(val)
}


ppmt <- function(x) {
  # Zscore 1.
  z.norms1 <- lapply(x, nscore)
  z.nscore1 <- do.call("cbind", lapply(z.norms1, function(o) o$nscore))

  # Projection pursuit.
  z.pp <- fastICA(z.nscore1, ncol(x), alg.typ="deflation", fun="logcosh", alpha=1, method="C",
                  row.norm=FALSE, maxit=200, tol=1e-6, verbose=FALSE)
  z.S <- as.data.frame(z.pp$S)
  colnames(z.S) <- colnames(z.nscore1)

  # Zscore 2.
  z.norms2 <- lapply(z.S, nscore)
  z.nscore2 <- do.call("cbind", lapply(z.norms2, function(o) o$nscore))

  y <- list(z.norms1=z.norms1, z.pp=z.pp, z.norms2=z.norms2, Y=z.nscore2)

  return(y)
}

ppmt_inv <- function(x, ppmt) {
  # Zscore 2 inverse.
  y <- sapply(1:ncol(x), function(i) nsbacktr(x[,i], ppmt$z.norms2[[i]], tails="none"))

  # Projection pursuit inverse.
  y <- y %*% ppmt$z.pp$A

  # Zscore 1 inverse.
  y <- sapply(1:ncol(y), function(i) nsbacktr(y[,i],  ppmt$z.norms1[[i]], tails="none"))
}

