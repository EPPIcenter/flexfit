vtoA <- function(ys, v, b, xfix, ifix) {
  bv <- b*(1 - v^(xfix[2] - xfix[1]))/(1 - v^(xfix[3] - xfix[1]))
  #  w <- (bv - 1)/(v^xfix[2]*(1 - bv*v^(xfix[3] - xfix[2])))
  w <- (bv - 1)/(v^xfix[2] - bv*v^xfix[3])
  A <- ys[ifix[2]] * (1 + w*v^xfix[1])*(1 + w*v^xfix[2]) /
    (w*(v^xfix[1] - v^xfix[2]))
  return(list(A = A, w = w))  # list append names later
}

getStart4par <- function(x, y, ifix = NULL, nv = 10, info = "") {
  y <- tapply(y, x, mean, names = NULL)  # tapply output is ordered by x
  names(y) <- NULL
  x <- sort(unique(x))
  n <- length(x)
  if (is.null(ifix)) ifix <- c(1, ceiling(n/2), n)
  yshift <- y - y[ifix[1]]  # first subtracted
  b <- yshift[ifix[3]]/yshift[ifix[2]]
  if (!is.finite(b)) {
    return(NA)
  }
  ivar <- setdiff(1:n, ifix)
  xfix <- x[ifix]

  ss <- numeric(nv - 1)
  vval <- seq(0, 1, length = nv + 1)[-c(1, nv + 1)]
  for(i in 1:(nv - 1)) {
    v <- vval[i]
    r <- vtoA(yshift, v, b, xfix, ifix)
    ss[i] <- sum((yshift[ivar] - r$A*(1/(1 + r$w*v^(x[ivar])) -
                                          1/(1 + r$w*v^(xfix[1]))))^2)
  }
  imin <- which.min(ss)
  # check fit: SS (sum of squares) too large
  if (ss[imin] >= (n - 3)*(mean(abs(diff(y)))/2)^2) {
    return(NA)
  }
  if (imin %in% c(1, nv - 1)) {
    vmin <- vval[imin]                                       # not a good fit
  } else {
    aa <- (ss[imin + 1] + ss[imin - 1] - 2*ss[imin]) * 5000  # a in quad eqn
    vmin <- vval[imin] - (ss[imin + 1] - ss[imin - 1])*25/aa
  }
  rmin <- vtoA(yshift, vmin, b, xfix, ifix)

  #***=========== optional update to provide fit even if A, w <= 0 =========***#
  if (rmin$w <= 0) {    # then A <= 0
    warning(paste(info, "Potential fit problem, manual check recommended"))
    high <- imin >= nv/2
    while(rmin$w <= 0) {
      imin <- high*(imin - 1) + (!high)*(imin + 1)  # increment/decrement imin
      vmin <- vval[imin]
      rmin <- vtoA(yshift, vmin, b, xfix, ifix)
    }
  }
  #============================================================================#

  if (rmin$w <= 0) {  #*** possibly unnecessary with the update?
    return(NA)
  }
  A0 <- y[ifix[1]] - rmin$A/(1 + rmin$w*vmin^(xfix[1]))
  # optional add mean error to A0
  #  A0 <- A0 + mean(y - (A0 + rmin$A/(1 + rmin$w*vmin^x)))
  s <- -1/log(vmin)
  X0 <- s*log(rmin$w)
  return(c(Aup = A0 + rmin$A, Alow = A0, Xmid = X0, Scale = s))
}

# Equal distance between three pts in X required
getvwEql <- function(yshift, c21, c31, x3) {
  v = (c31 - c21)/(c31*(c21 - 1))
  w <- (c21 - 1)/(1 - c21*v)
  s <- -(x3[2] - x3[1])/log(v)
  return(list(v = v, w = w, s = s))
}

# General: does not require equal distance in X (uses Newton's method for root)
g0 <- function(v, c21, c31, x21, x31) {
  (c21*c31 - c31)*v^x31 - (c21*c31 - c21)*v^x21 + c31 - c21
}
# first derivative
g1 <- function(v, c21, c31, x21, x31) {
  x31*(c21*c31 - c31)*v^(x31 - 1) - x21*(c21*c31 - c21)*v^(x21 - 1)
}
getvwGen <- function(yshift, c21, c31, x3, nv = 10, niter = 3) {
  x21 <- x3[2] - x3[1]
  x31 <- x3[3] - x3[1]
  rr <- numeric(nv - 1)
  vval <- seq(0, 1, length = nv + 1)[-c(1, nv + 1)]
  for (i in 1:(nv - 1)) {
    v <- vval[i]
    rr[i] <- g0(v, c21, c31, x21, x31)/(1 - v)
  }
  vmin <- vval[which.min(abs(rr))]
  # Newton's method
  for (k in 1:niter) {
    vmin <- vmin - g0(vmin, c21, c31, x21, x31)/g1(vmin, c21, c31, x21, x31)
  }
  w <- (c21 - 1)/(1 - c21*vmin^x21)
  s <- -1/log(vmin)
  return(list(v = vmin, w = w, s = s))
}

# optionally: when called, use ifix = c(1, floor(n/2) + 1 , n)
getStart3par <- function(x, y, A0, a = 1, ifix = NULL, nv = 10, niter = 3) {
  y <- tapply(y, x, mean)  # tapply output is ordered by x
  names(y) <- NULL
  x <- sort(unique(x))
  n <- length(x)
  equal <- FALSE
  if (!is.null(ifix)) {
    if (x[ifix[2]] - x[ifix[1]] == x[ifix[3]] - x[ifix[2]]) {
      equal <- TRUE
    }
  } else {
    mid <- ceiling(n/2)
    ifix <- c(1, mid, mid*2 - 1)
    if (x[ifix[2]] - x[ifix[1]] == x[ifix[3]] - x[ifix[2]]) {
      equal <- TRUE
    } else if (ifix[3] != n) {                # n is even, thus another option
      ifix <- ifix + 1
      if (x[ifix[2]] - x[ifix[1]] == x[ifix[3]] - x[ifix[2]]) {
        equal <- TRUE
      } else {
        ifix[1] <- 1                          # to include endpoints
      }
    }
  }
  x3 <- x[ifix]
  y3 <- y[ifix]
  yshift <- y3 - A0
  c21 <- (yshift[2]/yshift[1])^(1/a)
  c31 <- (yshift[3]/yshift[1])^(1/a)
  if (!all(diff(c(1, c21, c31)) > 0)) {
    return(NA)
  }
  if (equal) {
    r <- getvwEql(yshift, c21, c31, x3)
  } else {
    r <- getvwGen(yshift, c21, c31, x3, nv, niter)
  }
  if (r$w <= 0) {
    return(NA)
  }
  X0 <- r$s*log(r$w) + x3[1]
  A <- yshift[1]*(1 + r$w)^a
  return(c(Aup = A0 + A, Scale = r$s, Xmid = X0))
}
