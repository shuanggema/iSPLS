#=============================================================================================
# Program: iSPLS_Functions.R
# Note   : code submission for paper "Integrative Sparse Partial Least Squares".
# In this file:
# ispls.dv,  c_value_homo, c_value_hetero,ro, ro_d1st define the functions for ispls
#

#=============================================================================================
#=============================     ispls         =============================================  
#=============================================================================================
# Function: ispls 
# Source  : ispls.dv,  c_value_homo, c_value_hetero,ro, ro_d1st
# Input   :
#                x : matrix of explanatory variables for L datasets
#                y : matrix of dependent variables for L datasets
#                L : numeric,number of datasets;
#              mu1 : numeric, L1 penalty parameter; 
#              mu3 : numeric, contrasted penalty parameter
#            kappa : 0<kappa<0.5;
#            type1 : character, "homo" or "hetero" type of structure
#            type2 : character, "m" or "s" type of contrasted penalty
#          scale.x : character, "TRUE" or "FALSE", scale x or not
#          scale.y : character, "TRUE" or "FALSE", scale y ot not
#          maxstep : numeric, maximum iteration steps
# Output  : A list of results including the estimated loadings and coefficients

ispls <- function(x, y, L, mu1, mu3, kappa, type1, type2, scale.x = TRUE, scale.y = FALSE, maxstep = 50, trace = FALSE) {
 
  # initialization
  
  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(x)
  p <- ncol(x)/L
  q <- ncol(y)/L
  ip <- c(1:p)
  one <- matrix(1, 1, n)
  
  # center & scale x & y
  
  meany <- drop(one %*% y / n)
  y <- scale(y, meany, FALSE)
  meanx <- drop(one %*% x) / n
  x <- scale(x, meanx, FALSE)
  
  if (scale.x) {
    normx <- sqrt(drop(one %*% (x^2)) / (n - 1))
    if (any(normx < .Machine$double.eps)) {
      stop("Some of the columns of the predictor matrix have zero variance.")
    }
    x <- scale(x, FALSE, normx)
  } else {
    normx <- rep(1, p)
  }
  
  
  if (scale.y) {
    normy <- sqrt(drop(one %*% (y^2)) / (n - 1))
    if (any(normy < .Machine$double.eps)) {
      stop("Some of the columns of the response matrix have zero variance.")
    }
    y <- scale(y, FALSE, normy)
  } else {
    normy <- rep(1, q)
  }
  
  # initilize objects
  
  what <- matrix(0, p, L)
  
  # main iteration
  
  if (is.null(colnames(x))) {
    xnames <- rep(c(1:p), L)
  } else {
    xnames <- colnames(x)
  }
  
  if (trace) {
    cat("The variables that join the set of selected variables at each step:\n")
  }
  
  # define Z
  
  fun.1 <- function(l) {
    Z_l <- t(x[, ((l - 1) * p + 1):(l * p)]) %*% y[, ((l - 1) * q + 1):(l * q)]
    Z_l <- Z_l / n
  }
  Z <- matrix(mapply(fun.1, c(1:L)), nr = p)
  
  # fit direction vector
  
  what <- ispls.dv(Z, p, q, n,  L, mu1, mu3, kappa, type1, type2, maxstep)
  what_cut <- ifelse(abs(what) > 10^(-5), what, 0)
  
  # normalization
  
  what_cut_norm <- sqrt(colSums(what_cut^2)) + 0.0001
  what_dir <- t(t(what_cut) / what_cut_norm)
  what <- ifelse(abs(what_dir) > 10^(-5), what_dir, 0)
  
  # selected variables
  
  new2A <- list()
  new2A <- mapply(function(l) {
    A <- ip[what[, l] != 0]
    return(A)
  }, c(1:L), SIMPLIFY = FALSE)
  
  # fit y with component t=xw
  
  betahat <- matrix(0, nr = p, nc = q * L)
  
  fun.fit <- function(l) {
    x_l <- x[, ((l - 1) * p + 1):(l * p), drop = FALSE]
    w_l <- what[, l]
    t_l <- x_l %*% w_l
    
    if (sum(w_l == 0) != p) {
      y_l <- y[, ((l - 1) * q + 1):(l * q), drop = FALSE]
      fit_l <- lm(y_l ~ t_l - 1)
      betahat_l <- matrix(w_l %*% coef(fit_l), nr = p, nc = q)
    }
    else {
      betahat_l <- matrix(0, nr = p, nc = q)
    }
    return(betahat_l)
  }
  
  betahat <- matrix(mapply(fun.fit, c(1:L)), nr = p, nc = q * L)
  
  if (!is.null(colnames(x))) {
    rownames(betahat) <- c(1:p)
  }
  if (q > 1 & !is.null(colnames(y))) {
    colnames(betahat) <- colnames(y)
  }
  
  # print out variables that join the active set
  
  if (trace) {
    for (l in 1:L)
    {
      new2A_l <- new2A[[l]]
      if (length(new2A_l) <= 10) {
        cat(paste("DataSet ", l, "- ", ":\n", sep = ""))
        cat(paste("X", new2A_l, ", ", sep = " "))
        cat("\n")
      } else {
        cat(paste("DataSet ", l, "- ", ":\n", sep = ""))
        nlines <- ceiling(length(new2A_l) / 10)
        for (i in 0:(nlines - 2))
        {
          cat(paste("X", new2A_l[(10 * i + 1):(10 * (i + 1))], ", ", sep = " "))
          cat("\n")
        }
        cat(paste("X", new2A_l[(10 * (nlines - 1) + 1):length(new2A_l)], ", ", sep = " "))
        cat("\n")
      }
    }
  }
  
  # return objects
  object <- list(
    x = x, y = y, betahat = betahat, loadings = what, new2A = new2A,
    type1 = type1, type2 = type2, mu1 = mu1, mu3 = mu3, kappa = kappa,
    meanx = meanx, normx = normx, meany = meany, normy = normy
  )
  class(object) <- "ispls"
  return(object)
}

# =============================================================================================
# =============================      ispls.dv      ============================================
# =============================================================================================
# Function: ispls.dv.  
#           Fit iSPLS direction vectors
# Source  : c_value_homo, c_value_hetero,ro, ro_d1st
# Input   :
#                Z : matrix 
#                p : numeric, number of independent variables
#                q : numeric, number of dependent variables
#                n : numeric, sample size
#                L : numeric, number of datasets
#              mu1 : numeric, L1 penalty parameter
#              mu3 : numeric, contrasted penalty parameter
#            kappa : 0<kappa<0.5
#            type1 : character, "homo" or "hetero" type of structure
#            type2 : character, "m" or "s" type of contrasted penalty
#          maxstep : numeric, maximum iteration steps
# Output  : A vector of estimated values

"ispls.dv" <- function(Z, p, q, n, L, mu1, mu3, kappa, type1, type2, maxstep) {
  # homo-magnitude
  
  # define M
  
  fun.2 <- function(l) M_l <- Z[, ((l - 1) * q + 1):(l * q)] %*% t(Z[, ((l - 1) * q + 1):(l * q)]) / q
  M <- matrix(mapply(fun.2, c(1:L)), nr = p)
  dis <- 10
  i <- 1
 
  # main iteration: optimize c and a iteratively
  
  kappa2 <- (1 - kappa) / (1 - 2 * kappa)
  
  # initial value for c(l) (outside the unit circle)
  
  c <- matrix(1 / p, nr = p, nc = L)
  a <- matrix(1 / p, nr = p, nc = L)
  
  while (dis > 1e-4 & i <= maxstep) {
    # optimize a(l) for fixed c(l)
    c.old <- c
    fun.3 <- function(l) {
      h <- function(lambda) {
        alpha <- solve(M[, ((l - 1) * p + 1):(l * p)] + lambda * diag(p)) %*% M[, ((l - 1) * p + 1):(l * p)] %*% c[, l]
        obj <- t(alpha) %*% alpha - 1 / kappa2^2
        return(obj)
      }
      
      # control size of M_l & c_l if too small
      
      while (h(1e-4) * h(1e+12) > 0) { # while( h(1e-4) <= 1e+5 )
        {
          M[, ((l - 1) * p + 1):(l * p)] <- 2 * M[, ((l - 1) * p + 1):(l * p)]
          c[, l] <- 2 * c[, l]
        }
      }
      
      # optimization
      
      lambdas <- uniroot(h, c(1e-4, 1e+12))$root
      a_l <- kappa2 * solve(M[, ((l - 1) * p + 1):(l * p)] + lambdas * diag(p)) %*% M[, ((l - 1) * p + 1):(l * p)] %*% c[, l]
      return(a_l)
    }
    
    a <- mapply(fun.3, c(1:L))
    
    # optimize c(l) for fixed a(l)
    
    if (type1 == "homo") c <- c_value_homo(Z, a, c, p, q, n, L, mu1, mu3, type2)
    if (type1 == "hetero") c <- c_value_hetero(Z, a, c, p, q, n, L, mu1, mu3, type2)
    
    # calculate discrepancy between a & c
    c_norm <- sqrt(colSums(c^2)) + 0.0001
    c <- t(t(c) / c_norm)
    dis <- max(abs(c - c.old))

    i <- i + 1
    if (sum(apply(c, 2, function(x) sum(x <= 0.0001) == p)) > 0) {
      cat("The value of mu1 is too large");
      break} # exists an l such that c(l)=0
  }
  c_norm <- sqrt(colSums(c^2)) + 0.0001
  c <- t(t(c) / c_norm)
  
  return(c)
}

# =============================================================================================
# =============================      cmat_value      ==========================================
# =============================================================================================
# Function: c_value_homo, c_value_hetero.  
#           Optimize c for fixed a under the homo/hetero structure
# Source  : ro, ro_d1st
# Input   :
#                Z : matrix 
#                a : matrix
#                c : matrix
#                p : numeric, number of independent variables
#                q : numeric, number of dependent variables
#                n : numeric, sample size
#                L : numeric, number of datasets
#              mu1 : numeric, L1 penalty parameter
#              mu3 : numeric, contrasted penalty parameter
#            type2 : character, "m" or "s" type of contrasted penalty
# Output  : A vector of estimated values

c_value_homo <- function(Z, a, c,  p, q, n, L, mu1, mu3, type2) {
    # compute s[j,l]
    fun.s <- function(j, l) {
      Z_l <- Z[, ((l - 1) * q + 1):(l * q)]
      a_l <- matrix(a[, l], nc = 1)
      c_j <- c[j, ]
      s1 <- t(a_l) %*% Z_l %*% t(matrix(Z_l[j, ], nr = 1)) / (n^2)
      if (type2 == "m") s2 <- mu3 * sum(c_j[-l])
      if (type2 == "s") {
        s2 <- mu3 * sum(mapply(function(x) {
          r <- x / sqrt(x^2 + 0.5)
          return(r)
        }, c_j[-l])) / sqrt(c_j[l]^2 + 0.5)
      }
      s <- 2 * (s1 + s2)
      return(s)
    }
    result.s <- mapply(fun.s, rep(c(1:p), times = L), rep(c(1:L), each = p))
    s <- matrix(result.s, nr = p, nc = L)
    
    # compute ro'(||c_j||,mu1,a)
    
    norm_c_j <- apply(c, 1, function(x) {
      return(sqrt(sum(x^2)))
    })
    ro_d <- ro_d1st(norm_c_j, mu1, 6)
    
    # compute c[j,l]
    fun.c <- function(j, l) {
      s_norm <- sqrt(sum(s[j, ]^2))
      if (type2 == "m") c <- (s_norm > ro_d[j]) * s[j, l] * (s_norm - ro_d[j]) / (2 * (1 / q + mu3 * (L - 1)) * s_norm)
      
      if (type2 == "s") c <- (s_norm > ro_d[j]) * s[j, l] * (s_norm - ro_d[j]) / (2 * (1 / q + mu3 * (L - 1) / (c[j, l]^2 + 0.5)) * s_norm)
      
      return(c)
    }
    c <- matrix(mapply(fun.c, rep(c(1:p), times = L), rep(c(1:L), each = p)), nr = p, nc = L)
    return(c)
}

c_value_hetero <- function(Z, a, c, p, q, n, L, mu1, mu3, type2) {
    # compute mu[j,l]
    fun.mu <- function(j, l) {
      c_j <- c[j, ]
      ro_j <- mapply(ro, abs(c_j), mu1, 6)
      s_ro <- sum(as.data.frame(ro_j[1, ]))
      mu_jl <- ro_d1st(s_ro, 1, 1 / 2 * L * 6 * mu1^2) * ro_d1st(abs(c_j[l]), mu1, 6)
      return(mu_jl)
    }
    result.mu <- mapply(fun.mu, rep(c(1:p), times = L), rep(c(1:L), each = p))
    mu <- matrix(result.mu, nr = p, nc = L)
    
    # compute s[j,l]
    
    fun.s <- function(j, l) {
      Z_l <- Z[, ((l - 1) * q + 1):(l * q)]
      a_l <- matrix(a[, l], nc = 1)
      c_j <- c[j, ]
      s1 <- t(a_l) %*% Z_l %*% t(matrix(Z_l[j, ], nr = 1)) / n
      if (type2 == "m") s2 <- mu3 * sum(c_j[-l])
      if (type2 == "s") {
        s2 <- mu3 * sum(mapply(function(x) {
          r <- x / sqrt(x^2 + 0.5)
          return(r)
        }, c_j[-l])) / sqrt(c_j[l]^2 + 0.5)
      }
      s <- 2 * (s1 + s2)
      return(s)
    }
    result.s <- mapply(fun.s, rep(c(1:p), times = L), rep(c(1:L), each = p))
    s <- matrix(result.s, nr = p, nc = L)
    
    # compute c[j,l]
    
    fun.c <- function(j, l) {
      if (type2 == "m") c <- sign(s[j, l]) * (abs(s[j, l]) > mu[j, l]) * (abs(s[j, l]) - mu[j, l]) / (2 * (1 / q + mu3 * (L - 1)))
      
      if (type2 == "s") c <- sign(s[j, l]) * (abs(s[j, l]) > mu[j, l]) * (abs(s[j, l]) - mu[j, l]) / (2 * (1 / q + mu3 * (L - 1) / (c[j, l]^2 + 0.5)))
      
      return(c)
    }
    c <- matrix(mapply(fun.c, rep(c(1:p), times = L), rep(c(1:L), each = p)), nr = p, nc = L)
    return(c)
}


# =============================================================================================
# =============================      MCP & 1st derivative      =================================
# =============================================================================================
# Function: ro, ro_d1st
#           MCP and its first derivative
# Input   :
#                x : vector 
#               mu : numeric, MCP penalty parameter
#            alpha : numeric, MCP penalty parameter
# Output  : Vector of values
ro <- function(x, mu, alpha) {
  f <- function(x) mu * (1 > x / (mu * alpha)) * (1 - x / (mu * alpha))
  r <- integrate(f, 0, x)
  return(r)
}

ro_d1st <- function(x, mu, alpha) {
  r <- mu * (1 > x / (mu * alpha)) * (1 - x / (mu * alpha))
  return(r)
}

