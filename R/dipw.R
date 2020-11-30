# set the mean of each column to zero
center_colmeans <- function(x) {
  xcenter <- colMeans(x)
  x - rep(xcenter, rep.int(nrow(x), ncol(x)))
}

# Quadratic program
# This function is used to solve the quadratic program in Step 3b, Algorithm 1.
# The function is constructed based on the mosek optimization solver. The
# code is modified from the balanceHD package available at

# https://github.com/swager/balanceHD

# and licensed under GPL-3

# This implements a method from the paper:
# Athey, S., Imbens, G. W., & Wager, S. (2018). Approximate residual balancing:
# debiased inference of average treatment effects in high dimensions. Journal of
# the Royal Statistical Society Series B, 80(4), 597-623.

#' @import methods
#' @import Matrix
quad.prog = function(M,
                             target,
                             kappa = 0.5,
                             verbose = FALSE) {

  # The primal problem is:
  # Minimize 1/2 x' diag(qvec) x
  # subject to Ax <= b,
  # where the constraints indexed by "equality" are required to be equalities
  #
  # Here, the vector x is interpreted as (max.imbalance, mu)
  qvec.primal = 2 * c(kappa, rep(1 - kappa, nrow(M)))

  tM = Matrix::t(M)
  nvar = length(qvec.primal)

  # A.primal.list = list(
  #   matrix(c(0, rep(1, nrow(M))), 1, nvar),
  #   cbind(rep(-1, ncol(M)), tM),
  #   cbind(rep(-1, ncol(M)), -tM)
  # )
  # A.primal = Reduce(rbind, A.primal.list)
  A.primal <- rbind(
    matrix(c(0, rep(1, nrow(M))), 1, nvar),
    cbind(rep(-1, ncol(M)), tM),
    cbind(rep(-1, ncol(M)), -tM)
  )

  bvec = c(1,
           target,
           -target)

  equality.primal = c(TRUE, rep(FALSE, 2*ncol(M)))

  # The dual problem is then
  # Minimize 1/2 lambda A diag(1/qvec) A' lambda + b * lambda
  # Subject to lambda >= 0 for the lambdas corresponding to inequality constraints

  A.dual = 1/sqrt(qvec.primal) * Matrix::t(A.primal)

  # Next we turn this into a conic program:
  # Minimize t
  # Subject to bvec * lambda + q - t = 0
  #	A.dual * lambda - mu = 0
  #	lambda >= 0 (for the inequality constraints)
  #	2 * q >= ||mu||_2^2
  # Where we note that the last constraint is a rotated cone.
  #
  # Below, the solution vector is (lambda, mu, q, t, ONE), where ONE
  # is just a variable constrained to be 1.

  A.conic = rbind(c(bvec, rep(0, nvar), 1, -1, 0),
                          cbind(A.dual, Matrix::diag(-1, nvar, nvar), 0, 0, 0),
                          c(rep(0, length(bvec) + nvar + 2), 1))
  rhs.conic = c(rep(0, 1 + nvar), 1)

  blx.conic = rep(-Inf, ncol(A.conic))
  blx.conic[which(!equality.primal)] = 0
  bux.conic = rep(Inf, ncol(A.conic))

  obj.conic = c(rep(0, length(bvec) + nvar + 1), 1, 0)

  mosek.problem <- list()
  mosek.problem$sense <- "min"
  mosek.problem$c <- obj.conic
  mosek.problem$bx <- rbind(blx = blx.conic, bux = bux.conic)
  mosek.problem$A <- as(A.conic, "CsparseMatrix")
  mosek.problem$bc <- rbind(blc = rhs.conic, buc = rhs.conic)
  mosek.problem$cones <- cbind(list("RQUAD", c(length(bvec) + nvar + 1, length(bvec) + nvar + 3, length(bvec) + 1:nvar)))

  if (verbose) {
    mosek.out = Rmosek::mosek(mosek.problem)
  } else {
    mosek.out = Rmosek::mosek(mosek.problem, opts=list(verbose=0))
  }

  primal = -1/qvec.primal * (t(A.primal) %*% mosek.out$sol$itr$xx[1:nrow(A.primal)])
  delta = primal[1]
  mu = primal[1 + 1:nrow(M)]
  mu
}

# algorithm with B = 1 and deterministic split
dripw.qp.ate <- function(X, Y, W, prop.raw, kappa=0.5, ...){

  # start estimating the mu
  nfolds <- 2 #set up the number of folds
  folds <- cut(seq(1,nrow(X)),breaks=nfolds,labels=FALSE)
  mu <- rep(0, nrow(X))
  yy <- W * Y * (1 - prop.raw) / prop.raw + (1 - W) * Y * prop.raw / (1 - prop.raw)
  beta.fit <- cv.glmnet(X, yy, ...)
  # fit the residual with lasso
  residual <- yy - predict(beta.fit, newx=X)
  for(i in 1:nfolds){
    idx <- which(folds==i,arr.ind=TRUE)
    Xy <- t(center_colmeans(X[-idx,])) %*% (residual[-idx] - mean(residual[-idx])) / (nrow(X) - length(idx))
    mu[idx] <- quad.prog(center_colmeans(X[idx,]), Xy, kappa=kappa)
    mu[idx] <- mu[idx] * length(idx) + mean(residual[idx])
    #mu[idx] <- mu[idx] * length(idx)
  }

  tau <- sum((Y[W==1] - predict(beta.fit, newx=X[W==1,]) - mu[W==1]) / prop.raw[W==1]) / sum(1 / prop.raw[W==1])
  tau <- tau - sum((Y[W==0] - predict(beta.fit, newx=X[W==0,]) - mu[W==0]) / (1 - prop.raw[W==0])) / sum(1 / (1 - prop.raw[W==0]))
  tau
}

#algorithm with random splits
dripw.random.ate <- function(X, Y, W, prop.raw, kappa=0.5, B=1, ...){

  #use lasso to remove the linear signal
  yy <- W * Y * (1 - prop.raw) / prop.raw + (1 - W) * Y * prop.raw / (1 - prop.raw)
  beta.fit <- cv.glmnet(X, yy, ...)
  residual <- yy - predict(beta.fit, newx=X)

  #start estimating the mu
  nfolds <- 2 #set up the number of folds
  tau <- 0

  #start random splitting
  for (b in 1:B) {
    folds <- rep(0, nrow(X))
    idx <- sample.int(nrow(X), size=as.integer(nrow(X) / 2))
    folds[idx] <- 1
    folds[folds == 0] <- 2
    mu <- rep(0, nrow(X))
    for(i in 1:nfolds){
      idx <- which(folds==i,arr.ind=TRUE)
      Xy <- t(center_colmeans(X[-idx,])) %*% (residual[-idx] - mean(residual[-idx])) / (nrow(X) - length(idx))
      mu[idx] <- quad.prog(center_colmeans(X[idx,]), Xy, kappa=kappa)
      mu[idx] <- mu[idx] * length(idx) + mean(residual[idx])
      #mu[idx] <- mu[idx] * length(idx)
    }

    tau <- tau + sum((Y[W==1] - predict(beta.fit, newx=X[W==1,]) - mu[W==1]) / prop.raw[W==1]) / sum(1 / prop.raw[W==1])
    tau <- tau - sum((Y[W==0] - predict(beta.fit, newx=X[W==0,]) - mu[W==0]) / (1 - prop.raw[W==0])) / sum(1 / (1 - prop.raw[W==0]))
  }
  tau / B
}

# algorithm with B = 3 and deterministic split
dripw.split.ate <- function(X, Y, W, prop.raw, kappa=0.5, ...){

  # use lasso to remove the linear signal
  yy <- W * Y * (1 - prop.raw) / prop.raw + (1 - W) * Y * prop.raw / (1 - prop.raw)
  beta.fit <- cv.glmnet(X, yy, ...)
  residual <- yy - predict(beta.fit, newx=X)

  # start estimating the mu using multiple folds
  nfolds <- 4 #set up the number of folds
  sfolds <- cut(seq(1,nrow(X)),breaks=nfolds,labels=FALSE)
  taus <- 0
  for(foldid in c(2,3,4)) {
    mu <- rep(0, nrow(X)); folds <- rep(0, nrow(X))
    folds[(sfolds == 1) | (sfolds == foldid)] <- 1
    folds[!((sfolds == 1) | (sfolds == foldid))] <- 2
    #fit the residual with lasso
    for(i in 1:2){
      idx <- which(folds==i,arr.ind=TRUE)
      Xy <- t(center_colmeans(X[-idx,])) %*% (residual[-idx] - mean(residual[-idx])) / (nrow(X) - length(idx))
      mu[idx] <- quad.prog(center_colmeans(X[idx,]), Xy, kappa=kappa)
      mu[idx] <- mu[idx] * length(idx) + mean(residual[idx])
      #mu[idx] <- mu[idx] * length(idx)
    }

    tau <- sum((Y[W==1] - predict(beta.fit, newx=X[W==1,]) - mu[W==1]) / prop.raw[W==1]) / sum(1 / prop.raw[W==1])
    tau <- tau - sum((Y[W==0] - predict(beta.fit, newx=X[W==0,]) - mu[W==0]) / (1 - prop.raw[W==0])) / sum(1 / (1 - prop.raw[W==0]))
    taus <- taus + tau
  }
  taus / 3
}

#' Estimate the Average treatment effect E[Y(1) - Y(0)] from observational data
#'
#' @param X the n by p input covariance matrix
#' @param Y the n dimensional observed response
#' @param W the n dimensional binary vector indicating treatment assignment
#' @param r1 optional n dimensional vector of an initial estimate of E[Y(1) |
#'   X_i] for i = 1, ..., n. The default is NULL
#' @param r0 optional n dimensional vector of an initial estimate of E[Y(0) |
#'   X_i] for i = 1, ..., n. The default is NULL
#' @param kappa the weight parameter for quadratic programming. Default is 0.5
#' @param splitting the options for splitting. "1" means B = 1 split, "3" means
#'   B = 3 splits, "random" means random splits.
#' @param B the number of iterations for random splits, the default is 1. Only
#'   used when splitting is set to "random".
#' @param \dots additional arguments that can be passed to
#'   \code{\link[glmnet]{cv.glmnet}}
#'
#' @examples
#' \dontrun{
#' # Estimating average treatment effect with a toy data
#' # Notice that the external optimisation software \code{MOSEK}
#' # must be installed separately before running the example code.
#' # Without \code{MOSEK}, the example code is not executable.
#' # For how to install \code{MOSEK}, see documentation of \code{\link[Rmosek]{Rmosek}}.
#' set.seed(1)
#' n <- 100; p <- 200
#' X <- scale(matrix(rnorm(n*p), n, p))
#' W <- rbinom(n, 1, 1 / (1 + exp(-X[, 1])))
#' Y <- X[,1] + W * X[,2] + rnorm(n)
#' # Getting an estimate of average treatment effect
#' (est <- dipw.ate(X, Y, W))
#' }
#' @references Wang, Y., Shah, R. D. (2020) \emph{Debiased inverse propensity
#'   score weighting for estimation of average treatment effects with
#'   high-dimensional confounders} \url{https://arxiv.org/abs/2011.08661}
#'
#' @import glmnet
#' @importFrom stats predict
#'
#' @return tau the estimated average treatment effect
#'
#' @export dipw.ate

dipw.ate <- function(X, Y, W, r1=NULL, r0=NULL, kappa=0.5, splitting=c("1", "3", "random"), B=1, ...){

  splitting = match.arg(splitting)

  if(!(splitting %in% c("1", "3", "random"))){
    stop("invalid splitting argument")
  }

  # Check if X, Y, W, r1 and r0 are in correct form
  if (!is.matrix(X)) stop("X should be a matrix with at least one column")
  np <- dim(X)
  if (is.null(np) | (np[2] < 1L))
    stop("X should be a matrix with at least one column")
  n <- as.integer(np[1])
  p <- as.integer(np[2])
  Y <- as.numeric(Y)
  if (length(Y) != n) stop("Y must have nrow(X) components")
  if (length(W) != n) stop("W must have nrow(X) components")
  if (!all(W %in% 0:1)) stop("W must be all binary")

  if(!is.null(r1)) {
    r1 <- as.numeric(r1)
    if (length(r1) != n) stop("r1 must be either null or have nrow(X) components")
  }

  if(!is.null(r0)) {
    r0 <- as.numeric(r0)
    if (length(r0) != n) stop("r0 must be either null or have nrow(X) components")
  }

	#estimate the propensity weights
	prop.fit = cv.glmnet(X, W, family = "binomial", ...)
	prop.raw = predict(prop.fit, newx = X, type = "response")

	residual <- Y
	if(!is.null(r1))
	  residual <- residual - W * r1
	if(!is.null(r0))
	  residual <- residual - (1 - W) * r0

	if(splitting == "1")
	  tau <- dripw.qp.ate(X, residual, W, prop.raw, kappa=kappa, ...)
	else if (splitting == "3")
	  tau <- dripw.split.ate(X, residual, W, prop.raw, kappa=kappa, ...)
	else
	  tau <- dripw.random.ate(X, residual, W, prop.raw, kappa=kappa, B=B, ...)

	if(!is.null(r1)) {
	  tau <- tau + mean(r1)
	}

	if(!is.null(r0)) {
	  tau <- tau - mean(r0)
	}

	tau
}

#' Estimation of E[Y(1)] or E[Y(0)] from observational data
#'
#' @param X the n by p input covariance matrix
#' @param Y the n dimensional observed response
#' @param W the n dimensional binary vector indicating treatment assignment
#' @param Treated \code{TRUE} if we seek to estimate E[Y(1)], \code{FALSE} if we instead wish
#' to estimate E[Y(0)]. The default is TRUE
#' @param r optional n dimensional vector containing initial estimates of
#' E[Y(\code{Treated}) | X_i] for i = 1, ..., n. The default is NULL
#' @param kappa the weight parameter for quadratic programming. Default is 0.5
#' @param splitting the options for splitting. "1" means B = 1 split, "3"
#' means B = 3 splits, "random" means random splits.
#' @param B the number of iterations for random splits, the default is 1.
#' Only valid when splitting is set to "random".
#' @param \dots additional arguments that can be passed to \code{\link[glmnet]{cv.glmnet}}
#'
#' @examples
#' \dontrun{
#' # Estimating mean of the potential outcome with a toy data
#' # Notice that the external optimisation software \code{MOSEK}
#' # must be installed separately before running the example code.
#' # Without \code{MOSEK}, the example code is not executable.
#' # For how to install \code{MOSEK}, see documentation of \code{\link[Rmosek]{Rmosek}}.
#' set.seed(1)
#' n <- 100; p <- 200
#' X <- scale(matrix(rnorm(n*p), n, p))
#' W <- rbinom(n, 1, 1 / (1 + exp(-X[, 1])))
#' Y <- X[,1] + W * X[,2] + rnorm(n)
#' # Getting an estimate of potential outcome mean
#' (est <- dipw.mean(X, Y, W, Treated=TRUE))
#' }
#' @references Wang, Y., Shah, R. D. (2020) \emph{Debiased inverse propensity
#' score weighting for estimation of average treatment effects with high-dimensional
#' confounders} \url{https://arxiv.org/abs/2011.08661}
#'
#' @import glmnet
#' @importFrom stats predict
#'
#' @return the expectation E[Y(1)] or E[Y(0)]
#'
#' @export dipw.mean

dipw.mean <- function(X, Y, W, Treated=TRUE, r=NULL, kappa=0.5, splitting=c("1", "3", "random"), B=1, ...){

  splitting = match.arg(splitting)

  if(!(splitting %in% c("1", "3", "random"))){
    stop('splitting must be chosen from "1", "3", "random"')
  }

  # Check if X, Y, W, r1 and r0 are in correct form
  if (!is.matrix(X)) stop("X should be a matrix with at least one column")
  np <- dim(X)
  if (is.null(np) | (np[2] < 1L))
    stop("X should be a matrix with at least one column")
  n <- as.integer(np[1])
  p <- as.integer(np[2])
  Y <- as.numeric(Y)
  if (length(Y) != n) stop("Y must have nrow(X) components")
  if (length(W) != n) stop("W must have nrow(X) components")
  if (!all(W %in% 0:1)) stop("W must be all binary")

  if(!is.null(r)){
    r <- as.numeric(r)
    if (length(r) != n) stop("r1 must be either null or have nrow(X) components")
  }

  #Start the main program
  if(!Treated) {
    W = !W
  }
  Y[W==0] <- 0
  dipw.ate(X, Y, W, r1=r, kappa=kappa, splitting=splitting, B=B, ...)
}
