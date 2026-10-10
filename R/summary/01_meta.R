# meta.R
# Random-effects meta-regression using metafor.
# Port of meta.py: REML estimation and Knapp–Hartung inference.
#
# Dependency: install.packages("metafor")

rma <- function(y, v, X) {
  if (!requireNamespace("metafor", quietly = TRUE)) {
    stop("Package 'metafor' is required.", call. = FALSE)
  }
  
  y <- as.numeric(y)
  v <- as.numeric(v)
  X <- as.matrix(X)
  
  k <- nrow(X)
  p <- ncol(X)
  
  stopifnot(
    length(y) == k,
    length(v) == k,
    all(is.finite(y)),
    all(is.finite(v)),
    all(is.finite(X)),
    all(v > 0),
    k > p,
    qr(X)$rank == p
  )
  
  # X already contains the complete design matrix.
  model <- metafor::rma(
    yi = y,
    vi = v,
    mods = X,
    intercept = FALSE,
    method = "REML",
    test = "knha",
    level = 95
  )
  
  list(
    b = as.numeric(stats::coef(model)),
    cov = unname(stats::vcov(model)),
    df = k - p,
    tau2 = as.numeric(model$tau2),
    w = 1 / (v + model$tau2),
    k = k,
    p = p,
    model = model
  )
}


contrast <- function(fit, c) {
  # Linear combination c'b, with a two-sided t-test and 95% CI.
  stopifnot(length(c) == fit$p, all(is.finite(c)))
  
  est <- sum(c * fit$b)
  se <- sqrt(drop(t(c) %*% fit$cov %*% c))
  tq <- stats::qt(0.975, df = fit$df)
  
  c(
    estimate = est,
    ci_low = est - tq * se,
    ci_high = est + tq * se,
    p = 2 * stats::pt(
      abs(est / se),
      df = fit$df,
      lower.tail = FALSE
    )
  )
}


prediction_interval <- function(fit, c) {
  # 95% prediction interval for a new species.
  # Retains the original script's t distribution with k - p df.
  stopifnot(length(c) == fit$p, all(is.finite(c)))
  
  est <- sum(c * fit$b)
  se2 <- drop(t(c) %*% fit$cov %*% c)
  h <- stats::qt(0.975, df = fit$df) *
    sqrt(se2 + fit$tau2)
  
  c(
    pi_low = est - h,
    pi_high = est + h
  )
}


ftest <- function(fit, L) {
  # Omnibus Wald-type F-test of L b = 0.
  if (is.null(dim(L))) {
    L <- matrix(L, nrow = 1L)
  } else {
    L <- as.matrix(L)
  }
  
  stopifnot(
    ncol(L) == fit$p,
    nrow(L) >= 1L,
    all(is.finite(L)),
    qr(t(L))$rank == nrow(L)
  )
  
  b <- L %*% fit$b
  V <- L %*% fit$cov %*% t(L)
  df1 <- nrow(L)
  
  F_value <- drop(t(b) %*% solve(V, b)) / df1
  
  c(
    F = F_value,
    df1 = df1,
    df2 = fit$df,
    p = stats::pf(
      F_value,
      df1 = df1,
      df2 = fit$df,
      lower.tail = FALSE
    )
  )
}