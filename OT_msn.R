# Multivariate Skew-Normal (MSN) Transformation & Transport
# Based on Azzalini and Capitanio (1999, 2014)
# Author: Elena Bortolato
# Date: 29 04 2025

# -------------------------------------
# Libraries
# -------------------------------------
library(sn)
library(MASS)
# install.packages("CensMFM")  # Optional alternative sampler

# -------------------------------------
# Clear Environment
# -------------------------------------
rm(list = ls())

# -------------------------------------
# Core Function: Transport map from the d dimensional MSN to d dimensional multivariate standard normal

# n = number of points to sample for visualisation, 
# xi = location parameter of the multivariate skew normal
# Omega = scale parameter
# alpha = skewness
# psi = specific vector (of length d) for which it is of interest to obtain the transformation
# -------------------------------------


ot_msn <- function(n, xi, Omega, alpha, psi) {
  d <- length(alpha) 

  # Generate samples from multivariate skew normal
  samples <- rmsn(n, xi, Omega = Omega, alpha = alpha)
  
  # fix the first value as the value of interest
  samples[1, ] <- psi

  # Delta and transformation matrix Q
  delta <- Omega %*% (alpha / sqrt(1 + c(t(alpha) %*% Omega %*% alpha)))
  alphanorm <- alpha / sqrt(sum(alpha^2))
  Q <- qr.Q(qr(cbind(alphanorm, diag(d)[, -1])))
  if (det(Q) < 0) Q[, 1] <- -Q[, 1]

  # Centered samples with MAP point injected
  centered <- t(t(samples) - xi)
  Z <- centered %*% Q

  # Step 1: Symmetrize
  omega1 <- sqrt((t(Q) %*% Omega %*% Q)[1, 1])
  uniforms <- psn(Z[, 1], alpha = (t(Q) %*% alpha)[1], omega = omega1)
  Z_symm <- Z
  Z_symm[, 1] <- qnorm(uniforms)

  # Step 2: Normalize
  mu <- c(t(Q) %*% delta * sqrt(2 / pi))
  mu[1] <- 0
  V <- t(Q) %*% (Omega - 2 / pi * delta %*% t(delta)) %*% Q
  V[1, -1] <- V[-1, 1] <- (t(diag(c(1, rep(1 / sqrt(V[1, 1]), d - 1)))) %*% V %*%
    t(diag(c(1, rep(1 / sqrt(V[1, 1]), d - 1))))) [-1, 1]
  V[1, 1] <- 1

  # Final normalized samples
  U <- t(t(Z_symm) - mu) %*% t(chol(solve(V)))

  list(
    Z_samples = Z,
    Zprime_samples = Z_symm,
    Normal_samples = U,
    Q = Q,
    delta = delta,
    samples = samples,
    omega1 = omega1
  )
}

# -------------------------------------
# Simulation & Visualization
# -------------------------------------
set.seed(4219)
n <- 10000
xi <- c(0.6, -5)
Omega <- matrix(c(1, 0.85, 0.85, 1), 2, 2)
alpha <- c(0.4, 7)
psi <- c(2, -4)

# Generate and transform samples
result <- ot_msn(n, xi, Omega, alpha, psi)
X <- result$samples
Z <- result$Z_samples
Z_symm <- result$Zprime_samples
U <- result$Normal_samples
Q <- result$Q
omega1 <- result$omega1

# Step: Compute quantile levels
p_levels <- c(0.25, 0.5, 0.75, 0.95)
alpha_rot <- sqrt(sum(alpha^2))
qvals <- qsn(p_levels, xi = 0, omega = omega1, alpha = alpha_rot)

# Quantile lines back-transformed to X-space
lines_in_X <- lapply(qvals, function(q) {
  z_line <- cbind(rep(q, 100), seq(-3, 3, length.out = 100))
  t(Q %*% t(z_line)) + matrix(rep(xi, each = 100), ncol = 2, byrow = FALSE)
})

# -------------------------------------
# Plots
# -------------------------------------
par(mfrow = c(1, 3))

plot(X, col = rgb(0.015, 0, 0.015, 0.015), pch = 16,
     main = "SN", xlab = "X1", ylab = "X2", xlim = c(-3, 3), ylim = c(-3, 4.2))
points(X[1, 1], X[1, 2], col = 2, pch = "x")

plot(Z, col = rgb(0.015, 0, 0.015, 0.015), pch = 16,
     main = "Rotated SN", xlab = "Z1", ylab = "Z2", xlim = c(-3, 3), ylim = c(-3, 4.2))
points(Z[1, 1], Z[1, 2], col = 2, pch = "x")
abline(v = qvals, col = "red", lty = 2, lwd = 0.5)

plot(U, col = rgb(0.015, 0, 0.015, 0.015), pch = 16,
     main = "Normal", xlab = "U1", ylab = "U2")
points(U[1, 1], U[1, 2], col = 2, pch = "x")
abline(h = 0, v = 0)

# -------------------------------------
# Sanity Checks
# -------------------------------------
cat("Mean of U:\n")
print(colMeans(U))
cat("Variance of U:\n")
print(var(U))

# -------------------------------------
# MSN Mean & Variance Calculation
# -------------------------------------
set.seed(11)
X <- rmsn(200000, xi, Omega, alpha = alpha)
delta <- Omega %*% (alpha / sqrt(1 + c(t(alpha) %*% Omega %*% alpha)))
mu_theoretical <- delta * sqrt(2 / pi)

cat("Empirical mean:\n")
print(colMeans(X))
cat("Theoretical mean:\n")
print(mu_theoretical)

cat("Empirical var:\n")
print(var(X))
cat("Theoretical var:\n")
print(Omega - tcrossprod(mu_theoretical))

