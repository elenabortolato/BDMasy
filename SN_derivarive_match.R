# Multivariate Skew-Normal (MSN) Derivarive matching 
# Based on Zhou, Grazian, Ormerod (2024)
# Author: Elena Bortolato
# Date: 29 04 2025

# -------------------------------------
# Libraries
# -------------------------------------
library(Matrix)

# Define zeta functions
zeta1 <- function(kappa) { dnorm(kappa) / pnorm(kappa) }
zeta2 <- function(kappa) {
  z1 <- zeta1(kappa)
  -z1 * (kappa + z1)
}
zeta3 <- function(kappa) {
  z1 <- zeta1(kappa)
  z2 <- zeta2(kappa)
  -z2 * (kappa + 2 * z1)
}

# Line search equation for kappa
solve_kappa <- function(R) {
  f <- function(kappa) {
    z1 <- zeta1(kappa)
    z2 <- zeta2(kappa)
    z3 <- zeta3(kappa)
    (kappa * z3^(2/3) - R * z1) * (z3^(2/3) + R * z2) + R^2 * z1 * z2
  }
  uniroot(f, lower = -10, upper = 10)$root
}

# Derivative Matching (DM) main function
derivative_matching <- function(m, J, t) {
  p <- length(m)
  
  # Step 1: Compute R
  t_1_3 <- sign(t) * abs(t)^(1/3)
  R <- sum(t_1_3 * solve(J, t_1_3))  # equivalent to (t^(1/3))' J^{-1} (t^(1/3))
  
  # Step 2: Solve for kappa
  kappa_star <- solve_kappa(R)
  
  # Step 3: Compute d*, Sigma*, mu*
  d_star <- sign(t) * abs(t / zeta3(kappa_star))^(1/3)
  Sigma_star <- solve(J + zeta2(kappa_star) * (d_star %*% t(d_star)))
  mu_star <- m - zeta1(kappa_star) * Sigma_star %*% d_star
  
  list(mu = mu_star, Sigma = Sigma_star, d = d_star)
}
