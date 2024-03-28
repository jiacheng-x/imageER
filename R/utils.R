#####################################################################
###################### Kernel Functions###############################
#####################################################################

ker_xy <- function(x, y) {
  z <- 9 / 16 * (1 - x^2) * (1 - y^2) * (abs(x) < 1) * (abs(y) < 1)
  return(z)
}

ker_x <- function(x) {
  z <- 3 / 4 * (1 - x^2) * (abs(x) < 1)
  return(z)
}

######################################################################
######## Function used to estimate least square coefficients##########
######################################################################


# x and y are coordinates, n is sample size, z is the surface
LSEstimate <- function(x, y, z, n, h) {
  beta1 <- matrix(rep(0, n * n), nrow = n) # define empty matrix to store beta values
  beta2 <- matrix(rep(0, n * n), nrow = n) # define empty matrix to store beta values

  for (i in 1:n) {
    for (j in 1:n) {
      x_wt_s <- sum(ker_x((x - x[i]) / h))
      x_wt_s1 <- sum(ker_x((x - x[i]) / h) * (x - x[i]))
      x_wt_s2 <- sum(ker_x((x - x[i]) / h) * (x - x[i])^2)
      x_wt_s3 <- sum(ker_x((x - x[i]) / h) * (x - x[i])^3)
      x_wt_s4 <- sum(ker_x((x - x[i]) / h) * (x - x[i])^4)
      y_wt_s <- sum(ker_x((y - y[j]) / h))
      y_wt_s1 <- sum(ker_x((y - y[j]) / h) * (y - y[j]))
      y_wt_s2 <- sum(ker_x((y - y[j]) / h) * (y - y[j])^2)
      y_wt_s3 <- sum(ker_x((y - y[j]) / h) * (y - y[j])^3)
      y_wt_s4 <- sum(ker_x((y - y[j]) / h) * (y - y[j])^4)
      X.temp <- matrix(
        c(
          x_wt_s * y_wt_s, y_wt_s * x_wt_s1, x_wt_s * y_wt_s1, y_wt_s * x_wt_s2, x_wt_s * y_wt_s2, 2 * x_wt_s1 * y_wt_s1,
          y_wt_s * x_wt_s1, y_wt_s * x_wt_s2, x_wt_s1 * y_wt_s1, y_wt_s * x_wt_s3, x_wt_s1 * y_wt_s2, 2 * x_wt_s2 * y_wt_s1,
          x_wt_s * y_wt_s1, x_wt_s1 * y_wt_s1, x_wt_s * y_wt_s2, y_wt_s1 * x_wt_s2, x_wt_s * y_wt_s3, 2 * x_wt_s1 * y_wt_s2,
          y_wt_s * x_wt_s2, y_wt_s * x_wt_s3, y_wt_s1 * x_wt_s2, y_wt_s * x_wt_s4, x_wt_s2 * y_wt_s2, 2 * x_wt_s3 * y_wt_s1,
          x_wt_s * y_wt_s2, x_wt_s1 * y_wt_s2, x_wt_s * y_wt_s3, x_wt_s2 * y_wt_s2, x_wt_s * y_wt_s4, 2 * x_wt_s1 * y_wt_s3,
          2 * x_wt_s1 * y_wt_s1, 2 * x_wt_s2 * y_wt_s1, 2 * x_wt_s1 * y_wt_s2, 2 * x_wt_s3 * y_wt_s1, 2 * x_wt_s1 * y_wt_s3, 4 * x_wt_s2 * y_wt_s2
        ),
        byrow = T, ncol = 6
      )
      Z.temp <- z
      Z.temp2 <- matrix(rep(0, n * n), nrow = n)
      Z.temp3 <- matrix(rep(0, n * n), nrow = n)
      Z.temp4 <- matrix(rep(0, n * n), nrow = n)
      Z.temp5 <- matrix(rep(0, n * n), nrow = n)
      Z.temp6 <- matrix(rep(0, n * n), nrow = n)
      wt_mt <- ker_x((x - x[i]) / h) %o% ker_x((y - y[j]) / h)
      x_mt <- matrix(rep(x - x[i], n), nrow = n)
      y_mt <- matrix(rep(y - y[j], n), ncol = n, byrow = T)
      Z.temp <- Z.temp * wt_mt
      Z.temp2 <- Z.temp * x_mt
      Z.temp3 <- Z.temp * y_mt
      Z.temp4 <- Z.temp * x_mt^2
      Z.temp5 <- Z.temp * y_mt^2
      Z.temp6 <- 2 * Z.temp * x_mt * y_mt
      Z.m <- matrix(c(sum(Z.temp), sum(Z.temp2), sum(Z.temp3), sum(Z.temp4), sum(Z.temp5), sum(Z.temp6)), ncol = 1)
      beta.temp <- solve(X.temp) %*% Z.m
      beta1[i, j] <- beta.temp[2]
      beta2[i, j] <- beta.temp[3]
    }
  }
  beta_values <- list("Beta_1" = beta1, "Beta_2" = beta2)
  return(beta_values)
}
