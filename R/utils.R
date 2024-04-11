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

ker <- function(x, y) {
  z <- 144 / 121 * (1 - x^2) * (1 - y^2) * (abs(x) <= 1 / 2) * (abs(y) <=
    1 / 2)
  return(z)
}

ker.1 <- function(x) {
  z <- 12 / 11 * (1 - x^2) * ((-1 / 2 <= x) & (x <= 0))
  return(z)
}

ker.2 <- function(x) {
  z <- 12 / 11 * (1 - x^2) * ((0 <= x) & (x <= 1 / 2))
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

#####################################################
### Function used to have symmetric boundary padding##
#####################################################
rotate <- function(x) {
  x[] <- rev(x)
  x
}

BoundaryPadding <- function(z, k) {
  n <- dim(z)[1]
  z_pad <- z
  kk <- (k + 1) / 2
  n_pad <- n + 2 * (kk + k)
  x_pad <- (1:n_pad) / n_pad
  y_pad <- (1:n_pad) / n_pad

  # symmetric padding
  emp_row <- matrix(rep(0, (kk + k) * n), nrow = (kk + k))
  z_pad <- rbind(emp_row, z_pad, emp_row)
  emp_col <- matrix(rep(0, (n + (kk + k) * 2) * (kk + k)), nrow = (n + (kk + k) * 2))
  z_pad <- cbind(emp_col, z_pad, emp_col)

  # Left boundary
  z_pad[((kk + k) + 1):(n_pad - (kk + k)), (kk + k):1] <- z[1:n, 1:(kk + k)]

  # Right boundary
  z_pad[((kk + k) + 1):(n_pad - (kk + k)), (n_pad - (kk + k) + 1):n_pad] <-
    z[1:n, n:(n - (kk + k) + 1)]

  # Top boundary
  z_pad[1:(kk + k), ((kk + k) + 1):(n_pad - (kk + k))] <- z[(kk + k):1, 1:n]

  # Bottom boundary
  z_pad[(n_pad - (kk + k) + 1):n_pad, ((kk + k) + 1):(n_pad - (kk + k))] <-
    z[n:(n - (kk + k) + 1), 1:n]

  # Top-Left corner
  z_pad[1:(kk + k), 1:(kk + k)] <- rotate(z[1:(kk + k), 1:(kk + k)])

  # Top-Right corner
  z_pad[1:(kk + k), (n_pad - (kk + k) + 1):n_pad] <- rotate(z[1:(kk + k), (n - (kk + k) + 1):n])

  # Bottom-Left corner
  z_pad[(n_pad - (kk + k) + 1):n_pad, 1:(kk + k)] <- rotate(z[(n - (kk + k) + 1):n, 1:(kk + k)])

  # Bottom-Right corner
  z_pad[(n_pad - (kk + k) + 1):n_pad, (n_pad - (kk + k) + 1):n_pad] <-
    rotate(z[(n - (kk + k) + 1):n, (n - (kk + k) + 1):n])
  return(z_pad)
}

##########################################################
#### function used to fill non-square image to square#####
##########################################################

SquareFill <- function(z) {
  width <- dim(z)[1]
  height <- dim(z)[2]
  diff <- width - height
  if (diff > 0) {
    cat("Longer width\n")
    zero_matrix <- matrix(rep(0, width * abs(diff)), nrow = width)
    z_fill <- cbind(z, zero_matrix)
  } else if (diff < 0) {
    cat("Longer height\n")
    zero_matrix <- matrix(rep(0, abs(diff) * height), nrow = abs(diff))
    z_fill <- rbind(z, zero_matrix)
  } else {
    cat("Square image, no need to fill\n")
    return(z)
  }
  return(z_fill)
}


###############################################################
#### Side functions for creating Jump preserving estimator#####
###############################################################

B_s1s2_r1r2 <- function(n, x, y, x.i, y.j, s1, s2, r1, r2, h_n) {
  temp0 <- 0
  if (s1 == 1 & s2 == 1) {
    X.temp <- ((x - x.i)^r1 * ker.1((x - x.i) / h_n))
    dim(X.temp) <- c(length(X.temp), 1)
    Y.temp <- ((y - y.j)^r2 * ker.1((y - y.j) / h_n))
    dim(Y.temp) <- c(1, length(Y.temp))
    temp0 <- sum(X.temp %*% Y.temp)
  } else if (s1 == 1 & s2 == 2) {
    X.temp <- ((x - x.i)^r1 * ker.1((x - x.i) / h_n))
    dim(X.temp) <- c(length(X.temp), 1)
    Y.temp <- ((y - y.j)^r2 * ker.2((y - y.j) / h_n))
    dim(Y.temp) <- c(1, length(Y.temp))
    temp0 <- sum(X.temp %*% Y.temp)
  } else if (s1 == 2 & s2 == 1) {
    X.temp <- ((x - x.i)^r1 * ker.2((x - x.i) / h_n))
    dim(X.temp) <- c(length(X.temp), 1)
    Y.temp <- ((y - y.j)^r2 * ker.1((y - y.j) / h_n))
    dim(Y.temp) <- c(1, length(Y.temp))
    temp0 <- sum(X.temp %*% Y.temp)
  } else if (s1 == 2 & s2 == 2) {
    X.temp <- ((x - x.i)^r1 * ker.2((x - x.i) / h_n))
    dim(X.temp) <- c(length(X.temp), 1)
    Y.temp <- ((y - y.j)^r2 * ker.2((y - y.j) / h_n))
    dim(Y.temp) <- c(1, length(Y.temp))
    temp0 <- sum(X.temp %*% Y.temp)
  }
  return(temp0)
}

A_s1s2_p <- function(p, n, x, y, x.i, y.j, s1, s2, h_n) {
  if (p == 1) {
    return(
      B_s1s2_r1r2(n, x, y, x.i, y.j, s1, s2, 2, 0, h_n) * B_s1s2_r1r2(n, x, y, x.i, y.j, s1, s2, 0, 2, h_n) -
        B_s1s2_r1r2(n, x, y, x.i, y.j, s1, s2, 1, 1, h_n) * B_s1s2_r1r2(n, x, y, x.i, y.j, s1, s2, 1, 1, h_n)
    )
  } else if (p == 2) {
    return(
      B_s1s2_r1r2(n, x, y, x.i, y.j, s1, s2, 0, 1, h_n) * B_s1s2_r1r2(n, x, y, x.i, y.j, s1, s2, 1, 1, h_n) -
        B_s1s2_r1r2(n, x, y, x.i, y.j, s1, s2, 1, 0, h_n) * B_s1s2_r1r2(n, x, y, x.i, y.j, s1, s2, 0, 2, h_n)
    )
  } else if (p == 3) {
    return(
      B_s1s2_r1r2(n, x, y, x.i, y.j, s1, s2, 1, 0, h_n) * B_s1s2_r1r2(n, x, y, x.i, y.j, s1, s2, 1, 1, h_n) -
        B_s1s2_r1r2(n, x, y, x.i, y.j, s1, s2, 0, 1, h_n) * B_s1s2_r1r2(n, x, y, x.i, y.j, s1, s2, 2, 0, h_n)
    )
  } else if (p == 4) {
    return(
      B_s1s2_r1r2(n, x, y, x.i, y.j, s1, s2, 0, 0, h_n) * B_s1s2_r1r2(n, x, y, x.i, y.j, s1, s2, 0, 2, h_n) -
        B_s1s2_r1r2(n, x, y, x.i, y.j, s1, s2, 0, 1, h_n) * B_s1s2_r1r2(n, x, y, x.i, y.j, s1, s2, 0, 1, h_n)
    )
  } else if (p == 5) {
    return(
      B_s1s2_r1r2(n, x, y, x.i, y.j, s1, s2, 0, 1, h_n) * B_s1s2_r1r2(n, x, y, x.i, y.j, s1, s2, 1, 0, h_n) -
        B_s1s2_r1r2(n, x, y, x.i, y.j, s1, s2, 0, 0, h_n) * B_s1s2_r1r2(n, x, y, x.i, y.j, s1, s2, 1, 1, h_n)
    )
  } else if (p == 6) {
    return(
      B_s1s2_r1r2(n, x, y, x.i, y.j, s1, s2, 0, 0, h_n) * B_s1s2_r1r2(n, x, y, x.i, y.j, s1, s2, 2, 0, h_n) -
        B_s1s2_r1r2(n, x, y, x.i, y.j, s1, s2, 1, 0, h_n) * B_s1s2_r1r2(n, x, y, x.i, y.j, s1, s2, 1, 0, h_n)
    )
  }
}

delta_s1s2 <- function(n, x, y, x.i, y.j, s1, s2, h_n) {
  delta_temp <- matrix(
    c(
      B_s1s2_r1r2(n, x, y, x.i, y.j, s1, s2, 0, 0, h_n),
      B_s1s2_r1r2(n, x, y, x.i, y.j, s1, s2, 1, 0, h_n),
      B_s1s2_r1r2(n, x, y, x.i, y.j, s1, s2, 0, 1, h_n),
      B_s1s2_r1r2(n, x, y, x.i, y.j, s1, s2, 1, 0, h_n),
      B_s1s2_r1r2(n, x, y, x.i, y.j, s1, s2, 2, 0, h_n),
      B_s1s2_r1r2(n, x, y, x.i, y.j, s1, s2, 1, 1, h_n),
      B_s1s2_r1r2(n, x, y, x.i, y.j, s1, s2, 0, 1, h_n),
      B_s1s2_r1r2(n, x, y, x.i, y.j, s1, s2, 1, 1, h_n),
      B_s1s2_r1r2(n, x, y, x.i, y.j, s1, s2, 0, 2, h_n)
    ),
    nrow = 3,
    byrow = TRUE
  )
  return(det(delta_temp))
}

a_hat_s1s2 <- function(n, x, y, x.i, y.j, s1, s2, h_n, z) {
  A_s1s2_1 <- A_s1s2_p(1, n, x, y, x.i, y.j, s1, s2, h_n)
  A_s1s2_2 <- A_s1s2_p(2, n, x, y, x.i, y.j, s1, s2, h_n)
  A_s1s2_3 <- A_s1s2_p(3, n, x, y, x.i, y.j, s1, s2, h_n)
  delta_temp <- delta_s1s2(n, x, y, x.i, y.j, s1, s2, h_n)
  A_1_mt <- matrix(rep(A_s1s2_1, n * n), nrow = n)
  x_mt <- matrix(rep(x - x.i, n), nrow = n)
  y_mt <- matrix(rep(y - y.j, n), ncol = n, byrow = T)
  if (s1 == 1 & s2 == 1) {
    wt_mt <- ker.1((x - x.i) / h_n) %o% ker.1((y - y.j) / h_n)
    a_sum <- sum((A_1_mt + A_s1s2_2 * x_mt + A_s1s2_3 * y_mt) * wt_mt * z)
  } else if (s1 == 1 & s2 == 2) {
    wt_mt <- ker.1((x - x.i) / h_n) %o% ker.2((y - y.j) / h_n)
    a_sum <- sum((A_1_mt + A_s1s2_2 * x_mt + A_s1s2_3 * y_mt) * wt_mt * z)
  } else if (s1 == 2 & s2 == 1) {
    wt_mt <- ker.2((x - x.i) / h_n) %o% ker.1((y - y.j) / h_n)
    a_sum <- sum((A_1_mt + A_s1s2_2 * x_mt + A_s1s2_3 * y_mt) * wt_mt * z)
  } else if (s1 == 2 & s2 == 2) {
    wt_mt <- ker.2((x - x.i) / h_n) %o% ker.2((y - y.j) / h_n)
    a_sum <- sum((A_1_mt + A_s1s2_2 * x_mt + A_s1s2_3 * y_mt) * wt_mt * z)
  }
  return(a_sum / delta_temp)
}

b_hat_s1s2 <- function(n, x, y, x.i, y.j, s1, s2, h_n, z) {
  temp0 <- 0
  A_s1s2_2 <- A_s1s2_p(2, n, x, y, x.i, y.j, s1, s2, h_n)
  A_s1s2_4 <- A_s1s2_p(4, n, x, y, x.i, y.j, s1, s2, h_n)
  A_s1s2_5 <- A_s1s2_p(5, n, x, y, x.i, y.j, s1, s2, h_n)
  delta_temp <- delta_s1s2(n, x, y, x.i, y.j, s1, s2, h_n)
  A_2_mt <- matrix(rep(A_s1s2_2, n * n), nrow = n)
  x_mt <- matrix(rep(x - x.i, n), nrow = n)
  y_mt <- matrix(rep(y - y.j, n), ncol = n, byrow = T)
  if (s1 == 1 & s2 == 1) {
    wt_mt <- ker.1((x - x.i) / h_n) %o% ker.1((y - y.j) / h_n)
    b_sum <- sum((A_2_mt + A_s1s2_4 * x_mt + A_s1s2_5 * y_mt) * wt_mt * z)
  } else if (s1 == 1 & s2 == 2) {
    wt_mt <- ker.1((x - x.i) / h_n) %o% ker.2((y - y.j) / h_n)
    b_sum <- sum((A_2_mt + A_s1s2_4 * x_mt + A_s1s2_5 * y_mt) * wt_mt * z)
  } else if (s1 == 2 & s2 == 1) {
    wt_mt <- ker.2((x - x.i) / h_n) %o% ker.1((y - y.j) / h_n)
    b_sum <- sum((A_2_mt + A_s1s2_4 * x_mt + A_s1s2_5 * y_mt) * wt_mt * z)
  } else if (s1 == 2 & s2 == 2) {
    wt_mt <- ker.2((x - x.i) / h_n) %o% ker.2((y - y.j) / h_n)
    b_sum <- sum((A_2_mt + A_s1s2_4 * x_mt + A_s1s2_5 * y_mt) * wt_mt * z)
  }
  return(b_sum / delta_temp)
}

c_hat_s1s2 <- function(n, x, y, x.i, y.j, s1, s2, h_n, z) {
  temp0 <- 0
  A_s1s2_3 <- A_s1s2_p(3, n, x, y, x.i, y.j, s1, s2, h_n)
  A_s1s2_5 <- A_s1s2_p(5, n, x, y, x.i, y.j, s1, s2, h_n)
  A_s1s2_6 <- A_s1s2_p(6, n, x, y, x.i, y.j, s1, s2, h_n)
  delta_temp <- delta_s1s2(n, x, y, x.i, y.j, s1, s2, h_n)
  A_3_mt <- matrix(rep(A_s1s2_3, n * n), nrow = n)
  x_mt <- matrix(rep(x - x.i, n), nrow = n)
  y_mt <- matrix(rep(y - y.j, n), ncol = n, byrow = T)
  if (s1 == 1 & s2 == 1) {
    wt_mt <- ker.1((x - x.i) / h_n) %o% ker.1((y - y.j) / h_n)
    c_sum <- sum((A_3_mt + A_s1s2_5 * x_mt + A_s1s2_6 * y_mt) * wt_mt)
  } else if (s1 == 1 & s2 == 2) {
    wt_mt <- ker.1((x - x.i) / h_n) %o% ker.2((y - y.j) / h_n)
    c_sum <- sum((A_3_mt + A_s1s2_5 * x_mt + A_s1s2_6 * y_mt) * wt_mt)
  } else if (s1 == 2 & s2 == 1) {
    wt_mt <- ker.2((x - x.i) / h_n) %o% ker.1((y - y.j) / h_n)
    c_sum <- sum((A_3_mt + A_s1s2_5 * x_mt + A_s1s2_6 * y_mt) * wt_mt)
  } else if (s1 == 2 & s2 == 2) {
    wt_mt <- ker.2((x - x.i) / h_n) %o% ker.2((y - y.j) / h_n)
    c_sum <- sum((A_3_mt + A_s1s2_5 * x_mt + A_s1s2_6 * y_mt) * wt_mt)
  }
  return(c_sum / delta_temp)
}

RSS_a_hat_s1s2 <- function(n, x, y, x.i, y.j, s1, s2, h_n, z) {
  a_hat_temp <- a_hat_s1s2(n, x, y, x.i, y.j, s1, s2, h_n, z)
  b_hat_temp <- b_hat_s1s2(n, x, y, x.i, y.j, s1, s2, h_n, z)
  c_hat_temp <- c_hat_s1s2(n, x, y, x.i, y.j, s1, s2, h_n, z)
  a_hat_mt <- matrix(rep(a_hat_temp, n * n), nrow = n)
  x_mt <- matrix(rep(x - x.i, n), nrow = n)
  y_mt <- matrix(rep(y - y.j, n), ncol = n, byrow = T)
  if (s1 == 1 & s2 == 1) {
    wt_mt <- ker.1((x - x.i) / h_n) %o% ker.1((y - y.j) / h_n)
    rss_sum <-
      sum((z - (a_hat_mt + b_hat_temp * x_mt + c_hat_temp * y_mt))^2 * wt_mt)
  } else if (s1 == 1 & s2 == 2) {
    wt_mt <- ker.1((x - x.i) / h_n) %o% ker.2((y - y.j) / h_n)
    rss_sum <-
      sum((z - (a_hat_mt + b_hat_temp * x_mt + c_hat_temp * y_mt))^2 * wt_mt)
  } else if (s1 == 2 & s2 == 1) {
    wt_mt <- ker.2((x - x.i) / h_n) %o% ker.1((y - y.j) / h_n)
    rss_sum <-
      sum((z - (a_hat_mt + b_hat_temp * x_mt + c_hat_temp * y_mt))^2 * wt_mt)
  } else if (s1 == 2 & s2 == 2) {
    wt_mt <- ker.2((x - x.i) / h_n) %o% ker.2((y - y.j) / h_n)
    rss_sum <-
      sum((z - (a_hat_mt + b_hat_temp * x_mt + c_hat_temp * y_mt))^2 * wt_mt)
  }
  return(rss_sum)
}


# threshold u_n

beta_s1s2_r1r2 <- function(s1, s2, r1, r2) {
  temp0 <- 0
  s <- seq(-1 / 2, 1 / 2, by = 0.01)
  if (s1 == 1 & s2 == 1) {
    for (i in s) {
      temp0 <- temp0 + 0.01^2 * ((i >= -1 / 2) & (i < 0)) * i^r1 * ker.1(i) *
        sum(s[(s >= -1 / 2) & (s < 0)]^r2 * ker.1(s[(s >= -1 / 2) &
          (s < 0)]))
    }
  } else if (s1 == 1 & s2 == 2) {
    for (i in s) {
      temp0 <- temp0 + 0.01^2 * ((i >= -1 / 2) & (i < 0)) * i^r1 * ker.1(i) *
        sum(s[(s > 0) & (s <= 1 / 2)]^r2 * ker.1(s[(s > 0) &
          (s <= 1 / 2)]))
    }
  } else if (s1 == 2 & s2 == 1) {
    for (i in s) {
      temp0 <- temp0 + 0.01^2 * ((i > 0) & (i <= 1 / 2)) * i^r1 * ker.1(i) *
        sum(s[(s >= -1 / 2) & (s < 0)]^r2 * ker.1(s[(s >= -1 / 2) &
          (s < 0)]))
    }
  } else if (s1 == 2 & s2 == 2) {
    for (i in s) {
      temp0 <- temp0 + 0.01^2 * ((i > 0) & (i <= 1 / 2)) * i^r1 * ker.1(i) *
        sum(s[(s > 0) & (s <= 1 / 2)]^r2 * ker.1(s[(s > 0) &
          (s <= 1 / 2)]))
    }
  }
  return(temp0)
}

gamma_s1s2_p <- function(p, s1, s2) {
  if (p == 1) {
    return(
      beta_s1s2_r1r2(s1, s2, 2, 0) * beta_s1s2_r1r2(s1, s2, 0, 2) -
        beta_s1s2_r1r2(s1, s2, 1, 1) * beta_s1s2_r1r2(s1, s2, 1, 1)
    )
  } else if (p == 2) {
    return(
      beta_s1s2_r1r2(s1, s2, 0, 1) * beta_s1s2_r1r2(s1, s2, 1, 1) -
        beta_s1s2_r1r2(s1, s2, 1, 0) * beta_s1s2_r1r2(s1, s2, 0, 2)
    )
  } else if (p == 3) {
    return(
      beta_s1s2_r1r2(s1, s2, 1, 0) * beta_s1s2_r1r2(s1, s2, 1, 1) -
        beta_s1s2_r1r2(s1, s2, 0, 1) * beta_s1s2_r1r2(s1, s2, 2, 0)
    )
  } else if (p == 4) {
    return(
      beta_s1s2_r1r2(s1, s2, 0, 0) * beta_s1s2_r1r2(s1, s2, 0, 2) -
        beta_s1s2_r1r2(s1, s2, 0, 1) * beta_s1s2_r1r2(s1, s2, 0, 1)
    )
  } else if (p == 5) {
    return(
      beta_s1s2_r1r2(s1, s2, 0, 1) * beta_s1s2_r1r2(s1, s2, 1, 0) -
        beta_s1s2_r1r2(s1, s2, 0, 0) * beta_s1s2_r1r2(s1, s2, 1, 1)
    )
  } else if (p == 6) {
    return(
      beta_s1s2_r1r2(s1, s2, 0, 0) * beta_s1s2_r1r2(s1, s2, 2, 0) -
        beta_s1s2_r1r2(s1, s2, 1, 0) * beta_s1s2_r1r2(s1, s2, 1, 0)
    )
  }
}

delta_hat_s1s2 <- function(s1, s2) {
  delta_temp <- matrix(
    c(
      beta_s1s2_r1r2(s1, s2, 0, 0),
      beta_s1s2_r1r2(s1, s2, 1, 0),
      beta_s1s2_r1r2(s1, s2, 0, 1),
      beta_s1s2_r1r2(s1, s2, 1, 0),
      beta_s1s2_r1r2(s1, s2, 2, 0),
      beta_s1s2_r1r2(s1, s2, 1, 1),
      beta_s1s2_r1r2(s1, s2, 0, 1),
      beta_s1s2_r1r2(s1, s2, 1, 1),
      beta_s1s2_r1r2(s1, s2, 0, 2)
    ),
    nrow = 3,
    byrow = TRUE
  )
  return(det(delta_temp))
}

############### This calculation is quite slow########################

# phi <- 0
# s <- seq(-1/2,1/2,by=0.01)
# for (i in s) {
#  phi_temp <- 0
#  for (j in s) {
#    phi_temp <- phi_temp + 1/delta_hat_s1s2(1,1)^2*0.01*0.01*
#      (gamma_s1s2_p(1,1,1)+gamma_s1s2_p(2,1,1)*i+gamma_s1s2_p(3,1,1)*j)^2*
#      ((i>=-1/2) & (i<0))*((j>=-1/2) & (j<0))*ker(i,j)^2
#  }
#  phi <- phi + phi_temp
# }
############### \phi^2=30.59158#########################


u_n <- function(sigma_hat, alpha, n, h_n, phi = 5.327194697) {
  R_alpha <-
    data.frame(
      "alpha" = c(
        0.5,
        0.4,
        0.3,
        0.2,
        0.15,
        0.1,
        0.075,
        0.05,
        0.025,
        0.01,
        0.001,
        0.0001
      ),
      "R" = c(
        1.9746,
        2.2070,
        2.4596,
        2.7785,
        2.9747,
        3.2314,
        3.3959,
        3.6253,
        3.9647,
        4.3911,
        5.2724,
        5.8557
      )
    )
  u_n <- sigma_hat * phi * R_alpha$R[R_alpha$alpha == alpha] / (n * h_n)
  return(u_n)
}


# Generate jump-preserving estimator
rotate <- function(x) {
  x[] <- rev(x)
  x
}

# without symmetric padding
fHat <- function(x.i, y.j, n, x, y, h_n, z, u_n) {
  a_hat_11 <- a_hat_s1s2(n, x, y, x.i, y.j, 1, 1, h_n, z)
  a_hat_12 <- a_hat_s1s2(n, x, y, x.i, y.j, 1, 2, h_n, z)
  a_hat_21 <- a_hat_s1s2(n, x, y, x.i, y.j, 2, 1, h_n, z)
  a_hat_22 <- a_hat_s1s2(n, x, y, x.i, y.j, 2, 2, h_n, z)
  a_hat_all <- c(a_hat_11, a_hat_12, a_hat_21, a_hat_22)
  if ((max(a_hat_all) - min(a_hat_all)) <= u_n) {
    f_hat <- mean(a_hat_all)
  } else {
    rss_all <-
      c(
        RSS_a_hat_s1s2(n, x, y, x.i, y.j, 1, 1, h_n, z),
        RSS_a_hat_s1s2(n, x, y, x.i, y.j, 1, 2, h_n, z),
        RSS_a_hat_s1s2(n, x, y, x.i, y.j, 2, 1, h_n, z),
        RSS_a_hat_s1s2(n, x, y, x.i, y.j, 2, 2, h_n, z)
      )
    f_hat <- a_hat_all[which.min(rss_all)]
  }
  return(f_hat)
}
