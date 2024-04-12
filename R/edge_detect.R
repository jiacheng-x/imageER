#' Function to detect edges/jumps on SQUARE image.
#'
#' Detect edges/jumps of single channel square images (equal width and height).
#' @param z Input image, single channel only, square image only.
#' @param k Tuning parameter, size of neighborhood for edge calculation, should be an odd integer.
#' @param alpha_n Parameter to adjust detection sensitivity, default 0.05.
#' @param h Parameter to adjust size of neighborhood during least square calculation, should be an integer.
#'
#' @return An 0, 1 matrix of edge locations with 1 being edge location
#' @references
#' Qiu, Peihua. “Discontinuous Regression Surfaces Fitting.” The Annals of Statistics 26, no. 6 (1998): 2218–45.
#' @export
#' @importFrom stats qnorm
#' @importFrom stats var
#'
#' @examples
#' \dontrun{
#' #Generate a surface with a circular jump in the middle
#' n <- 64
#' x <- y <- (1:n)/n
#' f <- matrix(0, n , n)
#'
#' for (i in 1:n){
#'   for (j in 1:n){
#'     if ((x[i]-0.5)^2+(y[j]-0.5)^2 > 0.25^2){
#'       f[i,j] <- -2*(x[i]-0.5)^2-2*(y[j]-0.5)^2
#'     } else{
#'       f[i,j] <- -2*(x[i]-0.5)^2-2*(y[j]-0.5)^2+1
#'     }
#'   }
#' }
#'
#' #Adding normal noise
#' sigma <- 1/2 * sd(as.vector(f))
#' set.seed(1234)
#' noise <- matrix(rnorm(n*n,0,sigma), n, n)
#' z <- f + noise
#'
#' #Detecting Edge
#' edge_map <- edge_detect(z, k=5, h=4, alpha_n=0.05)
#' }

edge_detect <- function(z, k, h, alpha_n = 0.05) {
  if (k %% 2 == 0) {
    stop("k is not odd, please ensure k is an odd number.")
  }

  if (dim(z)[1] != dim(z)[2]) {
    stop("The image is not a square image, please use DetectionRecoverBlocking function instead.")
  }

  if (length(dim(z)) > 2) {
    stop("The image is not single channel, please use a single channel from color image or grayscale image.")
  }
  n <- dim(z)[1]
  x <- (1:n) / n
  y <- (1:n) / n
  h_n <- h/n
  Betas <- LSEstimate(x, y, z, n, h_n)
  beta1 <- Betas$Beta_1
  beta2 <- Betas$Beta_2
  k1 <- (k + 1) / 2
  neigh <- matrix(rep(0, 4), nrow = 2)
  delta <- matrix(rep(0, n * n), nrow = n)
  threshold <- matrix(rep(0, n * n), nrow = n)
  edge <- matrix(rep(0, n * n), nrow = n)
  for (i in (k + k1):(n - k - k1 + 1)) {
    for (j in (k + k1):(n - k - k1 + 1)) {
      deno <- beta1[i, j]^2 + beta2[i, j]^2
      if (deno == 0) {
        cat("No edge detected at (", i, ", ", j, ").\n")
        break
      } else {
        cond <- abs(beta2[i, j]) / sqrt(deno)
        if (cond > sqrt(2) / 2) {
          neigh[1, 1] <- trunc(i + k / beta2[i, j] * beta1[i, j])
          neigh[1, 2] <- j + k
          neigh[2, 1] <- trunc(i - k / beta2[i, j] * beta1[i, j])
          neigh[2, 2] <- j - k
        } else {
          neigh[1, 1] <- i + k
          neigh[1, 2] <- trunc(j + k * beta2[i, j] / beta1[i, j])
          neigh[2, 1] <- i - k
          neigh[2, 2] <- trunc(j - k * beta2[i, j] / beta1[i, j])
        }
        temp1 <- 0
        temp2 <- 0
        temp3 <- 0
        z_temp1 <- c()
        z_temp2 <- c()
        z_temp3 <- c()
        for (i1 in (max(neigh[1, 1] - k1 + 1, 1)):(min(neigh[1, 1] + k1 - 1, n))) {
          for (j1 in (max(neigh[1, 2] - k1 + 1, 1)):(min(neigh[1, 2] + k1 - 1, n))) {
            temp1 <- temp1 + z[i1, j1]
            z_temp1 <- c(z_temp1, z[i1, j1])
          }
        }
        for (i1 in (max(neigh[2, 1] - k1 + 1, 1)):(min(neigh[2, 1] + k1 - 1, n))) {
          for (j1 in (max(neigh[2, 2] - k1 + 1, 1)):(min(neigh[2, 2] + k1 - 1, n))) {
            temp2 <- temp2 + z[i1, j1]
            z_temp2 <- c(z_temp2, z[i1, j1])
          }
        }
        for (i1 in (i - k1 + 1):(i + k1 - 1)) {
          for (j1 in (j - k1 + 1):(j + k1 - 1)) {
            temp3 <- temp3 + z[i1, j1]
            z_temp3 <- c(z_temp3, z[i1, j1])
          }
        }

        delta[i, j] <- min(abs(temp1 - temp3), abs(temp2 - temp3)) / (k * k)
        var_hat <- (var(z_temp3) + var(z_temp1)) * (abs(temp1 - temp3) > abs(temp2 - temp3)) +
          (var(z_temp3) + var(z_temp2)) * (1 - (abs(temp1 - temp3) > abs(temp2 - temp3))) # use sample variance to estimate variance in threshold calculation
        threshold[i, j] <- qnorm(1 - alpha_n / 2) * sqrt(2 * var_hat) / k
      }
    }
  }

  for (i in (k + k1):(n - k - k1 + 1)) {
    for (j in (k + k1):(n - k - k1 + 1)) {
      if (delta[i, j] > threshold[i, j]) {
        edge[i, j] <- 1
      }
    }
  }

  ## modify four corners
  # up-left
  sum <- 0
  for (i in (k + k1):(2 * k + k1 - 1)) {
    for (j in (k + k1):(2 * k + k1 - 1)) {
      sum <- sum + edge[i, j]
    }
  }
  if (sum <= k1) {
    for (i in (k + k1):(2 * k)) {
      for (j in (k + k1):(2 * k)) {
        edge[i, j] <- 0
      }
    }
  }
  # up-right
  sum <- 0
  for (i in (n - 2 * k - k1 + 2):(n - k - k1 + 1)) {
    for (j in (k + k1):(2 * k + k1 - 1)) {
      sum <- sum + edge[i, j]
    }
  }
  if (sum <= k1) {
    for (i in (n - 2 * k + 1):(n - k - k1 + 1)) {
      for (j in (k + k1):(2 * k)) {
        edge[i, j] <- 0
      }
    }
  }
  # down-left
  sum <- 0
  for (i in (k + k1):(2 * k + k1 - 1)) {
    for (j in (n - 2 * k - k1 + 2):(n - k - k1 + 1)) {
      sum <- sum + edge[i, j]
    }
  }
  if (sum <= k1) {
    for (i in (k + k1):(2 * k)) {
      for (j in (n - 2 * k + 1):(n - k - k1 + 1)) {
        edge[i, j] <- 0
      }
    }
  }
  # down-right
  sum <- 0
  for (i in (n - 2 * k - k1 + 2):(n - k - k1 + 1)) {
    for (j in (n - 2 * k - k1 + 2):(n - k - k1 + 1)) {
      sum <- sum + edge[i, j]
    }
  }
  if (sum <= k1) {
    for (i in (n - 2 * k + 1):(n - k - k1 + 1)) {
      for (j in (n - 2 * k + 1):(n - k - k1 + 1)) {
        edge[i, j] <- 0
      }
    }
  }
  # upper boundary
  for (i in (2 * k + 1):(n - 2 * k)) {
    for (j in (k + k1):(2 * k)) {
      sum <- 0
      if (edge[i, j] == 1) {
        for (i1 in (i - k1 + 1):(i + k1 - 1)) {
          for (j1 in (k + k1):(2 * k + k1 - 1)) {
            sum <- sum + edge[i1, j1]
          }
        }
        if (sum <= k1) {
          edge[i, j] <- 0
        }
      }
    }
  }
  # down boundary
  for (i in (2 * k + 1):(n - 2 * k)) {
    for (j in (n - 2 * k + 1):(n - k - k1 + 1)) {
      sum <- 0
      if (edge[i, j] == 1) {
        for (i1 in (i - k1 + 1):(i + k1 - 1)) {
          for (j1 in (n - 2 * k - k1 + 2):(n - k - k1 + 1)) {
            sum <- sum + edge[i1, j1]
          }
        }
        if (sum <= k1) {
          edge[i, j] <- 0
        }
      }
    }
  }
  # left boundary
  for (j in (2 * k + 1):(n - 2 * k)) {
    for (i in (k + k1):(2 * k)) {
      sum <- 0
      if (edge[i, j] == 1) {
        for (j1 in (j - k1 + 1):(j + k1 - 1)) {
          for (i1 in (k + k1):(2 * k + k1 - 1)) {
            sum <- sum + edge[i1, j1]
          }
        }
        if (sum <= k1) {
          edge[i, j] <- 0
        }
      }
    }
  }
  # right boundary
  for (j in (2 * k + 1):(n - 2 * k)) {
    for (i in (n - 2 * k + 1):(n - k - k1 + 1)) {
      sum <- 0
      if (edge[i, j] == 1) {
        for (j1 in (j - k1 + 1):(j + k1 - 1)) {
          for (i1 in (n - 2 * k - k1 + 2):(n - k - k1 + 1)) {
            sum <- sum + edge[i1, j1]
          }
        }
        if (sum <= k1) {
          edge[i, j] <- 0
        }
      }
    }
  }
  # center boundary
  for (i in (2 * k + 1):(n - 2 * k)) {
    for (j in (2 * k + 1):(n - 2 * k)) {
      sum <- 0
      if (edge[i, j] == 1) {
        for (i1 in (i - k1 + 1):(i + k1 - 1)) {
          for (j1 in (j - k1 + 1):(j + k1 - 1)) {
            sum <- sum + edge[i1, j1]
          }
        }
        if (sum <= k1) {
          edge[i, j] <- 0
        }
      }
    }
  }
  return(edge)
}
