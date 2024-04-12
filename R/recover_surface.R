#' Function to recover noisy surface with jump detail preserved.
#'
#' Denoise single channel square images, jump/edge information will be retained instead of being blurred. Will need edge information (output of \code{edge_detect()}) as input.
#' @param z Input image, single channel only, square image only.
#' @param k1 Tuning parameter, size of neighborhood for local weighted average calculation.
#' @param edge A matrix of detected edges, should be ourput from \code{edge_detect}.
#'
#' @return A denoised image.
#' @references
#' Qiu, Peihua. “Discontinuous Regression Surfaces Fitting.” The Annals of Statistics 26, no. 6 (1998): 2218–45.
#' @export
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
#'
#' #Recover/Denoise Surface
#' recover_surface(z, k1=5, edge=edge_map)
#' }


recover_surface <- function(z, k1, edge) {
  if (dim(z)[1] != dim(z)[2]) {
    stop("The image is not a square image, please use DetectionRecoverBlocking function instead.")
  }

  if (length(dim(z)) > 2) {
    stop("The image is not single channel, please use a single channel from color image or grayscale image.")
  }
  n <- dim(z)[1]
  k11 <- trunc((k1 + 1) / 2)
  # Upper bar
  fit <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:(k1 + k11 - 1)) {
      fit[i, j] <- 0
      temp1 <- 0
      for (i1 in (i - (k1 - 1) / 2):(i + (k1 - 1) / 2)) {
        for (j1 in (j - (k1 - 1) / 2):(j + (k1 - 1) / 2)) {
          if (i1 >= 1 & i1 <= n & j1 >= 1) {
            fit[i, j] <- fit[i, j] + z[i1, j1] * ker_xy((i1 - i) * 2 / (k1 - 1), (j1 - j) * 2 / (k1 - 1))
            temp1 <- temp1 + ker_xy((i1 - i) * 2 / (k1 - 1), (j1 - j) * 2 / (k1 - 1))
          }
        }
      }
      fit[i, j] <- fit[i, j] / temp1
    }
  }

  # lower bar
  for (i in 1:n) {
    for (j in (n - (k1 + k11) + 2):n) {
      fit[i, j] <- 0
      temp1 <- 0
      for (i1 in (i - (k1 - 1) / 2):(i + (k1 - 1) / 2)) {
        for (j1 in (j - (k1 - 1) / 2):(j + (k1 - 1) / 2)) {
          if (i1 >= 1 & i1 <= n & j1 <= n) {
            fit[i, j] <- fit[i, j] + z[i1, j1] * ker_xy((i1 - i) * 2 / (k1 - 1), (j1 - j) * 2 / (k1 - 1))
            temp1 <- temp1 + ker_xy((i1 - i) * 2 / (k1 - 1), (j1 - j) * 2 / (k1 - 1))
          }
        }
      }
      fit[i, j] <- fit[i, j] / temp1
    }
  }

  # left bar
  for (i in 1:(k1 + k11 - 1)) {
    for (j in (k1 + k11):(n - (k1 + k11) + 1)) {
      fit[i, j] <- 0
      temp1 <- 0
      for (i1 in (i - (k1 - 1) / 2):(i + (k1 - 1) / 2)) {
        for (j1 in (j - (k1 - 1) / 2):(j + (k1 - 1) / 2)) {
          if (i1 >= 1) {
            fit[i, j] <- fit[i, j] + z[i1, j1] * ker_xy((i1 - i) * 2 / (k1 - 1), (j1 - j) * 2 / (k1 - 1))
            temp1 <- temp1 + ker_xy((i1 - i) * 2 / (k1 - 1), (j1 - j) * 2 / (k1 - 1))
          }
        }
      }
      fit[i, j] <- fit[i, j] / temp1
    }
  }

  # right bar
  for (i in (n - (k1 + k11) + 2):n) {
    for (j in (k1 + k11):(n - (k1 + k11) + 1)) {
      fit[i, j] <- 0
      temp1 <- 0
      for (i1 in (i - (k1 - 1) / 2):(i + (k1 - 1) / 2)) {
        for (j1 in (j - (k1 - 1) / 2):(j + (k1 - 1) / 2)) {
          if (i1 <= n) {
            fit[i, j] <- fit[i, j] + z[i1, j1] * ker_xy((i1 - i) * 2 / (k1 - 1), (j1 - j) * 2 / (k1 - 1))
            temp1 <- temp1 + ker_xy((i1 - i) * 2 / (k1 - 1), (j1 - j) * 2 / (k1 - 1))
          }
        }
      }
      fit[i, j] <- fit[i, j] / temp1
    }
  }

  # center
  for (i in (k1 + k11):(n - (k1 + k11) + 1)) {
    for (j in (k1 + k11):(n - (k1 + k11) + 1)) {
      njump <- 0
      fit[i, j] <- 0
      temp1 <- 0
      for (i1 in (i - (k1 - 1) / 2):(i + (k1 - 1) / 2)) {
        for (j1 in (j - (k1 - 1) / 2):(j + (k1 - 1) / 2)) {
          njump <- njump + edge[i1, j1]
          fit[i, j] <- fit[i, j] + z[i1, j1] * ker_xy((i1 - i) * 2 / (k1 - 1), (j1 - j) * 2 / (k1 - 1))
          temp1 <- temp1 + ker_xy((i1 - i) * 2 / (k1 - 1), (j1 - j) * 2 / (k1 - 1))
        }
      }
      if (njump <= k11) {
        fit[i, j] <- fit[i, j] / temp1
      } else {
        fit[i, j] <- 0
        temp1 <- 0
        sigmaxx <- 0
        sigmayy <- 0
        sigmaxy <- 0
        xbar <- 0
        ybar <- 0

        for (i1 in (i - (k1 - 1) / 2):(i + (k1 - 1) / 2)) {
          for (j1 in (j - (k1 - 1) / 2):(j + (k1 - 1) / 2)) {
            if (edge[i1, j1] == 1) {
              xbar <- xbar + i1
              ybar <- ybar + j1
              sigmaxx <- sigmaxx + i1 * i1
              sigmayy <- sigmayy + j1 * j1
              sigmaxy <- sigmaxy + i1 * j1
            }
          }
        }
        xbar <- xbar / njump
        ybar <- ybar / njump
        sigmaxx <- sigmaxx / njump - xbar^2
        sigmayy <- sigmayy / njump - ybar^2
        sigmaxy <- sigmaxy / njump - xbar * ybar
        del <- (sigmaxx - sigmayy)^2 + 4 * sigmaxy^2
        lambda1 <- (sigmayy - sigmaxx - sqrt(del)) / 2
        prin <- sigmaxy * (i - xbar) + lambda1 * (j - ybar)

        if (prin >= 0) {
          for (i1 in (i - (k1 - 1) / 2):(i + (k1 - 1) / 2)) {
            for (j1 in (j - (k1 - 1) / 2):(j + (k1 - 1) / 2)) {
              temp2 <- sigmaxy * (i1 - xbar) + lambda1 * (j1 - ybar)
              if (temp2 >= 0) {
                fit[i, j] <- fit[i, j] + z[i1, j1] * ker_xy((i1 - i) * 2 / (k1 - 1), (j1 - j) * 2 / (k1 - 1))
                temp1 <- temp1 + ker_xy((i1 - i) * 2 / (k1 - 1), (j1 - j) * 2 / (k1 - 1))
              }
            }
          }
          fit[i, j] <- fit[i, j] / temp1
        } else {
          for (i1 in (i - (k1 - 1) / 2):(i + (k1 - 1) / 2)) {
            for (j1 in (j - (k1 - 1) / 2):(j + (k1 - 1) / 2)) {
              temp2 <- sigmaxy * (i1 - xbar) + lambda1 * (j1 - ybar)
              if (temp2 < 0) {
                fit[i, j] <- fit[i, j] + z[i1, j1] * ker_xy((i1 - i) * 2 / (k1 - 1), (j1 - j) * 2 / (k1 - 1))
                temp1 <- temp1 + ker_xy((i1 - i) * 2 / (k1 - 1), (j1 - j) * 2 / (k1 - 1))
              }
            }
          }
          fit[i, j] <- fit[i, j] / temp1
        }
      }
    }
  }
  return(fit)
}
