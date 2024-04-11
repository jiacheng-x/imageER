#' Function to recover noisy surface with jump detail preserved.
#'
#' @param kk Tuning parameter, size of neighborhood for local weighted average calculation.
#' @param z Input image, single channel only.
#' @param edge Detected edges, should be ourput from \code{edge_detect}.
#'
#' @return A denoised image.
#' @references
#' Qiu, Peihua. “Discontinuous Regression Surfaces Fitting.” The Annals of Statistics 26, no. 6 (1998): 2218–45.
#' @export
#'
#' @examples
recover_surface <- function(z, kk, edge) {
  n <- dim(z)[1]
  kk1 <- trunc((kk + 1) / 2)
  # Upper bar
  fit <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:(kk + kk1 - 1)) {
      fit[i, j] <- 0
      temp1 <- 0
      for (i1 in (i - (kk - 1) / 2):(i + (kk - 1) / 2)) {
        for (j1 in (j - (kk - 1) / 2):(j + (kk - 1) / 2)) {
          if (i1 >= 1 & i1 <= n & j1 >= 1) {
            fit[i, j] <- fit[i, j] + z[i1, j1] * ker_xy((i1 - i) * 2 / (kk - 1), (j1 - j) * 2 / (kk - 1))
            temp1 <- temp1 + ker_xy((i1 - i) * 2 / (kk - 1), (j1 - j) * 2 / (kk - 1))
          }
        }
      }
      fit[i, j] <- fit[i, j] / temp1
    }
  }

  # lower bar
  for (i in 1:n) {
    for (j in (n - (kk + kk1) + 2):n) {
      fit[i, j] <- 0
      temp1 <- 0
      for (i1 in (i - (kk - 1) / 2):(i + (kk - 1) / 2)) {
        for (j1 in (j - (kk - 1) / 2):(j + (kk - 1) / 2)) {
          if (i1 >= 1 & i1 <= n & j1 <= n) {
            fit[i, j] <- fit[i, j] + z[i1, j1] * ker_xy((i1 - i) * 2 / (kk - 1), (j1 - j) * 2 / (kk - 1))
            temp1 <- temp1 + ker_xy((i1 - i) * 2 / (kk - 1), (j1 - j) * 2 / (kk - 1))
          }
        }
      }
      fit[i, j] <- fit[i, j] / temp1
    }
  }

  # left bar
  for (i in 1:(kk + kk1 - 1)) {
    for (j in (kk + kk1):(n - (kk + kk1) + 1)) {
      fit[i, j] <- 0
      temp1 <- 0
      for (i1 in (i - (kk - 1) / 2):(i + (kk - 1) / 2)) {
        for (j1 in (j - (kk - 1) / 2):(j + (kk - 1) / 2)) {
          if (i1 >= 1) {
            fit[i, j] <- fit[i, j] + z[i1, j1] * ker_xy((i1 - i) * 2 / (kk - 1), (j1 - j) * 2 / (kk - 1))
            temp1 <- temp1 + ker_xy((i1 - i) * 2 / (kk - 1), (j1 - j) * 2 / (kk - 1))
          }
        }
      }
      fit[i, j] <- fit[i, j] / temp1
    }
  }

  # right bar
  for (i in (n - (kk + kk1) + 2):n) {
    for (j in (kk + kk1):(n - (kk + kk1) + 1)) {
      fit[i, j] <- 0
      temp1 <- 0
      for (i1 in (i - (kk - 1) / 2):(i + (kk - 1) / 2)) {
        for (j1 in (j - (kk - 1) / 2):(j + (kk - 1) / 2)) {
          if (i1 <= n) {
            fit[i, j] <- fit[i, j] + z[i1, j1] * ker_xy((i1 - i) * 2 / (kk - 1), (j1 - j) * 2 / (kk - 1))
            temp1 <- temp1 + ker_xy((i1 - i) * 2 / (kk - 1), (j1 - j) * 2 / (kk - 1))
          }
        }
      }
      fit[i, j] <- fit[i, j] / temp1
    }
  }

  # center
  for (i in (kk + kk1):(n - (kk + kk1) + 1)) {
    for (j in (kk + kk1):(n - (kk + kk1) + 1)) {
      njump <- 0
      fit[i, j] <- 0
      temp1 <- 0
      for (i1 in (i - (kk - 1) / 2):(i + (kk - 1) / 2)) {
        for (j1 in (j - (kk - 1) / 2):(j + (kk - 1) / 2)) {
          njump <- njump + edge[i1, j1]
          fit[i, j] <- fit[i, j] + z[i1, j1] * ker_xy((i1 - i) * 2 / (kk - 1), (j1 - j) * 2 / (kk - 1))
          temp1 <- temp1 + ker_xy((i1 - i) * 2 / (kk - 1), (j1 - j) * 2 / (kk - 1))
        }
      }
      if (njump <= kk1) {
        fit[i, j] <- fit[i, j] / temp1
      } else {
        fit[i, j] <- 0
        temp1 <- 0
        sigmaxx <- 0
        sigmayy <- 0
        sigmaxy <- 0
        xbar <- 0
        ybar <- 0

        for (i1 in (i - (kk - 1) / 2):(i + (kk - 1) / 2)) {
          for (j1 in (j - (kk - 1) / 2):(j + (kk - 1) / 2)) {
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
          for (i1 in (i - (kk - 1) / 2):(i + (kk - 1) / 2)) {
            for (j1 in (j - (kk - 1) / 2):(j + (kk - 1) / 2)) {
              temp2 <- sigmaxy * (i1 - xbar) + lambda1 * (j1 - ybar)
              if (temp2 >= 0) {
                fit[i, j] <- fit[i, j] + z[i1, j1] * ker_xy((i1 - i) * 2 / (kk - 1), (j1 - j) * 2 / (kk - 1))
                temp1 <- temp1 + ker_xy((i1 - i) * 2 / (kk - 1), (j1 - j) * 2 / (kk - 1))
              }
            }
          }
          fit[i, j] <- fit[i, j] / temp1
        } else {
          for (i1 in (i - (kk - 1) / 2):(i + (kk - 1) / 2)) {
            for (j1 in (j - (kk - 1) / 2):(j + (kk - 1) / 2)) {
              temp2 <- sigmaxy * (i1 - xbar) + lambda1 * (j1 - ybar)
              if (temp2 < 0) {
                fit[i, j] <- fit[i, j] + z[i1, j1] * ker_xy((i1 - i) * 2 / (kk - 1), (j1 - j) * 2 / (kk - 1))
                temp1 <- temp1 + ker_xy((i1 - i) * 2 / (kk - 1), (j1 - j) * 2 / (kk - 1))
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
