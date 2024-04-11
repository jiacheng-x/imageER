#' Generate a jump-preserving surface estimator from noisy data.
#'
#' @param z Input surface, should be single channel image
#' @param sigma_hat Estimated noise level (noise variance) of input surface
#' @param alpha Significance level to control the reconstruction performance of jump, default 0.05
#' @param h Parameter to adjust size of neighborhood for local surface estimation, should be an integer.
#'
#' @return Estimated surface.
#' @references
#' Qiu, Peihua. “The Local Piecewisely Linear Kernel Smoothing Procedure for Fitting Jump Regression Surfaces.” Technometrics 46, no. 1 (February 1, 2004): 87–98. https://doi.org/10.1198/004017004000000149.
#'
#' Qiu, Peihua. “Jump-Preserving Surface Reconstruction from Noisy Data.” Annals of the Institute of Statistical Mathematics 61, no. 3 (September 1, 2009): 715–51. https://doi.org/10.1007/s10463-007-0166-9.
#' @export
#'
#' @examples
#'
JPEstimator <- function(z, h, sigma_hat, alpha=0.05, h_n) {
  n <- dim(z)[1]
  h_n <- h/n
  z_pad <- z
  h1 <- h_n * n
  h2 <- trunc((h1 + 1) / 2)
  n_pad <- n + 2 * h2
  x_pad <- (1:n_pad) / n_pad
  y_pad <- (1:n_pad) / n_pad

  #symmetric padding
  emp_row <- matrix(rep(0, h2 * n), nrow = h2)
  z_pad <- rbind(emp_row, z_pad, emp_row)
  emp_col <- matrix(rep(0, (n + h2 * 2) * h2), nrow = (n + h2 * 2))
  z_pad <- cbind(emp_col, z_pad, emp_col)

  #Left boundary
  z_pad[(h2 + 1):(n_pad - h2), h2:1] <- z[1:n, 1:h2]

  #Right boundary
  z_pad[(h2 + 1):(n_pad - h2), (n_pad - h2 + 1):n_pad] <-
    z[1:n, n:(n - h2 + 1)]

  #Top boundary
  z_pad[1:h2, (h2 + 1):(n_pad - h2)] <- z[h2:1, 1:n]

  #Bottom boundary
  z_pad[(n_pad - h2 + 1):n_pad, (h2 + 1):(n_pad - h2)] <-
    z[n:(n - h2 + 1), 1:n]

  #Top-Left corner
  z_pad[1:h2, 1:h2] <- rotate(z[1:h2, 1:h2])

  #Top-Right corner
  z_pad[1:h2, (n_pad - h2 + 1):n_pad] <- rotate(z[1:h2, (n - h2 + 1):n])

  #Bottom-Left corner
  z_pad[(n_pad - h2 + 1):n_pad, 1:h2] <- rotate(z[(n - h2 + 1):n, 1:h2])

  #Bottom-Right corner
  z_pad[(n_pad - h2 + 1):n_pad, (n_pad - h2 + 1):n_pad] <-
    rotate(z[(n - h2 + 1):n, (n - h2 + 1):n])

  #fitting the surface
  f_hat_pad <- z_pad

  h_n_pad = h_n*n/n_pad

  for (i in (h2 + 1):(n_pad - h2)) {
    for (j in (h2 + 1):(n_pad - h2)) {
      f_hat_pad[i, j] <-
        fHat(i / n_pad, j / n_pad, n_pad, x_pad, y_pad, h_n_pad, z_pad, u_n)
    }
  }
  f_hat <- f_hat_pad[(h2 + 1):(n_pad - h2), (h2 + 1):(n_pad - h2)]
  return(f_hat)
}
