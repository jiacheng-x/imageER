#' Title
#'
#' @param z
#' @param h_n
#' @param u_n
#'
#' @return
#' @export
#'
#' @examples
JPEstimator <- function(z, h_n, u_n) {
  n <- dim(z)[1]
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
