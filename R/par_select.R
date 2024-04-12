#' Function to select optimal tuning parameter for edge detection (k) and surface recovery (k1)
#'
#' Choose optimal tuning parameter for jump/edge detection and surface recovery/denoising. Omega value can be adjusted for better edge detection or better denoising.
#' @param k A vector of tuning parameter k, should all be odd numbers.
#' @param k1 A vector of tuning parameter k1.
#' @param z Input image, single channel only, square image only.
#' @param h Parameter to adjust size of neighborhood during least square calculation, should be an integer.
#' @param alpha_n Parameter to adjust detection sensitivity, default 0.05.
#' @param alpha_jp Significance level to control the reconstruction performance of jump, default 0.05. Can be chosen from (0.5, 0.4, 0.3, 0.2, 0.15, 0.1, 0.075, 0.05, 0.025, 0.01, 0.001, 0.0001)
#' @param omega Parameter to adjust for selection preference, default 0.5. If greater than 0.5, focus more on edge detection, if less than 0.5, focus more on surface recovery.
#'
#' @return A message of selection.
#' @export
#' @importFrom pracma hausdorff_dist
#'
#' @examples
#' \dontrun{
#' # Generate a surface with a circular jump in the middle
#' n <- 64
#' x <- y <- (1:n) / n
#' f <- matrix(0, n, n)
#'
#' for (i in 1:n) {
#'   for (j in 1:n) {
#'     if ((x[i] - 0.5)^2 + (y[j] - 0.5)^2 > 0.25^2) {
#'       f[i, j] <- -2 * (x[i] - 0.5)^2 - 2 * (y[j] - 0.5)^2
#'     } else {
#'       f[i, j] <- -2 * (x[i] - 0.5)^2 - 2 * (y[j] - 0.5)^2 + 1
#'     }
#'   }
#' }
#'
#' # Adding normal noise
#' sigma <- 1 / 2 * sd(as.vector(f))
#' set.seed(1234)
#' noise <- matrix(rnorm(n * n, 0, sigma), n, n)
#' z <- f + noise
#'
#' k <- c(3, 5, 7, 9, 11)
#' k1 <- c(5, 10, 15, 20, 25)
#'
#' par_select(k, k1, z, h = 3, alpha_n = 0.05, alpha_jp = 0.05, omega = 0.5)
#' }
par_select <- function(k, k1, z, h, alpha_n = 0.05, alpha_jp = 0.05, omega = 0.5) {
  if (dim(z)[1] != dim(z)[2]) {
    stop("The image is not a square image, please use DetectionRecoverBlocking function instead.")
  }

  if (length(dim(z)) > 2) {
    stop("The image is not single channel, please use a single channel from color image or grayscale image.")
  }

  n <- dim(z)[1]

  # Generate Jump-Presering Bootstrap samples for True Hausdroff Distance estimation
  f_hat_int <- jump_preserving_estimator(z, h, sigma_hat = 0, alpha = alpha_jp)
  sigma_est <- sqrt(sum((z - f_hat_int)^2) / n^2)
  f_hat_upd <- jump_preserving_estimator(z, h, sigma_hat = sigma_est, alpha = alpha_jp)

  # Residual set
  epsilon <- z - f_hat_upd

  # Generate bootstrap sample
  for (i in 1:10) {
    assign(paste0("f_hat_", i), f_hat_upd + matrix(sample(epsilon, n * n, replace = TRUE), nrow = n, ncol = n))
  }

  out <- data.frame(
    "k" = sort(rep(k, length(k1))),
    "k1" = rep(k1, length(k)),
    "SSEO" = rep(0, length(k) * length(k1)),
    "HD_est_tm" = rep(0, length(k) * length(k1))
  )
  for (i in k) {
    edge <- edge_detect(z, k = i, h, alpha_n)
    #### Hausdorff distance
    assign(paste0("jump_detect_", i), as.matrix(which(edge == 1, arr.ind = T) / n))
    assign(paste0("hd_temp_", i), rep(NA, 10))
    for (j in 1:10) {
      assign(paste0("edge_temp_", i, j), edge_detect(get(paste0("f_hat_", j)), k = i, h, alpha_n))
      assign(paste0("jump_detect_temp_", i, j), as.matrix(which(get(paste0("edge_temp_", i, j)) == 1, arr.ind = T) / n))
      hd_temp_value <- hausdorff_dist(get(paste0("jump_detect_temp_", i, j)), get(paste0("jump_detect_", i)))
      hd_tempvec_name <- paste0("hd_temp_", i, "[", j, "]")
      hd_tempval_name <- "hd_temp_value"
      eval(parse(text = paste(hd_tempvec_name, "<-", hd_tempval_name)))
    }
    out$HD_est_tm[out$k == i] <- mean(get(paste0("hd_temp_", i)), trim = 0.1)
    for (l in k1) {
      # Recover the surface with first-order approximation of the JLC
      assign(paste0("fit_", i, l), recover_surface(z, k1 = l, edge))

      #### Sum of squared deviations from observed responses
      out$SSEO[out$k == i & out$k1 == l] <- sum((z - get(paste0("fit_", i, l)))^2)
    }
  }
  optimal <- crit_func(out, omega)
  print(paste("The optimal k is", optimal$k, ", the optimal k1 is", optimal$k1))
}
