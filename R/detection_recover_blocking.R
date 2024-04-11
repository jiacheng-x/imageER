#' Title Function used to detect edge and denoise image with blocking capability
#'
#' Detect edges and recover images in one step, efficiently handling non-square images by dividing them into smaller blocks to speed up calculations.
#' @param z Input image, single channel only.
#' @param h Parameter to adjust size of neighborhood during least square calculation, should be an integer.
#' @param k Tuning parameter, size of neighborhood for edge calculation.
#' @param k1 Tuning parameter, size of neighborhood for denoising.
#' @param block_size Size of blocks you wish to use, default size 64x64.
#'
#' @return A list of detected edges (an 0, 1 matrix of edge locations with 1 being edge location) and denoised image.
#' @export
#'
#' @examples
DetectionRecoverBlocking <- function(z, h, k, k1, block_size = 64) {
  # Get dimensions
  width <- dim(z)[1]
  height <- dim(z)[2]

  # Calculate remainder and adjust overlap
  remainder_width <- block_size - width %% block_size
  remainder_height <- block_size - height %% block_size

  # Calculate the number of horizontal and vertical blocks
  n_blocks_width <- width %/% block_size
  n_blocks_height <- height %/% block_size

  # Calculate overlap size
  overlap_width <- remainder_width %/% n_blocks_width
  overlap_height <- remainder_height %/% n_blocks_height

  # Calculate last overlap size
  overlap_width_last <- remainder_width %% overlap_width + overlap_width
  overlap_height_last <- remainder_height %% overlap_height + overlap_height

  # Calculate non-overlap size
  non_overlap_width <- block_size - overlap_width
  non_overlap_height <- block_size - overlap_height

  # Print the number of horizontal and vertical blocks
  cat("Number of horizontal blocks:", n_blocks_width + 1, "\n")
  cat("Number of vertical blocks:", n_blocks_height + 1, "\n")
  cat("Horizontal overlapping size:", overlap_width, "\n")
  cat("Vertical overlapping size:", overlap_height, "\n")
  cat("Last horizontal overlapping size:", overlap_width_last, "\n")
  cat("Last vertical overlapping size:", overlap_height_last, "\n")
  cat("Horizontal non-overlapping size:", non_overlap_width, "\n")
  cat("Vertical non-overlapping size:", non_overlap_height, "\n")

  blocks <- list()
  result_matrix <- matrix(1, nrow = width, ncol = height)
  recover_matrix <- matrix(1, nrow = width, ncol = height)
  kk <- (k + 1) / 2
  n_pad <- block_size + 2 * (kk + k)
  x <- y <- (1:n_pad) / n_pad
  h_pad <- h * block_size / n_pad
  p <- 0
  # Loop through the image to extract blocks
  for (i in seq(1, 1 + (n_blocks_height - 1) * non_overlap_height, length.out = n_blocks_height)) {
    for (j in seq(1, 1 + (n_blocks_width - 1) * non_overlap_width, length.out = n_blocks_width)) {
      p <- p + 1
      # Special treat for first block
      if (i == 1 & j == 1) {
        cat(
          "Block # ", p, ":Block width from", j, "to", (j + block_size - 1),
          "height from", i, "to", (i + block_size - 1), "\n"
        )
        block <- z[j:(j + block_size - 1), i:(i + block_size - 1)]
        # blocks = append(blocks, list(block))
        block_pad <- BoundaryPadding(block, k)
        ls_estimator <- LSEstimate(x, y, block_pad, n_pad, h_pad)
        ls_beta1 <- ls_estimator$Beta_1
        ls_beta2 <- ls_estimator$Beta_2
        edge_temp <- EdgeDectect(block_pad, k, n_pad, ls_beta1, ls_beta2, 0.05)
        # record final edge matrix
        result_matrix[
          j:(j + block_size - 1 - overlap_width / 2),
          i:(i + block_size - 1 - overlap_height / 2)
        ] <-
          edge_temp[
            (1 + (kk + k)):(block_size + (kk + k)),
            (1 + (kk + k)):(block_size + (kk + k))
          ][
            1:(1 + block_size - 1 - overlap_width / 2),
            1:(1 + block_size - 1 - overlap_height / 2)
          ]
        # record final recovery matrix
        recover_matrix[
          j:(j + block_size - 1 - overlap_width / 2),
          i:(i + block_size - 1 - overlap_height / 2)
        ] <-
          RecoverSurface(n_pad, k1, block_pad, edge_temp)[
            (1 + (kk + k)):(block_size + (kk + k)),
            (1 + (kk + k)):(block_size + (kk + k))
          ][
            1:(1 + block_size - 1 - overlap_width / 2),
            1:(1 + block_size - 1 - overlap_height / 2)
          ]
      } else if (j == 1) {
        cat(
          "Block # ", p, ":Block width from", j, "to", (j + block_size - 1),
          "height from", i, "to", (i + block_size - 1), "\n"
        )
        block <- z[j:(j + block_size - 1), i:(i + block_size - 1)]
        # blocks = append(blocks, list(block))
        block_pad <- BoundaryPadding(block, k)
        ls_estimator <- LSEstimate(x, y, block_pad, n_pad, h_pad)
        ls_beta1 <- ls_estimator$Beta_1
        ls_beta2 <- ls_estimator$Beta_2
        edge_temp <- EdgeDectect(block_pad, k, n_pad, ls_beta1, ls_beta2, 0.05)
        # record final edge matrix
        result_matrix[
          j:(j + block_size - 1 - overlap_width / 2),
          (i + overlap_height / 2):(i + block_size - 1 - overlap_height / 2)
        ] <-
          edge_temp[
            (1 + (kk + k)):(block_size + (kk + k)),
            (1 + (kk + k)):(block_size + (kk + k))
          ][
            1:(1 + block_size - 1 - overlap_width / 2),
            (1 + overlap_height / 2):(1 + block_size - 1 - overlap_height / 2)
          ]
        # record final recovery matrix
        recover_matrix[
          j:(j + block_size - 1 - overlap_width / 2),
          (i + overlap_height / 2):(i + block_size - 1 - overlap_height / 2)
        ] <-
          RecoverSurface(n_pad, k1, block_pad, edge_temp)[
            (1 + (kk + k)):(block_size + (kk + k)),
            (1 + (kk + k)):(block_size + (kk + k))
          ][
            1:(1 + block_size - 1 - overlap_width / 2),
            (1 + overlap_height / 2):(1 + block_size - 1 - overlap_height / 2)
          ]
      } else if (i == 1) {
        cat(
          "Block # ", p, ":Block width from", j, "to", (j + block_size - 1),
          "height from", i, "to", (i + block_size - 1), "\n"
        )
        block <- z[j:(j + block_size - 1), i:(i + block_size - 1)]
        # blocks = append(blocks, list(block))
        block_pad <- BoundaryPadding(block, k)
        ls_estimator <- LSEstimate(x, y, block_pad, n_pad, h_pad)
        ls_beta1 <- ls_estimator$Beta_1
        ls_beta2 <- ls_estimator$Beta_2
        edge_temp <- EdgeDectect(block_pad, k, n_pad, ls_beta1, ls_beta2, 0.05)
        # record final edge matrix
        result_matrix[
          (j + overlap_width / 2):(j + block_size - 1 - overlap_width / 2),
          i:(i + block_size - 1 - overlap_height / 2)
        ] <-
          edge_temp[
            (1 + (kk + k)):(block_size + (kk + k)),
            (1 + (kk + k)):(block_size + (kk + k))
          ][
            (1 + overlap_width / 2):(1 + block_size - 1 - overlap_width / 2),
            1:(1 + block_size - 1 - overlap_height / 2)
          ]
        # record final recovery matrix
        recover_matrix[
          (j + overlap_width / 2):(j + block_size - 1 - overlap_width / 2),
          i:(i + block_size - 1 - overlap_height / 2)
        ] <-
          RecoverSurface(n_pad, k1, block_pad, edge_temp)[
            (1 + (kk + k)):(block_size + (kk + k)),
            (1 + (kk + k)):(block_size + (kk + k))
          ][
            (1 + overlap_width / 2):(1 + block_size - 1 - overlap_width / 2),
            1:(1 + block_size - 1 - overlap_height / 2)
          ]
      } else {
        # Extract block
        cat(
          "Block # ", p, ":Block width from", j, "to", (j + block_size - 1),
          "height from", i, "to", (i + block_size - 1), "\n"
        )
        block <- z[j:(j + block_size - 1), i:(i + block_size - 1)]
        # blocks = append(blocks, list(block))
        block_pad <- BoundaryPadding(block, k)
        ls_estimator <- LSEstimate(x, y, block_pad, n_pad, h_pad)
        ls_beta1 <- ls_estimator$Beta_1
        ls_beta2 <- ls_estimator$Beta_2
        edge_temp <- EdgeDectect(block_pad, k, n_pad, ls_beta1, ls_beta2, 0.05)
        # record final edge matrix
        result_matrix[
          (j + overlap_width / 2):(j + block_size - 1 - overlap_width / 2),
          (i + overlap_height / 2):(i + block_size - 1 - overlap_height / 2)
        ] <-
          edge_temp[
            (1 + (kk + k)):(block_size + (kk + k)),
            (1 + (kk + k)):(block_size + (kk + k))
          ][
            (1 + overlap_width / 2):(1 + block_size - 1 - overlap_width / 2),
            (1 + overlap_height / 2):(1 + block_size - 1 - overlap_height / 2)
          ]
        # record final recovery matrix
        recover_matrix[
          (j + overlap_width / 2):(j + block_size - 1 - overlap_width / 2),
          (i + overlap_height / 2):(i + block_size - 1 - overlap_height / 2)
        ] <-
          RecoverSurface(n_pad, k1, block_pad, edge_temp)[
            (1 + (kk + k)):(block_size + (kk + k)),
            (1 + (kk + k)):(block_size + (kk + k))
          ][
            (1 + overlap_width / 2):(1 + block_size - 1 - overlap_width / 2),
            (1 + overlap_height / 2):(1 + block_size - 1 - overlap_height / 2)
          ]
      }
    }
    p <- p + 1
    if (i == 1) {
      cat(
        "Block # ", p, ":Horizontally last block width from", width - block_size + 1, "to", width,
        "height from", i, "to", (i + block_size - 1), "\n"
      )
      block_last_h <- z[(width - block_size + 1):width, i:(i + block_size - 1)]
      # blocks = append(blocks, list(block_last_h))
      block_pad <- BoundaryPadding(block_last_h, k)
      ls_estimator <- LSEstimate(x, y, block_pad, n_pad, h_pad)
      ls_beta1 <- ls_estimator$Beta_1
      ls_beta2 <- ls_estimator$Beta_2
      edge_temp <- EdgeDectect(block_pad, k, n_pad, ls_beta1, ls_beta2, 0.05)
      # record final edge matrix
      result_matrix[
        (width - block_size + overlap_width / 2 + 1):width,
        i:(i + block_size - 1 - overlap_height / 2)
      ] <-
        edge_temp[
          (1 + (kk + k)):(block_size + (kk + k)),
          (1 + (kk + k)):(block_size + (kk + k))
        ][
          (1 + overlap_width / 2):block_size,
          i:(i + block_size - 1 - overlap_height / 2)
        ]
      # record final recovery matrix
      recover_matrix[
        (width - block_size + overlap_width / 2 + 1):width,
        i:(i + block_size - 1 - overlap_height / 2)
      ] <-
        RecoverSurface(n_pad, k1, block_pad, edge_temp)[
          (1 + (kk + k)):(block_size + (kk + k)),
          (1 + (kk + k)):(block_size + (kk + k))
        ][
          (1 + overlap_width / 2):block_size,
          i:(i + block_size - 1 - overlap_height / 2)
        ]
    } else {
      cat(
        "Block # ", p, ":Horizontally last block width from", width - block_size + 1, "to", width,
        "height from", i, "to", (i + block_size - 1), "\n"
      )
      block_last_h <- z[(width - block_size + 1):width, i:(i + block_size - 1)]
      # blocks = append(blocks, list(block_last_h))
      block_pad <- BoundaryPadding(block_last_h, k)
      ls_estimator <- LSEstimate(x, y, block_pad, n_pad, h_pad)
      ls_beta1 <- ls_estimator$Beta_1
      ls_beta2 <- ls_estimator$Beta_2
      edge_temp <- EdgeDectect(block_pad, k, n_pad, ls_beta1, ls_beta2, 0.05)
      # record final edge matrix
      result_matrix[
        (width - block_size + overlap_width / 2 + 1):width,
        (i + overlap_height / 2):(i + block_size - 1 - overlap_height / 2)
      ] <-
        edge_temp[
          (1 + (kk + k)):(block_size + (kk + k)),
          (1 + (kk + k)):(block_size + (kk + k))
        ][
          (1 + overlap_width / 2):block_size,
          (1 + overlap_height / 2):(1 + block_size - 1 - overlap_height / 2)
        ]
      # record final recovery matrix
      recover_matrix[
        (width - block_size + overlap_width / 2 + 1):width,
        (i + overlap_height / 2):(i + block_size - 1 - overlap_height / 2)
      ] <-
        RecoverSurface(n_pad, k1, block_pad, edge_temp)[
          (1 + (kk + k)):(block_size + (kk + k)),
          (1 + (kk + k)):(block_size + (kk + k))
        ][
          (1 + overlap_width / 2):block_size,
          (1 + overlap_height / 2):(1 + block_size - 1 - overlap_height / 2)
        ]
    }
  }
  for (j in seq(1, 1 + (n_blocks_width - 1) * non_overlap_width, length.out = n_blocks_width)) {
    p <- p + 1
    if (j == 1) {
      cat(
        "Block # ", p, ":Vertically last block width from", j, "to", j + block_size - 1,
        "height from", height - block_size + 1, "to", height, "\n"
      )
      block <- z[j:(j + block_size - 1), (height - block_size + 1):height]
      # blocks = append(blocks, list(block))
      block_pad <- BoundaryPadding(block, k)
      ls_estimator <- LSEstimate(x, y, block_pad, n_pad, h_pad)
      ls_beta1 <- ls_estimator$Beta_1
      ls_beta2 <- ls_estimator$Beta_2
      edge_temp <- EdgeDectect(block_pad, k, n_pad, ls_beta1, ls_beta2, 0.05)
      # record final edge matrix
      result_matrix[
        j:(j + block_size - 1 - overlap_width / 2),
        (height - block_size + overlap_height / 2 + 1):height
      ] <-
        edge_temp[
          (1 + (kk + k)):(block_size + (kk + k)),
          (1 + (kk + k)):(block_size + (kk + k))
        ][
          1:(1 + block_size - 1 - overlap_width / 2),
          (1 + overlap_height / 2):block_size
        ]
      # record final recovery matrix
      recover_matrix[
        j:(j + block_size - 1 - overlap_width / 2),
        (height - block_size + overlap_height / 2 + 1):height
      ] <-
        RecoverSurface(n_pad, k1, block_pad, edge_temp)[
          (1 + (kk + k)):(block_size + (kk + k)),
          (1 + (kk + k)):(block_size + (kk + k))
        ][
          1:(1 + block_size - 1 - overlap_width / 2),
          (1 + overlap_height / 2):block_size
        ]
    } else {
      cat(
        "Block # ", p, ":Vertically last block width from", j, "to", j + block_size - 1,
        "height from", height - block_size + 1, "to", height, "\n"
      )
      block <- z[j:(j + block_size - 1), (height - block_size + 1):height]
      # blocks = append(blocks, list(block))
      block_pad <- BoundaryPadding(block, k)
      ls_estimator <- LSEstimate(x, y, block_pad, n_pad, h_pad)
      ls_beta1 <- ls_estimator$Beta_1
      ls_beta2 <- ls_estimator$Beta_2
      edge_temp <- EdgeDectect(block_pad, k, n_pad, ls_beta1, ls_beta2, 0.05)
      # record final edge matrix
      result_matrix[
        (j + overlap_width / 2):(j + block_size - 1 - overlap_width / 2),
        (height - block_size + overlap_height / 2 + 1):height
      ] <-
        edge_temp[
          (1 + (kk + k)):(block_size + (kk + k)),
          (1 + (kk + k)):(block_size + (kk + k))
        ][
          (1 + overlap_width / 2):(1 + block_size - 1 - overlap_width / 2),
          (1 + overlap_height / 2):block_size
        ]
      # record final recovery matrix
      recover_matrix[
        (j + overlap_width / 2):(j + block_size - 1 - overlap_width / 2),
        (height - block_size + overlap_height / 2 + 1):height
      ] <-
        RecoverSurface(n_pad, k1, block_pad, edge_temp)[
          (1 + (kk + k)):(block_size + (kk + k)),
          (1 + (kk + k)):(block_size + (kk + k))
        ][
          (1 + overlap_width / 2):(1 + block_size - 1 - overlap_width / 2),
          (1 + overlap_height / 2):block_size
        ]
    }
  }
  # blocks = append(blocks, list(z[(width - block_size + 1):width, (height - block_size + 1):height]))
  block_pad <- BoundaryPadding(z[(width - block_size + 1):width, (height - block_size + 1):height], k)
  ls_estimator <- LSEstimate(x, y, block_pad, n_pad, h_pad)
  ls_beta1 <- ls_estimator$Beta_1
  edge_temp <- EdgeDectect(block_pad, k, n_pad, ls_beta1, ls_beta2, 0.05)
  # record final edge matrix
  result_matrix[
    (width - block_size + overlap_width / 2 + 1):width,
    (height - block_size + overlap_height / 2 + 1):height
  ] <-
    edge_temp[
      (1 + (kk + k)):(block_size + (kk + k)),
      (1 + (kk + k)):(block_size + (kk + k))
    ][
      (1 + overlap_width / 2):block_size,
      (1 + overlap_height / 2):block_size
    ]
  # record final recovery matrix
  recover_matrix[
    (width - block_size + overlap_width / 2 + 1):width,
    (height - block_size + overlap_height / 2 + 1):height
  ] <-
    RecoverSurface(n_pad, k1, block_pad, edge_temp)[
      (1 + (kk + k)):(block_size + (kk + k)),
      (1 + (kk + k)):(block_size + (kk + k))
    ][
      (1 + overlap_width / 2):block_size,
      (1 + overlap_height / 2):block_size
    ]
  return(list("Jumps" = result_matrix, "Recover" = recover_matrix))
  # return(list("Blocks" = blocks, "Jumps" = result_matrix, "Revocer" = recover_matrix))
}
