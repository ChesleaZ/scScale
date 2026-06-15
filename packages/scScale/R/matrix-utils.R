# Small matrix helpers, isolated so the preprocessing code stays readable.

as_dense_matrix <- function(x) {
  if (is.matrix(x)) return(x)
  as.matrix(x)
}

row_sums <- function(x) {
  if (inherits(x, "sparseMatrix") && requireNamespace("Matrix", quietly = TRUE)) {
    return(Matrix::rowSums(x))
  }
  rowSums(as_dense_matrix(x))
}

row_means <- function(x) {
  if (inherits(x, "sparseMatrix") && requireNamespace("Matrix", quietly = TRUE)) {
    return(Matrix::rowMeans(x))
  }
  rowMeans(as_dense_matrix(x))
}

row_means_log1p <- function(x) {
  if (inherits(x, "sparseMatrix")) {
    y <- x
    y@x <- log1p(y@x)
    return(row_means(y))
  }
  rowMeans(log1p(as_dense_matrix(x)))
}

row_means_log1p_square <- function(x) {
  if (inherits(x, "sparseMatrix")) {
    y <- x
    y@x <- log1p(y@x)^2
    return(row_means(y))
  }
  rowMeans(log1p(as_dense_matrix(x))^2)
}
