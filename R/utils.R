`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

clean_eigenvalues <- function(eigenvalues) {
  sort(eigenvalues[is.finite(eigenvalues) & eigenvalues > 0], decreasing = TRUE)
}
