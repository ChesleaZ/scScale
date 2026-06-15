read_eigenvalues_csv <- function(path) {
  df <- utils::read.csv(path, check.names = FALSE)
  value_col <- intersect(c("eigenvalue", "eigenvalues", "mu", "value"), names(df))
  if (!length(value_col)) {
    stop("Eigenvalue CSV must contain one of: eigenvalue, eigenvalues, mu, value.", call. = FALSE)
  }
  df[[value_col[[1]]]]
}

read_counts_csv <- function(path) {
  df <- utils::read.csv(path, check.names = FALSE)
  first_col <- df[[1]]
  has_feature_column <- !is.numeric(first_col) || anyDuplicated(names(df))
  if (has_feature_column) {
    feature_names <- as.character(first_col)
    mat <- as.matrix(df[, -1, drop = FALSE])
    rownames(mat) <- feature_names
  } else {
    mat <- as.matrix(df)
  }
  mode(mat) <- "numeric"
  mat
}
