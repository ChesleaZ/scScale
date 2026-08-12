# Scaling-law fits for external cell representations.

scscale_cell_subspace <- function(x, rank) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (nrow(x) < 3L || ncol(x) < 1L || any(!is.finite(x))) {
    stop("A representation must be a finite cell-by-dimension matrix with at least three cells.", call. = FALSE)
  }
  x <- sweep(x, 2L, colMeans(x), "-")
  keep <- colSums(x^2) > .Machine$double.eps
  x <- x[, keep, drop = FALSE]
  rank <- min(as.integer(rank), ncol(x), nrow(x) - 1L)
  if (!is.finite(rank) || rank < 1L) {
    stop("The requested representation rank is not available.", call. = FALSE)
  }
  decomposition <- eigen(crossprod(x), symmetric = TRUE)
  values <- pmax(decomposition$values, 0)
  positive <- values > max(values) * 1e-12
  if (sum(positive) < rank) {
    stop("The numerical representation rank is below the requested rank.", call. = FALSE)
  }
  scores <- sweep(
    x %*% decomposition$vectors[, seq_len(rank), drop = FALSE],
    2L,
    sqrt(values[seq_len(rank)]),
    "/"
  )
  scores <- qr.Q(qr(scores))[, seq_len(rank), drop = FALSE]
  rownames(scores) <- rownames(x)
  colnames(scores) <- paste0("z", seq_len(rank))
  scores
}

scscale_validate_cell_matrix <- function(x, name) {
  if (!is.matrix(x) && !inherits(x, "Matrix")) {
    stop(name, " must be a cell-by-dimension matrix.", call. = FALSE)
  }
  x <- as.matrix(x)
  if (is.null(rownames(x)) || anyNA(rownames(x)) || any(!nzchar(rownames(x))) ||
      anyDuplicated(rownames(x))) {
    stop(name, " must have unique, nonempty cell IDs as row names.", call. = FALSE)
  }
  x
}

#' Construct a scScale external representation
#'
#' @export
scscale_representation <- function(
  reference,
  umi = NULL,
  umi_values = NULL,
  name = "external",
  rank = ncol(reference),
  rmt_dimension = ncol(reference)
) {
  reference <- scscale_validate_cell_matrix(reference, "reference")
  rank <- as.integer(rank)
  rmt_dimension <- as.numeric(rmt_dimension)
  if (length(rank) != 1L || !is.finite(rank) || rank < 1L ||
      length(rmt_dimension) != 1L || !is.finite(rmt_dimension) || rmt_dimension <= 0) {
    stop("rank and rmt_dimension must be positive scalars.", call. = FALSE)
  }
  if (is.null(umi)) umi <- list()
  if (!is.list(umi)) stop("umi must be a list of cell-by-dimension matrices.", call. = FALSE)
  if (length(umi)) {
    if (is.null(umi_values)) {
      umi_values <- suppressWarnings(as.numeric(names(umi)))
    }
    if (length(umi_values) != length(umi) || any(!is.finite(umi_values) | umi_values < 0)) {
      stop("Provide one nonnegative umi_values entry per UMI representation, or use numeric list names.", call. = FALSE)
    }
  } else {
    umi_values <- numeric(0)
  }
  umi_subspaces <- lapply(seq_along(umi), function(i) {
    matrix_i <- scscale_validate_cell_matrix(umi[[i]], paste0("umi[[", i, "]]"))
    scscale_cell_subspace(matrix_i, rank)
  })
  reference_subspace <- scscale_cell_subspace(reference, rank)
  out <- list(
    name = as.character(name)[[1L]],
    rank = ncol(reference_subspace),
    native_dimension = ncol(reference),
    rmt_dimension = rmt_dimension,
    cells = rownames(reference_subspace),
    reference = reference_subspace,
    umi = umi_subspaces,
    umi_values = as.numeric(umi_values),
    source = "external"
  )
  class(out) <- "scscale_representation"
  out
}

#' Construct a PCA representation from count matrices
#'
#' @export
scscale_pca_representation <- function(
  reference_counts,
  umi_counts = NULL,
  umi_values = NULL,
  n_features = 2000,
  rank = 10,
  min_cells = 10,
  target_depth = 1e4,
  count_transform = c("log1p_cpm", "pearson_residual", "log1p"),
  center = TRUE,
  scale = FALSE
) {
  count_transform <- match.arg(count_transform)
  reference_counts <- counts_matrix(reference_counts)
  if (is.null(rownames(reference_counts)) || is.null(colnames(reference_counts))) {
    stop("reference_counts must have feature and cell names.", call. = FALSE)
  }
  features <- select_hvgs(reference_counts, n_features = n_features, min_cells = min_cells)
  make_embedding <- function(counts, label) {
    counts <- counts_matrix(counts)
    missing <- setdiff(features, rownames(counts))
    if (length(missing)) stop(label, " is missing ", length(missing), " selected features.", call. = FALSE)
    normalized <- scscale_normalize_counts(
      counts[features, , drop = FALSE],
      target_depth = target_depth,
      count_transform = count_transform,
      center = center,
      scale = scale
    )
    vectors <- right_singular_vectors(normalized, r = rank, use_irlba = FALSE)
    rownames(vectors) <- colnames(counts)
    vectors
  }
  reference <- make_embedding(reference_counts, "reference_counts")
  if (is.null(umi_counts)) umi_counts <- list()
  if (!is.list(umi_counts)) stop("umi_counts must be a list of count matrices.", call. = FALSE)
  umi <- lapply(seq_along(umi_counts), function(i) {
    make_embedding(umi_counts[[i]], paste0("umi_counts[[", i, "]]"))
  })
  out <- scscale_representation(
    reference = reference,
    umi = umi,
    umi_values = umi_values,
    name = "PCA",
    rank = rank,
    rmt_dimension = length(features)
  )
  out$source <- "PCA"
  out$features <- features
  out$preprocessing <- list(
    target_depth = target_depth,
    count_transform = count_transform,
    center = center,
    scale = scale,
    min_cells = min_cells
  )
  out
}

scscale_target_subspace <- function(
  target, cells, rank, input, target_depth, count_transform, center, scale
) {
  target_is_matrix <- is.matrix(target) || inherits(target, "Matrix")
  if (target_is_matrix) {
    target <- as.matrix(target)
    if (!is.null(rownames(target)) && all(cells %in% rownames(target))) {
      target <- target[cells, , drop = FALSE]
    } else if (!is.null(colnames(target)) && all(cells %in% colnames(target))) {
      target <- t(target[, cells, drop = FALSE])
    } else {
      stop("The target matrix must identify the representation cells by rows or columns.", call. = FALSE)
    }
  } else {
    target <- align_vector_to_cells(target, cells, "target")
    if (is.numeric(target) || is.integer(target)) {
      target <- matrix(as.numeric(target), ncol = 1L, dimnames = list(cells, "target"))
    } else {
      target <- t(label_target_matrix(target))
      rownames(target) <- cells
    }
  }
  if (input == "counts" && target_is_matrix) {
    normalized <- scscale_normalize_counts(
      t(target),
      target_depth = target_depth,
      count_transform = count_transform,
      center = center,
      scale = scale
    )
    output <- right_singular_vectors(normalized, r = rank, use_irlba = FALSE)
    rownames(output) <- rownames(target)
    return(output)
  }
  scscale_cell_subspace(target, rank)
}

#' Fit UMI scaling for a representation
#'
#' @export
scscale_representation_umi_fit <- function(
  representation,
  target,
  target_rank = 12,
  target_input = c("counts", "normalized", "representation"),
  target_depth = 1e4,
  count_transform = c("log1p_cpm", "pearson_residual", "log1p"),
  center = TRUE,
  scale = FALSE,
  weighting = c("inverse_variance", "equal"),
  n_starts = 1L,
  seed = 1,
  I_infinity_max_factor = 2,
  control = list(maxit = 5000)
) {
  if (!inherits(representation, "scscale_representation")) {
    stop("representation must be a scscale_representation.", call. = FALSE)
  }
  if (length(representation$umi) < 5L) {
    stop("At least five UMI representations are required.", call. = FALSE)
  }
  target_input <- match.arg(target_input)
  count_transform <- match.arg(count_transform)
  weighting <- match.arg(weighting)
  target_subspace <- scscale_target_subspace(
    target, representation$cells, target_rank,
    input = target_input,
    target_depth = target_depth,
    count_transform = count_transform,
    center = center,
    scale = scale
  )
  rows <- lapply(seq_along(representation$umi), function(i) {
    x <- representation$umi[[i]]
    common <- intersect(rownames(x), rownames(target_subspace))
    if (length(common) < 5L) stop("Too few common cells at UMI entry ", i, ".", call. = FALSE)
    value <- subspace_overlap_mi(
      x[common, , drop = FALSE],
      target_subspace[common, , drop = FALSE]
    )
    data.frame(
      U = representation$umi_values[[i]],
      I_empirical = value$mi,
      replicate = stats::ave(
        seq_along(representation$umi_values),
        representation$umi_values,
        FUN = seq_along
      )[[i]],
      n = length(common),
      x_rank = ncol(x),
      y_rank = ncol(target_subspace)
    )
  })
  points <- do.call(rbind, rows)
  n_reference <- nrow(representation$reference)
  fit <- scscale_effective_fit(
    points,
    d = representation$rmt_dimension,
    n_reference = n_reference,
    U_reference = max(points$U),
    weighting = weighting,
    n_starts = n_starts,
    seed = seed,
    I_infinity_max_factor = I_infinity_max_factor,
    control = control
  )
  fit$representation <- representation
  fit$target_subspace <- target_subspace
  fit$umi_points <- points
  fit$representation_name <- representation$name
  class(fit) <- c("scscale_representation_umi_fit", class(fit))
  fit
}

scscale_inverse_sqrt <- function(x, relative_floor = 1e-8) {
  decomposition <- eigen(x, symmetric = TRUE)
  floor <- max(decomposition$values) * relative_floor
  values <- pmax(decomposition$values, floor)
  decomposition$vectors %*% diag(1 / sqrt(values), nrow = length(values)) %*%
    t(decomposition$vectors)
}

scscale_fit_cca_directions <- function(x, y, rank, relative_floor = 1e-8) {
  mean_x <- colMeans(x)
  mean_y <- colMeans(y)
  xc <- sweep(x, 2L, mean_x, "-")
  yc <- sweep(y, 2L, mean_y, "-")
  denominator <- nrow(x) - 1
  wx <- scscale_inverse_sqrt(crossprod(xc) / denominator, relative_floor)
  wy <- scscale_inverse_sqrt(crossprod(yc) / denominator, relative_floor)
  decomposition <- svd(wx %*% (crossprod(xc, yc) / denominator) %*% wy,
    nu = rank, nv = rank)
  list(
    mean_x = mean_x,
    mean_y = mean_y,
    direction_x = wx %*% decomposition$u[, seq_len(rank), drop = FALSE],
    direction_y = wy %*% decomposition$v[, seq_len(rank), drop = FALSE]
  )
}

scscale_heldout_cca_mi <- function(fit, x, y, eps = 1e-10) {
  scores_x <- sweep(x, 2L, fit$mean_x, "-") %*% fit$direction_x
  scores_y <- sweep(y, 2L, fit$mean_y, "-") %*% fit$direction_y
  rho <- diag(stats::cor(scores_x, scores_y))
  rho2 <- pmin(pmax(rho^2, 0), 1 - eps)
  list(mi = -0.5 * sum(log1p(-rho2)), rho2 = rho2)
}

#' Fit held-out cell-number scaling for a representation
#'
#' @export
scscale_representation_cell_fit <- function(
  object,
  cell_numbers,
  replicates = 50,
  evaluation_cells = 5000,
  seed = 1,
  weighting = c("inverse_variance", "equal"),
  relative_floor = 1e-8
) {
  if (!inherits(object, "scscale_representation_umi_fit")) {
    stop("object must be a scscale_representation_umi_fit.", call. = FALSE)
  }
  weighting <- match.arg(weighting)
  x <- object$representation$reference
  y <- object$target_subspace
  common <- intersect(rownames(x), rownames(y))
  x <- x[common, , drop = FALSE]
  y <- y[common, , drop = FALSE]
  evaluation_cells <- as.integer(evaluation_cells)
  replicates <- as.integer(replicates)
  if (evaluation_cells < 5L || evaluation_cells >= length(common) || replicates < 1L) {
    stop("evaluation_cells must leave a nonempty training pool and replicates must be positive.", call. = FALSE)
  }
  set.seed(seed)
  split <- sample.int(length(common))
  evaluation_index <- split[seq_len(evaluation_cells)]
  training_pool <- split[-seq_len(evaluation_cells)]
  anchor_n <- length(training_pool)
  cell_numbers <- sort(unique(c(as.integer(cell_numbers), anchor_n)))
  if (any(!is.finite(cell_numbers) | cell_numbers < 3L | cell_numbers > anchor_n)) {
    stop("cell_numbers must lie between 3 and the training-pool size.", call. = FALSE)
  }
  rank <- min(ncol(x), ncol(y))
  if (min(cell_numbers) <= rank + 1L) {
    stop("Every cell-number level must exceed the fitted CCA rank plus one.", call. = FALSE)
  }
  permutations <- lapply(seq_len(replicates), function(i) sample(training_pool))
  rows <- vector("list", length(cell_numbers) * replicates)
  components <- vector("list", length(rows))
  index <- 0L
  for (replicate in seq_len(replicates)) {
    for (n in cell_numbers) {
      index <- index + 1L
      selected <- permutations[[replicate]][seq_len(n)]
      cca <- scscale_fit_cca_directions(
        x[selected, , drop = FALSE], y[selected, , drop = FALSE], rank, relative_floor
      )
      evaluated <- scscale_heldout_cca_mi(
        cca, x[evaluation_index, , drop = FALSE], y[evaluation_index, , drop = FALSE]
      )
      rows[[index]] <- data.frame(
        n = n,
        U = object$U_reference,
        I_empirical = evaluated$mi,
        replicate = replicate
      )
      components[[index]] <- data.frame(
        n = n,
        replicate = replicate,
        component = seq_len(rank),
        rho2 = evaluated$rho2
      )
    }
  }
  points <- do.call(rbind, rows)
  fit <- object
  umi_n_reference <- fit$n_reference
  fit$n_reference <- anchor_n
  fit <- scscale_effective_cell_calibrate_B(
    fit,
    points,
    weighting = weighting
  )
  fit$cell_points <- points
  fit$cell_components <- do.call(rbind, components)
  fit$cell_design <- list(
    evaluation_cells = common[evaluation_index],
    training_pool = common[training_pool],
    anchor_n = anchor_n,
    umi_n_reference = umi_n_reference,
    rank = rank,
    seed = seed,
    replicates = replicates
  )
  class(fit) <- c("scscale_representation_cell_fit", class(fit))
  fit
}
