################################################################################
# Real multimodal CITE-seq validation of two batch-effect scaling laws
#
# Biological information:
#   I_bio = -1/2 sum_k log(1 - gamma_k)
# where gamma_k are squared canonical correlations between RNA PCA and ADT PCA
# cell-side subspaces. The main fitting target is normalized MI, I_bio / r_eff.
#
# Scaling laws:
#   Fixed total cells, varying batch number:
#     I_bio(m) = I_inf - C log(1 - A / m)
#   Fixed batch number, varying cells per batch:
#     I_bio(s) = I_inf - C log(1 + A / s)
################################################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(hdf5r)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(irlba)
  library(tibble)
  library(purrr)
})

################################################################################
# 0. User settings
################################################################################

run_mode <- Sys.getenv("RUN_MODE", unset = "trend")  # "smoke", "trend", or "full"

data_source <- "bmmc_h5ad"  # "bmmc_h5ad" or "h5seurat"
h5_file <- "/home/ubuntu/data/GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad"
batch_col_preferred <- "batch"
adt_normalization <- "log_scale"  # "log_scale" or "clr"
feature_selection_method <- Sys.getenv("FEATURE_SELECTION", unset = "first")  # "variance" or "first"
out_dir <- if (run_mode == "trend") {
  "/home/ubuntu/bmmc_batch_size_tuned_results/fixed_m3_cell_per_batch_results"
} else {
  "bmmc_batch_scaling_results"
}
base_seed <- 1
clean_output_dir <- TRUE

if (run_mode == "smoke") {
  n_rep <- 2
  min_fit_points <- 3
  n_hvf_grid <- c(300)
  r_embed_grid <- c(5)
  r_target_grid <- c(5)
  n0_fixed_grid <- c(1500)
  m_grid_fixed_n <- c(2, 3, 5, 6, 10, 12, 15)
  m_fixed_grid <- c(6)
  s_grid_fixed_m <- c(100, 150, 200, 300, 400)
} else if (run_mode == "trend") {
  n_rep <- 5
  min_fit_points <- 5
  n_hvf_grid <- c(800)
  r_embed_grid <- c(10)
  r_target_grid <- c(10)
  n0_fixed_grid <- c(1500)
  m_grid_fixed_n <- c(2, 3, 5, 6, 10, 12, 15)
  m_fixed_grid <- c(3)
  s_grid_fixed_m <- c(
    100, 200, 300, 500, 800, 1000, 1500, 2000,
    3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000
  )
} else if (run_mode == "full") {
  n_rep <- 10
  min_fit_points <- 5
  n_hvf_grid <- c(500, 1000, 1500)
  r_embed_grid <- c(5, 10, 15, 20)
  r_target_grid <- c(5, 10, 15, 20)
  n0_fixed_grid <- c(2400, 3600, 4800)
  m_grid_fixed_n <- c(2, 3, 4, 6, 8, 12)
  m_fixed_grid <- c(4, 6)
  s_grid_fixed_m <- c(200, 300, 500, 800, 1000, 1500, 2000, 3000)
} else {
  stop("run_mode must be one of 'smoke', 'trend', or 'full'.")
}

set.seed(base_seed)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
if (clean_output_dir) {
  old_files <- list.files(out_dir, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  if (length(old_files) > 0) {
    unlink(old_files, recursive = TRUE, force = TRUE)
  }
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("Run mode: ", run_mode)
message("Data source: ", data_source)
message("Results directory: ", normalizePath(out_dir, mustWork = FALSE))

################################################################################
# 1. Load real multimodal CITE-seq object
################################################################################

choose_batch_col <- function(meta, preferred) {
  if (!is.null(preferred) && preferred %in% colnames(meta)) {
    return(preferred)
  }

  candidates <- grep(
    "donor|sample|batch|time|orig|library|lane|stim",
    colnames(meta),
    value = TRUE,
    ignore.case = TRUE
  )

  usable <- candidates[vapply(candidates, function(cc) {
    length(unique(na.omit(meta[[cc]]))) >= 2
  }, logical(1))]

  if (length(usable) < 1) {
    stop("No usable batch-related metadata column found.")
  }

  usable[[1]]
}

get_assay_matrix <- function(obj, assay, layer_use = c("counts", "data")) {
  layer_use <- match.arg(layer_use)

  out <- tryCatch(
    GetAssayData(obj, assay = assay, slot = layer_use),
    error = function(e_slot) {
      tryCatch(
        GetAssayData(obj, assay = assay, layer = layer_use),
        error = function(e_layer) NULL
      )
    }
  )

  out
}

extract_first_available <- function(obj, assays, layer_order, label) {
  for (assay in assays) {
    if (!assay %in% Assays(obj)) next
    for (layer_use in layer_order) {
      mat <- get_assay_matrix(obj, assay = assay, layer_use = layer_use)
      if (!is.null(mat) && nrow(mat) > 0 && ncol(mat) > 0) {
        return(list(matrix = mat, source = paste(assay, layer_use)))
      }
    }
  }
  stop("Could not extract ", label, " matrix from assays: ", paste(assays, collapse = ", "))
}

read_h5ad_string_vector <- function(group, field) {
  node <- group[[field]]
  if (inherits(node, "H5Group") && "codes" %in% names(node)) {
    return(as.character(node[["codes"]]$read()))
  }
  as.character(node$read())
}

load_bmmc_h5ad_counts <- function(h5ad_file) {
  message("Loading BMMC CITE-seq h5ad counts: ", h5ad_file)
  h5 <- H5File$new(path.expand(h5ad_file), mode = "r")
  on.exit(h5$close_all(), add = TRUE)

  shape <- h5[["layers/counts"]]$attr_open("shape")$read()
  n_cells <- shape[[1]]
  n_features <- shape[[2]]
  message("H5AD counts shape cells x features: ", n_cells, " x ", n_features)

  feature_codes <- as.integer(h5[["var/feature_types"]]$read())
  if (!all(feature_codes %in% c(0L, 1L))) {
    stop("Expected two h5ad feature type codes for ADT/GEX, found: ", paste(unique(feature_codes), collapse = ", "))
  }
  gex_code <- as.integer(names(which.max(table(feature_codes))))
  adt_code <- setdiff(sort(unique(feature_codes)), gex_code)
  if (length(adt_code) != 1) {
    stop("Could not infer ADT feature code.")
  }
  adt_code <- adt_code[[1]]

  feature_names <- make.unique(as.character(h5[["var/_index"]]$read()))
  cell_names <- make.unique(as.character(h5[["obs/_index"]]$read()))

  meta <- data.frame(row.names = cell_names, check.names = FALSE)
  obs_fields <- intersect(
    c("batch", "Samplename", "Site", "DonorNumber", "DonorID", "cell_type"),
    names(h5[["obs"]])
  )
  for (field in obs_fields) {
    meta[[field]] <- read_h5ad_string_vector(h5[["obs"]], field)
  }

  data <- h5[["layers/counts/data"]]$read()
  indices <- as.integer(h5[["layers/counts/indices"]]$read())
  indptr <- as.integer(h5[["layers/counts/indptr"]]$read())
  cell_index <- rep.int(seq_len(n_cells), diff(indptr))

  feature_code_by_nnz <- feature_codes[indices + 1L]

  make_feature_matrix <- function(code, assay_name) {
    keep <- feature_code_by_nnz == code
    feature_global <- indices[keep] + 1L
    features_use <- which(feature_codes == code)
    feature_local <- match(feature_global, features_use)

    mat <- sparseMatrix(
      i = feature_local,
      j = cell_index[keep],
      x = data[keep],
      dims = c(length(features_use), n_cells),
      dimnames = list(feature_names[features_use], cell_names)
    )
    message(assay_name, " features x cells: ", paste(dim(mat), collapse = " x "),
            "; nnz = ", length(mat@x))
    as(mat, "dgCMatrix")
  }

  list(
    X_rna_all = make_feature_matrix(gex_code, "GEX/RNA"),
    Y_adt_all = make_feature_matrix(adt_code, "ADT"),
    meta_all_use = meta,
    source = "GSE194122 Open Problems NeurIPS 2021 BMMC CITE-seq h5ad"
  )
}

load_h5seurat_counts <- function(h5seurat_file) {
  message("Loading h5Seurat object: ", h5seurat_file)
  obj <- LoadH5Seurat(h5seurat_file, images = FALSE)
  obj <- UpdateSeuratObject(obj)

  message("Object loaded.")
  print(obj)

  message("Assays: ", paste(Assays(obj), collapse = ", "))
  message("Metadata columns:")
  print(colnames(obj@meta.data))

  rna_extracted <- extract_first_available(
    obj,
    assays = c("SCT", "RNA"),
    layer_order = c("counts", "data"),
    label = "RNA"
  )

  adt_extracted <- extract_first_available(
    obj,
    assays = c("ADT"),
    layer_order = c("counts", "data"),
    label = "ADT"
  )

  list(
    X_rna_all = rna_extracted$matrix,
    Y_adt_all = adt_extracted$matrix,
    meta_all_use = obj@meta.data,
    source = paste(rna_extracted$source, "/", adt_extracted$source)
  )
}

loaded_data <- if (identical(data_source, "bmmc_h5ad")) {
  load_bmmc_h5ad_counts(h5_file)
} else if (identical(data_source, "h5seurat")) {
  load_h5seurat_counts(h5_file)
} else {
  stop("Unknown data_source: ", data_source)
}

X_rna_all <- loaded_data$X_rna_all
Y_adt_all <- loaded_data$Y_adt_all
meta_all_use <- loaded_data$meta_all_use

batch_col <- choose_batch_col(meta_all_use, batch_col_preferred)
message("Using batch column: ", batch_col)
print(table(meta_all_use[[batch_col]], useNA = "ifany"))

common_cells <- Reduce(
  intersect,
  list(colnames(X_rna_all), colnames(Y_adt_all), rownames(meta_all_use))
)

if (length(common_cells) < 100) {
  stop("Too few common cells across RNA, ADT, and metadata: ", length(common_cells))
}

X_rna_all <- X_rna_all[, common_cells, drop = FALSE]
Y_adt_all <- Y_adt_all[, common_cells, drop = FALSE]
meta_all_use <- meta_all_use[common_cells, , drop = FALSE]
meta_all_use[[batch_col]] <- as.character(meta_all_use[[batch_col]])

stopifnot(identical(colnames(X_rna_all), colnames(Y_adt_all)))
stopifnot(identical(colnames(X_rna_all), rownames(meta_all_use)))

message("Common cells: ", length(common_cells))
message("RNA features x cells: ", paste(dim(X_rna_all), collapse = " x "))
message("ADT features x cells: ", paste(dim(Y_adt_all), collapse = " x "))

batch_sizes <- sort(table(meta_all_use[[batch_col]]), decreasing = TRUE)
message("Batch sizes:")
print(batch_sizes)

write.csv(
  tibble(batch = names(batch_sizes), n_cells = as.integer(batch_sizes)),
  file = file.path(out_dir, "batch_sizes.csv"),
  row.names = FALSE
)

################################################################################
# 2. Matrix and MI helper functions
################################################################################

select_hvfs_sparse <- function(counts, n_hvf = 1500, min_cells = 20) {
  X <- as(counts, "dgCMatrix")

  detected <- Matrix::rowSums(X > 0)
  if (identical(feature_selection_method, "first")) {
    candidates <- rownames(X)[detected >= min_cells]
    return(candidates[seq_len(min(n_hvf, length(candidates)))])
  }

  X_log <- X
  X_log@x <- log1p(X_log@x)

  X2 <- X_log
  X2@x <- X2@x^2

  feature_mean <- Matrix::rowMeans(X_log)
  feature_second <- Matrix::rowMeans(X2)
  feature_var <- feature_second - feature_mean^2
  feature_var[!is.finite(feature_var)] <- 0
  names(feature_var) <- rownames(X)

  keep <- detected >= min_cells
  feature_var[!keep] <- 0

  out <- names(sort(feature_var, decreasing = TRUE))
  out <- out[feature_var[out] > 0]
  out[seq_len(min(n_hvf, length(out)))]
}

build_rna_matrix <- function(X_rna, features_use) {
  features_use <- intersect(features_use, rownames(X_rna))
  if (length(features_use) < 2) {
    stop("Fewer than two RNA features available after HVF intersection.")
  }

  X <- as.matrix(X_rna[features_use, , drop = FALSE])
  X <- log1p(X)
  X <- scale(t(X), center = TRUE, scale = TRUE)
  X[!is.finite(X)] <- 0
  X <- t(X)

  rownames(X) <- features_use
  colnames(X) <- colnames(X_rna)
  X
}

build_adt_target_matrix <- function(Y_adt, normalization = adt_normalization) {
  A <- as.matrix(Y_adt)

  if (identical(normalization, "clr")) {
    A_use <- apply(A, 2, function(x) {
      log1p(x / exp(mean(log1p(x))))
    })
  } else if (identical(normalization, "log_scale")) {
    A_use <- log1p(A)
  } else {
    stop("Unknown ADT normalization: ", normalization)
  }

  A_cells <- t(A_use)
  A_cells <- scale(A_cells, center = TRUE, scale = TRUE)
  A_cells[!is.finite(A_cells)] <- 0

  Y <- t(A_cells)
  rownames(Y) <- rownames(Y_adt)
  colnames(Y) <- colnames(Y_adt)
  Y
}

get_right_singular_vectors <- function(X, r = 20) {
  X <- as.matrix(X)
  r_eff <- min(r, nrow(X) - 1, ncol(X) - 1)
  if (r_eff < 1) return(NULL)

  sv <- irlba::irlba(
    X,
    nv = r_eff,
    nu = 0,
    maxit = 2000,
    work = max(3 * r_eff, 30)
  )

  Z <- sv$v[, seq_len(r_eff), drop = FALSE]
  rownames(Z) <- colnames(X)
  colnames(Z) <- paste0("z", seq_len(ncol(Z)))
  Z
}

orthonormalize_cell_matrix <- function(Z) {
  Z <- as.matrix(Z)
  Z <- scale(Z, center = TRUE, scale = FALSE)
  Z[!is.finite(Z)] <- 0

  keep_col <- apply(Z, 2, function(x) sd(x) > 1e-12)
  Z <- Z[, keep_col, drop = FALSE]
  if (ncol(Z) < 1) return(NULL)

  qrZ <- qr(Z)
  Q <- qr.Q(qrZ)
  Q <- Q[, seq_len(qrZ$rank), drop = FALSE]
  rownames(Q) <- rownames(Z)
  Q
}

compute_subspace_overlap_mi <- function(Z1, Z2, eps = 1e-12) {
  common_cells <- intersect(rownames(Z1), rownames(Z2))

  if (length(common_cells) < 5) {
    return(list(mi = NA_real_, mi_norm = NA_real_, gamma = numeric(0), r_eff = 0))
  }

  Z1 <- Z1[common_cells, , drop = FALSE]
  Z2 <- Z2[common_cells, , drop = FALSE]

  Q1 <- orthonormalize_cell_matrix(Z1)
  Q2 <- orthonormalize_cell_matrix(Z2)

  if (is.null(Q1) || is.null(Q2)) {
    return(list(mi = NA_real_, mi_norm = NA_real_, gamma = numeric(0), r_eff = 0))
  }

  r_eff <- min(ncol(Q1), ncol(Q2))
  if (r_eff < 1) {
    return(list(mi = NA_real_, mi_norm = NA_real_, gamma = numeric(0), r_eff = 0))
  }

  Cmat <- crossprod(Q1, Q2)
  gamma <- svd(Cmat, nu = 0, nv = 0)$d^2
  gamma <- gamma[seq_len(min(r_eff, length(gamma)))]
  gamma <- pmin(pmax(gamma, 0), 1 - eps)

  mi <- -0.5 * sum(log(1 - gamma))
  list(
    mi = mi,
    mi_norm = mi / length(gamma),
    gamma = gamma,
    r_eff = length(gamma)
  )
}

make_onehot_target <- function(labels, cells) {
  labels <- as.factor(labels)
  Y <- model.matrix(~ 0 + labels)
  rownames(Y) <- cells
  colnames(Y) <- gsub("^labels", "", colnames(Y))
  Y
}

sample_cells_from_batches <- function(meta, batch_col, m_batch, cells_per_batch, seed, batches_use = NULL) {
  set.seed(seed)

  batch_sizes_local <- table(meta[[batch_col]])
  eligible_batches <- names(batch_sizes_local)[batch_sizes_local >= cells_per_batch]

  if (is.null(batches_use)) {
    if (length(eligible_batches) < m_batch) {
      return(NULL)
    }
    batches_use <- sample(eligible_batches, m_batch, replace = FALSE)
  } else {
    batches_use <- as.character(batches_use)
    if (length(batches_use) != m_batch || any(!batches_use %in% eligible_batches)) {
      return(NULL)
    }
  }

  cells_use <- unlist(lapply(batches_use, function(bb) {
    cells_bb <- rownames(meta)[meta[[batch_col]] == bb]
    sample(cells_bb, cells_per_batch, replace = FALSE)
  }), use.names = FALSE)

  list(cells = cells_use, batches = batches_use)
}

################################################################################
# 3. Experiment execution
################################################################################

feature_cache <- new.env(parent = emptyenv())

get_features_for_n_hvf <- function(n_hvf) {
  key <- as.character(n_hvf)
  if (!exists(key, envir = feature_cache, inherits = FALSE)) {
    message("Selecting global HVFs: n_hvf = ", n_hvf)
    assign(
      key,
      select_hvfs_sparse(X_rna_all, n_hvf = n_hvf, min_cells = 20),
      envir = feature_cache
    )
  }
  get(key, envir = feature_cache, inherits = FALSE)
}

choose_fixed_batches_for_rep <- function(m_batch, max_cells_per_batch, rep_id, setting) {
  seed_setting_id <- if ("seed_setting_id" %in% names(setting)) {
    setting$seed_setting_id
  } else {
    setting$setting_id
  }

  seed <- base_seed + seed_setting_id * 100000L + rep_id * 1000L + m_batch * 10L + max_cells_per_batch
  set.seed(seed)

  eligible_batches <- names(batch_sizes)[batch_sizes >= max_cells_per_batch]
  if (length(eligible_batches) < m_batch) {
    return(NULL)
  }

  sample(eligible_batches, m_batch, replace = FALSE)
}

run_one_batch_design <- function(setting, m_batch, cells_per_batch, rep_id, experiment_label, batches_use = NULL) {
  seed_setting_id <- if ("seed_setting_id" %in% names(setting)) {
    setting$seed_setting_id
  } else {
    setting$setting_id
  }

  seed <- base_seed +
    seed_setting_id * 100000L +
    rep_id * 1000L +
    m_batch * 10L +
    cells_per_batch

  sampled <- sample_cells_from_batches(
    meta = meta_all_use,
    batch_col = batch_col,
    m_batch = m_batch,
    cells_per_batch = cells_per_batch,
    seed = seed,
    batches_use = batches_use
  )

  if (is.null(sampled)) {
    return(tibble(
      setting_id = setting$setting_id,
      experiment = experiment_label,
      rep = rep_id,
      m_batch = m_batch,
      cells_per_batch = cells_per_batch,
      n_cells = m_batch * cells_per_batch,
      batch_col = batch_col,
      batches_used = NA_character_,
      status = "skipped_infeasible_sampling",
      message = "Not enough batches with requested cells per batch.",
      I_bio = NA_real_,
      I_bio_norm = NA_real_,
      r_bio_eff = NA_integer_
    ))
  }

  tryCatch({
    cells_use <- sampled$cells
    X_rna <- X_rna_all[, cells_use, drop = FALSE]
    Y_adt <- Y_adt_all[, cells_use, drop = FALSE]
    meta_sub <- meta_all_use[cells_use, , drop = FALSE]

    features_global <- get_features_for_n_hvf(setting$n_hvf)
    X_mat <- build_rna_matrix(X_rna = X_rna, features_use = features_global)
    Y_mat <- build_adt_target_matrix(Y_adt)

    Z_rna <- get_right_singular_vectors(X = X_mat, r = setting$r_embed)
    Z_adt <- get_right_singular_vectors(X = Y_mat, r = setting$r_target)

    if (is.null(Z_rna) || is.null(Z_adt)) {
      stop("PCA/SVD returned NULL embedding.")
    }

    bio_out <- compute_subspace_overlap_mi(Z1 = Z_rna, Z2 = Z_adt)

    tibble(
      setting_id = setting$setting_id,
      experiment = experiment_label,
      rep = rep_id,
      m_batch = m_batch,
      cells_per_batch = cells_per_batch,
      n_cells = length(cells_use),
      batch_col = batch_col,
      batches_used = paste(sampled$batches, collapse = ","),
      status = "ok",
      message = NA_character_,
      I_bio = bio_out$mi,
      I_bio_norm = bio_out$mi_norm,
      r_bio_eff = bio_out$r_eff
    )
  }, error = function(e) {
    tibble(
      setting_id = setting$setting_id,
      experiment = experiment_label,
      rep = rep_id,
      m_batch = m_batch,
      cells_per_batch = cells_per_batch,
      n_cells = m_batch * cells_per_batch,
      batch_col = batch_col,
      batches_used = if (!is.null(sampled$batches)) paste(sampled$batches, collapse = ",") else NA_character_,
      status = "failed",
      message = conditionMessage(e),
      I_bio = NA_real_,
      I_bio_norm = NA_real_,
      r_bio_eff = NA_integer_
    )
  })
}

summarize_results <- function(res) {
  if (nrow(res) < 1) return(tibble())

  res %>%
    filter(status == "ok") %>%
    group_by(setting_id, experiment, m_batch, cells_per_batch, n_cells) %>%
    summarise(
      mean_I_bio = mean(I_bio, na.rm = TRUE),
      sd_I_bio = sd(I_bio, na.rm = TRUE),
      se_I_bio = sd_I_bio / sqrt(n()),
      mean_I_bio_norm = mean(I_bio_norm, na.rm = TRUE),
      sd_I_bio_norm = sd(I_bio_norm, na.rm = TRUE),
      se_I_bio_norm = sd_I_bio_norm / sqrt(n()),
      n_rep_observed = n(),
      .groups = "drop"
    ) %>%
    arrange(setting_id, experiment, m_batch, cells_per_batch)
}

feasible_fixed_n_design <- function(n0, m_grid) {
  rows <- lapply(m_grid, function(mm) {
    if (n0 %% mm != 0) {
      return(tibble(
        m_batch = mm,
        cells_per_batch = NA_integer_,
        n_cells = n0,
        feasible = FALSE,
        reason = "n0_not_divisible_by_m"
      ))
    }

    ss <- as.integer(n0 / mm)
    enough <- sum(batch_sizes >= ss) >= mm
    tibble(
      m_batch = mm,
      cells_per_batch = ss,
      n_cells = n0,
      feasible = enough,
      reason = ifelse(enough, "ok", "not_enough_batches_with_required_size")
    )
  })

  bind_rows(rows)
}

feasible_fixed_m_design <- function(m_fixed, s_grid) {
  rows <- lapply(s_grid, function(ss) {
    enough <- sum(batch_sizes >= ss) >= m_fixed
    tibble(
      m_batch = m_fixed,
      cells_per_batch = as.integer(ss),
      n_cells = as.integer(m_fixed * ss),
      feasible = enough,
      reason = ifelse(enough, "ok", "not_enough_batches_with_required_size")
    )
  })

  bind_rows(rows)
}

make_settings <- function() {
  expand.grid(
    n_hvf = n_hvf_grid,
    r_embed = r_embed_grid,
    r_target = r_target_grid,
    n0_fixed = n0_fixed_grid,
    m_fixed = m_fixed_grid,
    KEEP.OUT.ATTRS = FALSE
  ) %>%
    as_tibble() %>%
    mutate(setting_id = row_number()) %>%
    select(setting_id, everything()) %>%
    mutate(
      seed_setting_id = ifelse(run_mode == "trend", 24L, setting_id),
      n_rep_requested = n_rep,
      run_mode = run_mode,
      batch_col = batch_col,
      adt_normalization = adt_normalization,
      feature_selection_method = feature_selection_method,
      m_grid_fixed_n = paste(m_grid_fixed_n, collapse = ";"),
      s_grid_fixed_m = paste(s_grid_fixed_m, collapse = ";")
    )
}

settings <- make_settings()
write.csv(settings, file.path(out_dir, "settings_grid.csv"), row.names = FALSE)
message("Number of candidate settings: ", nrow(settings))

################################################################################
# 4. Fitting helpers
################################################################################

metric_for_target <- function(dat, target = c("norm", "raw")) {
  target <- match.arg(target)
  if (target == "norm") {
    dat %>% mutate(I_fit = mean_I_bio_norm, se_fit = se_I_bio_norm)
  } else {
    dat %>% mutate(I_fit = mean_I_bio, se_fit = se_I_bio)
  }
}

calc_fit_stats <- function(y, pred) {
  resid <- y - pred
  denom <- sum((y - mean(y))^2)
  R2 <- ifelse(denom > 0, 1 - sum(resid^2) / denom, NA_real_)
  tibble(
    R2 = R2,
    RMSE = sqrt(mean(resid^2)),
    MAE = mean(abs(resid))
  )
}

try_nls <- function(expr, dat, starts, lower, upper) {
  best_fit <- NULL
  best_rmse <- Inf
  errors <- character(0)

  for (ii in seq_along(starts)) {
    fit <- tryCatch(
      nls(
        expr,
        data = dat,
        start = starts[[ii]],
        algorithm = "port",
        lower = lower,
        upper = upper,
        control = nls.control(maxiter = 1000, warnOnly = TRUE)
      ),
      error = function(e) {
        errors <<- c(errors, conditionMessage(e))
        NULL
      }
    )

    if (!is.null(fit)) {
      pred <- as.numeric(predict(fit, newdata = dat))
      rmse <- sqrt(mean((dat$I_fit - pred)^2))
      if (is.finite(rmse) && rmse < best_rmse) {
        best_fit <- fit
        best_rmse <- rmse
      }
    }
  }

  if (is.null(best_fit)) {
    return(list(fit = NULL, error = paste(unique(errors), collapse = " | ")))
  }

  list(fit = best_fit, error = NA_character_)
}

fit_batch_number_law <- function(summary_df, target = c("norm", "raw")) {
  target <- match.arg(target)
  dat <- summary_df %>%
    filter(experiment == "fixed_n_vary_m") %>%
    metric_for_target(target) %>%
    filter(is.finite(I_fit), m_batch > 0) %>%
    arrange(m_batch)

  n_points <- nrow(dat)
  if (n_points < min_fit_points) {
    return(list(
      ok = FALSE,
      data = dat,
      pred = tibble(),
      message = paste0("Too few points for fit: ", n_points, " < ", min_fit_points)
    ))
  }

  min_m <- min(dat$m_batch)
  max_abs_a <- 0.95 * min_m
  y_rng <- diff(range(dat$I_fit))
  c_start <- max(abs(y_rng), stats::sd(dat$I_fit), 1e-3)

  starts <- list(
    list(I_inf = mean(dat$I_fit), C = c_start, A = 0.25 * min_m),
    list(I_inf = min(dat$I_fit), C = c_start, A = -0.25 * min_m),
    list(I_inf = max(dat$I_fit), C = 0.5 * c_start, A = 0.5 * min_m),
    list(I_inf = mean(dat$I_fit), C = 2 * c_start, A = -0.5 * min_m)
  )

  fit_out <- try_nls(
    I_fit ~ I_inf - C * log(1 - A / m_batch),
    dat = dat,
    starts = starts,
    lower = c(I_inf = -Inf, C = 0, A = -max_abs_a),
    upper = c(I_inf = Inf, C = Inf, A = max_abs_a)
  )

  if (is.null(fit_out$fit)) {
    return(list(ok = FALSE, data = dat, pred = tibble(), message = fit_out$error))
  }

  fit <- fit_out$fit
  co <- coef(fit)
  dat <- dat %>%
    mutate(I_pred = as.numeric(predict(fit, newdata = .)), resid = I_fit - I_pred)

  pred <- tibble(
    m_batch = seq(min(dat$m_batch), max(dat$m_batch), length.out = 200)
  ) %>%
    mutate(I_pred = co[["I_inf"]] - co[["C"]] * log(1 - co[["A"]] / m_batch))

  list(
    ok = TRUE,
    fit = fit,
    data = dat,
    pred = pred,
    stats = calc_fit_stats(dat$I_fit, dat$I_pred),
    coef = co,
    message = "ok"
  )
}

fit_cells_per_batch_law <- function(summary_df, target = c("norm", "raw")) {
  target <- match.arg(target)
  dat <- summary_df %>%
    filter(experiment == "fixed_m_vary_s") %>%
    metric_for_target(target) %>%
    filter(is.finite(I_fit), cells_per_batch > 0) %>%
    arrange(cells_per_batch)

  n_points <- nrow(dat)
  if (n_points < min_fit_points) {
    return(list(
      ok = FALSE,
      data = dat,
      pred = tibble(),
      message = paste0("Too few points for fit: ", n_points, " < ", min_fit_points)
    ))
  }

  y_rng <- diff(range(dat$I_fit))
  c_start <- max(abs(y_rng), stats::sd(dat$I_fit), 1e-3)
  s_med <- median(dat$cells_per_batch)
  s_max <- max(dat$cells_per_batch)

  starts <- list(
    list(I_inf = max(dat$I_fit), C = c_start, A = s_med),
    list(I_inf = mean(dat$I_fit), C = 0.5 * c_start, A = 0.5 * s_med),
    list(I_inf = max(dat$I_fit) + c_start, C = 2 * c_start, A = 2 * s_med),
    list(I_inf = min(dat$I_fit), C = c_start, A = 0.1 * s_med),
    list(I_inf = max(dat$I_fit), C = 0.25 * c_start, A = 5 * s_med),
    list(I_inf = max(dat$I_fit) + 0.5 * c_start, C = 4 * c_start, A = 0.25 * s_med)
  )

  fit_out <- try_nls(
    I_fit ~ I_inf - C * log(1 + A / cells_per_batch),
    dat = dat,
    starts = starts,
    lower = c(I_inf = -Inf, C = 0, A = 0),
    upper = c(I_inf = Inf, C = Inf, A = 20 * s_max)
  )

  if (is.null(fit_out$fit)) {
    return(list(ok = FALSE, data = dat, pred = tibble(), message = fit_out$error))
  }

  fit <- fit_out$fit
  co <- coef(fit)
  dat <- dat %>%
    mutate(I_pred = as.numeric(predict(fit, newdata = .)), resid = I_fit - I_pred)

  pred <- tibble(
    cells_per_batch = seq(
      min(dat$cells_per_batch),
      max(dat$cells_per_batch),
      length.out = 200
    )
  ) %>%
    mutate(I_pred = co[["I_inf"]] - co[["C"]] * log(1 + co[["A"]] / cells_per_batch))

  list(
    ok = TRUE,
    fit = fit,
    data = dat,
    pred = pred,
    stats = calc_fit_stats(dat$I_fit, dat$I_pred),
    coef = co,
    message = "ok"
  )
}

fit_row_from_result <- function(fit_out, setting, law, target, summary_df) {
  n_points <- nrow(fit_out$data)
  mean_se <- if (n_points > 0) mean(fit_out$data$se_fit, na.rm = TRUE) else NA_real_
  smooth_score <- if (n_points > 2 && all(is.finite(fit_out$data$I_fit))) {
    stats::cor(seq_len(n_points), fit_out$data$I_fit, method = "spearman")
  } else {
    NA_real_
  }

  base <- tibble(
    setting_id = setting$setting_id,
    law = law,
    target = target,
    fit_status = ifelse(fit_out$ok, "ok", "failed"),
    fit_message = fit_out$message,
    n_points = n_points,
    n_rep_requested = n_rep,
    mean_empirical_se = mean_se,
    smooth_spearman_by_x = smooth_score
  )

  if (!fit_out$ok) {
    return(base %>%
      mutate(I_inf = NA_real_, C = NA_real_, A = NA_real_, R2 = NA_real_, RMSE = NA_real_, MAE = NA_real_))
  }

  base %>%
    mutate(
      I_inf = unname(fit_out$coef[["I_inf"]]),
      C = unname(fit_out$coef[["C"]]),
      A = unname(fit_out$coef[["A"]]),
      R2 = fit_out$stats$R2,
      RMSE = fit_out$stats$RMSE,
      MAE = fit_out$stats$MAE
    )
}

################################################################################
# 5. Plotting helpers
################################################################################

format_fit_title <- function(prefix, fit_row) {
  paste0(
    "R2 = ", round(fit_row$R2, 3)
  )
}

plot_bio_fit_m <- function(fit_out, fit_row, target = c("norm", "raw")) {
  target <- match.arg(target)
  y_lab <- ifelse(target == "norm", "Normalized biological information", "Biological information")
  mean_col <- ifelse(target == "norm", "mean_I_bio_norm", "mean_I_bio")

  ggplot(fit_out$data, aes(x = m_batch, y = .data[[mean_col]])) +
    geom_point(aes(color = "Empirical MI"), size = 2.4, alpha = 0.95) +
    geom_line(data = fit_out$pred, aes(x = m_batch, y = I_pred, color = "RMT scaling law"), linewidth = 1.0) +
    scale_x_continuous(breaks = sort(unique(fit_out$data$m_batch))) +
    scale_color_manual(
      name = NULL,
      values = c("Empirical MI" = "black", "RMT scaling law" = "#1f77b4")
    ) +
    labs(
      x = "Number of batches (m)",
      y = y_lab,
      title = format_fit_title("Batch-number scaling", fit_row)
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 10),
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = scales::alpha("white", 0.85), color = "grey70", linewidth = 0.25),
      legend.key = element_blank()
    )
}

plot_bio_fit_s <- function(fit_out, fit_row, target = c("norm", "raw")) {
  target <- match.arg(target)
  y_lab <- ifelse(target == "norm", "Normalized biological information", "Biological information")
  mean_col <- ifelse(target == "norm", "mean_I_bio_norm", "mean_I_bio")

  ggplot(fit_out$data, aes(x = cells_per_batch, y = .data[[mean_col]])) +
    geom_point(aes(color = "Empirical MI"), size = 2.4, alpha = 0.95) +
    geom_line(data = fit_out$pred, aes(x = cells_per_batch, y = I_pred, color = "RMT scaling law"), linewidth = 1.0) +
    scale_x_log10(
      breaks = c(100, 300, 1000, 3000, 6000, 10000),
      labels = scales::comma_format()
    ) +
    scale_color_manual(
      name = NULL,
      values = c("Empirical MI" = "black", "RMT scaling law" = "#1f77b4")
    ) +
    labs(
      x = "Cells per batch (s)",
      y = y_lab,
      title = format_fit_title("Batch-size scaling", fit_row)
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
      plot.title = element_text(size = 10),
      legend.position = c(0.98, 0.02),
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = scales::alpha("white", 0.85), color = "grey70", linewidth = 0.25),
      legend.key = element_blank()
    )
}

save_plot_pair <- function(filename_base, plot_obj, width = 5.5, height = 3.0) {
  ggsave(file.path(out_dir, paste0(filename_base, ".png")), plot_obj, width = width, height = height, dpi = 300)
  ggsave(file.path(out_dir, paste0(filename_base, ".pdf")), plot_obj, width = width, height = height)
}

################################################################################
# 6. Grid search
################################################################################

all_raw_results <- list()
all_summaries <- list()
all_fit_rows <- list()
all_feasibility <- list()
setting_status_rows <- list()

for (ii in seq_len(nrow(settings))) {
  setting <- settings[ii, ]
  message("\n================================================================")
  message("Setting ", setting$setting_id, " / ", nrow(settings),
          ": n_hvf=", setting$n_hvf,
          ", r_embed=", setting$r_embed,
          ", r_target=", setting$r_target,
          ", n0_fixed=", setting$n0_fixed,
          ", m_fixed=", setting$m_fixed)
  message("================================================================")

  fixed_n_design <- feasible_fixed_n_design(setting$n0_fixed, m_grid_fixed_n) %>%
    mutate(setting_id = setting$setting_id, experiment = "fixed_n_vary_m")
  fixed_m_design <- feasible_fixed_m_design(setting$m_fixed, s_grid_fixed_m) %>%
    mutate(setting_id = setting$setting_id, experiment = "fixed_m_vary_s")

  design_all <- bind_rows(fixed_n_design, fixed_m_design) %>%
    select(setting_id, experiment, everything())
  all_feasibility[[ii]] <- design_all

  message("Fixed-n feasible x-axis points: ", sum(fixed_n_design$feasible), " / ", nrow(fixed_n_design))
  if (any(!fixed_n_design$feasible)) {
    print(fixed_n_design %>% filter(!feasible) %>% select(m_batch, cells_per_batch, reason))
  }
  message("Fixed-m feasible x-axis points: ", sum(fixed_m_design$feasible), " / ", nrow(fixed_m_design))
  if (any(!fixed_m_design$feasible)) {
    print(fixed_m_design %>% filter(!feasible) %>% select(cells_per_batch, reason))
  }

  fixed_n_ok <- fixed_n_design %>% filter(feasible)
  fixed_m_ok <- fixed_m_design %>% filter(feasible)

  if (nrow(fixed_n_ok) < min_fit_points || nrow(fixed_m_ok) < min_fit_points) {
    msg <- paste0(
      "Skipped setting because feasible points are below minimum: fixed_n=",
      nrow(fixed_n_ok), ", fixed_m=", nrow(fixed_m_ok), ", required=", min_fit_points
    )
    message(msg)

    setting_status_rows[[ii]] <- tibble(
      setting_id = setting$setting_id,
      setting_status = "skipped_infeasible_grid",
      setting_message = msg
    )

    all_fit_rows[[ii]] <- bind_rows(
      tibble(setting_id = setting$setting_id, law = "fixed_n_vary_m", target = c("norm", "raw"), fit_status = "failed", fit_message = msg, n_points = nrow(fixed_n_ok), n_rep_requested = n_rep, mean_empirical_se = NA_real_, smooth_spearman_by_x = NA_real_, I_inf = NA_real_, C = NA_real_, A = NA_real_, R2 = NA_real_, RMSE = NA_real_, MAE = NA_real_),
      tibble(setting_id = setting$setting_id, law = "fixed_m_vary_s", target = c("norm", "raw"), fit_status = "failed", fit_message = msg, n_points = nrow(fixed_m_ok), n_rep_requested = n_rep, mean_empirical_se = NA_real_, smooth_spearman_by_x = NA_real_, I_inf = NA_real_, C = NA_real_, A = NA_real_, R2 = NA_real_, RMSE = NA_real_, MAE = NA_real_)
    )
    next
  }

  setting_status_rows[[ii]] <- tibble(
    setting_id = setting$setting_id,
    setting_status = "attempted",
    setting_message = "ok"
  )

  raw_rows <- list()
  row_idx <- 1L

  for (rep_id in seq_len(n_rep)) {
    message("  Replicate ", rep_id, " / ", n_rep)
    fixed_m_batches <- choose_fixed_batches_for_rep(
      m_batch = setting$m_fixed,
      max_cells_per_batch = max(fixed_m_ok$cells_per_batch),
      rep_id = rep_id,
      setting = setting
    )

    for (jj in seq_len(nrow(fixed_n_ok))) {
      mm <- fixed_n_ok$m_batch[[jj]]
      ss <- fixed_n_ok$cells_per_batch[[jj]]
      message("    fixed_n_vary_m: m=", mm, ", s=", ss, ", n=", setting$n0_fixed)
      raw_rows[[row_idx]] <- run_one_batch_design(setting, mm, ss, rep_id, "fixed_n_vary_m")
      row_idx <- row_idx + 1L
      gc(verbose = FALSE)
    }

    for (jj in seq_len(nrow(fixed_m_ok))) {
      mm <- fixed_m_ok$m_batch[[jj]]
      ss <- fixed_m_ok$cells_per_batch[[jj]]
      message("    fixed_m_vary_s: m=", mm, ", s=", ss, ", n=", mm * ss)
      raw_rows[[row_idx]] <- run_one_batch_design(
        setting,
        mm,
        ss,
        rep_id,
        "fixed_m_vary_s",
        batches_use = fixed_m_batches
      )
      row_idx <- row_idx + 1L
      gc(verbose = FALSE)
    }
  }

  raw_df <- bind_rows(raw_rows)
  summary_df <- summarize_results(raw_df)

  all_raw_results[[ii]] <- raw_df
  all_summaries[[ii]] <- summary_df

  write.csv(raw_df, file.path(out_dir, sprintf("setting_%04d_raw_results.csv", setting$setting_id)), row.names = FALSE)
  write.csv(summary_df, file.path(out_dir, sprintf("setting_%04d_summary.csv", setting$setting_id)), row.names = FALSE)

  fit_m_norm <- fit_batch_number_law(summary_df, target = "norm")
  fit_s_norm <- fit_cells_per_batch_law(summary_df, target = "norm")
  fit_m_raw <- fit_batch_number_law(summary_df, target = "raw")
  fit_s_raw <- fit_cells_per_batch_law(summary_df, target = "raw")

  all_fit_rows[[ii]] <- bind_rows(
    fit_row_from_result(fit_m_norm, setting, "fixed_n_vary_m", "norm", summary_df),
    fit_row_from_result(fit_s_norm, setting, "fixed_m_vary_s", "norm", summary_df),
    fit_row_from_result(fit_m_raw, setting, "fixed_n_vary_m", "raw", summary_df),
    fit_row_from_result(fit_s_raw, setting, "fixed_m_vary_s", "raw", summary_df)
  )

  write.csv(bind_rows(all_fit_rows[[ii]]), file.path(out_dir, sprintf("setting_%04d_fit_results.csv", setting$setting_id)), row.names = FALSE)
}

raw_results_all <- bind_rows(all_raw_results)
summary_all <- bind_rows(all_summaries)
fit_results_all <- bind_rows(all_fit_rows) %>%
  left_join(settings %>% select(-n_rep_requested), by = "setting_id")
feasibility_all <- bind_rows(all_feasibility) %>%
  left_join(settings, by = "setting_id")
setting_status <- bind_rows(setting_status_rows) %>%
  left_join(settings, by = "setting_id")

write.csv(raw_results_all, file.path(out_dir, "all_replicate_raw_results.csv"), row.names = FALSE)
write.csv(summary_all, file.path(out_dir, "all_empirical_summaries.csv"), row.names = FALSE)
write.csv(fit_results_all, file.path(out_dir, "all_fit_results.csv"), row.names = FALSE)
write.csv(feasibility_all, file.path(out_dir, "all_feasibility_results.csv"), row.names = FALSE)
write.csv(setting_status, file.path(out_dir, "setting_status.csv"), row.names = FALSE)

################################################################################
# 7. Ranking and best-setting plots
################################################################################

norm_fit_wide <- fit_results_all %>%
  filter(target == "norm") %>%
  select(
    setting_id,
    law,
    fit_status,
    fit_message,
    n_points,
    mean_empirical_se,
    R2,
    RMSE,
    MAE,
    I_inf,
    C,
    A
  ) %>%
  pivot_wider(
    names_from = law,
    values_from = c(fit_status, fit_message, n_points, mean_empirical_se, R2, RMSE, MAE, I_inf, C, A),
    names_sep = "__"
  )

ranked_settings <- settings %>%
  left_join(setting_status %>% select(setting_id, setting_status, setting_message), by = "setting_id") %>%
  left_join(norm_fit_wide, by = "setting_id") %>%
  mutate(
    both_fits_ok = fit_status__fixed_n_vary_m == "ok" & fit_status__fixed_m_vary_s == "ok",
    mean_R2 = rowMeans(cbind(R2__fixed_n_vary_m, R2__fixed_m_vary_s), na.rm = TRUE),
    mean_RMSE = rowMeans(cbind(RMSE__fixed_n_vary_m, RMSE__fixed_m_vary_s), na.rm = TRUE),
    mean_MAE = rowMeans(cbind(MAE__fixed_n_vary_m, MAE__fixed_m_vary_s), na.rm = TRUE),
    mean_SE = rowMeans(cbind(mean_empirical_se__fixed_n_vary_m, mean_empirical_se__fixed_m_vary_s), na.rm = TRUE),
    min_points = pmin(n_points__fixed_n_vary_m, n_points__fixed_m_vary_s, na.rm = TRUE),
    rank_score = ifelse(
      both_fits_ok,
      mean_R2 - 0.1 * mean_RMSE - 0.05 * mean_MAE - 0.05 * mean_SE + 0.01 * pmin(min_points, 10),
      -Inf
    )
  ) %>%
  arrange(desc(both_fits_ok), desc(rank_score), desc(mean_R2), mean_RMSE, mean_SE)

write.csv(ranked_settings, file.path(out_dir, "ranked_setting_results.csv"), row.names = FALSE)

if (nrow(ranked_settings) > 0 && any(ranked_settings$both_fits_ok, na.rm = TRUE)) {
  best_setting_id <- ranked_settings$setting_id[which(ranked_settings$both_fits_ok)[1]]
  best_setting <- settings %>% filter(setting_id == best_setting_id)
  best_summary <- summary_all %>% filter(setting_id == best_setting_id)
  best_fits <- fit_results_all %>% filter(setting_id == best_setting_id)

  write.csv(best_setting, file.path(out_dir, "best_setting_summary.csv"), row.names = FALSE)
  write.csv(best_summary, file.path(out_dir, "best_setting_empirical_summary.csv"), row.names = FALSE)
  write.csv(best_fits, file.path(out_dir, "best_setting_fit_results.csv"), row.names = FALSE)

  message("Best setting ID: ", best_setting_id)
  print(best_setting)

  fit_m_norm <- fit_batch_number_law(best_summary, target = "norm")
  fit_s_norm <- fit_cells_per_batch_law(best_summary, target = "norm")
  fit_m_raw <- fit_batch_number_law(best_summary, target = "raw")
  fit_s_raw <- fit_cells_per_batch_law(best_summary, target = "raw")

  fit_row_m_norm <- best_fits %>% filter(law == "fixed_n_vary_m", target == "norm")
  fit_row_s_norm <- best_fits %>% filter(law == "fixed_m_vary_s", target == "norm")
  fit_row_m_raw <- best_fits %>% filter(law == "fixed_n_vary_m", target == "raw")
  fit_row_s_raw <- best_fits %>% filter(law == "fixed_m_vary_s", target == "raw")

  if (fit_m_norm$ok) {
    save_plot_pair("best_batch_number_scaling_norm", plot_bio_fit_m(fit_m_norm, fit_row_m_norm, "norm"))
  }
  if (fit_s_norm$ok) {
    save_plot_pair("best_batch_size_scaling_norm", plot_bio_fit_s(fit_s_norm, fit_row_s_norm, "norm"), width = 7.2, height = 3.2)
  }
  if (fit_m_raw$ok) {
    save_plot_pair("best_batch_number_scaling_raw", plot_bio_fit_m(fit_m_raw, fit_row_m_raw, "raw"))
  }
  if (fit_s_raw$ok) {
    save_plot_pair("best_batch_size_scaling_raw", plot_bio_fit_s(fit_s_raw, fit_row_s_raw, "raw"), width = 7.2, height = 3.2)
  }
} else {
  warning("No setting had successful normalized fits for both scaling laws. See ranked_setting_results.csv and all_fit_results.csv.")
  write.csv(tibble(message = "No successful best setting found."), file.path(out_dir, "best_setting_summary.csv"), row.names = FALSE)
}

message("\nCompleted. Results saved to: ", normalizePath(out_dir, mustWork = FALSE))
