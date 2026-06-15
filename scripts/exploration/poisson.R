compute_empirical_mi <- function(X1, X2, r = 5) {
  X1 <- log1p(X1)
  X2 <- log1p(X2)
  
  X1 <- sweep(X1, 1, rowMeans(X1), "-")
  X2 <- sweep(X2, 1, rowMeans(X2), "-")
  
  V1 <- svd(X1, nu = 0, nv = r)$v
  V2 <- svd(X2, nu = 0, nv = r)$v
  
  gamma <- svd(t(V1) %*% V2, nu = 0, nv = 0)$d^2
  gamma <- pmin(pmax(gamma, 0), 1 - 1e-12)
  
  -0.5 * sum(log(1 - gamma))
}


set.seed(1)

p <- 1000
n <- 2000
total_umi <- 1000
shape <- 1
K <- 5
r <- K

# ----------------------------
# 1. Cell clusters
# ----------------------------
cell_cluster <- rep(seq_len(K), length.out = n)
cell_cluster <- sample(cell_cluster)

# ----------------------------
# 2. Ground-truth mean matrix
# genes × clusters
# ----------------------------
GeneMean_K <- matrix(0, nrow = p, ncol = K)

for (k in seq_len(K)) {
  gm <- rlnorm(p, meanlog = -1, sdlog = 1)
  gm <- gm / sum(gm) * total_umi
  GeneMean_K[, k] <- gm
}

# Ground-truth genes × cells mean matrix
Lambda_base <- GeneMean_K[, cell_cluster]

# ----------------------------
# 3. Two independent Gamma biological replicates
# ----------------------------
Lambda1 <- matrix(
  rgamma(p * n, shape = shape, scale = as.vector(Lambda_base) / shape),
  nrow = p,
  ncol = n
)

Lambda2 <- matrix(
  rgamma(p * n, shape = shape, scale = as.vector(Lambda_base) / shape),
  nrow = p,
  ncol = n
)

# ----------------------------
# 4. Split high / low expression genes
# ----------------------------
gene_mean <- rowMeans(Lambda_base)

ord <- order(gene_mean, decreasing = TRUE)
gene_mean_sorted <- gene_mean[ord]

cum_counts <- cumsum(gene_mean_sorted)
total_counts <- sum(gene_mean_sorted)

k_half <- which(cum_counts >= total_counts / 2)[1]

high_idx <- ord[1:k_half]
low_idx  <- ord[(k_half + 1):length(ord)]

# ----------------------------
# 5. MI comparisons
# ----------------------------
library(ggplot2)

sampling_rates <- c(0.1, 0.25, 0.5, 1, 2, 5, 10, 100)

res_list <- list()
pca_list <- list()

for (s in sampling_rates) {
  print(paste0("sampling rate: ", s))
  X1_s <- matrix(
    rpois(p * n, lambda = as.vector(s * Lambda1)),
    nrow = p,
    ncol = n
  )
  
  X2_s <- matrix(
    rpois(p * n, lambda = as.vector(s * Lambda2)),
    nrow = p,
    ncol = n
  )
  
  avg_umi_full <- mean(colSums(X1_s))
  avg_umi_high <- mean(colSums(X1_s[high_idx, , drop = FALSE]))
  avg_umi_low  <- mean(colSums(X1_s[low_idx, , drop = FALSE]))
  
  res_s <- data.frame(
    sampling_rate = s,
    avg_umi_per_cell = c(avg_umi_full, avg_umi_high, avg_umi_low),
    group = c("full", "high", "low"),
    
    mi_gamma_gamma = c(
      compute_empirical_mi(Lambda1, Lambda2, r = r),
      compute_empirical_mi(Lambda1[high_idx, ], Lambda2[high_idx, ], r = r),
      compute_empirical_mi(Lambda1[low_idx, ], Lambda2[low_idx, ], r = r)
    ),
    
    mi_poisson_poisson = c(
      compute_empirical_mi(X1_s, X2_s, r = r),
      compute_empirical_mi(X1_s[high_idx, ], X2_s[high_idx, ], r = r),
      compute_empirical_mi(X1_s[low_idx, ], X2_s[low_idx, ], r = r)
    ),
    
    n_genes = c(
      p,
      length(high_idx),
      length(low_idx)
    )
  )
  
  res_s$mi_loss_due_to_poisson <- 
    res_s$mi_gamma_gamma - res_s$mi_poisson_poisson
  
  res_s$retention <- 
    res_s$mi_poisson_poisson / res_s$mi_gamma_gamma
  
  res_list[[as.character(s)]] <- res_s
  
  # ----------------------------
  # PCA visualization data
  # ----------------------------
  X_log <- log1p(X1_s)
  X_log <- sweep(X_log, 1, rowMeans(X_log), "-")
  
  pca <- prcomp(t(X_log), center = FALSE, scale. = FALSE)
  
  pca_df <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    cluster = factor(cell_cluster),
    sampling_rate = s,
    avg_umi_per_cell = avg_umi_full
  )
  
  pca_list[[as.character(s)]] <- pca_df
}

res_depth <- do.call(rbind, res_list)
rownames(res_depth) <- NULL

pca_df_all <- do.call(rbind, pca_list)
rownames(pca_df_all) <- NULL


res_full <- subset(res_depth, group == "full")

bound_val <- unique(res_full$mi_gamma_gamma)

ggplot(res_full, aes(
  x = avg_umi_per_cell,
  y = mi_poisson_poisson
)) +
  geom_line(linewidth = 1.2, color = "black") +
  geom_point(size = 2.5, color = "black") +
  
  # theoretical upper bound
  geom_hline(
    yintercept = bound_val,
    linetype = "dashed",
    linewidth = 1
  ) +
  
  scale_x_log10() +
  theme_classic(base_size = 14) +
  labs(
    x = "Average UMI per cell",
    y = "Poisson–Poisson empirical MI",
    title = "Information retained vs sequencing depth",
    subtitle = "Dashed line = Gamma–Gamma MI (upper bound)"
  )


ggplot(pca_df_all, aes(PC1, PC2, color = cluster)) +
  geom_point(size = 0.6, alpha = 0.6) +
  facet_wrap(~ avg_umi_per_cell, scales = "free") +
  theme_classic(base_size = 14) +
  labs(
    title = "PCA visualization across sequencing depths",
    subtitle = "Each panel uses Poisson-sampled counts at a different UMI depth",
    color = "Cluster"
  )


res_hl <- subset(res_depth, group %in% c("high", "low"))

ggplot(res_hl, aes(
  x = factor(sampling_rate),
  y = mi_poisson_poisson,
  fill = group
)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  
  theme_classic(base_size = 14) +
  
  labs(
    x = "Sampling rate (sequencing depth multiplier)",
    y = "Poisson–Poisson empirical MI",
    title = "High vs Low expression: information retention",
    fill = "Gene group"
  )









set.seed(1)

# ----------------------------
# Fixed parameters
# ----------------------------
p <- 1000
total_umi <- 500
shape <- 1
K <- 5
r <- K
sampling_rate <- 0.5

n_list <- c(100, 200, 500, 1000, 2000, 5000, 50000)

res_list <- list()
pca_list <- list()

for (n in n_list) {
  print(paste0("n cells: ", n))
  # ----------------------------
  # 1. Cell clusters
  # ----------------------------
  cell_cluster <- rep(seq_len(K), length.out = n)
  cell_cluster <- sample(cell_cluster)
  
  # ----------------------------
  # 2. Ground-truth mean (fixed distribution)
  # ----------------------------
  GeneMean_K <- matrix(0, nrow = p, ncol = K)
  
  for (k in seq_len(K)) {
    gm <- rlnorm(p, meanlog = -1, sdlog = 0.4)
    gm <- gm / sum(gm) * total_umi
    GeneMean_K[, k] <- gm
  }
  
  Lambda_base <- GeneMean_K[, cell_cluster]
  
  # ----------------------------
  # 3. Gamma biological noise
  # ----------------------------
  Lambda1 <- matrix(
    rgamma(p * n, shape = shape, scale = as.vector(Lambda_base) / shape),
    nrow = p
  )
  
  Lambda2 <- matrix(
    rgamma(p * n, shape = shape, scale = as.vector(Lambda_base) / shape),
    nrow = p
  )
  
  # ----------------------------
  # 4. Poisson technical noise
  # ----------------------------
  X1 <- matrix(
    rpois(p * n, lambda = as.vector(sampling_rate * Lambda1)),
    nrow = p
  )
  
  X2 <- matrix(
    rpois(p * n, lambda = as.vector(sampling_rate * Lambda2)),
    nrow = p
  )
  
  # ----------------------------
  # 5. MI
  # ----------------------------
  mi_gamma <- compute_empirical_mi(Lambda1, Lambda2, r = r)
  mi_pois  <- compute_empirical_mi(X1, X2, r = r)
  
  res_list[[as.character(n)]] <- data.frame(
    n_cells = n,
    mi_gamma_gamma = mi_gamma,
    mi_poisson_poisson = mi_pois,
    retention = mi_pois / mi_gamma
  )
  
  # ----------------------------
  # 6. PCA (use X1)
  # ----------------------------
  X_log <- log1p(X1)
  X_log <- sweep(X_log, 1, rowMeans(X_log), "-")
  
  pca <- prcomp(t(X_log), center = FALSE, scale. = FALSE)
  
  pca_df <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    cluster = factor(cell_cluster),
    n_cells = n
  )
  
  pca_list[[as.character(n)]] <- pca_df
}

res_n <- do.call(rbind, res_list)
pca_n <- do.call(rbind, pca_list)


ggplot(res_n, aes(x = n_cells, y = mi_poisson_poisson)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  
  geom_hline(
    yintercept = unique(res_n$mi_gamma_gamma),
    linetype = "dashed"
  ) +
  
  scale_x_log10() +
  theme_classic(base_size = 14) +
  labs(
    x = "Number of cells",
    y = "Poisson–Poisson empirical MI",
    title = "Effect of sample size on MI estimation",
    subtitle = "Dashed = Gamma–Gamma upper bound"
  )

ggplot(pca_n, aes(PC1, PC2, color = cluster)) +
  geom_point(size = 0.6, alpha = 0.6) +
  facet_wrap(~ n_cells, scales = "free") +
  theme_classic(base_size = 14) +
  labs(
    title = "PCA improves with more cells",
    subtitle = "Same noise level, better estimation with larger n",
    color = "Cluster"
  )





set.seed(1)

library(ggplot2)

# ----------------------------
# Fixed parameters
# ----------------------------
p <- 500
n <- 1000
total_umi <- 1000
K <- 5
r <- K
sampling_rate <- 1

shape_list <- c(0.25, 0.5, 1, 2, 5, 10, 50, 100)

# ----------------------------
# Fixed ground truth mean matrix
# ----------------------------
cell_cluster <- rep(seq_len(K), length.out = n)
cell_cluster <- sample(cell_cluster)

GeneMean_K <- matrix(0, nrow = p, ncol = K)

for (k in seq_len(K)) {
  gm <- rlnorm(p, meanlog = -1, sdlog = 0.4)
  gm <- gm / sum(gm) * total_umi
  GeneMean_K[, k] <- gm
}

Lambda_base <- GeneMean_K[, cell_cluster]

# ----------------------------
# Run simulation across Gamma variance
# ----------------------------
res_list <- list()
pca_list <- list()

for (shape in shape_list) {
  
  Lambda1 <- matrix(
    rgamma(p * n, shape = shape, scale = as.vector(Lambda_base) / shape),
    nrow = p,
    ncol = n
  )
  
  Lambda2 <- matrix(
    rgamma(p * n, shape = shape, scale = as.vector(Lambda_base) / shape),
    nrow = p,
    ncol = n
  )
  
  X1 <- matrix(
    rpois(p * n, lambda = as.vector(sampling_rate * Lambda1)),
    nrow = p,
    ncol = n
  )
  
  X2 <- matrix(
    rpois(p * n, lambda = as.vector(sampling_rate * Lambda2)),
    nrow = p,
    ncol = n
  )
  
  mi_gamma <- compute_empirical_mi(Lambda1, Lambda2, r = r)
  mi_pois  <- compute_empirical_mi(X1, X2, r = r)
  
  res_list[[as.character(shape)]] <- data.frame(
    shape = shape,
    gamma_cv = 1 / sqrt(shape),
    gamma_var_factor = 1 / shape,
    avg_umi_per_cell = mean(colSums(X1)),
    mi_gamma_gamma = mi_gamma,
    mi_poisson_poisson = mi_pois,
    retention = mi_pois / mi_gamma
  )
  
  # PCA visualization
  X_log <- log1p(X1)
  X_log <- sweep(X_log, 1, rowMeans(X_log), "-")
  
  pca <- prcomp(t(X_log), center = FALSE, scale. = FALSE)
  
  pca_list[[as.character(shape)]] <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    cluster = factor(cell_cluster),
    shape = shape,
    gamma_cv = 1 / sqrt(shape)
  )
}

res_gamma_var <- do.call(rbind, res_list)
pca_gamma_var <- do.call(rbind, pca_list)

ggplot(res_gamma_var, aes(x = gamma_cv)) +
  geom_line(aes(y = mi_gamma_gamma), linewidth = 1.2, linetype = "dashed") +
  geom_point(aes(y = mi_gamma_gamma), size = 2.5) +
  geom_line(aes(y = mi_poisson_poisson), linewidth = 1.2) +
  geom_point(aes(y = mi_poisson_poisson), size = 2.5) +
  scale_x_log10() +
  theme_classic(base_size = 14) +
  labs(
    x = "Gamma noise CV = 1 / sqrt(shape)",
    y = "Empirical MI",
    title = "Gamma variability sets the information ceiling",
    subtitle = "Dashed = Gamma–Gamma MI, solid = Poisson–Poisson MI"
  )

ggplot(pca_gamma_var, aes(PC1, PC2, color = cluster)) +
  geom_point(size = 0.6, alpha = 0.6) +
  facet_wrap(~ gamma_cv, scales = "free") +
  theme_classic(base_size = 14) +
  labs(
    title = "PCA across Gamma biological variability",
    subtitle = "Higher Gamma CV means stronger cell-to-cell biological noise",
    color = "Cluster"
  )











set.seed(1)
library(ggplot2)

# ----------------------------
# Parameters
# ----------------------------
p <- 500
n <- 1000
K <- 5
r <- K

umi_list <- c(200, 500, 1000, 2000, 5000, 50000)
shape <- 0.25   # strong biological noise (CV = 2)

# ----------------------------
# Cell clusters
# ----------------------------
cell_cluster <- rep(seq_len(K), length.out = n)
cell_cluster <- sample(cell_cluster)

# ----------------------------
# Base gene mean
# ----------------------------
base_mean <- rlnorm(p, meanlog = -1, sdlog = 0.4)

# ----------------------------
# Fix one ground-truth mean (use any reference UMI for normalization)
# ----------------------------
umi_ref <- 1000
base_scaled <- base_mean / sum(base_mean) * umi_ref

# cluster means (fixed)
GeneMean_K <- matrix(0, nrow = p, ncol = K)
for (k in seq_len(K)) {
  gm <- base_scaled + rnorm(p, sd = 0.3 * base_scaled)
  gm <- pmax(gm, 1e-6)
  gm <- gm / sum(gm) * umi_ref
  GeneMean_K[, k] <- gm
}

Lambda_base_ref <- GeneMean_K[, cell_cluster]

# ----------------------------
# Sample Gamma ONCE (fixed biology)
# ----------------------------
Lambda1_ref <- matrix(
  rgamma(p * n, shape = shape,
         scale = as.vector(Lambda_base_ref) / shape),
  nrow = p
)

Lambda2_ref <- matrix(
  rgamma(p * n, shape = shape,
         scale = as.vector(Lambda_base_ref) / shape),
  nrow = p
)

# Gamma–Gamma bound (fixed)
mi_gamma <- compute_empirical_mi(Lambda1_ref, Lambda2_ref, r = r)

# ----------------------------
# Sweep UMI via scaling
# ----------------------------
res_list <- list()
pca_list <- list()

for (umi in umi_list) {
  
  scale_factor <- umi / umi_ref
  
  Lambda1 <- Lambda1_ref * scale_factor
  Lambda2 <- Lambda2_ref * scale_factor
  
  # Poisson sampling
  X1 <- matrix(
    rpois(p * n, lambda = as.vector(Lambda1)),
    nrow = p
  )
  
  X2 <- matrix(
    rpois(p * n, lambda = as.vector(Lambda2)),
    nrow = p
  )
  
  mi_pois <- compute_empirical_mi(X1, X2, r = r)
  
  res_list[[as.character(umi)]] <- data.frame(
    umi = umi,
    mi_gamma_gamma = mi_gamma,
    mi_poisson_poisson = mi_pois,
    retention = mi_pois / mi_gamma
  )
  
  # PCA
  X_log <- log1p(X1)
  X_log <- sweep(X_log, 1, rowMeans(X_log), "-")
  
  pca <- prcomp(t(X_log), center = FALSE, scale. = FALSE)
  
  pca_list[[as.character(umi)]] <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    cluster = factor(cell_cluster),
    umi = umi
  )
}

res_umi <- do.call(rbind, res_list)
pca_umi <- do.call(rbind, pca_list)

ggplot(res_umi, aes(x = umi, y = mi_poisson_poisson)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = unique(res_umi$mi_gamma_gamma),
             linetype = "dashed") +
  scale_x_log10() +
  theme_classic(base_size = 14) +
  labs(
    x = "Total UMI per cell",
    y = "Empirical MI",
    title = "MI approaches Gamma–Gamma bound as sequencing depth increases",
    subtitle = "Dashed = irreducible biological limit"
  )

ggplot(pca_umi, aes(PC1, PC2, color = cluster)) +
  geom_point(size = 0.6, alpha = 0.6) +
  facet_wrap(~ umi, scales = "free") +
  theme_classic(base_size = 14) +
  labs(
    title = "PCA with fixed biology and increasing sequencing depth",
    subtitle = "Even extreme depth (50k UMI) cannot beat biological variability",
    color = "Cluster"
  )
