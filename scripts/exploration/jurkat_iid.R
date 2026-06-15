library(Seurat)
library(Matrix)
library(irlba)
library(ggplot2)

# ============================================================
# Load data
# ============================================================

# replace with your path
data_dir <- "data/Jurkat/sample_filtered_feature_bc_matrix"

counts <- Read10X(data.dir = data_dir)

seu <- CreateSeuratObject(counts)

# ============================================================
# Basic preprocessing
# ============================================================

seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 1e4)

seu <- FindVariableFeatures(
  seu,
  selection.method = "vst",
  nfeatures = 2000
)

genes_use <- VariableFeatures(seu)

X <- GetAssayData(seu, slot = "data")[genes_use, ]

# ============================================================
# Convert to dense matrix
# ============================================================

X <- as.matrix(X)

# center genes
X <- X - rowMeans(X)

# dimensions
p <- nrow(X)
n <- ncol(X)

cat("p =", p, "\n")
cat("n =", n, "\n")
cat("p/n =", p/n, "\n")

# ============================================================
# PCA spectrum
# ============================================================

k <- min(100, min(dim(X)) - 1)

svd_res <- irlba(X, nv = k)

eigvals <- (svd_res$d^2) / n

df <- data.frame(
  PC = 1:length(eigvals),
  Eigenvalue = eigvals
)

# ============================================================
# Marchenko-Pastur upper edge
# ============================================================

# estimate noise variance from median tail eigenvalues
tail_start <- floor(length(eigvals) * 0.5)

sigma2_hat <- median(eigvals[tail_start:length(eigvals)]) /
  (1 + sqrt(p/n))^2

mp_upper <- sigma2_hat * (1 + sqrt(p/n))^2

cat("Estimated MP upper edge =", mp_upper, "\n")

# ============================================================
# Plot spectrum
# ============================================================

ggplot(df, aes(x = PC, y = Eigenvalue)) +
  geom_point(size = 2) +
  geom_line() +
  geom_hline(
    yintercept = mp_upper,
    color = "red",
    linetype = "dashed",
    linewidth = 1
  ) +
  theme_bw(base_size = 16) +
  ggtitle("Jurkat Eigenvalue Spectrum") +
  ylab("Covariance Eigenvalue")

# ============================================================
# Plot log spectrum
# ============================================================

ggplot(df, aes(x = PC, y = Eigenvalue)) +
  geom_point(size = 2) +
  geom_line() +
  geom_hline(
    yintercept = mp_upper,
    color = "red",
    linetype = "dashed",
    linewidth = 1
  ) +
  scale_y_log10() +
  theme_bw(base_size = 16) +
  ggtitle("Jurkat Eigenvalue Spectrum (log scale)") +
  ylab("Covariance Eigenvalue")

# ============================================================
# Optional: visualize PCs
# ============================================================

emb <- svd_res$v[, 1:2]

plot_df <- data.frame(
  PC1 = emb[,1],
  PC2 = emb[,2]
)

ggplot(plot_df, aes(PC1, PC2)) +
  geom_point(size = 0.3, alpha = 0.3) +
  theme_bw(base_size = 16) +
  ggtitle("Jurkat PC space")
