source("R/00_utils.R")
suppressPackageStartupMessages(library(Rfast)) # for fast distance if available
set.seed(611)

# ------------------------------
# Shell generator
# generate_shell_clusters(n_shells, k_per_shell, max_radius, noise_sd = 0.1)
# Radii are equally spaced in (r_min, max_radius]; r_min > 0 to avoid origin collapse.

generate_shell_clusters <- function(n_shells, k_per_shell, max_radius, noise_sd = 0.1, r_min = 1e-3) {
  if (max_radius <= 0) stop("max_radius must be positive")
  radii <- seq(r_min, max_radius, length.out = n_shells)
  make_point <- function(r) {
    # Sample direction uniformly on S^2
    v <- rnorm(3)
    v <- v / sqrt(sum(v^2))
    # Radial noise (thickness)
    rr <- abs(r + rnorm(1, sd = noise_sd))
    rr * v
  }
  pts <- lapply(radii, function(r) t(replicate(k_per_shell, make_point(r))))
  X <- do.call(rbind, pts)
  y <- rep(seq_along(radii), each = k_per_shell)
  list(X = X, y = y, radii = radii)
}

# ------------------------------
# Spectral clustering wrapper compatible with clusGap
# We compute adjacency A by thresholding pairwise Euclidean distance with d_threshold.
# Then compute normalized Laplacian L_sym and take the k eigenvectors for embedding.
# Finally run k-means on rows of the embedding and return the kmeans object.

spectral_gap_fun <- function(d_threshold = 1.0) {
  function(x, k) {
    # Pairwise distances
    D <- if (requireNamespace("Rfast", quietly = TRUE)) Rfast::Dist(x) else as.matrix(dist(x))
    A <- (D < d_threshold) * 1.0
    diag(A) <- 0
    # Degree and normalized Laplacian
    deg <- rowSums(A)
    # Guard against isolated points by adding tiny ridge to degrees
    deg[deg == 0] <- 1e-8
    Dm12 <- Diagonal(x = 1 / sqrt(deg))
    L <- Diagonal(n = nrow(A), x = deg) - A
    Lsym <- Dm12 %*% L %*% Dm12
    
    # Smallest k eigenvectors
    es <- eigen(Lsym, symmetric = TRUE)
    idx <- order(es$values, decreasing = FALSE)[seq_len(k)]
    U <- es$vectors[, idx, drop = FALSE]
    # Row-normalize U for stability
    U <- U / sqrt(rowSums(U^2) + 1e-12)
    
    # K-means on U
    stats::kmeans(U, centers = k, nstart = 20, iter.max = 50)
  }
}

# Preview shells with plotly
preview <- generate_shell_clusters(n_shells = 4, k_per_shell = 150, max_radius = 6, noise_sd = 0.1)
plt <- plot_ly(x = preview$X[,1], y = preview$X[,2], z = preview$X[,3],
               type = "scatter3d", mode = "markers",
               marker = list(size = 2), color = factor(preview$y))
htmlwidgets::saveWidget(plt, file = "output/figures/task2_shells_preview.html", selfcontained = FALSE)
cat("Saved: output/figures/task2_shells_preview.html\n")

# ------------------------------
# Simulation: shrink max_radius from 10 to 0; estimate K* with clusGap using spectral wrapper

n_shells <- 4
k_per_shell <- 100
noise_sd <- 0.1
radii_grid <- seq(10, 0, by = -1)
B <- 20

est_k <- sapply(radii_grid, function(Rmax) {
  if (Rmax <= 0) return(1L)
  dat <- generate_shell_clusters(n_shells, k_per_shell, Rmax, noise_sd)
  X <- dat$X
  gap <- cluster::clusGap(X, FUNcluster = spectral_gap_fun(d_threshold = 1.0), K.max = n_shells + 3, B = B)
  estimate_kstar(gap)
})

shell_df <- data.frame(max_radius = radii_grid, est_k = est_k)

p2 <- ggplot(shell_df, aes(x = max_radius, y = est_k)) +
  geom_line() + geom_point() +
  geom_hline(yintercept = n_shells, linetype = 2) +
  labs(title = "Spectral clustering (clusGap) on concentric shells",
       subtitle = "Adjacency threshold d_threshold = 1.0; K* vs. max_radius",
       x = "max_radius (decreasing compresses shells)",
       y = "Estimated clusters (K*)") +
  theme_bw(base_size = 13)

save_png(p2, "output/figures/task2_gap_vs_radius.png")

# Heuristic failure point: first radius where K* < n_shells
fail_idx <- which(est_k < n_shells)
first_fail <- if (length(fail_idx)) min(radii_grid[fail_idx]) else NA_real_
write.csv(data.frame(first_drop_max_radius = first_fail), "output/tables/task2_gap_thresholds.csv", row.names = FALSE)
cat("Saved: output/figures/task2_gap_vs_radius.png\n")
cat("Saved: output/tables/task2_gap_thresholds.csv\n")

cat(sprintf("First failure radius (K* < %d): %s\n", n_shells, as.character(first_fail)))