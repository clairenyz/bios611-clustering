suppressPackageStartupMessages({
  library(cluster)
  library(ggplot2)
  library(plotly)
  library(Matrix)
})

# ------------------------------
# Task 1 generator: hypercube clusters
# ------------------------------
# generate_hypercube_clusters(n, k, side_length, noise_sd = 1.0)
# n: dimensions and also number of clusters
# k: points per cluster
# side_length: L (distance of cluster centers along axes)
# noise_sd: Gaussian noise SD around each center

generate_hypercube_clusters <- function(n, k, side_length, noise_sd = 1.0) {
  # Centers are the positive corners: (L,0,...,0), (0,L,0,...,0), ..., (0,...,0,L)
  centers <- diag(side_length, nrow = n, ncol = n)
  # For each center, generate k points with N(center, noise_sd^2 I)
  X <- lapply(1:n, function(i) {
    matrix(rnorm(k * n, mean = 0, sd = noise_sd), nrow = k, ncol = n) +
      matrix(rep(centers[i, ], each = k), nrow = k)
  })
  X <- do.call(rbind, X)
  # True labels (1..n repeated k times)
  y <- rep(1:n, each = k)
  list(X = X, y = y)
}

# K-means wrapper with robust initialization for clusGap
kmeans_gap <- function(x, k) {
  set.seed(NULL)
  stats::kmeans(x, centers = k, nstart = 20, iter.max = 50)
}

# Compute a tidy dataframe from clusGap result
extract_gap_df <- function(gap_obj) {
  data.frame(
    k = gap_obj$Tab[, "k"],
    logW = gap_obj$Tab[, "logW"],
    ElogW = gap_obj$Tab[, "E.logW"],
    gap = gap_obj$Tab[, "gap"],
    se_sim = gap_obj$Tab[, "SE.sim"]
  )
}

# Helper: find modal K* estimate over repeated sims (using clusGap's firstSEmax rule)
estimate_kstar <- function(gap_obj) {
  cluster::maxSE(f = gap_obj$Tab[, "gap"], SE.f = gap_obj$Tab[, "SE.sim"], method = "firstSEmax")
}

# Save plot PNG utility
save_png <- function(plot, file, width = 1400, height = 900, res = 150) {
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  png(file, width = width, height = height, res = res)
  print(plot)
  dev.off()
}