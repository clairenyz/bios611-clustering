source("R/00_utils.R")
set.seed(611)

# Simulation parameters
ns <- c(6, 5, 4, 3, 2)                 # dimensions & number of clusters
side_lengths_by_n <- lapply(ns, function(n) seq(10, 1, by = -1))
k_per_cluster <- 100
noise_sd <- 1.0
B <- 20                                # clusGap Monte Carlo replicates (adjust for speed/accuracy)

results <- list()
thresholds <- data.frame()

for (idx in seq_along(ns)) {
  n <- ns[idx]
  side_grid <- side_lengths_by_n[[idx]]
  
  est_k_by_L <- sapply(side_grid, function(L) {
    dat <- generate_hypercube_clusters(n = n, k = k_per_cluster, side_length = L, noise_sd = noise_sd)
    X <- dat$X
    
    # clusGap over k = 1..(n+2) to be safe
    gap <- cluster::clusGap(X, FUNcluster = kmeans_gap, K.max = n + 2, B = B, verbose = FALSE)
    est_k <- estimate_kstar(gap)
    est_k
  })
  
  results[[as.character(n)]] <- data.frame(n = n, side_length = side_grid, est_k = est_k_by_L)
  
  # Find earliest side_length where estimate drops below true n consistently
  # Consistency heuristic: first L where est_k < n and remains <= n for the rest of the grid
  drop_idx <- which(est_k_by_L < n)
  drop_L <- if (length(drop_idx)) min(side_grid[drop_idx]) else NA_real_
  
  thresholds <- rbind(thresholds, data.frame(n = n, first_drop_side_length = drop_L))
}

plot_df <- do.call(rbind, results)

p <- ggplot(plot_df, aes(x = side_length, y = est_k, group = factor(n))) +
  geom_line() + geom_point() +
  geom_hline(aes(yintercept = n, color = factor(n)), linetype = 2) +
  scale_x_reverse(breaks = unique(plot_df$side_length)) +
  labs(title = "Estimated number of clusters vs side_length by dimension",
       subtitle = "Hypercube centers with Gaussian noise; clusGap + kmeans (nstart=20, iter.max=50)",
       x = "side_length (L) — decreasing means centers move closer",
       y = "Estimated clusters (K*)",
       color = "Dimension n") +
  theme_bw(base_size = 13)

save_png(p, "output/figures/task1_gap_by_dim.png")

write.csv(thresholds, "output/tables/task1_gap_thresholds.csv", row.names = FALSE)
cat("Saved: output/figures/task1_gap_by_dim.png\n")
cat("Saved: output/tables/task1_gap_thresholds.csv\n")