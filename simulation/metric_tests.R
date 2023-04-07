library(tidyverse)
library(ggthemes)
library(scales)

library(kBET)
library(cluster)
library(lisi)

theme_set(theme_few())

set.seed(1)

stdev <- 0.5
n_cluster <- 200

cluster1 <- data.frame(
  X1 = rnorm(n_cluster, mean = 0, sd = stdev),
  X2 = rnorm(n_cluster, mean = 0, sd = stdev),
  true_cluster = "Cluster 1"
  ) %>%
  mutate(
    bad_mixing = if_else(X1 - X2 > 0, "Batch 1", "Batch 2"),
    good_mixing = if_else(row_number() %% 2 == 0, "Batch 1", "Batch 2")
  )

cluster1

cluster2 <- data.frame(
  X1 = rnorm(n_cluster, mean = 0, sd = stdev),
  X2 = rnorm(n_cluster, mean = 0, sd = stdev),
  true_cluster = "Cluster 2"
  ) %>%
  mutate(
    # bad_mixing = if_else(X1 - X2 < 0, "Batch 1", "Batch 2"),
    bad_mixing = if_else(X1^2 + X2^2 <= 0.5, "Batch 1", "Batch 2"),
    good_mixing = if_else(row_number() %% 2 == 0, "Batch 1", "Batch 2"),
    X1 = X1 + 3,
    X2 = X2 + 3
  )

cluster3 <- data.frame(
  X1 = rnorm(n_cluster, mean = 0, sd = stdev),
  X2 = rnorm(n_cluster, mean = 0, sd = stdev),
  true_cluster = "Cluster 3"
  ) %>%
  mutate(
    bad_mixing = if_else(abs(X2) + abs(X1) > 0.75, "Batch 1", "Batch 2"),
    good_mixing = if_else(row_number() %% 2 == 0, "Batch 1", "Batch 2"),
    X1 = X1 + 5,
    X2 = X2 - 2.5
  )

all_cluster <- bind_rows(cluster1, cluster2)

# Visualization of cluster labels
all_cluster %>%
  pivot_longer(cols = c(bad_mixing, good_mixing), names_to = "mixtype", values_to = "batch") %>%
  mutate(mixtype = if_else(mixtype == "bad_mixing", "Bad Mixing", "Good Mixing")) %>%
  ggplot(aes(x = X1, y = X2, color = batch, shape = true_cluster)) + 
  geom_point() + facet_wrap(mixtype ~ .) +
  labs(color = "Batch Number", shape = "Cell Cluster") +
  scale_color_tableau() +
  guides(shape = guide_legend(order = 2), col = guide_legend(order = 1)) +
  theme(legend.position = "bottom")

###########################################
# kBET
###########################################
sample_size <- floor(seq(0.05, 0.25, 0.05) * nrow(all_cluster))

# bad mixing
bad_mixing <- sample_size %>%
  lapply(FUN = function(x) {
    kBET(all_cluster[, 1:2], all_cluster$bad_mixing, k0 = x, plot = FALSE)
    }
  )

bad_mixing_df <- bad_mixing %>% lapply(FUN = function(x) {
  bind_cols(k0 = x[["params"]][["k0"]], as.data.frame(x[["stats"]]))
    }
  ) %>%
  list_rbind()

bad_mixing_plot <- bad_mixing_df %>%
  select(-kBET.signif) %>%
  pivot_longer(cols = c(kBET.expected, kBET.observed)) %>%
  mutate(name = if_else(name == "kBET.expected", "Expected", "Observed")) %>%
  ggplot(aes(x = factor(k0), y = value, color = name)) +
  geom_boxplot() +
  facet_wrap(name ~ ., nrow = 1) +
  ylim(0, 1) +
  labs(color = "") +
  scale_color_tableau() +
  guides(shape = guide_legend(order = 2), col = guide_legend(order = 1)) +
  theme(legend.position = "bottom")

bad_mixing_plot

# good mixing
good_mixing <- sample_size %>%
  lapply(FUN = function(x) {
    kBET(all_cluster[, 1:2], all_cluster$good_mixing, k0 = x, plot = FALSE)
    }
  )

good_mixing_df <- good_mixing %>% lapply(FUN = function(x) {
  bind_cols(k0 = x[["params"]][["k0"]], as.data.frame(x[["stats"]]))
    }
  ) %>%
  list_rbind()

good_mixing_plot <- good_mixing_df %>%
  select(-kBET.signif) %>%
  pivot_longer(cols = c(kBET.expected, kBET.observed)) %>%
  mutate(name = if_else(name == "kBET.expected", "Expected", "Observed")) %>%
  ggplot(aes(x = factor(k0), y = value, color = name)) +
  geom_boxplot() +
  facet_wrap(name ~ ., nrow = 1) +
  ylim(0, 1) +
  labs(color = "") +
  scale_color_tableau() +
  guides(shape = guide_legend(order = 2), col = guide_legend(order = 1)) +
  theme(legend.position = "bottom")

good_mixing_plot

###########################################
# ASW
###########################################

dist_matrix <- dist(all_cluster[, 1:2])

bad_mixing_sw <- silhouette(as.numeric(factor(all_cluster$bad_mixing)), dist_matrix)
good_mixing_sw <- silhouette(as.numeric(factor(all_cluster$good_mixing)), dist_matrix)
cluster_sw <- silhouette(as.numeric(factor(all_cluster$true_cluster)), dist_matrix)


all_cluster <- all_cluster %>%
  bind_cols(
    bad_sil_width = bad_mixing_sw[, "sil_width"],
    good_sil_width = good_mixing_sw[, "sil_width"],
    cluster_sil_width = cluster_sw[, "sil_width"]
  )

all_cluster %>%
  pivot_longer(cols = c(bad_mixing, good_mixing), names_to = "mixtype", values_to = "batch") %>%
  mutate(
    mixtype = if_else(mixtype == "bad_mixing", "Bad Mixing", "Good Mixing"),
    sil_width = if_else(mixtype == "Bad Mixing", bad_sil_width, good_sil_width)
  ) %>%
  ggplot(aes(x = X1, y = X2, color = sil_width, shape = batch)) +
  geom_point() +
  facet_wrap(mixtype ~ .) +
  scale_color_viridis_c(end = 0.9)
