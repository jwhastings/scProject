library(tidyverse)
library(ggthemes)
library(scales)
library(cowplot)

library(Rtsne)
library(umap)
library(kBET)
library(cluster)
library(lisi)

theme_set(theme_few())

set.seed(1)

####################################################################################################
# Qualitative Evaluation Metrics
####################################################################################################

library(dslabs) # MNIST dataset

# mnist <- read_mnist(download = TRUE, destdir = getwd())
mnist <- read_mnist(getwd())

pixels <- mnist[["train"]][["images"]]
dim(pixels)
# remove columns with no information
pixels <- pixels[, colSums(pixels) > 0]
dim(pixels)

labels <- factor(mnist[["train"]][["labels"]])

mnist_df <- bind_cols(data.frame(digit = labels), as.data.frame(pixels))

# t-SNE
digits_tsne <- readRDS("digits_tsne.rds")

digits_tsne_df <- data.frame(
  X1 = digits_tsne[["Y"]][, 1],
  X2 = digits_tsne[["Y"]][, 2],
  digit = labels
)

digits_tsne_label <- digits_tsne_df %>%
  group_by(digit) %>%
  summarize(
    X1_median = median(X1),
    X2_median = median(X2)
  )

digits_tsne_plot <- ggplot() +
  geom_point(
    data = digits_tsne_df,
    aes(x = X1, y = X2, color = digit),
    alpha = 0.2,
    show.legend = FALSE
  ) +
  geom_label(
    data = digits_tsne_label,
    aes(x = X1_median, y = X2_median, label = digit)
  ) +
  labs(x = "t-SNE 1", y = "t-SNE 2")

ggsave("tsne_plot.png", digits_tsne_plot)

# UMAP
digits_umap <- readRDS("digits_umap.rds")

digits_umap_df <- data.frame(
  X1 = digits_umap[, 1],
  X2 = digits_umap[, 2],
  digit = labels
)

digits_umap_label <- digits_umap_df %>%
  group_by(digit) %>%
  summarize(
    X1_median = median(X1),
    X2_median = median(X2)
  )

digits_umap_plot <- ggplot() +
  geom_point(
    data = digits_umap_df,
    aes(x = X1, y = X2, color = digit),
    alpha = 0.2,
    show.legend = F
  ) +
  geom_label(
    data = digits_umap_label,
    aes(x = X1_median, y = X2_median, label = digit)
  ) +
  labs(x = "UMAP 1", y = "UMAP 2")

ggsave("umap_plot.png", digits_umap_plot)

####################################################################################################
# Quantitative Evaluation Metrics
####################################################################################################

# Generate data

set.seed(1)
stdev <- 0.75
n_cluster <- 500
bax_mix_shift <- 0.75

cluster1 <- data.frame(
  X1 = rnorm(n_cluster, mean = 0, sd = stdev),
  X2 = rnorm(n_cluster, mean = 0, sd = stdev),
  true_cluster = "Cluster 1"
  ) %>%
  mutate(
    good_mixing = if_else(row_number() %% 5 == 0, "Batch 1", "Batch 2"),
    bad_mixing = good_mixing
  ) %>%
  pivot_longer(
    cols = c(bad_mixing, good_mixing),
    names_to = "mixtype",
    values_to = "batch"
  ) %>%
  mutate(
    X1 = if_else(mixtype == "bad_mixing" & batch == "Batch 1", X1 + bax_mix_shift, X1),
    X2 = if_else(mixtype == "bad_mixing" & batch == "Batch 1", X2 + bax_mix_shift, X2),
    X1 = if_else(mixtype == "bad_mixing" & batch == "Batch 2", X1 - bax_mix_shift, X1),
    X2 = if_else(mixtype == "bad_mixing" & batch == "Batch 2", X2 - bax_mix_shift, X2)
  )

# cluster1 %>%
#   ggplot(aes(x = X1, y = X2, color = batch, shape = true_cluster)) +
#   geom_point() +
#   facet_wrap(mixtype ~ .)

cluster2 <- data.frame(
  X1 = rnorm(n_cluster, mean = 0, sd = stdev),
  X2 = rnorm(n_cluster, mean = 0, sd = stdev),
  true_cluster = "Cluster 2"
  ) %>%
  mutate(
    good_mixing = if_else(row_number() %% 2 == 0, "Batch 1", "Batch 2"),
    bad_mixing = good_mixing,
    X1 = X1 + 10
  ) %>%
  pivot_longer(
    cols = c(bad_mixing, good_mixing),
    names_to = "mixtype",
    values_to = "batch"
  ) %>%
  mutate(
    X1 = if_else(mixtype == "bad_mixing" & batch == "Batch 2", X1 +  bax_mix_shift, X1)
  )

all_cluster <- bind_rows(
  cluster1
  # , cluster2
)

# Visualization of cluster labels
cells_metrics_plot <- all_cluster %>%
  mutate(mixtype = if_else(mixtype == "bad_mixing", "Batch Effect", "No Batch Effect")) %>%
  ggplot(aes(x = X1, y = X2, color = batch)) + 
  geom_point() +
  facet_wrap(mixtype ~ .) +
  labs(color = "Batch Number") +
  scale_color_tableau() +
  guides(shape = guide_legend(order = 2), col = guide_legend(order = 1)) +
  theme(legend.position = "bottom")

ggsave("cells_metrics_plot.png", cells_metrics_plot)

####################################################################################################
# kBET
####################################################################################################

sample_size <- floor(seq(0.05, 0.25, 0.05) * nrow(all_cluster))

bad_mixing <- filter(all_cluster, mixtype == "bad_mixing")
good_mixing <- filter(all_cluster, mixtype == "good_mixing")

# bad mixing
bad_mixing_kBET <- sample_size %>%
  lapply(FUN = function(x) {
    kBET(bad_mixing[, 1:2], bad_mixing$batch, k0 = x, plot = FALSE)
    }
  )

bad_mixing_kBET_df <- bad_mixing_kBET %>% lapply(FUN = function(x) {
  bind_cols(k0 = x[["params"]][["k0"]], as.data.frame(x[["stats"]]))
    }
  ) %>%
  list_rbind()

bad_mixing_kBET_plot <- bad_mixing_kBET_df %>%
  select(-kBET.signif) %>%
  pivot_longer(cols = c(kBET.expected, kBET.observed)) %>%
  mutate(name = if_else(name == "kBET.expected", "Expected", "Observed")) %>%
  ggplot(aes(x = factor(k0), y = value, color = name)) +
  geom_boxplot() +
  facet_wrap(name ~ ., nrow = 1) +
  ylim(0, 1) +
  labs(x = "Neighbors", y = "Rejection Rate", color = "") +
  scale_color_tableau() +
  guides(shape = guide_legend(order = 2), col = guide_legend(order = 1)) +
  theme(legend.position = "bottom")

# good mixing
good_mixing_kBET <- sample_size %>%
  lapply(FUN = function(x) {
    kBET(good_mixing[, 1:2], good_mixing$batch, k0 = x, plot = FALSE)
    }
  )

good_mixing_kBET_df <- good_mixing_kBET %>% lapply(FUN = function(x) {
  bind_cols(k0 = x[["params"]][["k0"]], as.data.frame(x[["stats"]]))
    }
  ) %>%
  list_rbind()

good_mixing_kBET_plot <- good_mixing_kBET_df %>%
  select(-kBET.signif) %>%
  pivot_longer(cols = c(kBET.expected, kBET.observed)) %>%
  mutate(name = if_else(name == "kBET.expected", "Expected", "Observed")) %>%
  ggplot(aes(x = factor(k0), y = value, color = name)) +
  geom_boxplot() +
  facet_wrap(name ~ ., nrow = 1) +
  ylim(0, 1) +
  labs(x = "Neighbors", y = "Rejection Rate", color = "") +
  scale_color_tableau() +
  guides(shape = guide_legend(order = 2), col = guide_legend(order = 1)) +
  theme(legend.position = "bottom")

# bad vs. good mixing
kbet_metric_plot <- bind_rows(
  select(bad_mixing_kBET_df, k0, kBET.observed) %>% mutate(mixtype = "Batch Effect"),
  select(good_mixing_kBET_df, k0, kBET.observed) %>% mutate(mixtype = "No Batch Effect")
  ) %>%
  ggplot(aes(x = factor(k0), y = kBET.observed, color = mixtype)) +
  geom_boxplot() +
  facet_wrap(mixtype ~ ., nrow = 1) +
  ylim(0, 1) +
  labs(x = "Number of Neighbors Considered", y = "Rejection Rate", color = "") +
  scale_color_brewer(palette = "Set1") +
  guides(shape = guide_legend(order = 2), col = guide_legend(order = 1)) +
  theme(legend.position = "bottom")

ggsave("kbet_metric_plot.png", kbet_metric_plot)

####################################################################################################
# ASW
####################################################################################################

bad_mixing_sw <- silhouette(as.numeric(factor(bad_mixing$batch)), dist(bad_mixing[, 1:2])) %>%
  as.data.frame() %>%
  select(batch_sil_width = sil_width)
good_mixing_sw <- silhouette(as.numeric(factor(good_mixing$batch)), dist(good_mixing[, 1:2])) %>%
  as.data.frame() %>%
  select(batch_sil_width = sil_width)

# bad_cluster_sw <- silhouette(as.numeric(factor(bad_mixing$true_cluster)), dist(bad_mixing[, 1:2])) %>%
#   as.data.frame() %>%
#   select(cluster_sil_width = sil_width)
# good_cluster_sw <- silhouette(as.numeric(factor(good_mixing$true_cluster)), dist(good_mixing[, 1:2])) %>%
#   as.data.frame() %>%
#   select(cluster_sil_width = sil_width)

bad_mixing_sw_df <- bind_cols(
  bad_mixing
  , bad_mixing_sw
  # , bad_cluster_sw
)
good_mixing_sw_df <- bind_cols(
  good_mixing
  , good_mixing_sw
  # , good_cluster_sw
)

###

all_mixing_sw_df <- bind_rows(
  bad_mixing_sw_df,
  good_mixing_sw_df
  ) %>%
  mutate(mixtype = if_else(mixtype == "bad_mixing", "Batch Effect", "No Batch Effect"))
  

all_mixing_sw_df_stat <- all_mixing_sw_df %>%
  group_by(mixtype) %>%
  summarize(
    mean = mean(batch_sil_width),
    median = median(batch_sil_width)
  )

asw_plot <- all_mixing_sw_df %>%
  ggplot(aes(x = X1, y = X2, color = batch_sil_width, shape = batch)) +
  geom_point() +
  facet_wrap(mixtype ~ .) +
  scale_color_viridis_c(end = 0.90) +
  labs(color = "Silhouette Width", shape = "Batch Number")

asw_density <- ggplot() +
  geom_density(data = all_mixing_sw_df, aes(x = batch_sil_width, fill = mixtype), alpha = 0.5) +
  geom_vline(data = all_mixing_sw_df_stat, aes(xintercept = mean, color = mixtype), linewidth = 1, linetype = "dashed") +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Silhouette Width", y = "Density", color = "", fill = "")

asw_plot_both <- plot_grid(
  asw_plot,
  asw_density,
  ncol = 1,
  align = "hv",
  axis = "tlbr",
  rel_heights = c(1, 0.5)
)

ggsave("asw_plot_both.png", asw_plot_both)

####################################################################################################
# LISI
####################################################################################################

perp <- 40

bad_lisi <- compute_lisi(
  bad_mixing[, 1:2],
  data.frame(
    LISI_Batch = bad_mixing$batch
  ),
  c("LISI_Batch"),
  perplexity = perp
)

good_lisi <- compute_lisi(
  good_mixing[, 1:2],
  data.frame(
    LISI_Batch = bad_mixing$batch
  ),
  c("LISI_Batch"),
  perplexity = perp
)

bad_lisi_df <- bind_cols(bad_mixing, bad_lisi)
good_lisi_df <- bind_cols(good_mixing, good_lisi)

all_lisi_df <- bind_rows(
  bad_lisi_df,
  good_lisi_df
  ) %>%
  mutate(mixtype = if_else(mixtype == "bad_mixing", "Batch Effect", "No Batch Effect"))

all_lisi_df_stat <- all_lisi_df %>%
  group_by(mixtype) %>%
  summarize(
    mean = mean(LISI_Batch),
    median = median(LISI_Batch)
  )

lisi_plot <- all_lisi_df %>%
  ggplot(aes(x = X1, y = X2, color = LISI_Batch, shape = batch)) +
  geom_point() +
  facet_wrap(mixtype ~ .) +
  scale_color_viridis_c(end = 0.9) +
  labs(color = "LISI", shape = "Batch Number")

lisi_density <- ggplot() +
  geom_density(data = all_lisi_df, aes(x = LISI_Batch, fill = mixtype), alpha = 0.5) +
  geom_vline(data = all_lisi_df_stat, aes(xintercept = median, color = mixtype), linewidth = 1, linetype = "dashed") +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Local Inverse Simpson's Index", y = "Density", color = "", fill = "")

lisi_plot_both <- plot_grid(
  lisi_plot,
  lisi_density,
  ncol = 1,
  align = "hv",
  axis = "tlbr",
  rel_heights = c(1, 0.5)
)

ggsave("lisi_plot_both.png", lisi_plot_both)

### SI

si <- tribble(
  ~Type, ~n,
  "A", 80,
  "B", 125,
  "C", 95
) %>% 
  mutate(pi = n / sum(n))

1 / sum(si$pi^2)
