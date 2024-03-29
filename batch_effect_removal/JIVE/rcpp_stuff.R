library(microbenchmark)
library(Rcpp)
library(RcppEigen)
library(RcppArmadillo)
library(tidyverse)
library(RSpectra)
library(r.jive)
library(cowplot)
source("jive_speedup.R")

# sourceCpp("matrix_multiplication.cpp")

set.seed(1)

########################
# jive() vs. jive_v2() #
########################

plot_data <- function(matrix) {
  matrix %>%
    as.data.frame() %>%
    # Data wrangling
    as_tibble() %>%
    rowid_to_column(var="X") %>%
    gather(key="Y", value="Z", -1) %>%
    # Change Y to numeric
    mutate(Y=as.numeric(gsub("V","",Y))) %>%
    # Viz
    ggplot(aes(X, Y, fill= Z)) + 
    geom_tile() +
    theme(legend.position="none") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    coord_flip() +
    theme_nothing()
}

# 50 x 100, two datasets
data(SimData)

tic("old")
jive_object <- jive(SimData, rankJ = 1, rankA = c(1, 1), method = "given")
toc()

summary(jive_object)

tic("new")
jive_v2_object <- jive_v2(SimData, rankJ = 1, rankA = c(1, 1), method = "given")
toc()

simulation_bench1 <- microbenchmark(
  jive(SimData, method = "given"),
  jive_v2(SimData, method = "given"),
  times = 100
  ) %>%
  as.data.frame() %>%
  mutate(time_ms = time * 1e-6)

# 200 x 1000, two datasets
N <- 200
k <- 1000

set.seed(1)

J <- t(sample(-2:2, k, replace = TRUE))
J <- J[rep(seq_len(nrow(J)), each = N / 2), ]
J <- rbind(matrix(0, N / 2, k), J)

A1 <- t(sort(rep(-2:2, k / 5)))
A1 <- A1[rep(seq_len(nrow(A1)), each = N), ]

A2 <- t(sample(rep(-2:2, k / 5)))
A2 <- A2[rep(seq_len(nrow(A2)), each = N), ]

E1 <- matrix(rnorm(N * k), N, k)
E2 <- matrix(rnorm(N * k), N, k)

D1 <- J + A1 + E1
D2 <- J + A2 + E2

orig_plot <- plot_grid(
  nrow = 2,
  plot_data(D1),
  plot_data(J),
  plot_data(A1),
  plot_data(E1),
  
  plot_data(D2),
  plot_data(J),
  plot_data(A2),
  plot_data(E2)
)

ggsave("output/simdata2_orig.png", orig_plot, width = 9, height = 4)

###

SimData2 <- list()
SimData2[["Data1"]] <- D1
SimData2[["Data2"]] <- D2

jive_object <- jive(SimData2, rankJ = 1, rankA = c(1, 1), method = "given")
print(jive_object)

jive_v2_object <- jive_v2(SimData2, rankJ = 1, rankA = c(1, 1), method = "given")
print(jive_v2_object)

est_plot <- plot_grid(
  nrow = 2,
  plot_data(jive_v2_object[["data"]][[1]]),
  plot_data(jive_v2_object[["joint"]][[1]]),
  plot_data(jive_v2_object[["individual"]][[1]]),
  plot_data(jive_v2_object[["data"]][[1]] - jive_v2_object[["joint"]][[1]] - jive_v2_object[["individual"]][[1]]),
  
  plot_data(jive_v2_object[["data"]][[2]]),
  plot_data(jive_v2_object[["joint"]][[2]]),
  plot_data(jive_v2_object[["individual"]][[2]]),
  plot_data(jive_v2_object[["data"]][[2]] - jive_v2_object[["joint"]][[2]] - jive_v2_object[["individual"]][[2]])
)

ggsave("output/simdata2_jive.png", est_plot, width = 9, height = 4)

###

simulation_bench2 <- microbenchmark(
  jive(SimData2, method = "given"),
  jive_v2(SimData2, method = "given"),
  times = 20
  ) %>%
  as.data.frame() %>%
  mutate(time_ms = time * 1e-6)

############
# Full SVD #
############

# Benchmark 1: N = k = 100
N <- 100
k <- 100

X <- matrix(rnorm(N*k), N, k)

test <- eigenBDCSVD(X, n_cores = 4)

full_svd_bench1 <- microbenchmark(
  svd(X),
  eigenBDCSVD(X, n_cores = 1),
  eigenBDCSVD(X, n_cores = 4),
  eigenBDCSVD(X, n_cores = 8),
  eigenBDCSVD(X, n_cores = 16),
  times = 100
  ) %>%
  as.data.frame() %>%
  mutate(time_ms = time * 1e-6)

# Benchmark 2: N = k = 1000 (square matrix)
N <- 1000
k <- 1000

X <- matrix(rnorm(N*k), N, k)

full_svd_bench2 <- microbenchmark(
  svd(X),
  eigenBDCSVD(X, n_cores = 1),
  eigenBDCSVD(X, n_cores = 4),
  eigenBDCSVD(X, n_cores = 8),
  eigenBDCSVD(X, n_cores = 16),
  times = 100
) %>%
  as.data.frame() %>%
  mutate(time_ms = time * 1e-6)

# Benchmark 3: N = 400, k = 2500 (fat matrix)
N <- 400
k <- 2500

X <- matrix(rnorm(N*k), N, k)

full_svd_bench3 <- microbenchmark(
  svd(X),
  eigenBDCSVD(X, n_cores = 1),
  eigenBDCSVD(X, n_cores = 4),
  eigenBDCSVD(X, n_cores = 8),
  eigenBDCSVD(X, n_cores = 16),
  times = 100
) %>%
  as.data.frame() %>%
  mutate(time_ms = time * 1e-6)

##################################################
# Truncated/Partial Singular Value Decomposition #
##################################################

# Benchmark 1: N = k = 100
N <- 100
k <- 100

X <- matrix(rnorm(N*k), N, k)

partial_svd_bench1 <- microbenchmark(
  svd(X, nu =  1, nv =  1),
  svd(X, nu =  5, nv =  5),
  svd(X, nu = 10, nv = 10),
  svds(X, k =  1),
  svds(X, k =  5),
  svds(X, k = 10),
  times = 100
  ) %>%
  as.data.frame() %>%
  mutate(time_ms = time * 1e-6)

# Benchmark 2: N = k = 1000
N <- 1000
k <- 1000

X <- matrix(rnorm(N*k), N, k)

partial_svd_bench2 <- microbenchmark(
  svd(X, nu =  1, nv =  1),
  svd(X, nu =  5, nv =  5),
  svd(X, nu = 10, nv = 10),
  svds(X, k =  1),
  svds(X, k =  5),
  svds(X, k = 10),
  times = 100
  ) %>%
  as.data.frame() %>%
  mutate(time_ms = time * 1e-6)

#########################
# Matrix Multiplication #
#########################

# Benchmark 1: N = k = 100
N <- 100
k <- 100

A <- matrix(rnorm(N*k), N, k)
B <- matrix(rnorm(N*k), k, N)

mm_bench1 <- microbenchmark(
  A%*%B, 
  armaMatMult(A, B),
  eigenMapMatMult2(A, B, n_cores = 1),
  eigenMapMatMult2(A, B, n_cores = 2),
  eigenMapMatMult2(A, B, n_cores = 4), 
  eigenMapMatMult2(A, B, n_cores = 8), 
  times = 100
  ) %>%
  as.data.frame() %>%
  mutate(time_ms = time * 1e-6)

# Benchmark 2: N = k = 1000
N <- 1000
k <- 1000

A <- matrix(rnorm(N*k), N, k)
B <- matrix(rnorm(N*k), k, N)

mm_bench2 <- microbenchmark(
  A%*%B, 
  armaMatMult(A, B),
  eigenMapMatMult2(A, B, n_cores = 1),
  eigenMapMatMult2(A, B, n_cores = 2),
  eigenMapMatMult2(A, B, n_cores = 4), 
  eigenMapMatMult2(A, B, n_cores = 8), 
  times = 100
  ) %>%
  as.data.frame() %>%
  mutate(time_ms = time * 1e-6)

################
# Save Results #
################

#####
theme_set(ggthemes::theme_few())

####

simulation_dat <- bind_rows(
  bind_cols(benchmark = "SimData", simulation_bench1),
  bind_cols(benchmark = "SimData2", simulation_bench2)
)

write_csv(simulation_dat, file = "output/simulation_benchmark.csv")

simulation_dat <- read_csv("output/simulation_benchmark.csv")

simulation_bench1_plot <- simulation_dat %>%
  filter(benchmark == "SimData") %>%
  ggplot(aes(x = expr, y = time_ms, color = expr)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "Method", y = "Time (milliseconds)",
       title = "Simulated Data Benchmark") +
  scale_colour_Publication() +
  theme_Publication()

simulation_dat %>% filter(benchmark == "SimData") %>% group_by(expr) %>% summarize(mean = mean(time_ms))

simulation_bench1_plot
ggsave("output/jive_v2_simdata_benchmark.png", simulation_bench1_plot, width = 9, height = 4)

simulation_bench2_plot <- simulation_dat %>%
  filter(benchmark == "SimData2") %>%
  ggplot(aes(x = expr, y = time_ms, color = expr)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "", y = "Time (milliseconds)",
       title = "Simulated Data", color = "Function") +
  scale_color_tableau()

simulation_dat %>% filter(benchmark == "SimData2") %>% group_by(expr) %>% summarize(mean = mean(time_ms))

simulation_bench2_plot
ggsave("output/jive_v2_simdata2_benchmark.png", simulation_bench2_plot, width = 9, height = 4)

###

full_svd_dat <- bind_rows(
  bind_cols(benchmark = "FullSVD100", full_svd_bench1),
  bind_cols(benchmark = "FullSVD1000", full_svd_bench2),
  bind_cols(benchmark = "FullSVDfat", full_svd_bench3)
)

write_csv(full_svd_dat, file = "output/full_svd_benchmark.csv")

full_svd_dat <- read_csv("output/full_svd_benchmark.csv") %>%
  mutate(
    expr = factor(
      expr,
      levels = c("svd(X)", "eigenBDCSVD(X, n_cores = 1)", "eigenBDCSVD(X, n_cores = 4)", "eigenBDCSVD(X, n_cores = 8)", "eigenBDCSVD(X, n_cores = 16)")
    )
  )

full_svd_dat %>%
  filter(benchmark == "FullSVD100") %>%
  ggplot(aes(x = expr, y = time_ms, color = expr)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "Method", y = "Time (milliseconds)",
       title = "Full SVD") +
  scale_colour_Publication() +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

full_svd_dat %>%
  filter(benchmark == "FullSVD1000") %>%
  ggplot(aes(x = expr, y = time_ms, color = expr)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "Method", y = "Time (milliseconds)",
       title = "Full SVD") +
  scale_colour_Publication() +
  theme_Publication()

full_svd_dat %>%
  filter(benchmark == "FullSVDfat") %>%
  ggplot(aes(x = expr, y = time_ms)) +
  geom_boxplot() +
  labs(x = "Method", y = "Time (milliseconds)",
       title = "Full SVD") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

###

partial_svd_dat <- bind_rows(
  bind_cols(benchmark = "PartialSVD100", partial_svd_bench1),
  bind_cols(benchmark = "PartialSVD1000", partial_svd_bench2)
)

write_csv(partial_svd_dat, file = "output/partial_svd_benchmark.csv")

partial_svd_dat <- read_csv("output/partial_svd_benchmark.csv") %>%
  mutate(
    expr = factor(
      expr,
      levels = c(
        "svd(X, nu = 1, nv = 1)", "svds(X, k = 1)",
        "svd(X, nu = 5, nv = 5)", "svds(X, k = 5)",
        "svd(X, nu = 10, nv = 10)", "svds(X, k = 10)"
      )
    ),
    library = case_when(
      expr %in% c("svd(X, nu = 1, nv = 1)", "svd(X, nu = 5, nv = 5)", "svd(X, nu = 10, nv = 10)") ~ "Base R",
      TRUE ~ "RSpectra"
    ),
    library = factor(
      library,
      levels = c("Base R", "RSpectra")
    )
  )

# ggplot(data = partial_svd_bench1, aes(x = expr, y = time_ms)) +
#   geom_boxplot() +
#   labs(x = "Method", y = "Time (milliseconds)",
#        title = "Partial SVD", subtitle = "(100x100 Matrix)") +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

partial_svd_dat %>% filter(benchmark == "PartialSVD1000") %>% group_by(expr) %>% summarize(mean = mean(time_ms))

partial_svd_plot <- partial_svd_dat %>%
  filter(benchmark == "PartialSVD1000") %>%
  ggplot(aes(x = expr, y = time_ms, color = library)) +
  geom_boxplot() +
  labs(x = "", y = "Time (milliseconds)",
       title = "Partial SVD", color = "") +
  scale_x_discrete(labels = c("svd 1", "svds 1", "svd 5", "svds 5", "svd 10", "svds 10")) +
  scale_color_tableau() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "bottom")

partial_svd_plot

ggsave("output/partial_svd_benchmark.png", partial_svd_plot, width = 9, height = 4)
ggsave("output/partial_svd_benchmark_thin.png", partial_svd_plot, width = 5, height = 4)

###

mm_dat <- bind_rows(
  bind_cols(benchmark = "MM100", mm_bench1),
  bind_cols(benchmark = "MM1000", mm_bench2)
)

write_csv(mm_dat, file = "output/mm_benchmark.csv")

mm_dat <- read_csv("output/mm_benchmark.csv") %>%
  mutate(
    expr = factor(
      expr,
      levels = c(
        "A %*% B", "armaMatMult(A, B)", "eigenMapMatMult2(A, B, n_cores = 1)",
        "eigenMapMatMult2(A, B, n_cores = 2)", "eigenMapMatMult2(A, B, n_cores = 4)",
        "eigenMapMatMult2(A, B, n_cores = 8)"
      )
    ),
    library = case_when(
      expr == "A %*% B" ~ "Base R",
      expr == "armaMatMult(A, B)" ~ "Armadillo",
      TRUE ~ "Eigen"
    ),
    library = factor(
      library,
      levels = c("Base R", "Armadillo", "Eigen")
    )
  )

mm_dat %>% filter(benchmark == "MM1000") %>% group_by(expr) %>% summarize(mean = mean(time_ms))

ggplot(data = mm_bench1, aes(x = expr, y = time_ms)) +
  geom_boxplot() +
  labs(x = "Method", y = "Time (milliseconds)",
       title = "Matrix Multiplication", subtitle = "(100x100 Matrix) * (100x100 Matrix)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


mm_plot <- mm_dat %>%
  filter(benchmark == "MM1000") %>%
  ggplot(aes(x = expr, y = time_ms, color = library)) +
  geom_boxplot() +
  labs(x = "", y = "Time (milliseconds)",
       title = "Matrix Multiplication", color = "") +
  scale_x_discrete(labels = c("Base R", "Armadillo", "Eigen 1", "Eigen 2", "Eigen 4", "Eigen 8")) +
  scale_color_tableau() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "bottom")

mm_plot

ggsave("output/mm_benchmark.png", mm_plot, width = 9, height = 4)
ggsave("output/mm_benchmark_thin.png", mm_plot, width = 5, height = 4)


### both?

ggsave("output/both_benchmark.png", plot_grid(partial_svd_plot, mm_plot, nrow = 1, labels = c("A", "B")), width = 9, height = 4)
