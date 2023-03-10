library(microbenchmark)
library(Rcpp)
library(RcppEigen)
library(RcppArmadillo)
library(tidyverse)
library(RSpectra)

sourceCpp("matrix_multiplication.cpp")

set.seed(1)

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

ggplot(data = full_svd_bench1, aes(x = expr, y = time_ms)) +
  geom_boxplot() +
  labs(x = "Method", y = "Time (milliseconds)",
       title = "Full SVD", subtitle = "(100x100 Matrix)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

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

ggplot(data = full_svd_bench2, aes(x = expr, y = time_ms)) +
  geom_boxplot() +
  labs(x = "Method", y = "Time (milliseconds)",
       title = "Full SVD", subtitle = paste0("(", N, "x", k, " Matrix)")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Benchmark 3: N = 400, k = 1000 (fat matrix)
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

ggplot(data = full_svd_bench3, aes(x = expr, y = time_ms)) +
  geom_boxplot() +
  labs(x = "Method", y = "Time (milliseconds)",
       title = "Full SVD", subtitle = paste0("(", N, "x", k, " Matrix)")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Benchmark 4: N = 2500, k = 400 (thin matrix)
N <- 2500
k <- 400

X <- matrix(rnorm(N*k), N, k)

full_svd_bench4 <- microbenchmark(
  svd(X),
  eigenBDCSVD(X, n_cores = 1),
  eigenBDCSVD(X, n_cores = 4),
  eigenBDCSVD(X, n_cores = 8),
  eigenBDCSVD(X, n_cores = 16),
  eigenBDCSVD_sv(X, n_cores = 1),
  eigenBDCSVD_sv(X, n_cores = 4),
  eigenBDCSVD_sv(X, n_cores = 8),
  eigenBDCSVD_sv(X, n_cores = 16),
  times = 100
) %>%
  as.data.frame() %>%
  mutate(time_ms = time * 1e-6)

ggplot(data = full_svd_bench4, aes(x = expr, y = time_ms)) +
  geom_boxplot() +
  labs(x = "Method", y = "Time (milliseconds)",
       title = "Full SVD", subtitle = paste0("(", N, "x", k, " Matrix)")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

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

ggplot(data = partial_svd_bench1, aes(x = expr, y = time_ms)) +
  geom_boxplot() +
  labs(x = "Method", y = "Time (milliseconds)",
       title = "Partial SVD", subtitle = "(100x100 Matrix)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

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

ggplot(data = partial_svd_bench2, aes(x = expr, y = time_ms)) +
  geom_boxplot() +
  labs(x = "Method", y = "Time (milliseconds)",
       title = "Partial SVD", subtitle = "(1000x1000 Matrix)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

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
  eigenMatMult(A, B, n_cores = 1),
  eigenMatMult(A, B, n_cores = 2),
  eigenMatMult(A, B, n_cores = 4),
  eigenMatMult(A, B, n_cores = 8),
  eigenMapMatMult2(A, B, n_cores = 1),
  eigenMapMatMult2(A, B, n_cores = 2),
  eigenMapMatMult2(A, B, n_cores = 4), 
  eigenMapMatMult2(A, B, n_cores = 8), 
  times = 100
  ) %>%
  as.data.frame() %>%
  mutate(time_ms = time * 1e-6)

ggplot(data = mm_bench1, aes(x = expr, y = time_ms)) +
  geom_boxplot() +
  labs(x = "Method", y = "Time (milliseconds)",
       title = "Matrix Multiplication", subtitle = "(100x100 Matrix) * (100x100 Matrix)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Benchmark 2: N = k = 1000
N <- 1000
k <- 1000

A <- matrix(rnorm(N*k), N, k)
B <- matrix(rnorm(N*k), k, N)

mm_bench2 <- microbenchmark(
  A%*%B, 
  armaMatMult(A, B),
  eigenMatMult(A, B, n_cores = 1),
  eigenMatMult(A, B, n_cores = 2),
  eigenMatMult(A, B, n_cores = 4),
  eigenMatMult(A, B, n_cores = 8),
  eigenMapMatMult2(A, B, n_cores = 1),
  eigenMapMatMult2(A, B, n_cores = 2),
  eigenMapMatMult2(A, B, n_cores = 4), 
  eigenMapMatMult2(A, B, n_cores = 8), 
  times = 100
  ) %>%
  as.data.frame() %>%
  mutate(time_ms = time * 1e-6)

ggplot(data = mm_bench2, aes(x = expr, y = time_ms)) +
  geom_boxplot() +
  labs(x = "Method", y = "Time (milliseconds)",
       title = "Matrix Multiplication", subtitle = "(1000x1000 Matrix) * (1000x1000 Matrix)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

################
# Save Results #
################

full_svd_dat <- bind_rows(
  bind_cols(benchmark = "FullSVD100", partial_svd_bench1),
  bind_cols(benchmark = "FullSVD1000", partial_svd_bench2),
  bind_cols(benchmark = "FullSVDfat", partial_svd_bench3),
  bind_cols(benchmark = "FullSVDthin", partial_svd_bench4),
)

write_csv(svd_dat, file = "output/partial_svd_benchmark.csv")

###

partial_svd_dat <- bind_rows(
  bind_cols(benchmark = "PartialSVD100", partial_svd_bench1),
  bind_cols(benchmark = "PartialSVD1000", partial_svd_bench2)
)

write_csv(svd_dat, file = "output/partial_svd_benchmark.csv")

###

mm_dat <- bind_rows(
  bind_cols(benchmark = "MM100", mm_bench1),
  bind_cols(benchmark = "MM1000", mm_bench2)
)

write_csv(mm_dat, file = "output/mm_benchmark.csv")
