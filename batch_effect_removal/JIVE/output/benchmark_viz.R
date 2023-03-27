library(tidyverse)


truncated_svd_dat <- read_csv("svd_benchmark.csv")

truncated_svd_plot <- truncated_svd_dat %>%
  filter(benchmark == "SVD1000") %>%
  ggplot(aes(x = expr, y = time_ms)) +
  geom_boxplot() +
  labs(x = "Method", y = "Time (milliseconds)",
       title = "Partial SVD", subtitle = "(1000x1000 Matrix)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave("truncated_svd_plot.png", truncated_svd_plot)

###

mm_dat <- read_csv("mm_benchmark.csv")

mm_plot <- mm_dat %>%
  filter(benchmark == "MM1000") %>%
  ggplot(aes(x = expr, y = time_ms)) +
  geom_boxplot() +
  labs(x = "Method", y = "Time (milliseconds)",
       title = "Matrix Multiplication", subtitle = "(1000x1000 Matrix) * (1000x1000 Matrix)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave("mm_plot.png", mm_plot)
