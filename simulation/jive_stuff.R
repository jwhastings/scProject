plot_thing <- function(x) {
  df <- data.frame(t(x)) %>%
    bind_cols(Group = grp_b1) %>%
    pivot_longer(cols = starts_with("Gene"))
  
  ggplot(data = df) +
    geom_density(aes(x = value, fill = Group)) +
    facet_wrap(~ Group)
}

temp <- SVDmiss(b1, ncomp = min(ncol(b1), nrow(b1)))[[1]]
new_data <- temp$u %*% diag(x = temp$d) %*% t(temp$v)

centerValues <- apply(new_data, 1, mean, na.rm = T)

center_subtraction <- matrix(rep(centerValues, ncol(new_data)), nrow = nrow(new_data))

centered_data <- new_data - centerValues

scale_value <- norm(centered_data, type = "f") * sqrt(sum(2 * nrow(centered_data) * ncol(centered_data)))

final_data <- centered_data / scale_value

####

jive_data1_centered_scaled <- JIVE_results_default$data$Batch1
dim(jive_data1_centered_scaled)

jive_data1_scale_factor <- JIVE_results_default$scale$`Scale Values`[[1]]
jive_data1_center_values <- JIVE_results_default$scale$`Center Values`[[1]]
length(jive_data1_center_values)

jive_data1_centered <- jive_data1_centered_scaled * jive_data1_scale_factor

jive_data1 <- jive_data1_centered + jive_data1_center_values

####

test <- final_data - jive_data1
median(test)
