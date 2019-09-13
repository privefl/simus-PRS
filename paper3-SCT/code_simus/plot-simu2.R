library(tidyverse)

res_simu <- list.files("res_simu2", full.names = TRUE) %>%
  map_dfr(~ {
    simu <- gsubfn::strapply(.x, "^res_simu2/(.+)_[0-9]+\\.rds$")[[1]]
    res <- as_tibble(readRDS(.x))
    bind_cols(simu = simu, res)
  })

auc_simu <- res_simu %>%
  select(1, 2, 4, 6, 8) %>%
  group_by(simu) %>%
  summarise_all(~ {
    boot <- replicate(1e4, mean(sample(.x, replace = TRUE)))
    q <- quantile(boot, c(0.025, 0.975))
    res <- c(mean = mean(boot), inf = q[[1]], sup = q[[2]])
    list(res)
  }) %>%
  gather("Method", "AUC", -simu) %>%
  mutate(AUC = map(AUC, ~as_tibble(as.list(.x))),
         Method = ordered(
           Method, levels = c("auc_std_prs", "auc_max_prs", "auc_SCT", "auc_lassosum"),
           labels = c("stdCT", "maxCT", "SCT", "lassosum")
         )) %>%
  unnest()

ggplot(auc_simu, aes(simu, mean, fill = Method)) +
  bigstatsr::theme_bigstatsr() +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_hline(yintercept = 0.875, linetype = 3, color = "blue") +
  geom_col(position = position_dodge(), alpha = 0.5, color = "black", size = 1) +
  geom_errorbar(aes(ymin = inf, ymax = sup),
                position = position_dodge(width = 0.9), color = "black", width = 0.2, size = 1) +
  scale_y_continuous(limits = c(0.5, 0.9), minor_breaks = 0:50 / 50,
                     oob = scales::rescale_none) +
  labs(x = "Simulation", y = "AUC") +
  theme(legend.position = c(0.37, 0.85)) +
  scale_fill_manual(values = viridis::viridis(5)[1:4])

ggsave("figures/AUC-simu2.pdf", width = 870, height = 600, scale = 1 / 100)

auc_simu %>%
  rename(Scenario = simu) %>%
  mutate_at(3:5, ~ round(. * 100, 1)) %>%
  mutate(AUC = sprintf("%.1f [%.1f-%.1f]", mean, inf, sup)) %>%
  select(-(3:5)) %>%
  tidyr::spread(Method, AUC) %>%
  print() %>%
  xtable::xtable(align = "|c|l|c|c|c|c|") %>%
  print(include.rownames = FALSE)
#   Scenario stdCT            maxCT            SCT              lassosum
# 1 100      77.4 [76.0-78.8] 83.9 [83.4-84.4] 83.1 [82.6-83.6] 80.1 [79.5-80.8]
# 2 10K      69.4 [68.4-70.5] 73.0 [72.5-73.4] 72.9 [72.5-73.3] 71.2 [70.6-71.7]
# 3 1M       64.0 [63.6-64.4] 64.0 [63.6-64.4] 62.7 [62.3-63.0] 64.1 [63.3-64.8]
# 4 2chr     70.0 [68.8-71.2] 74.4 [73.6-75.2] 78.5 [77.9-79.1] 73.2 [72.5-73.8]
# 5 err      67.0 [66.0-68.1] 68.6 [67.7-69.5] 69.5 [68.9-70.1] 65.6 [64.9-66.3]
# 6 HLA      74.8 [72.9-76.3] 75.3 [73.5-76.9] 76.4 [74.5-78.0] 75.8 [74.2-77.2]
