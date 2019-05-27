library(tidyverse)

res_simu <- list.files("res_simu", full.names = TRUE) %>%
  map_dfr(~ {
    simu <- gsubfn::strapply(.x, "^res_simu/(.+)_[0-9]+\\.rds$")[[1]]
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

ggplot(auc_simu, aes(simu, mean, fill = Method, color = Method)) +
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
  scale_fill_manual(values = c("#440154FF", "#2A788EFF", "#7AD151FF", "#FDE725FF")) +
  scale_color_manual(values = c("#440154FF", "#2A788EFF", "#7AD151FF", "#FDE725FF"))

ggsave("figures/AUC-simus.pdf", width = 870, height = 600, scale = 1 / 100)

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
# 1 100      79.8 [77.0-82.0] 86.9 [86.6-87.3] 86.3 [85.8-86.8] 83.2 [81.8-84.2]
# 2 10K      72.5 [71.8-73.3] 75.1 [74.7-75.5] 76.0 [75.5-76.6] 74.9 [74.3-75.6]
# 3 1M       68.9 [68.3-69.4] 69.5 [68.8-70.0] 69.0 [68.5-69.6] 70.4 [70.0-70.9]
# 4 2chr     77.2 [76.7-77.7] 78.6 [78.0-79.2] 82.2 [81.8-82.7] 78.9 [78.4-79.4]
# 5 err      69.8 [68.9-70.7] 70.7 [70.1-71.2] 73.2 [72.5-73.9] 72.1 [71.5-72.8]
# 6 HLA      78.7 [78.0-79.5] 79.8 [79.1-80.4] 80.7 [80.2-81.3] 79.4 [78.7-80.2]
