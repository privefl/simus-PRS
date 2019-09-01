# system("scp privef@luxor:UKBiobank/res_simu/[^f]*.rds ~/Bureau/simus-PRS/res_simu_with_ldpred/")

library(tidyverse)

res_simu <- list.files("res_simu_with_ldpred", full.names = TRUE) %>%
  map_dfr(~ {
    simu <- gsubfn::strapply(.x, "^res_simu_with_ldpred/(.+)_[0-9]+\\.rds$")[[1]]
    res <- as_tibble(readRDS(.x))
    bind_cols(simu = simu, res)
  })

auc_simu <- res_simu %>%
  select(1, 2, 4, 6, 8, 10) %>%
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
           Method, levels = c("auc_std_prs", "auc_max_prs", "auc_SCT", "auc_lassosum", "auc_ldpred"),
           labels = c("stdCT", "maxCT", "SCT", "lassosum", "LDpred")
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
  theme(legend.position = c(0.37, 0.85))

ggsave("figures/AUC-simu1.pdf", width = 870, height = 600, scale = 1 / 100)

auc_simu %>%
  rename(Scenario = simu) %>%
  mutate_at(3:5, ~ round(. * 100, 1)) %>%
  mutate(AUC = sprintf("%.1f [%.1f-%.1f]", mean, inf, sup)) %>%
  select(-(3:5)) %>%
  tidyr::spread(Method, AUC) %>%
  print() %>%
  xtable::xtable(align = "|c|l|c|c|c|c|c|") %>%
  print(include.rownames = FALSE)
#   Scenario stdCT            maxCT            SCT              lassosum         LDpred
# 1 100      82.0 [81.6-82.5] 86.8 [86.4-87.2] 86.1 [85.7-86.5] 83.9 [83.5-84.4] 76.3 [75.7-76.9]
# 2 10K      72.5 [71.5-73.5] 75.6 [75.1-76.1] 76.4 [75.9-77.0] 75.0 [74.4-75.6] 71.2 [70.6-71.9]
# 3 1M       68.7 [68.1-69.2] 69.0 [68.3-69.6] 69.1 [68.6-69.6] 70.1 [69.4-70.6] 71.0 [70.3-71.6]
# 4 2chr     77.4 [76.9-78.0] 78.9 [78.3-79.4] 82.1 [81.6-82.7] 79.0 [78.3-79.7] 74.6 [73.6-75.8]
# 5 err      70.1 [69.3-70.9] 71.0 [70.6-71.5] 73.4 [72.8-74.1] 72.0 [71.3-72.8] 70.7 [70.1-71.4]
# 6 HLA      78.7 [78.0-79.5] 79.8 [79.1-80.3] 80.7 [80.2-81.3] 79.4 [78.7-80.2] 76.9 [76.3-77.4]
