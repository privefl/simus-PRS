# system("scp privef@krakenator:UKBiobank/res_small/[^f]*.rds ~/Bureau/simus-PRS/res_small/")
# system("scp privef@luxor:UKBiobank/res_small/[^f]*.rds ~/Bureau/simus-PRS/res_small/")

library(tidyverse)

res_small_mult <- list.files("res_small", full.names = TRUE) %>%
  map_dfr(~ {
    trait <- gsubfn::strapply(.x, "^res_small/(.+)_[0-9]+\\.rds$")[[1]]
    res <- as_tibble(readRDS(.x))
    bind_cols(Trait = tools::toTitleCase(trait), res)
  })

auc_small_mult <- res_small_mult %>%
  select(1, 2, 4, 6, 8, 10) %>%
  group_by(Trait) %>%
  summarise_all(~ {
    boot <- replicate(1e4, mean(sample(.x, replace = TRUE)))
    q <- quantile(boot, c(0.025, 0.975))
    res <- c(mean = mean(boot), inf = q[[1]], sup = q[[2]])
    list(res)
  }) %>%
  gather("Method", "AUC", -Trait) %>%
  mutate(AUC = map(AUC, ~as_tibble(as.list(.x))),
         Method = ordered(
           Method, levels = c("auc_std_prs", "auc_max_prs", "auc_SCT", "auc_lassosum", "auc_LDpred"),
           labels = c("stdCT", "maxCT", "SCT", "lassosum", "LDpred")
         ),
         Trait = ordered(Trait, levels = c("BRCA", "RA", "T1D", "T2D", "PRCA",
                                           "MDD", "CAD", "Asthma"))) %>%
  unnest() %>%
  print()

ggplot(auc_small_mult, aes(Trait, mean, fill = Method)) +
  bigstatsr::theme_bigstatsr() +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_col(position = position_dodge(), alpha = 0.5, color = "black", size = 1) +
  geom_errorbar(aes(ymin = inf, ymax = sup),
                position = position_dodge(width = 0.9), color = "black", width = 0.2, size = 1) +
  scale_y_continuous(limits = c(0.5, 0.8), minor_breaks = 0:50 / 50,
                     oob = scales::rescale_none) +
  labs(y = "AUC") +
  theme(legend.position = c(0.88, 0.80))

ggsave("figures/AUC-real-small-mult.pdf", width = 900, height = 550, scale = 1 / 100)

auc_small_mult %>%
  mutate_at(3:5, ~ round(. * 100, 1)) %>%
  mutate(AUC = sprintf("%.1f [%.1f-%.1f]", mean, inf, sup)) %>%
  select(-(3:5)) %>%
  tidyr::spread(Method, AUC) %>%
  print() %>%
  xtable::xtable(align = "|c|l|c|c|c|c|c|") %>%
  print(include.rownames = FALSE)
#   Trait  stdCT            maxCT            SCT              lassosum         LDpred
# 1 BRCA   61.9 [61.8-62.0] 62.9 [62.6-63.2] 62.7 [62.4-63.1] 57.8 [57.1-58.5] 62.0 [62.0-62.1]
# 2 RA     59.1 [59.0-59.2] 59.5 [59.3-59.7] 59.6 [59.3-59.8] 58.8 [58.2-59.3] 59.8 [59.7-59.8]
# 3 T1D    75.4 [74.6-76.2] 76.0 [75.1-77.0] 78.4 [77.7-79.1] 75.5 [74.7-76.2] 75.5 [74.8-76.3]
# 4 T2D    59.6 [59.5-59.7] 60.4 [60.2-60.7] 60.9 [60.6-61.1] 63.2 [62.9-63.5] 61.0 [60.6-61.3]
# 5 PRCA   67.0 [66.8-67.2] 68.5 [68.3-68.7] 69.0 [68.5-69.4] 57.1 [56.3-58.0] 65.6 [65.5-65.7]
# 6 MDD    54.1 [53.7-54.4] 58.7 [58.3-59.0] 54.7 [54.4-55.0] 52.8 [52.0-53.7] 60.0 [60.0-60.1]
# 7 CAD    59.7 [59.5-59.9] 61.1 [60.6-61.5] 60.8 [60.6-61.1] 62.8 [62.6-62.9] 61.0 [60.9-61.1]
# 8 Asthma 56.3 [56.2-56.5] 56.9 [56.8-56.9] 56.3 [55.6-56.9] 57.7 [57.3-58.0] 56.5 [56.3-56.7]
