# system("scp privef@krakenator:UKBiobank/res_small/[^file]*.rds ~/Bureau/simus-PRS/res_small/")

library(tidyverse)

res_small_mult <- list.files("res_small", full.names = TRUE) %>%
  map_dfr(~ {
    trait <- gsubfn::strapply(.x, "^res_small/(.+)_[0-9]+\\.rds$")[[1]]
    res <- as_tibble(readRDS(.x))
    bind_cols(Trait = trait, res)
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
         )) %>%
  unnest() %>%
  print()

ggplot(auc_small_mult, aes(Trait, mean, fill = Method, color = Method)) +
  bigstatsr::theme_bigstatsr() +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_col(position = position_dodge(), alpha = 0.5, color = "black", size = 1) +
  geom_errorbar(aes(ymin = inf, ymax = sup),
                position = position_dodge(width = 0.9), color = "black", width = 0.2, size = 1) +
  scale_y_continuous(limits = c(0.5, 0.9), minor_breaks = 0:50 / 50,
                     oob = scales::rescale_none) +
  labs(y = "AUC") +
  theme(legend.position = c(0.37, 0.85))

ggsave("figures/AUC-real-small-mult.pdf", width = 870, height = 600, scale = 1 / 100)

auc_small_mult %>%
  mutate_at(3:5, ~ round(. * 100, 1)) %>%
  mutate(AUC = sprintf("%.1f [%.1f-%.1f]", mean, inf, sup)) %>%
  select(-(3:5)) %>%
  tidyr::spread(Method, AUC) %>%
  print() %>%
  xtable::xtable(align = "|c|l|c|c|c|c|c|") %>%
  print(include.rownames = FALSE)
