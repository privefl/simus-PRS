library(tidyverse)

res_simu <- list.files("res_simu", full.names = TRUE) %>%
  map_dfr(~ {
    simu <- gsubfn::strapply(.x, "^res_simu/(.+)_[0-9]+\\.rds$")[[1]]
    res <- as_tibble(readRDS(.x))
    bind_cols(simu = simu, res)
  })

auc_simu <- res_simu %>%
  select(1, 2, 4, 6) %>%
  group_by(simu) %>%
  summarise_all(~ {
    boot <- replicate(1e3, mean(sample(.x, replace = TRUE)))
    q <- quantile(boot, c(0.025, 0.975))
    res <- c(mean = mean(boot), inf = q[[1]], sup = q[[2]])
    list(res)
  }) %>%
  gather("Method", "AUC", -simu) %>%
  mutate(AUC = map(AUC, ~as_tibble(as.list(.x))),
         Method = ordered(Method, levels = c("auc_std_prs", "auc_max_prs", "auc_SCT"),
                          labels = c("std_PRS", "max_PRS", "SCT"))) %>%
  unnest()

ggplot(auc_simu, aes(simu, mean, fill = Method, color = Method)) +
  bigstatsr::theme_bigstatsr() +
  geom_col(position = position_dodge(), alpha = 0.5, color = "black", size = 1) +
  geom_errorbar(aes(ymin = inf, ymax = sup),
                position = position_dodge(width = 0.9), color = "black", width = 0.2, size = 1) +
  scale_y_continuous(limits = c(0.5, NA), minor_breaks = 0:20 / 20,
                     oob = scales::rescale_none)
