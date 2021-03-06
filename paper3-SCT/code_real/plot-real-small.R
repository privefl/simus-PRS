library(tidyverse)

"Breast cancer (BRCA) & 62.2 [61.6-62.7] & 63.4 [62.8-63.9] & 62.9 [62.4-63.5] & 57.8 [57.3-58.4] \\
Rheumatoid arthritis (RA) & 59.2 [58.4-60.0] & 59.5 [58.7-60.3] & 59.5 [58.7-60.3] & 58.0 [57.1-58.8] \\
Type 1 diabetes (T1D) & 75.6 [72.4-78.7] & 76.7 [73.6-79.8] & 78.7 [75.5-81.8] & 75.5 [72.1-78.7] \\
Type 2 diabetes (T2D) & 59.8 [59.3-60.3] & 60.2 [59.7-60.7] & 61.0 [60.6-61.5] & 63.6 [63.1-64.1] \\
Prostate cancer (PRCA) & 67.1 [66.4-67.8] & 68.7 [68.0-69.3] & 69.3 [68.7-70.0] & 56.2 [55.4-56.9] \\
Depression (MDD) & 54.5 [54.1-54.9] & 58.4 [58.0-58.8] & 54.7 [54.3-55.1] & 51.6 [51.2-52.0] \\
Coronary artery disease (CAD) & 59.7 [59.2-60.3] & 60.0 [59.5-60.5] & 61.4 [60.8-61.9] & 62.3 [61.8-62.8] \\
Asthma & 56.2 [55.9-56.4] & 56.9 [56.7-57.2] & 57.2 [56.9-57.4] & 57.0 [56.7-57.3] \\" %>%
  scan(text = ., what = "", sep = "\\") %>%
  scan(text = ., what = "", sep = "&") %>%
  trimws() %>%
  matrix(ncol = 5, byrow = TRUE) %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  setNames(c("Trait", "stdCT", "maxCT", "SCT", "lassosum")) %>%
  mutate(Trait = gsubfn::strapply(Trait, "\\((.+)\\)", empty = "Asthma", simplify = c),
         Trait = ordered(Trait, levels = Trait)) %>%
  gather("Method", "AUC", -1) %>%
  separate(AUC, c("AUC", "inf", "sup"), sep = "[^[0-9\\.]]+", convert = TRUE) %>%
  mutate_if(is.numeric, ~ . / 100) %>%
  mutate(Method = ordered(Method, levels = c("stdCT", "maxCT", "SCT", "lassosum"))) %>%
  ggplot(aes(Trait, AUC, fill = Method, color = Method)) +
  bigstatsr::theme_bigstatsr() +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_col(position = position_dodge(), alpha = 0.5, color = "black", size = 1) +
  geom_errorbar(aes(ymin = inf, ymax = sup), position = position_dodge(width = 0.9),
                color = "black", width = 0.2, size = 1) +
  scale_y_continuous(limits = c(0.5, NA), oob = scales::rescale_none,
                     breaks = 0:10 / 10, minor_breaks = 0:50 / 50) +
  theme(legend.position = c(0.8, 0.8))

ggsave("figures/AUC-real-small.pdf", width = 850, height = 500, scale = 1 / 100)
# ggsave("figures/AUC-real-small.png", width = 850, height = 500, scale = 1 / 100)
