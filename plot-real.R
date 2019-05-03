library(tidyverse)

"Breast cancer (BRCA) & 62.1 [60.5-63.6] & 63.3 [61.7-64.8] & 65.9 [64.4-67.4] \\
 & 6256 & 2572 & 670,050 \\
Rheumatoid arthritis (RA) & 59.8 [57.7-61.8] & 60.3 [58.3-62.4] & 61.3 [59.1-63.4] \\
& 12,220 & 88,556 & 317,456 \\
Type 1 diabetes (T1D) & 75.4 [72.4-78.4] & 76.9 [73.9-79.7] & 78.7 [75.7-81.7] \\
& 1112 & 267 & 135,991 \\
Type 2 diabetes (T2D) & 59.5 [58.5-60.5] & 60.7 [59.8-61.6] & 63.8 [62.9-64.8] \\
& 252 & 33,238 & 535,785 \\
Prostate cancer (PRCA) & 68.0 [66.5-69.5] & 69.3 [67.8-70.8] & 71.7 [70.2-73.1] \\
& 1035 & 356 & 696,575 \\
Depression (MDD) & 55.7 [54.4-56.9] & 59.2 [58.0-60.4] & 59.5 [58.2-60.7] \\
& 165,584 & 222,912 & 524,099 \\
Coronary artery disease (CAD) & 59.9 [58.6-61.2] & 61.1 [59.9-62.4] & 63.9 [62.7-65.1] \\
& 1182 & 87,577 & 315,165 \\
Asthma & 56.8 [56.2-57.5] & 57.3 [56.7-58.0] & 60.7 [60.0-61.3] \\
& 3034 & 360 & 446,120 \\" %>%
  scan(text = ., what = "", sep = "\\") %>%
  scan(text = ., what = "", sep = "&") %>%
  trimws() %>%
  matrix(ncol = 8, byrow = TRUE) %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  select_at(1:4) %>%
  setNames(c("Trait", "stdCT", "maxCT", "SCT")) %>%
  mutate(Trait = gsubfn::strapply(Trait, "\\((.+)\\)", empty = "Asthma", simplify = c),
         Trait = ordered(Trait, levels = Trait)) %>%
  gather("Method", "AUC", -1) %>%
  separate(AUC, c("AUC", "inf", "sup"), sep = "[^[0-9\\.]]+", convert = TRUE) %>%
  mutate_if(is.numeric, ~ . / 100) %>%
  mutate(Method = ordered(Method, levels = c("stdCT", "maxCT", "SCT"))) %>%
  ggplot(aes(Trait, AUC, fill = Method, color = Method)) +
  bigstatsr::theme_bigstatsr() +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_col(position = position_dodge(), alpha = 0.5, color = "black", size = 1) +
  geom_errorbar(aes(ymin = inf, ymax = sup), position = position_dodge(width = 0.9),
                color = "black", width = 0.2, size = 1) +
  scale_y_continuous(limits = c(0.5, NA), oob = scales::rescale_none,
                     breaks = 0:10 / 10, minor_breaks = 0:50 / 50) +
  theme(legend.position = c(0.8, 0.8))

ggsave("figures/AUC-real.pdf", width = 750, height = 500, scale = 1 / 100)
