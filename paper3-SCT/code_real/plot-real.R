library(tidyverse)

"Breast cancer (BRCA) & 62.1 [60.5-63.6] & 63.3 [61.7-64.8] & 65.9 [64.4-67.4] & 57.9 [56.3-59.5] \\
 & 6256 & 2572 & 670,050 & 322,003 \\
Rheumatoid arthritis (RA) & 59.8 [57.7-61.8] & 60.3 [58.3-62.4] & 61.3 [59.1-63.4] & 59.5 [57.5-61.7] \\
& 12,220 & 88,556 & 317,456 & 672,922 \\
Type 1 diabetes (T1D) & 75.4 [72.4-78.4] & 76.9 [73.9-79.7] & 78.7 [75.7-81.7] & 75.3 [72.2-78.3] \\
& 1112 & 267 & 135,991 & 204,785 \\
Type 2 diabetes (T2D) & 59.1 [58.1-60.1] & 60.7 [59.8-61.7] & 63.8 [62.9-64.7] & 63.2 [62.3-64.1] \\
& 177 & 33,235 & 548,343 & 256,353 \\
Prostate cancer (PRCA) & 68.0 [66.5-69.5] & 69.3 [67.8-70.8] & 71.7 [70.2-73.1] & 58.7 [57.1-60.3] \\
& 1035 & 356 & 696,575 & 121,660 \\
Depression (MDD) & 55.7 [54.4-56.9] & 59.2 [58.0-60.4] & 59.5 [58.2-60.7] & 52.0 [50.8-53.3] \\
& 165,584 & 222,912 & 524,099 & 625,732 \\
Coronary artery disease (CAD) & 59.9 [58.6-61.2] & 61.1 [59.9-62.4] & 63.9 [62.7-65.1] & 63.0 [61.8-64.2] \\
& 1182 & 87,577 & 315,165 & 290,204 \\
Asthma & 56.8 [56.2-57.5] & 57.3 [56.7-58.0] & 60.7 [60.0-61.3] & 58.7 [58.1-59.4] \\
& 3034 & 360 & 446,120 & 75,965 \\" %>%
  scan(text = ., what = "", sep = "\\") %>%
  scan(text = ., what = "", sep = "&") %>%
  trimws() %>%
  matrix(ncol = 10, byrow = TRUE) %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  select_at(1:5) %>%
  setNames(c("Trait", "stdCT", "maxCT", "SCT", "lassosum")) %>%
  mutate(Trait = gsubfn::strapply(Trait, "\\((.+)\\)", empty = "Asthma", simplify = c),
         Trait = ordered(Trait, levels = Trait)) %>%
  gather("Method", "AUC", -1) %>%
  separate(AUC, c("AUC", "inf", "sup"), sep = "[^[0-9\\.]]+", convert = TRUE) %>%
  mutate_if(is.numeric, ~ . / 100) %>%
  mutate(Method = ordered(Method, levels = c("stdCT", "maxCT", "SCT", "lassosum"))) %>%
  ggplot(aes(Trait, AUC, fill = Method)) +
  bigstatsr::theme_bigstatsr() +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_col(position = position_dodge(), alpha = 0.5, color = "black", size = 1) +
  geom_errorbar(aes(ymin = inf, ymax = sup), position = position_dodge(width = 0.9),
                color = "black", width = 0.2, size = 1) +
  scale_y_continuous(limits = c(0.5, NA), oob = scales::rescale_none,
                     breaks = 0:10 / 10, minor_breaks = 0:50 / 50) +
  theme(legend.position = c(0.8, 0.8)) +
  scale_fill_manual(values = viridis::viridis(5)[1:4])

ggsave("figures/AUC-real.pdf", width = 850, height = 500, scale = 1 / 100)
# ggsave("figures/AUC-real.png", width = 850, height = 500, scale = 1 / 100)
