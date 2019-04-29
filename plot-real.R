library(tidyverse)

"Breast cancer (BRCA) & 62.1 [60.5-63.6] & 63.3 [61.7-64.8] & 65.9 [64.4-67.4] \\
 & 6256 & 2572 & 670,050 \\
Rheumatoid arthritis (RA) & 59.8 [57.7-61.8] & 60.3 [58.3-62.4] & 61.3 [59.1-63.4] \\
& 12,220 & 88,556 & 317,456 \\
Type 1 diabetes (T1D) & 75.4 [72.4-78.4] & 76.9 [73.9-79.7] & 78.7 [75.7-81.7] \\
& 1112 & 267 & 135,991 \\
Type 2 diabetes (T2D) & 59.5 [58.5-60.5] & 60.7 [59.8-61.6] & 63.8 [62.9-64.8] \\
& 252 & 33,238 & 535,785 \\
Prostate cancer (PRCA) & & & \\
& & & \\
Depression (MDD) & 53.9 [52.6-55.2] & 59.6 [58.3-60.8] & 59.8 [58.5-61.0] \\
& 170,505 & 205,096 & 473,333 \\
Coronary artery disease (CAD) & 59.9 [58.6-61.2] & 61.1 [59.9-62.4] & 63.9 [62.7-65.1] \\
& 1182 & 87,577 & 315,165 \\
Asthma & & & \\
& & & \\" %>%
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
  scale_y_continuous(limits = c(0.5, NA), minor_breaks = 0:50 / 50,
                     oob = scales::rescale_none) +
  theme(legend.position = c(0.7, 0.8))

ggsave("figures/AUC-real.pdf", width = 750, height = 500, scale = 1 / 100)
