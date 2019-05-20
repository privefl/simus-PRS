# download.file("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/MichailidouK_29059683_GCST004988/oncoarray_bcac_public_release_oct17%20(1).txt.gz",
#               destfile = "sumstats_BRCA.txt.gz")
# R.utils::gunzip("sumstats_BRCA.txt.gz")
library(bigreadr)
sumstats <- fread2("sumstats_BRCA.txt", na.strings = "NULL",
                   select = c("chr", "position_b37", "a0", "a1",
                              "bcac_onco_icogs_gwas_beta",
                              "bcac_onco_icogs_gwas_P1df"),
                   col.names = c("chr", "pos", "a0", "a1", "beta", "p"))
nrow(sumstats)  # 11,792,542

sumstats <- subset(sumstats, p < 0.1)
nrow(sumstats)  # 1,768,354
hist(sumstats$beta)
hist(sumstats$p)

info_snp_UKBB <- rbind_df(lapply(1:22, function(chr) {
  file <- paste0("data/ukb_imp_mfi/ukb_mfi_chr", chr, "_v3.txt")
  df <- fread2(file, select = c(3:5, 8), col.names = c("pos", "a0", "a1", "info"))
  cbind.data.frame(chr = chr, df)
}))
info_snp <- bigsnpr::snp_match(sumstats, info_snp_UKBB)
# 1,768,354 variants in summary statistics.
# 244,598 ambiguous SNPs have been removed.
# 1,461,181 variants have been matched; 0 were flipped and 327 were reversed.
info_snp <- bigsnpr::snp_match(sumstats, info_snp_UKBB, strand_flip = FALSE)
# 1,698,342 variants have been matched; 0 were flipped and 327 were reversed.
info_snp <- subset(na.omit(info_snp), info > 0.3)

list_snp_id <- with(info_snp, split(paste(chr, pos, a0, a1, sep = "_"),
                                    factor(chr, levels = 1:22)))
beta <- info_snp$beta
lpval <- -log10(info_snp$p)
info <- info_snp$info


# subset samples
sample <- fread2("ukb25589_imp_chr1_v3_s487327.sample")
str(sample)
(N <- readBin("data/ukb_imp_chr1_v3.bgen", what = 1L, size = 4, n = 4)[4])
sample <- sample[-1, ]
stopifnot(nrow(sample) == N)

csv <- "ukb22544.csv"
df0 <- fread2(
  csv,
  select = c("eid", "22001-0.0", "22006-0.0"),
  col.names = c("eid", "sex", "is_caucasian")
)
ind.indiv <- match(df0$eid, sample$ID_2)
# related samples
df0$is_rel2 <- df0$eid %in% fread2("ukb25589_rel_s488346.dat")$ID2

df_cancer0 <- fread2(csv, select = paste0("2453-", 0:2, ".0"))
df_cancer1 <- fread2(csv, select = c(paste0("20001-0.", 0:5),
                                     paste0("20001-1.", 0:5),
                                     paste0("20001-2.", 0:5)))
df_cancer2 <- fread2(csv, select = c(paste0("40001-", 0:2, ".0"),
                                     paste0("40002-0.", 0:13),
                                     paste0("40002-1.", 0:13),
                                     paste0("40002-2.", 0:13),
                                     paste0("40006-", 0:31, ".0"),
                                     paste0("41202-0.", 0:379),
                                     paste0("41204-0.", 0:434)))
ind_BRCA <- sort(unique(unlist(c(
  lapply(df_cancer1,  function(x) which(x == 1002)),
  lapply(df_cancer2, function(x) which(substr(x, 1, 3) %in% c("C50", "D05")))
))))
table(df0$sex[ind_BRCA])
#     0     1
# 16681   116

y <- rep(NA, nrow(df0))
y[rowMeans(df_cancer0 == 0, na.rm = TRUE) == 1] <- 0
y[ind_BRCA] <- 1

sub <- which(!is.na(ind.indiv) & !df0$is_rel2 & df0$is_caucasian &
               df0$sex == 0 & !is.na(y))
length(sub)  # 169,969
table(y.sub <- y[sub])
#      0      1
# 158391  11578

NCORES <- 12
system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKBB_BRCA",
    ind_row = ind.indiv[sub],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = NCORES
  )
) # 4.8H


library(bigsnpr)
ukbb <- snp_attach("data/UKBB_BRCA.rds")
G <- ukbb$genotypes
file.size(G$backingfile) / 1024^3  # 269 GB
CHR <- as.integer(ukbb$map$chromosome)
POS <- ukbb$map$physical.pos

set.seed(1); ind.train <- sort(sample(length(sub), 150e3))
# set.seed(2); ind.train <- c(sample(which(y.sub == 0), 2000), sample(which(y.sub == 1), 500))
ind.test <- setdiff(seq_along(sub), ind.train)

system.time(
  all_keep <- snp_grid_clumping(
    G, CHR, POS, lpS = lpval, ind.row = ind.train, infos.imp = info,
    grid.thr.imp = c(0.3, 0.6, 0.9, 0.95),
    grid.thr.r2 = c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.95),
    grid.base.size = c(50, 100, 200, 500),
    ncores = 8
  )
) # 4.3H

system.time(
  multi_PRS <- snp_grid_PRS(
    G, all_keep, betas = beta, lpS = lpval, ind.row = ind.train,
    n_thr_lpS = 50, backingfile = "data/UKBB_BRCA_scores", ncores = NCORES
  )
) # 1.5H

system.time(
  final_mod <- snp_grid_stacking(
    multi_PRS, y.sub[ind.train], ncores = NCORES,
  )
) # 3H
mod <- final_mod$mod
plot(mod)
summary(mod)
# # A tibble: 3 x 6
#    alpha validation_loss intercept beta           nb_var message
#    <dbl>           <dbl>     <dbl> <list>          <int> <list>
# 1 0.0001           0.239     -1.98 <dbl [80,416]>  33073 <chr [10]>
# 2 0.01             0.239     -1.90 <dbl [80,416]>   6027 <chr [10]>
# 3 1                0.239     -1.88 <dbl [80,416]>   4631 <chr [10]>

# Save
save(all_keep, multi_PRS, final_mod, file = "data/res_BRCA.RData")

new_beta <- final_mod$beta.G

length(ind <- which(new_beta != 0))  # 670,050
summary(new_beta)
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -1.688e-01  0.000e+00  0.000e+00 -8.950e-06  0.000e+00  2.388e-01
summary(new_beta[which(sign(new_beta * beta) < 0)])
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -3.075e-02 -1.648e-05  4.400e-07  2.633e-05  3.295e-05  3.568e-02

pred <- final_mod$intercept +
  big_prodVec(G, new_beta[ind], ind.row = ind.test, ind.col = ind)

AUCBoot(pred, y.sub[ind.test])  # 65.9 [64.4-67.4] / 62.9 [62.4-63.5]


library(tidyverse)

ind2 <- sort(sample(ind, 10e3))
ggplot(data.frame(y = new_beta, x = beta)[ind, ]) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  geom_abline(slope = 0, intercept = 0, color = "blue") +
  geom_point(aes(x, y), size = 0.6) +
  theme_bigstatsr() +
  labs(x = "Effect sizes from GWAS", y = "Non-zero effect sizes from SCT")

grid2 <- attr(all_keep, "grid") %>%
  mutate(thr.lp = list(attr(multi_PRS, "grid.lpS.thr")), num = row_number()) %>%
  unnest()
s <- nrow(grid2)
grid2$auc <- big_apply(multi_PRS, a.FUN = function(X, ind, s, y.train) {
  # Sum over all chromosomes, for the same C+T parameters
  single_PRS <- rowSums(X[, ind + s * (0:21)])
  bigstatsr::AUC(single_PRS, y.train)
}, ind = 1:s, s = s, y.train = y.sub[ind.train],
a.combine = 'c', block.size = 1, ncores = NCORES)

std_prs <- grid2 %>%
  filter(thr.imp == 0.3, thr.r2 == 0.2, size == 500, thr.lp <= 8) %>%
  arrange(desc(auc)) %>%
  slice(1) %>%
  print()
#   size thr.r2 thr.imp num   thr.lp       auc
# 1  500    0.2     0.3  14 3.248484 0.6208515

ind.keep <- unlist(map(all_keep, std_prs$num))
# Verif on training set
AUC(
  snp_PRS(G, beta[ind.keep], ind.test = ind.train, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = std_prs$thr.lp),
  y.sub[ind.train]
) # 0.6208515
# Eval on test set
AUCBoot(
  snp_PRS(G, beta[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = std_prs$thr.lp),
  y.sub[ind.test]
) # 62.1 [60.5-63.6] / 62.2 [61.6-62.7]
sum(lpval[ind.keep] > std_prs$thr.lp)  # 6256

max_prs <- grid2 %>% arrange(desc(auc)) %>% slice(1:10) %>% print() %>% slice(1)
#    size thr.r2 thr.imp num   thr.lp       auc
# 1  2500    0.2    0.95 100 3.654696 0.6351843
# 2  5000    0.1    0.95  96 3.654696 0.6349136
# 3  1000    0.2    0.95  99 3.654696 0.6348574
# 4  2000    0.1    0.95  95 3.654696 0.6348562
# 5  2500    0.2    0.95 100 3.248484 0.6347538
# 6  2500    0.2    0.90  72 3.248484 0.6345844
# 7  2000    0.1    0.95  95 4.111703 0.6345656
# 8  5000    0.1    0.95  96 4.111703 0.6345364
# 9  5000    0.1    0.90  68 4.111703 0.6344550
# 10 2500    0.2    0.90  72 3.654696 0.6342989

ind.keep <- unlist(map(all_keep, max_prs$num))
AUCBoot(
  snp_PRS(G, beta[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = max_prs$thr.lp),
  y.sub[ind.test]
) # 63.3 [61.7-64.8] / 63.4 [62.8-63.9]
sum(lpval[ind.keep] > max_prs$thr.lp)  # 2572

ggplot(grid2) +
  geom_point(aes(thr.lp, auc)) +
  facet_grid(thr.imp ~ thr.r2 + size) +
  scale_x_log10(limits = c(1, 8), breaks = c(1, 2, 5), minor_breaks = 1:10) +
  ylim(0.61, NA) +
  theme_bigstatsr(size.rel = 0.7) +
  labs(x = "-log10(p-value) threshold (log scale)", y = "AUC")

# Check large effects
ind3 <- intersect(which(abs(beta) > 2), ind)
gwas <- big_univLogReg(G, y.sub[ind.train], ind.train, ind3)
cbind(beta_GWAS = beta[ind3], pval_GWAS = 10^(-lpval[ind3]),
      af = big_scale()(G, ind.row = ind.train, ind.col = ind3)$center / 2,
      gwas, pval = predict(gwas, log10 = FALSE))
#   beta_GWAS pval_GWAS           af       estim    std.err niter      score      pval
# 1    -2.6335  0.004775 3.185000e-04  -1.0766894  0.7337089     6 -1.4674613 0.1422506
# 2    -2.9353  0.007172 5.003333e-05  -3.8039739  6.0073820     9 -0.6332166 0.5265922
# 3    -2.2807  0.011160 8.803333e-05  -2.9472414  3.0811107     8 -0.9565516 0.3387936
# 4    -2.1291  0.011220 3.220333e-04  -0.1486095  0.4739562     4 -0.3135511 0.7538620
# 5    -2.1569  0.004768 8.073333e-05  -0.3729499  1.2722939     5 -0.2931319 0.7694213
# 6    -2.0555  0.004252 1.413333e-05 -34.2383336 32.7139349     9 -1.0465978 0.2952851
# 7    -2.5219  0.003705 1.543333e-04  -1.9786611  1.6408402     6 -1.2058828 0.2278627
# 8    -2.0048  0.004386 6.113000e-04  -0.3215922  0.3879102     5 -0.8290377 0.4070831
# 9    -2.3273  0.004739 5.470000e-05  -5.8596615  7.2622512     8 -0.8068657 0.4197439
# 10   -2.9494  0.010600 8.666667e-05  -2.9112447  2.9569211     7 -0.9845527 0.3248438
# 11   -2.5205  0.021100 1.607000e-04  -3.3157859  2.3529122     7 -1.4092264 0.1587682
