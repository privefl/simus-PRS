# download.file("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/OkadaY_24390342_GCST002318/RA_GWASmeta_TransEthnic_v2.txt.gz",
#               destfile = "sumstats_RA.txt.gz")
# R.utils::gunzip("sumstats_RA.txt.gz")
library(bigreadr)
sumstats <- fread2("sumstats_RA.txt", select = 2:7,
                   col.names = c("chr", "pos", "a1", "a0", "or", "p"))
nrow(sumstats)  # 9,739,303

library(tidyverse)
sumstats <- sumstats %>%
  filter(p < 0.1) %>%
  mutate(beta = log(or), or = NULL)
hist(sumstats$beta)
hist(sumstats$p)

info_snp_UKBB <- rbind_df(lapply(1:22, function(chr) {
  file <- paste0("data/ukb_imp_mfi/ukb_mfi_chr", chr, "_v3.txt")
  df <- fread2(file, select = c(3:5, 8), col.names = c("pos", "a0", "a1", "info"))
  cbind.data.frame(chr = chr, df)
}))
info_snp <- bigsnpr::snp_match(sumstats, info_snp_UKBB)
# 1,055,461 variants in summary statistics.
# 161,799 ambiguous SNPs have been removed.
# 857,327 variants have been matched; 0 were flipped and 383,746 were reversed.
info_snp <- bigsnpr::snp_match(sumstats, info_snp_UKBB, strand_flip = FALSE)
# 1,011,519 variants have been matched; 0 were flipped and 461,016 were reversed.
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
df0 <- fread2(csv, select = c("eid", "22006-0.0"),
              col.names = c("eid", "is_caucasian"))
ind.indiv <- match(df0$eid, sample$ID_2)
# related samples
df0$is_rel2 <- df0$eid %in% fread2("ukb25589_rel_s488346.dat")$ID2

df_illness <- fread2(csv, select = c(paste0("20002-0.", 0:28),
                                     paste0("20002-1.", 0:28),
                                     paste0("20002-2.", 0:28)))
df_ICD10 <- fread2(csv, select = c(paste0("40001-", 0:2, ".0"),
                                   paste0("40002-0.", 0:13),
                                   paste0("40002-1.", 0:13),
                                   paste0("40002-2.", 0:13),
                                   paste0("41202-0.", 0:379),
                                   paste0("41204-0.", 0:434)))
ind_RA <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x == 1464)),
  lapply(df_ICD10,   function(x) which(substr(x, 1, 3) %in% c("M05", "M06")))
))))
ind_muscu <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x %in% c(1295, 1464:1467, 1477, 1538))),
  lapply(df_ICD10,   function(x) which(substr(x, 1, 1) == "M"))
))))
y <- rep(0, nrow(df0)); y[ind_muscu] <- NA; y[ind_RA] <- 1

sub <- which(!is.na(ind.indiv) & !df0$is_rel2 & df0$is_caucasian & !is.na(y))
length(sub)  # 231,942
table(y.sub <- y[sub])
#      0      1
# 226327   5615


NCORES <- 12
system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKBB_RA",
    ind_row = ind.indiv[sub],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = NCORES
  )
) # 2.8H


library(bigsnpr)
ukbb <- snp_attach("data/UKBB_RA.rds")
G <- ukbb$genotypes
file.size(G$backingfile) / 1024^3  # 218 GB
CHR <- as.integer(ukbb$map$chromosome)
POS <- ukbb$map$physical.pos

set.seed(1)
ind.train <- sort(sample(length(sub), 200e3))
ind.test <- setdiff(seq_along(sub), ind.train)

system.time(
  all_keep <- snp_grid_clumping(
    G, CHR, POS, lpS = lpval, ind.row = ind.train, infos.imp = info,
    grid.thr.imp = c(0.3, 0.6, 0.9, 0.95),
    grid.thr.r2 = c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.95),
    grid.base.size = c(50, 100, 200, 500),
    ncores = 8
  )
) # 3H -> swapping if too many cores

system.time(
  multi_PRS <- snp_grid_PRS(
    G, all_keep, betas = beta, lpS = lpval, ind.row = ind.train,
    n_thr_lpS = 50, backingfile = "data/UKBB_RA_scores", ncores = NCORES
  )
) # 1H

system.time(
  final_mod <- snp_grid_stacking(
    multi_PRS, y.sub[ind.train], ncores = NCORES, n.abort = 5
  )
) # 1.4H
mod <- final_mod$mod
plot(mod)
summary(mod)
# # A tibble: 3 x 6
#    alpha validation_loss intercept beta           nb_var message
#    <dbl>           <dbl>     <dbl> <list>          <int> <list>
# 1 0.0001           0.112     -4.05 <dbl [63,504]>  13451 <chr [10]>
# 2 0.01             0.112     -4.05 <dbl [63,504]>   2996 <chr [10]>
# 3 1                0.112     -4.05 <dbl [63,504]>   2352 <chr [10]>

new_beta <- final_mod$beta.G

length(ind <- which(new_beta != 0))  # 317,456
summary(new_beta)
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -1.346e-01  0.000e+00  0.000e+00  2.485e-05  0.000e+00  1.002e-01
summary(new_beta[which(sign(new_beta * beta) < 0)])
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -7.088e-03 -4.590e-05 -4.762e-06 -3.725e-05  1.515e-05  6.365e-03

pred <- final_mod$intercept +
  big_prodVec(G, new_beta[ind], ind.row = ind.test, ind.col = ind)

AUCBoot(pred, y.sub[ind.test])  # 61.3 [59.1-63.4]


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
#   size thr.r2 thr.imp num  thr.lp       auc
# 1  500    0.2     0.3  14 2.20872 0.5917841

ind.keep <- unlist(map(all_keep, std_prs$num))
# Verif on training set
AUC(
  snp_PRS(G, beta[ind.keep], ind.test = ind.train, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = std_prs$thr.lp),
  y.sub[ind.train]
) # 0.5917841
# Eval on test set
AUCBoot(
  snp_PRS(G, beta[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = std_prs$thr.lp),
  y.sub[ind.test]
) # 59.8 [57.7-61.8]
sum(lpval[ind.keep] > std_prs$thr.lp)  # 12,220

max_prs <- grid2 %>% arrange(desc(auc)) %>% slice(1:10) %>% print() %>% slice(1)
#    size thr.r2 thr.imp num   thr.lp       auc
# 1   200    0.5    0.95 102 1.123950 0.5984288
# 2  1000    0.5    0.95 104 1.123950 0.5982656
# 3   400    0.5    0.95 103 1.123950 0.5980790
# 4   200    0.5    0.95 102 1.004264 0.5980635
# 5   200    0.5    0.95 102 1.257900 0.5980216
# 6  1000    0.5    0.95 104 1.257900 0.5979347
# 7   100    0.5    0.95 101 1.123950 0.5978981
# 8   100    0.5    0.95 101 1.004264 0.5977654
# 9   400    0.5    0.95 103 1.257900 0.5977316
# 10 1000    0.5    0.95 104 1.004264 0.5976001

ind.keep <- unlist(map(all_keep, max_prs$num))
AUCBoot(
  snp_PRS(G, beta[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = max_prs$thr.lp),
  y.sub[ind.test]
) # 60.3 [58.3-62.4]
sum(lpval[ind.keep] > max_prs$thr.lp)  # 88,556

ggplot(grid2) +
  geom_point(aes(thr.lp, auc)) +
  facet_grid(thr.imp ~ thr.r2 + size) +
  scale_x_log10(limits = c(1, 8), breaks = c(1, 2, 5), minor_breaks = 1:10) +
  ylim(0.58, NA) +
  theme_bigstatsr(size.rel = 0.7) +
  labs(x = "-log10(p-value) threshold (log scale)", y = "AUC")
