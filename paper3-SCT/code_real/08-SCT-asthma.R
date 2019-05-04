download.file("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/DemenaisF_29273806_GCST006862/TAGC_meta-analyses_results_for_asthma_risk.zip",
              destfile = "sumstats_asthma.zip")
unzip("sumstats_asthma.zip")
library(bigreadr)
sumstats <- fread2(
  "TAGC meta-analyses results for asthma risk/TAGC_Multiancestry_and_European-Ancestry_Meta-analyses_Results.tsv",
  select = c(1, 3:5, 18, 20),
  col.names = c("chr", "pos", "a0", "a1", "beta", "p")
)
nrow(sumstats)  # 2,001,280

sumstats <- subset(sumstats, p < 0.5)
nrow(sumstats)  # 1,013,806
hist(sumstats$beta)
hist(sumstats$p)

info_snp_UKBB <- rbind_df(lapply(1:22, function(chr) {
  file <- paste0("data/ukb_imp_mfi/ukb_mfi_chr", chr, "_v3.txt")
  df <- fread2(file, select = c(3:5, 8), col.names = c("pos", "a0", "a1", "info"))
  cbind.data.frame(chr = chr, df)
}))
info_snp <- bigsnpr::snp_match(sumstats, info_snp_UKBB)
# 1,013,806 variants in summary statistics.
# 0 ambiguous SNPs have been removed.
# 1,013,119 variants have been matched; 0 were flipped and 1 were reversed.
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

ind_asthma <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x == 1111)),
  lapply(df_ICD10,   function(x) which(substr(x, 1, 3) == "J45"))
))))
ind_respi <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x %in% 1111:1125)),
  lapply(df_ICD10,   function(x) which(substr(x, 1, 1) == "J"))
))))
y <- rep(0, nrow(df0)); y[ind_respi] <- NA; y[ind_asthma] <- 1

sub <- which(!is.na(ind.indiv) & !df0$is_rel2 & df0$is_caucasian & !is.na(y))
length(sub)  # 305,772
table(y.sub <- y[sub])
#      0      1
# 261985  43787


NCORES <- 12
system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKBB_asthma",
    ind_row = ind.indiv[sub],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = NCORES
  )
) # 3.8H


library(bigsnpr)
ukbb <- snp_attach("data/UKBB_asthma.rds")
G <- ukbb$genotypes
file.size(G$backingfile) / 1024^3  # 288 GB
CHR <- as.integer(ukbb$map$chromosome)
POS <- ukbb$map$physical.pos

set.seed(1)
ind.train <- sort(sample(length(sub), 250e3))
ind.test <- setdiff(seq_along(sub), ind.train)

system.time(
  all_keep <- snp_grid_clumping(
    G, CHR, POS, lpS = lpval, ind.row = ind.train, infos.imp = info,
    grid.thr.imp = c(0.3, 0.6, 0.9, 0.95),
    grid.thr.r2 = c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.95),
    grid.base.size = c(50, 100, 200, 500),
    ncores = NCORES
  )
) # 4.3H

system.time(
  multi_PRS <- snp_grid_PRS(
    G, all_keep, betas = beta, lpS = lpval, ind.row = ind.train,
    n_thr_lpS = 50, backingfile = "data/UKBB_asthma_scores", ncores = NCORES
  )
) # 1.7H

system.time(
  final_mod <- snp_grid_stacking(
    multi_PRS, y.sub[ind.train], ncores = NCORES, n.abort = 3
  )
) # 12H
mod <- final_mod$mod
plot(mod)
summary(mod)
# # A tibble: 3 x 6
#    alpha validation_loss intercept beta           nb_var message
#    <dbl>           <dbl>     <dbl> <list>          <int> <list>
# 1 0.0001           0.403     -1.09 <dbl [89,012]>  61609 <chr [10]>
# 2 0.01             0.403     -1.16 <dbl [89,012]>  20564 <chr [10]>
# 3 1                0.403     -1.10 <dbl [89,012]>  17094 <chr [10]>

new_beta <- final_mod$beta.G

length(ind <- which(new_beta != 0))  # 446,120
summary(new_beta)
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -8.224e-02  0.000e+00  0.000e+00  4.550e-06  0.000e+00  1.163e-01
summary(new_beta[which(sign(new_beta * beta) < 0)])
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -2.551e-02 -8.117e-05 -4.700e-07 -2.134e-05  5.690e-05  4.051e-02

pred <- final_mod$intercept +
  big_prodVec(G, new_beta[ind], ind.row = ind.test, ind.col = ind)

AUCBoot(pred, y.sub[ind.test])  # 60.7 [60.0-61.3]


library(tidyverse)

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
# 1  500    0.2     0.3  14 2.30809 0.5667112

ind.keep <- unlist(map(all_keep, std_prs$num))
# Verif on training set
AUC(
  snp_PRS(G, beta[ind.keep], ind.test = ind.train, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = std_prs$thr.lp),
  y.sub[ind.train]
) # 0.5667112
# Eval on test set
AUCBoot(
  snp_PRS(G, beta[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = std_prs$thr.lp),
  y.sub[ind.test]
) # 56.8 [56.2-57.5]
sum(lpval[ind.keep] > std_prs$thr.lp)  # 3034

max_prs <- grid2 %>% arrange(desc(auc)) %>% slice(1:10) %>% print() %>% slice(1)
#    size thr.r2 thr.imp num   thr.lp       auc
# 1  2500    0.2    0.90  72 3.667058 0.5699586
# 2  2500    0.2    0.30  16 3.667058 0.5699415
# 3  2500    0.2    0.60  44 3.667058 0.5699415
# 4  2500    0.2    0.95 100 3.667058 0.5698932
# 5  1000    0.5    0.90  76 1.748291 0.5698451
# 6  1000    0.5    0.95 104 1.748291 0.5697478
# 7  5000    0.1    0.95  96 4.022821 0.5696572
# 8  2500    0.2    0.95 100 4.022821 0.5696496
# 9  1000    0.5    0.60  48 1.748291 0.5696363
# 10 1000    0.5    0.30  20 1.748291 0.5696320

ind.keep <- unlist(map(all_keep, max_prs$num))
AUCBoot(
  snp_PRS(G, beta[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = max_prs$thr.lp),
  y.sub[ind.test]
) # 57.3 [56.7-58.0]
sum(lpval[ind.keep] > max_prs$thr.lp)  # 360

ggplot(grid2) +
  geom_point(aes(thr.lp, auc)) +
  facet_grid(thr.imp ~ thr.r2 + size) +
  scale_x_log10(limits = c(0.3, 10), breaks = c(1, 2, 5), minor_breaks = 1:10) +
  ylim(0.54, NA) +
  theme_bigstatsr(size.rel = 0.7) +
  labs(x = "-log10(p-value) threshold (log scale)", y = "AUC")
