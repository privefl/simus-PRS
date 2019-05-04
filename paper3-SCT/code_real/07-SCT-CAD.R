# download.file("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/NikpayM_26343387_GCST003116/cad.add.160614.website.txt",
#               destfile = "sumstats_CAD.txt")
library(bigreadr)
sumstats <- fread2("sumstats_CAD.txt",
                   select = c("chr", "bp_hg19", "noneffect_allele",
                              "effect_allele", "beta", "p_dgc"),
                   col.names = c("chr", "pos", "a0", "a1", "beta", "p"))
nrow(sumstats)  # 9,455,778

sumstats <- subset(sumstats, p < 0.1)
nrow(sumstats)  # 1,000,725
hist(sumstats$beta)
hist(sumstats$p)

info_snp_UKBB <- rbind_df(lapply(1:22, function(chr) {
  file <- paste0("data/ukb_imp_mfi/ukb_mfi_chr", chr, "_v3.txt")
  df <- fread2(file, select = c(3:5, 8), col.names = c("pos", "a0", "a1", "info"))
  cbind.data.frame(chr = chr, df)
}))
info_snp <- bigsnpr::snp_match(sumstats, info_snp_UKBB)
# 1,000,725 variants in summary statistics.
# 140,247 ambiguous SNPs have been removed.
# 761,308 variants have been matched; 0 were flipped and 577,018 were reversed.
info_snp <- bigsnpr::snp_match(sumstats, info_snp_UKBB, strand_flip = FALSE)
# 898,499 variants have been matched; 0 were flipped and 681,118 were reversed.
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


df_heart <- fread2(csv, select = c(paste0("6150-0.", 0:3),
                                   paste0("6150-1.", 0:3),
                                   paste0("6150-2.", 0:3)))
df_illness <- fread2(csv, select = c(paste0("20002-0.", 0:28),
                                     paste0("20002-1.", 0:28),
                                     paste0("20002-2.", 0:28)))
df_ICD10 <- fread2(csv, select = c(paste0("40001-", 0:2, ".0"),
                                   paste0("40002-0.", 0:13),
                                   paste0("40002-1.", 0:13),
                                   paste0("40002-2.", 0:13),
                                   paste0("41202-0.", 0:379),
                                   paste0("41204-0.", 0:434)))
ind_CAD <- sort(unique(unlist(c(
  lapply(df_heart,   function(x) which(x == 1)),
  lapply(df_illness, function(x) which(x == 1075)),
  lapply(df_ICD10,   function(x) which(substr(x, 1, 3) %in% paste0("I", 21:23))),
  lapply(df_ICD10,   function(x) which(x == "I252"))
))))
ind_heart <- sort(unique(unlist(c(
  lapply(df_heart,   function(x) which(x %in% 1:3)),
  lapply(df_illness, function(x) which(x %in% 1074:1080)),
  lapply(df_ICD10,   function(x) which(substr(x, 1, 1) == "I"))
))))
y <- rep(0, nrow(df0)); y[ind_heart] <- NA; y[ind_CAD] <- 1

sub <- which(!is.na(ind.indiv) & !df0$is_rel2 & df0$is_caucasian & !is.na(y) &
               !is.na(df0$sex))
length(sub)  # 238,190
table(y.sub <- y[sub])
#      0      1
# 225927  12263

NCORES <- 12
system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKBB_CAD",
    ind_row = ind.indiv[sub],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = NCORES
  )
) # 2.7H


library(bigsnpr)
ukbb <- snp_attach("data/UKBB_CAD.rds")
G <- ukbb$genotypes
file.size(G$backingfile) / 1024^3  # 199 GB
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
) # 1.8H
plot(lengths(all_keep[[1]]), pch = 20)

system.time(
  multi_PRS <- snp_grid_PRS(
    G, all_keep, betas = beta, lpS = lpval, ind.row = ind.train,
    n_thr_lpS = 50, backingfile = "data/UKBB_CAD_scores", ncores = NCORES
  )
) # 1H

system.time(
  final_mod <- snp_grid_stacking(
    multi_PRS, y.sub[ind.train], ncores = NCORES, n.abort = 3
  )
) # 3.1H
mod <- final_mod$mod
plot(mod)
summary(mod)
# # A tibble: 3 x 6
#    alpha validation_loss intercept beta           nb_var message
#    <dbl>           <dbl>     <dbl> <list>          <int> <list>
# 1 0.0001           0.197     -3.45 <dbl [65,940]>  26518 <chr [10]>
# 2 0.01             0.197     -3.42 <dbl [65,940]>   5816 <chr [10]>
# 3 1                0.197     -3.40 <dbl [65,940]>   3886 <chr [10]>

# Save
save(all_keep, multi_PRS, final_mod, file = "data/res_CAD.RData")

new_beta <- final_mod$beta.G

length(ind <- which(new_beta != 0))  # 315,165
summary(new_beta)
#      Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -0.1325315  0.0000000  0.0000000  0.0000162  0.0000000  0.3180476
summary(new_beta[which(sign(new_beta * beta) < 0)])
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -3.581e-02 -7.685e-05 -8.000e-07 -9.660e-06  6.376e-05  2.581e-02

pred <- final_mod$intercept +
  big_prodVec(G, new_beta[ind], ind.row = ind.test, ind.col = ind)

AUCBoot(pred, y.sub[ind.test])  # 63.9 [62.7-65.1]


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
# 1  500    0.2     0.3  14 3.371456 0.6036057

ind.keep <- unlist(map(all_keep, std_prs$num))
# Verif on training set
AUC(
  snp_PRS(G, beta[ind.keep], ind.test = ind.train, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = std_prs$thr.lp),
  y.sub[ind.train]
) # 0.6036057
# Eval on test set
AUCBoot(
  snp_PRS(G, beta[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = std_prs$thr.lp),
  y.sub[ind.test]
) # 59.9 [58.6-61.2]
sum(lpval[ind.keep] > std_prs$thr.lp)  # 1182

max_prs <- grid2 %>% arrange(desc(auc)) %>% slice(1:10) %>% print() %>% slice(1)
#    size thr.r2 thr.imp num   thr.lp       auc
# 1   526   0.95    0.95 112 1.453368 0.6190804
# 2   210   0.95    0.95 111 1.453368 0.6190598
# 3   526   0.95    0.95 112 1.323643 0.6188899
# 4   105   0.95    0.95 110 1.453368 0.6188860
# 5   210   0.95    0.95 111 1.323643 0.6188684
# 6   105   0.95    0.95 110 1.323643 0.6187670
# 7   250   0.80    0.95 107 1.453368 0.6185269
# 8   625   0.80    0.95 108 1.453368 0.6184453
# 9   625   0.80    0.95 108 1.595806 0.6184363
# 10  526   0.95    0.95 112 1.205497 0.6184120

ind.keep <- unlist(map(all_keep, max_prs$num))
AUCBoot(
  snp_PRS(G, beta[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = max_prs$thr.lp),
  y.sub[ind.test]
) # 61.1 [59.9-62.4]
sum(lpval[ind.keep] > max_prs$thr.lp)  # 87,577

ggplot(grid2) +
  geom_point(aes(thr.lp, auc)) +
  facet_grid(thr.imp ~ thr.r2 + size) +
  scale_x_log10(limits = c(1, 7), breaks = c(1, 2, 5), minor_breaks = 1:10) +
  ylim(0.58, NA) +
  theme_bigstatsr(size.rel = 0.7) +
  labs(x = "-log10(p-value) threshold (log scale)", y = "AUC")


## The 'effect' of sex
system.time(
  final_mod2 <- snp_grid_stacking(
    multi_PRS, y.sub[ind.train], ncores = NCORES, n.abort = 3,
    covar.train = as.matrix(df0$sex[sub[ind.train]]), pf.covar = 0
  )
) # 2.7H
mod2 <- final_mod2$mod
plot(mod2)
summary(mod2)
# # A tibble: 3 x 6
#    alpha validation_loss intercept beta           nb_var message
#    <dbl>           <dbl>     <dbl> <list>          <int> <list>
# 1 0.0001           0.184     -4.43 <dbl [65,941]>  29110 <chr [10]>
# 2 0.01             0.184     -4.38 <dbl [65,941]>   6660 <chr [10]>
# 3 1                0.184     -4.35 <dbl [65,941]>   6350 <chr [10]>

new_beta2 <- final_mod2$beta.G

length(ind2 <- which(new_beta2 != 0))  # 316,140
pred2 <- final_mod2$intercept +
  big_prodVec(G, new_beta2[ind2], ind.row = ind.test, ind.col = ind2) +
  final_mod2$beta.covar * df0$sex[sub[ind.test]]

AUCBoot(pred2, y.sub[ind.test])  # 74.4 [73.4-75.5]

ggplot(data.frame(
  pred = 1 / (1 + exp(-c(pred, pred2))),
  method = rep(c("Without variable 'sex'", "Using variable 'sex'"), each = length(pred)),
  pheno = factor(y.sub[ind.test], levels = 0:1, labels = c("Control", "Case")),
  sex = factor(df0$sex[sub[ind.test]], levels = 0:1, labels = c("Woman", "Man"))
)) +
  theme_bigstatsr() +
  geom_density(aes(pred, fill = pheno), alpha = 0.3) +
  facet_grid(sex ~ method) +
  scale_x_continuous(limits = c(0, 0.2)) +
  labs(x = "Probability", y = "Density", fill = "Phenotype") +
  theme(legend.position = c(0.85, 0.85))


is.men <- (df0$sex[sub[ind.train]] == 1)
is.men2 <- (df0$sex[sub[ind.test]] == 1)

AUCBoot(pred[ is.men2], y.sub[ind.test[ is.men2]])   # 64.9 [63.5-66.3]
AUCBoot(pred[!is.men2], y.sub[ind.test[!is.men2]])   # 62.5 [59.7-65.2]

AUCBoot(pred2[ is.men2], y.sub[ind.test[ is.men2]])  # 64.9 [63.5-66.3]
AUCBoot(pred2[!is.men2], y.sub[ind.test[!is.men2]])  # 62.5 [59.8-65.2]

plot(pred, pred2, pch = 20, col = df0$sex[sub[ind.test]] + 1); abline(0, 1, col = "blue")


system.time(
  final_mod_men <- snp_grid_stacking(
    multi_PRS, y.sub[ind.train[is.men]], ind.train = which(is.men),
    ncores = NCORES, n.abort = 3
  )
) # 61 min

system.time(
  final_mod_women <- snp_grid_stacking(
    multi_PRS, y.sub[ind.train[!is.men]], ind.train = which(!is.men),
    ncores = NCORES, n.abort = 3
  )
) # 19 min

new_beta3 <- final_mod_men$beta.G
length(ind3 <- which(new_beta3 != 0))  # 271,946
pred3 <- final_mod_men$intercept +
  big_prodVec(G, new_beta3[ind3], ind.row = ind.test, ind.col = ind3)

AUCBoot(pred3[ is.men2], y.sub[ind.test[ is.men2]])  # 64.7 [63.3-66.2]
AUCBoot(pred3[!is.men2], y.sub[ind.test[!is.men2]])  # 62.1 [59.3-64.9]

new_beta4 <- final_mod_women$beta.G
length(ind4 <- which(new_beta4 != 0))  # 288,163
pred4 <- final_mod_women$intercept +
  big_prodVec(G, new_beta4[ind4], ind.row = ind.test, ind.col = ind4)

AUCBoot(pred4[ is.men2], y.sub[ind.test[ is.men2]])  # 63.3 [61.9-64.8]
AUCBoot(pred4[!is.men2], y.sub[ind.test[!is.men2]])  # 61.4 [58.7-64.2]
