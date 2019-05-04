# download.file("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/WrayNR_29700475_GCST005839/MDD2018_ex23andMe.gz",
#               destfile = "sumstats_MDD.txt.gz")
# R.utils::gunzip("sumstats_MDD.txt.gz")
library(bigreadr)
sumstats <- fread2("sumstats_MDD.txt", fill = TRUE,
                   select = c("CHR", "BP", "A1", "A2", "OR", "P"),
                   col.names = c("chr", "pos", "a1", "a0", "or", "p"))
nrow(sumstats)  # 13,554,550
library(dplyr)
sumstats <- sumstats %>%
  filter(p < 0.1) %>%
  mutate(beta = log(or), or = NULL, chr = as.integer(chr))
nrow(sumstats)  # 1,610,995
hist(sumstats$beta)
hist(sumstats$p)


info_snp_UKBB <- rbind_df(lapply(1:22, function(chr) {
  file <- paste0("data/ukb_imp_mfi/ukb_mfi_chr", chr, "_v3.txt")
  df <- fread2(file, select = c(3:5, 8), col.names = c("pos", "a0", "a1", "info"))
  cbind.data.frame(chr = chr, df)
}))
info_snp <- bigsnpr::snp_match(sumstats, info_snp_UKBB)
# 1,610,995 variants in summary statistics.
# 224,231 ambiguous SNPs have been removed.
# 1,201,894 variants have been matched; 102 were flipped and 521,969 were reversed.
info_snp <- bigsnpr::snp_match(sumstats, info_snp_UKBB, strand_flip = FALSE)
# 1,420,221 variants have been matched; 0 were flipped and 629,580 were reversed.
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
df0 <- fread2(csv, select = c("eid", "22000-0.0", "22006-0.0"),
              col.names = c("eid", "batch", "is_caucasian"))
sum(df0$batch < 0, na.rm = TRUE)  # 49,946
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
ind_MDD <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x == 1286)),
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) %in% c("F32", "F33")))
))))
ind_psy <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x %in% 1286:1291)),
  lapply(df_ICD10, function(x) which(substr(x, 1, 1) == "F"))
))))
y <- rep(0, nrow(df0)); y[ind_psy] <- NA; y[ind_MDD] <- 1

sub <- which(!is.na(ind.indiv) & !df0$is_rel2 & df0$is_caucasian & !is.na(y) &
               df0$batch > 0)
length(sub)  # 277,604
table(y.sub <- y[sub])
#      0      1
# 255317  22287


NCORES <- 12
system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKBB_MDD",
    ind_row = ind.indiv[sub],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = NCORES
  )
) # 4H


library(bigsnpr)
ukbb <- snp_attach("data/UKBB_MDD.rds")
G <- ukbb$genotypes
file.size(G$backingfile) / 1024^3  # 365 GB
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
    ncores = 6
  )
) # 13.4H

system.time(
  multi_PRS <- snp_grid_PRS(
    G, all_keep, betas = beta, lpS = lpval, ind.row = ind.train,
    n_thr_lpS = 50, backingfile = "data/UKBB_MDD_scores", ncores = NCORES
  )
) # 4.5H

system.time(
  final_mod <- snp_grid_stacking(
    multi_PRS, y.sub[ind.train], ncores = NCORES, n.abort = 3
  )
) # 3.6H
mod <- final_mod$mod
plot(mod)
summary(mod)

# Save
save(all_keep, multi_PRS, final_mod, file = "data/res_MDD.RData")

new_beta <- final_mod$beta.G

length(ind <- which(new_beta != 0))  # 524,099
summary(new_beta)
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -2.681e-01  0.000e+00  0.000e+00  4.934e-05  0.000e+00  1.898e-01
summary(new_beta[which(sign(new_beta * beta) < 0)])
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -7.465e-02 -2.126e-04 -1.156e-05 -1.477e-04  4.437e-05  5.539e-02

pred <- final_mod$intercept +
  big_prodVec(G, new_beta[ind], ind.row = ind.test, ind.col = ind)

AUCBoot(pred, y.sub[ind.test])  # 59.5 [58.2-60.7]


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
#   size thr.r2 thr.imp num   thr.lp       auc
# 1  500    0.2     0.3  14 1.272754 0.5456847

ind.keep <- unlist(map(all_keep, std_prs$num))
# Verif on training set
AUC(
  snp_PRS(G, beta[ind.keep], ind.test = ind.train, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = std_prs$thr.lp),
  y.sub[ind.train]
) # 0.5456847
# Eval on test set
AUCBoot(
  snp_PRS(G, beta[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = std_prs$thr.lp),
  y.sub[ind.test]
) # 55.7 [54.4-56.9]
sum(lpval[ind.keep] > std_prs$thr.lp)  # 165,584

max_prs <- grid2 %>% arrange(desc(auc)) %>% slice(1:10) %>% print() %>% slice(1)
#    size thr.r2 thr.imp num    thr.lp       auc
# 1   625    0.8    0.95 108 0.9999434 0.5918091
# 2   250    0.8    0.95 107 0.9999434 0.5913374
# 3   625    0.8    0.95 108 1.0493715 0.5912436
# 4   250    0.8    0.95 107 1.0493715 0.5907580
# 5   625    0.8    0.95 108 1.1012427 0.5902952
# 6   125    0.8    0.95 106 0.9999434 0.5900350
# 7  1000    0.5    0.95 104 0.9999434 0.5898991
# 8   250    0.8    0.95 107 1.1012427 0.5898025
# 9  1000    0.5    0.95 104 1.0493715 0.5896351
# 10  400    0.5    0.95 103 0.9999434 0.5896332

ind.keep <- unlist(map(all_keep, max_prs$num))
AUCBoot(
  snp_PRS(G, beta[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = max_prs$thr.lp),
  y.sub[ind.test]
) # 59.2 [58.0-60.4]
sum(lpval[ind.keep] > max_prs$thr.lp)  # 222,912

ggplot(grid2) +
  geom_point(aes(thr.lp, auc)) +
  facet_grid(thr.imp ~ thr.r2 + size) +
  scale_x_log10(limits = c(1, 4), breaks = c(1, 2, 4), minor_breaks = 1:10) +
  ylim(0.52, NA) +
  theme_bigstatsr(size.rel = 0.7) +
  labs(x = "-log10(p-value) threshold (log scale)", y = "AUC")


ind3 <- intersect(which(abs(beta) > 8), ind)
gwas <- big_univLogReg(G, y.sub[ind.train], ind.train, ind3)
cbind(beta_GWAS = beta[ind3], pval_GWAS = 10^(-lpval[ind3]),
      af = big_scale()(G, ind.row = ind.train, ind.col = ind3)$center / 2,
      gwas, pval = predict(gwas, log10 = FALSE))
#    beta_GWAS pval_GWAS         af         estim      std.err niter      score       pval
# 1  -8.750700  0.088670 0.00031440   -0.44859668 4.433835e-01     5 -1.0117577 0.31165394
# 2  -8.613200  0.034370 0.00092246   -0.23701529 2.218916e-01     5 -1.0681580 0.28544923
# 3  -8.111728  0.059710 0.00000778  -14.32590356 1.473188e+01     6 -0.9724424 0.33083049
# 4  -8.111728  0.014470 0.00000900    2.06582809 1.441504e+00     6  1.4331062 0.15182749
# 5  -8.383900  0.048800 0.00000388 -865.96733183 4.497506e+03    NA -0.1925439 0.84731619
# 6  -8.207900  0.014250 0.00054366    0.06429999 2.435848e-01     4  0.2639737 0.79180016
# 7   9.210340  0.084280 0.00438248    0.09832244 7.732358e-02     4  1.2715711 0.20352554
# 8  -8.517193  0.080990 0.00004366   -1.52485520 2.165943e+00     5 -0.7040143 0.48142385
# 9   8.765100  0.043380 0.00146470   -0.16032271 2.035233e-01     4 -0.7877363 0.43085096
# 10 -8.628700  0.041030 0.00395916   -0.08737906 9.762673e-02     4 -0.8950321 0.37076990
# 11  9.054000  0.092760 0.00001072    1.43283142 4.022248e+00     5  0.3562265 0.72167091
# 12 -8.866900  0.021650 0.00296426   -0.11194072 1.043994e-01     4 -1.0722351 0.28361445
# 13 -9.210340  0.019190 0.00000790 -813.83196680 4.896593e+03    NA -0.1662037 0.86799663
# 14  8.262000  0.026560 0.00064012    0.45151789 2.607551e-01     5  1.7315781 0.08334871
# 15  8.136900  0.053150 0.00072848    0.06021477 2.303720e-01     4  0.2613806 0.79379903
# 16 -8.111728  0.032100 0.00000062 -765.49217238 3.848189e+03    NA -0.1989227 0.84232320
# 17  8.121100  0.046650 0.00307300   -0.17136069 1.075134e-01     4 -1.5938547 0.11096857
# 18 -8.087900  0.071530 0.00004466   -0.26093054 1.181459e+00     5 -0.2208545 0.82520576
# 19  8.517193  0.050130 0.00268724    0.01236054 1.083695e-01     4  0.1140592 0.90919086
# 20 -8.517193  0.028230 0.00003914    0.37660199 9.616836e-01     5  0.3916069 0.69534866
# 21  8.257000  0.049410 0.00015768    0.57221809 3.588226e-01     5  1.5947105 0.11077699
# 22 -8.414900  0.018160 0.00208406    0.03599946 1.177297e-01     4  0.3057807 0.75977161
# 23  8.517193  0.078270 0.00118494   -0.07220227 1.809734e-01     4 -0.3989661 0.68991817
# 24 -8.517193  0.006523 0.00435750    0.01928752 8.799206e-02     4  0.2191962 0.82649723
# 25 -9.210340  0.076300 0.00008886    0.16226111 5.696757e-01     4  0.2848307 0.77577389
# 26 -8.456700  0.058450 0.00566904    0.03675442 7.077464e-02     4  0.5193162 0.60354026
# 27 -8.517193  0.064840 0.00035564   -0.29049609 3.404072e-01     5 -0.8533781 0.39344964
# 28  8.111728  0.005863 0.99995000    4.45390334 2.963645e+00     7  1.5028463 0.13287868
# 29  8.895500  0.084360 0.00000470    3.38686750 3.810110e+00     6  0.8889160 0.37404823
# 30 -9.210340  0.049990 0.00002176   -7.34929367 6.696494e+00     7 -1.0974839 0.27242994
# 31 -8.558400  0.079470 0.00026840    0.45306138 3.412574e-01     5  1.3276237 0.18430245
# 32 -8.517193  0.085690 0.00213094   -0.05227573 1.645161e-01     4 -0.3177546 0.75067112
# 33 -8.111728  0.099890 0.00016460   -0.39249888 5.321898e-01     5 -0.7375167 0.46080818
# 34 -8.111728  0.007289 0.00019558    0.30148880 3.631188e-01     5  0.8302759 0.40638282
# 35 -8.027300  0.096170 0.00367082   -0.03978635 1.008721e-01     4 -0.3944239 0.69326808
# 36 -8.036600  0.081220 0.00010924   -0.57166721 6.691692e-01     5 -0.8542940 0.39294211
# 37  8.654100  0.003695 0.00280600    0.22256560 9.660467e-02     5  2.3038804 0.02122936
# 38  9.054400  0.035150 0.00005486    0.27778991 7.690223e-01     5  0.3612248 0.71793140
