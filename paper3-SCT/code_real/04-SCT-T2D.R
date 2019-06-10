# Download from http://diagram-consortium.org/downloads.html
# DIAGRAM 1000G GWAS meta-analysis Stage 1 Summary statistics
# Published in Scott et al (2017)
# unzip("METAANALYSIS_DIAGRAM_SE1.zip")
library(bigreadr)
sumstats <- fread2("METAANALYSIS_DIAGRAM_SE1.txt", select = c(1:4, 6))
sumstats <- tidyr::separate(sumstats, "Chr:Position", c("chr", "pos"), convert = TRUE)
names(sumstats) <- c("chr", "pos", "a1", "a0", "beta", "p")
nrow(sumstats)  # 12,056,346
hist(sumstats$beta)
hist(sumstats$p)
sumstats <- subset(sumstats, p < 0.1)

info_snp_UKBB <- rbind_df(lapply(1:22, function(chr) {
  file <- paste0("data/ukb_imp_mfi/ukb_mfi_chr", chr, "_v3.txt")
  df <- fread2(file, select = c(3:5, 8), col.names = c("pos", "a0", "a1", "info"))
  cbind.data.frame(chr = chr, df)
}))
info_snp <- bigsnpr::snp_match(sumstats, info_snp_UKBB)
# 1,408,672 variants in summary statistics.
# 215,821 ambiguous SNPs have been removed.
# 1,145,260 variants have been matched; 38 were flipped and 499,125 were reversed.
info_snp <- bigsnpr::snp_match(sumstats, info_snp_UKBB, strand_flip = FALSE)
# 1,350,844 variants have been matched; 0 were flipped and 602,001 were reversed.
info_snp <- subset(na.omit(info_snp), info > 0.3)

list_snp_id <- with(info_snp, split(paste(chr, pos, a0, a1, sep = "_"),
                                    factor(chr, levels = 1:22)))
beta <- info_snp$beta
lpval <- -log10(info_snp$p)
info <- info_snp$info

# subset samples
library(bigreadr)
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
ind_diabetes <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x %in% 1220:1223)),
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) %in% paste0("E", 10:14)))
))))
ind_TD1 <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x == 1222)),
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "E10"))
))))
ind_TD2 <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x == 1223)),
  lapply(df_ICD10, function(x) which(substr(x, 1, 3) == "E11"))
))))
y <- rep(0, nrow(df_illness))
y[ind_diabetes] <- NA
y[ind_TD2] <- 1
y[ind_TD1] <- NA

sub <- which(!is.na(ind.indiv) & !df0$is_rel2 & df0$is_caucasian & !is.na(y))
length(sub)  # 328,723
table(y.sub <- y[sub])
#      0      1
# 314547  14176

NCORES <- 12
system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKBB_T2D",
    ind_row = ind.indiv[sub],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = NCORES
  )
) # 4H

library(bigsnpr)
ukbb <- snp_attach("data/UKBB_T2D.rds")
G <- ukbb$genotypes
file.size(G$backingfile) / 1024^3  # 410 GB
CHR <- as.integer(ukbb$map$chromosome)
POS <- ukbb$map$physical.pos

set.seed(1); ind.train <- sort(sample(length(sub), 250e3))
# set.seed(2); ind.train <- c(sample(which(y.sub == 0), 2000), sample(which(y.sub == 1), 500))
ind.test <- setdiff(seq_along(sub), ind.train)

system.time(
  all_keep <- snp_grid_clumping(
    G, CHR, POS, lpS = lpval, ind.row = ind.train, infos.imp = info,
    grid.thr.imp = c(0.3, 0.6, 0.9, 0.95),
    grid.thr.r2 = c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.95),
    grid.base.size = c(50, 100, 200, 500),
    ncores = 2  # use less cores because of swapping if not enough memory
  )
) # 16H -> 3.7H

system.time(
  multi_PRS <- snp_grid_PRS(
    G, all_keep, betas = beta, lpS = lpval, ind.row = ind.train,
    n_thr_lpS = 50, backingfile = "data/UKBB_T2D_scores", ncores = NCORES
  )
) # 2.7H

system.time(
  final_mod <- snp_grid_stacking(
    multi_PRS, y.sub[ind.train], ncores = NCORES
  )
) # 6.2H
mod <- final_mod$mod
plot(mod)
summary(mod)
# # A tibble: 3 x 6
#    alpha validation_loss intercept beta           nb_var message
#    <dbl>           <dbl>     <dbl> <list>          <int> <list>
# 1 0.0001           0.172     -2.21 <dbl [65,520]>  23548 <chr [10]>
# 2 0.01             0.172     -2.20 <dbl [65,520]>   6778 <chr [10]>
# 3 1                0.172     -2.21 <dbl [65,520]>   4924 <chr [10]>

new_beta <- final_mod$beta.G

length(ind <- which(new_beta != 0))  # 548,343
summary(new_beta)
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -0.3556050  0.0000000  0.0000000  0.0000178  0.0000000  0.2404328
summary(new_beta[which(sign(new_beta * beta) < 0)])
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -4.446e-02 -3.097e-04 -6.010e-06 -8.926e-05  1.902e-04  3.728e-02

pred <- final_mod$intercept +
  big_prodVec(G, new_beta[ind], ind.row = ind.test, ind.col = ind)

AUCBoot(pred, y.sub[ind.test])  # 63.8 [62.9-64.7] / 61.0 [60.6-61.5]

# Save
save(all_keep, multi_PRS, final_mod, file = "data/res_T2D.RData")

library(tidyverse)

ind2 <- sample(ind, 10e3)
ggplot(data.frame(y = new_beta, x = beta)[ind, ]) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  geom_abline(slope = 0, intercept = 0, color = "blue") +
  geom_point(aes(x, y), size = 0.8) +
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
# 1  500    0.2     0.3  14 5.580693 0.5989295

ind.keep <- unlist(map(all_keep, std_prs$num))
# Verif on training set
AUC(
  snp_PRS(G, beta[ind.keep], ind.test = ind.train, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = std_prs$thr.lp),
  y.sub[ind.train]
) # 0.5989295
# Eval on test set
AUCBoot(
  snp_PRS(G, beta[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = std_prs$thr.lp),
  y.sub[ind.test]
) # 59.1 [58.1-60.1] / 59.8 [59.3-60.3]
sum(lpval[ind.keep] > std_prs$thr.lp)  # 177

max_prs <- grid2 %>% arrange(desc(auc)) %>% slice(1:10) %>% print() %>% slice(1)
#    size thr.r2 thr.imp num   thr.lp       auc
# 1   625    0.8    0.95 108 1.956627 0.6141710
# 2   250    0.8    0.95 107 1.956627 0.6140094
# 3   125    0.8    0.95 106 1.956627 0.6139443
# 4   625    0.8    0.95 108 1.778804 0.6129134
# 5   125    0.8    0.95 106 1.778804 0.6128044
# 6   250    0.8    0.95 107 1.778804 0.6127813
# 7   625    0.8    0.95 108 2.152227 0.6127131
# 8   625    0.8    0.95 108 2.604043 0.6126175
# 9   250    0.8    0.95 107 2.152227 0.6125586
# 10  250    0.8    0.95 107 2.604043 0.6125395

ind.keep <- unlist(map(all_keep, max_prs$num))
AUCBoot(
  snp_PRS(G, beta[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = max_prs$thr.lp),
  y.sub[ind.test]
) # 60.7 [59.8-61.7] / 60.2 [59.7-60.7]
sum(lpval[ind.keep] > max_prs$thr.lp)  # 33,235

ggplot(grid2) +
  geom_point(aes(thr.lp, auc)) +
  facet_grid(thr.imp ~ thr.r2 + size) +
  scale_x_log10(limits = c(2, 10), breaks = c(1, 2, 5), minor_breaks = 1:10) +
  ylim(0.58, NA) +
  theme_bigstatsr(size.rel = 0.7) +
  labs(x = "-log10(p-value) threshold (log scale)", y = "AUC")


# lassosum
library(lassosum)
ukbb$map <- dplyr::mutate(ukbb$map, genetic.dist = 0, rsid = NULL,
                          chromosome = as.integer(chromosome))
ukbb$fam <- snp_fake(nrow(G), 1)$fam
bed <- snp_writeBed(ukbb, bedfile = "data/UKBB_T2D_lassosum.bed",
                    ind.row = ind.train)
library(doParallel)
registerDoParallel(cl <- makeCluster(NCORES))
system.time(
  out <- lassosum.pipeline(
    cor = p2cor(p = 10^-lpval, n = 26676 + 132532, sign = beta),
    snp = ukbb$map$marker.ID,
    A1 = ukbb$map$allele1,
    test.bfile = "data/UKBB_T2D_lassosum",
    LDblocks = "EUR.hg19",
    cluster = cl,
    sample = 20e3,
    exclude.ambiguous = FALSE
  )
) # 9H
stopCluster(cl)

v <- validate(out, pheno = y.sub[ind.train], validate.function = AUC)
length(ind <- which(v$best.beta != 0))  # 256,353
pred_lassosum <- big_prodVec(G, v$best.beta[ind], ind.row = ind.test, ind.col = ind)
AUCBoot(pred_lassosum, y.sub[ind.test])  # 63.2 [62.3-64.1] / 63.6 [63.1-64.1]
