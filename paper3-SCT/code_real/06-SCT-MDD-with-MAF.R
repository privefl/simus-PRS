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
  df <- fread2(file, select = c(3:6, 8), col.names = c("pos", "a0", "a1", "maf", "info"))
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
maf <- info_snp$maf
hist(log10(maf))


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


NCORES <- 16
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

set.seed(1); ind.train <- sort(sample(length(sub), 250e3))
ind.test <- setdiff(seq_along(sub), ind.train)

MAF <- 10^(-(6:1))
groups <- lapply(MAF, function(maf.thr) {
  which(maf > maf.thr)
})

system.time(
  all_keep <- snp_grid_clumping(
    G, CHR, POS, lpS = lpval, ind.row = ind.train, infos.imp = info,
    grid.thr.imp = c(0.3, 0.6, 0.9, 0.95),
    grid.thr.r2 = c(0.5, 0.8, 0.95),
    grid.base.size = c(100, 500),
    groups = groups,
    ncores = 4  # use less cores to prevent swapping if not enough RAM
  )
) # 13.4H -> 4H

system.time(
  multi_PRS <- snp_grid_PRS(
    G, all_keep, betas = beta, lpS = lpval, ind.row = ind.train,
    n_thr_lpS = 50, backingfile = "data/UKBB_MDD_scores2", ncores = NCORES
  )
) # 2H

system.time(
  final_mod <- snp_grid_stacking(
    multi_PRS, y.sub[ind.train], ncores = NCORES
  )
) # 3.6H
mod <- final_mod$mod
plot(mod)
summary(mod)

new_beta <- final_mod$beta.G

length(ind <- which(new_beta != 0))  # 518,809

pred <- final_mod$intercept +
  big_prodVec(G, new_beta[ind], ind.row = ind.test, ind.col = ind)

AUCBoot(pred, y.sub[ind.test])  # 60.7 [59.5-61.9] instead of 59.5 [58.2-60.7]


library(tidyverse)

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

max_prs <- grid2 %>% arrange(desc(auc)) %>% slice(1:10) %>% print() %>% slice(1)
#    size thr.r2 grp.num thr.imp num    thr.lp       auc
# 1  1000    0.5       5    0.95 134 0.9999434 0.6033700
# 2  1000    0.5       5    0.95 134 1.0493715 0.6029239
# 3  1000    0.5       5    0.90  98 0.9999434 0.6025065
# 4  1000    0.5       5    0.95 134 1.1012427 0.6021977
# 5  1000    0.5       5    0.90  98 1.0493715 0.6020040
# 6   200    0.5       5    0.95 133 0.9999434 0.6018330
# 7   200    0.5       5    0.95 133 1.0493715 0.6014175
# 8   200    0.5       5    0.90  97 0.9999434 0.6013103
# 9  1000    0.5       5    0.90  98 1.1012427 0.6012757
# 10 1000    0.5       5    0.95 134 1.1556781 0.6009329

ind.keep <- unlist(map(all_keep, max_prs$num))
sum(lpval[ind.keep] > max_prs$thr.lp)  # 135,850
AUCBoot(
  snp_PRS(G, beta[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = max_prs$thr.lp),
  y.sub[ind.test]
) # 60.7 [59.5-62.0] instead of 59.2 [58.0-60.4]

grid2 %>%
  filter(thr.r2 == 0.5, size == 1000) %>%
  mutate(thr.maf = MAF[grp.num]) %>%
  ggplot() +
  geom_point(aes(thr.lp, auc)) +
  facet_grid(thr.imp ~ thr.maf) +
  scale_x_log10(limits = c(1, 4), breaks = c(1, 2, 4), minor_breaks = 1:10) +
  ylim(0.53, NA) +
  theme_bigstatsr(size.rel = 0.7) +
  labs(x = "-log10(p-value) threshold (log scale)", y = "AUC")


# lassosum
library(lassosum)
ind_common <- which(maf > 0.01)
ukbb$map <- dplyr::mutate(ukbb$map, genetic.dist = 0, rsid = NULL,
                          chromosome = as.integer(chromosome))
ukbb$fam <- snp_fake(nrow(G), 1)$fam
bed <- snp_writeBed(ukbb, bedfile = "data/UKBB_MDD_lassosum2.bed",
                    ind.row = ind.train, ind.col = ind_common)
library(doParallel)
registerDoParallel(cl <- makeCluster(NCORES))
system.time(
  out <- lassosum.pipeline(
    cor = p2cor(10^-lpval, n = 59851 + 113154, sign = beta),
    snp = ukbb$map$marker.ID,
    A1 = ukbb$map$allele1,
    test.bfile = "data/UKBB_MDD_lassosum2",
    LDblocks = "EUR.hg19",
    cluster = cl,
    sample = 20e3,
    exclude.ambiguous = FALSE
  )
) # 6.4H
stopCluster(cl)

v <- validate(out, pheno = y.sub[ind.train], validate.function = AUC)
length(ind <- which(v$best.beta != 0))  # 387,817
pred_lassosum <- big_prodVec(G, v$best.beta[ind], ind.row = ind.test,
                             ind.col = ind_common[ind])
AUCBoot(pred_lassosum, y.sub[ind.test]) # 52.0 [50.7-53.2]
