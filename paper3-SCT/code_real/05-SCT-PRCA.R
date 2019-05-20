# download.file("http://practical.icr.ac.uk/blog/wp-content/uploads/uploadedfiles/oncoarray/MetaSummaryData/meta_v3_onco_euro_overall_ChrAll_1_release.zip",
#               destfile = "sumstats_PRCA.zip")
# unzip("sumstats_PRCA.zip")
library(bigreadr)
sumstats <- fread2(
  "meta_v3_onco_euro_overall_ChrAll_1_release.txt",
  select = c("Chr", "position", "Allele1", "Allele2", "Effect", "Pvalue"),
  col.names = c("chr", "pos", "a1", "a0", "beta", "p")
)
nrow(sumstats)  # 20,370,946
library(dplyr)
sumstats <- sumstats %>%
  filter(p < 0.05) %>%
  mutate(a0 = toupper(a0), a1 = toupper(a1))
nrow(sumstats)  # 1,416,082

hist(sumstats$beta)
hist(sumstats$p)

info_snp_UKBB <- rbind_df(lapply(1:22, function(chr) {
  file <- paste0("data/ukb_imp_mfi/ukb_mfi_chr", chr, "_v3.txt")
  df <- fread2(file, select = c(3:5, 8), col.names = c("pos", "a0", "a1", "info"))
  cbind.data.frame(chr = chr, df)
}))
info_snp <- bigsnpr::snp_match(sumstats, info_snp_UKBB)
# 1,416,082 variants in summary statistics.
# 198,181 ambiguous SNPs have been removed.
# 1,156,517 variants have been matched; 51 were flipped and 502,204 were reversed.
info_snp <- bigsnpr::snp_match(sumstats, info_snp_UKBB, strand_flip = FALSE)
# 1,345,888 variants have been matched; 0 were flipped and 597,133 were reversed.
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
ind_PRCA <- sort(unique(unlist(c(
  lapply(df_cancer1, function(x) which(x == 1044)),
  lapply(df_cancer2, function(x) which(x %in% c("C61", "D075")))
))))
table(df0$sex[ind_PRCA])
# 0    1
# 9 9282

y <- rep(NA, nrow(df0))
y[rowMeans(df_cancer0 == 0, na.rm = TRUE) == 1] <- 0
y[ind_PRCA] <- 1

sub <- which(!is.na(ind.indiv) & !df0$is_rel2 & df0$is_caucasian &
               df0$sex == 1 & !is.na(y))
length(sub)  # 147,964
table(y.sub <- y[sub])
#      0      1
# 141321   6643

NCORES <- 12
system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKBB_PRCA",
    ind_row = ind.indiv[sub],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = NCORES
  )
) # 3.3H


library(bigsnpr)
ukbb <- snp_attach("data/UKBB_PRCA.rds")
G <- ukbb$genotypes
file.size(G$backingfile) / 1024^3  # 181 GB
CHR <- as.integer(ukbb$map$chromosome)
POS <- ukbb$map$physical.pos

set.seed(1); ind.train <- sort(sample(length(sub), 120e3))
# set.seed(2); ind.train <- c(sample(which(y.sub == 0), 2000), sample(which(y.sub == 1), 500))
ind.test <- setdiff(seq_along(sub), ind.train)

system.time(
  all_keep <- snp_grid_clumping(
    G, CHR, POS, lpS = lpval, ind.row = ind.train, infos.imp = info,
    grid.thr.imp = c(0.3, 0.6, 0.9, 0.95),
    grid.thr.r2 = c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.95),
    grid.base.size = c(50, 100, 200, 500),
    ncores = 6  # use less cores to prevent swapping if not enough RAM
  )
) # 9H -> 3H

system.time(
  multi_PRS <- snp_grid_PRS(
    G, all_keep, betas = beta, lpS = lpval, ind.row = ind.train,
    n_thr_lpS = 50, backingfile = "data/UKBB_PRCA_scores", ncores = NCORES
  )
) # 1H

system.time(
  final_mod <- snp_grid_stacking(multi_PRS, y.sub[ind.train], ncores = NCORES)
) # 2H
mod <- final_mod$mod
plot(mod)
summary(mod)
# # A tibble: 3 x 6
#    alpha validation_loss intercept beta           nb_var message
#    <dbl>           <dbl>     <dbl> <list>          <int> <list>
# 1 0.0001           0.171     0.778 <dbl [83,804]>  28828 <chr [10]>
# 2 0.01             0.171     1.00  <dbl [83,804]>   6353 <chr [10]>
# 3 1                0.171     1.05  <dbl [83,804]>   4501 <chr [10]>

# Save
save(all_keep, multi_PRS, final_mod, file = "data/res_PRCA.RData")


new_beta <- final_mod$beta.G

length(ind <- which(new_beta != 0))  # 696,575
summary(new_beta)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# -1.857336 -0.000003  0.000000 -0.000072  0.000003  3.749961
summary(new_beta[which(sign(new_beta * beta) < 0)])
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -0.3079258 -0.0001564  0.0000003  0.0000647  0.0001823  0.3774212

pred <- final_mod$intercept +
  big_prodVec(G, new_beta[ind], ind.row = ind.test, ind.col = ind)

AUCBoot(pred, y.sub[ind.test])  # 71.7 [70.2-73.1] / 69.3 [68.7-70.0]


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
# 1  500    0.2     0.3  14 5.95265 0.6695709

ind.keep <- unlist(map(all_keep, std_prs$num))
# Verif on training set
AUC(
  snp_PRS(G, beta[ind.keep], ind.test = ind.train, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = std_prs$thr.lp),
  y.sub[ind.train]
) # 0.6695709
# Eval on test set
AUCBoot(
  snp_PRS(G, beta[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = std_prs$thr.lp),
  y.sub[ind.test]
) # 68.0 [66.5-69.5] / 67.1 [66.4-67.8]
sum(lpval[ind.keep] > std_prs$thr.lp)  # 1035

max_prs <- grid2 %>% arrange(desc(auc)) %>% slice(1:10) %>% print() %>% slice(1)
#     size thr.r2 thr.imp num   thr.lp       auc
# 1  10000   0.01     0.9  58 5.378748 0.6865972
# 2  20000   0.01     0.9  59 5.378748 0.6865972
# 3  50000   0.01     0.9  60 5.378748 0.6865972
# 4   5000   0.01     0.9  57 5.378748 0.6865827
# 5  10000   0.01     0.9  58 5.952650 0.6857405
# 6  20000   0.01     0.9  59 5.952650 0.6857405
# 7  50000   0.01     0.9  60 5.952650 0.6857405
# 8   5000   0.01     0.9  57 5.952650 0.6857095
# 9   5000   0.01     0.6  29 5.952650 0.6850979
# 10  5000   0.01     0.9  57 4.860176 0.6850812

ind.keep <- unlist(map(all_keep, max_prs$num))
AUCBoot(
  snp_PRS(G, beta[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = max_prs$thr.lp),
  y.sub[ind.test]
) # 69.3 [67.8-70.8] / 68.7 [68.0-69.3]
sum(lpval[ind.keep] > max_prs$thr.lp)  # 356

ggplot(grid2) +
  geom_point(aes(thr.lp, auc)) +
  facet_grid(thr.imp ~ thr.r2 + size) +
  scale_x_log10(limits = c(1.5, 10), breaks = c(1, 2, 5), minor_breaks = 1:10) +
  ylim(0.62, NA) +
  theme_bigstatsr(size.rel = 0.7) +
  labs(x = "-log10(p-value) threshold (log scale)", y = "AUC")

sum(abs(beta) > 10) # 24088
