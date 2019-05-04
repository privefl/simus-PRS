# urls <- gsubfn::strapply(
#   readLines("https://datadryad.org//resource/doi:10.5061/dryad.ns8q3"),
#   "<a href=\"(/bitstream/handle/10255/dryad\\.[0-9]+/meta_chr_[0-9]+\\?sequence=1)\">",
#   simplify = 'c')
#
# sumstats <- purrr::map_dfr(urls, ~ {
#   download.file(paste0("https://datadryad.org", .x),
#                 destfile = (tmp <- tempfile(fileext = ".txt")))
#   sumstats <- bigreadr::fread2(
#     tmp, select = c("chromosome", "position", "a0", "a1", "beta.meta", "p.meta"),
#     col.names = c("chr", "pos", "a0", "a1", "beta", "p"))
#   na.omit(sumstats)
# })
#
# saveRDS(sumstats, "sumstats_T1D.rds")

sumstats <- readRDS("sumstats_T1D.rds")
nrow(sumstats)  # 8,996,866
hist(sumstats$p)
hist(sumstats$beta)

sumstats <- subset(sumstats, p < 0.1)

library(bigreadr)
info_snp_UKBB <- rbind_df(lapply(1:22, function(chr) {
  file <- paste0("data/ukb_imp_mfi/ukb_mfi_chr", chr, "_v3.txt")
  df <- fread2(file, select = c(3:5, 8), col.names = c("pos", "a0", "a1", "info"))
  cbind.data.frame(chr = chr, df)
}))
info_snp <- bigsnpr::snp_match(sumstats, info_snp_UKBB)
# 1,241,172 variants in summary statistics.
# 172,374 ambiguous SNPs have been removed.
# 1,051,668 variants have been matched; 76 were flipped and 522 were reversed.
info_snp <- bigsnpr::snp_match(sumstats, info_snp_UKBB, strand_flip = FALSE)
# 1,222,801 variants have been matched; 0 were flipped and 471 were reversed.
info_snp <- subset(na.omit(info_snp), info > 0.3)

list_snp_id <- with(info_snp, split(paste(chr, pos, a0, a1, sep = "_"),
                                    factor(chr, levels = 1:22)))
beta <- info_snp$beta
lpval <- -log10(info_snp$p)
info <- info_snp$info

# Infinite values because of large effects on chromosome 6
lpval <- pmin(lpval, -log10(.Machine$double.xmin) + abs(beta))

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
y[ind_TD1] <- 1
y[ind_TD2] <- NA

sub <- which(!is.na(ind.indiv) & !df0$is_rel2 & df0$is_caucasian & !is.na(y))
length(sub)  # 315,318
table(y.sub <- y[sub])
#      0      1
# 314547    771

NCORES <- 12
system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKBB_T1D",
    ind_row = ind.indiv[sub],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = NCORES
  )
) # 4H

library(bigsnpr)
ukbb <- snp_attach("data/UKBB_T1D.rds")
G <- ukbb$genotypes
file.size(G$backingfile) / 1024^3  # 359 GB
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
    ncores = 6  # use less cores to prevent swapping if not enough RAM
  )
) # 4.4H

system.time(
  multi_PRS <- snp_grid_PRS(
    G, all_keep, betas = beta, lpS = lpval, ind.row = ind.train,
    n_thr_lpS = 50, backingfile = "data/UKBB_T1D_scores", ncores = NCORES
  )
) # 1.8H

system.time(
  final_mod <- snp_grid_stacking(
    multi_PRS, y.sub[ind.train], ncores = NCORES, K = 5
  )
) # 30 min
mod <- final_mod$mod
plot(mod)
summary(mod)
# # A tibble: 3 x 6
#    alpha validation_loss intercept beta           nb_var message
#    <dbl>           <dbl>     <dbl> <list>          <int> <list>
# 1 0.0001          0.0153     -7.33 <dbl [55,300]>   5646 <chr [5]>
# 2 0.01            0.0153     -7.33 <dbl [55,300]>   1775 <chr [5]>
# 3 1               0.0153     -7.32 <dbl [55,300]>   1394 <chr [5]>

new_beta <- final_mod$beta.G

length(ind <- which(new_beta != 0))  # 135,991
summary(new_beta)
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -0.2350732  0.0000000  0.0000000  0.0000032  0.0000000  0.5000515
summary(new_beta[which(sign(new_beta * beta) < 0)])
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -2.277e-02 -7.396e-05 -3.300e-08 -2.162e-05  6.201e-05  1.647e-02

pred <- final_mod$intercept +
  big_prodVec(G, new_beta[ind], ind.row = ind.test, ind.col = ind)

AUCBoot(pred, y.sub[ind.test])  # 78.7 [75.7-81.7]

# Save
save(all_keep, multi_PRS, final_mod, file = "data/res_T1D.RData")

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
# 1  500    0.2     0.3  14 4.580883 0.7489942

ind.keep <- unlist(map(all_keep, std_prs$num))
# Verif on training set
AUC(
  snp_PRS(G, beta[ind.keep], ind.test = ind.train, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = std_prs$thr.lp),
  y.sub[ind.train]
) # 0.7489942
# Eval on test set
AUCBoot(
  snp_PRS(G, beta[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = std_prs$thr.lp),
  y.sub[ind.test]
) # 75.4 [72.4-78.4]
sum(lpval[ind.keep] > std_prs$thr.lp)  # 1112

max_prs <- grid2 %>% arrange(desc(auc)) %>% slice(1:10) %>% print() %>% slice(1)
#     size thr.r2 thr.imp num   thr.lp       auc
# 1  10000   0.01     0.9  58 4.580883 0.7730035
# 2  20000   0.01     0.9  59 4.580883 0.7730035
# 3  50000   0.01     0.9  60 4.580883 0.7730035
# 4   5000   0.01     0.9  57 4.580883 0.7722518
# 5   5000   0.01     0.9  57 5.149850 0.7716009
# 6  10000   0.01     0.9  58 5.149850 0.7713099
# 7  20000   0.01     0.9  59 5.149850 0.7713099
# 8  50000   0.01     0.9  60 5.149850 0.7713099
# 9   5000   0.01     0.9  57 7.316962 0.7691121
# 10  5000   0.01     0.9  57 5.789486 0.7686678

ind.keep <- unlist(map(all_keep, max_prs$num))
AUCBoot(
  snp_PRS(G, beta[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
          lpS.keep = lpval[ind.keep], thr.list = max_prs$thr.lp),
  y.sub[ind.test]
) # 76.9 [73.9-79.7]
sum(lpval[ind.keep] > max_prs$thr.lp)  # 267

ggplot(grid2) +
  geom_point(aes(thr.lp, auc)) +
  facet_grid(thr.imp ~ thr.r2 + size) +
  scale_x_log10(limits = c(1, 420), breaks = c(1, 10, 200)) +
  ylim(0.72, NA) +
  theme_bigstatsr(size.rel = 0.7) +
  labs(x = "-log10(p-value) threshold (log scale)", y = "AUC")


# Check large effects
ind3 <- intersect(which(abs(beta) > 2), ind)
gwas <- big_univLogReg(G, y.sub[ind.train], ind.train, ind3)
cbind(beta_GWAS = beta[ind3], pval_GWAS = 10^(-lpval[ind3]),
      af = big_scale()(G, ind.row = ind.train, ind.col = ind3)$center / 2,
      gwas, pval = predict(gwas, log10 = FALSE))
#    beta_GWAS     pval_GWAS          af      estim    std.err niter      score         pval
# 1     -2.001  2.251049e-02 0.999744425  1.7967667 5.15478608     7  0.3485628 7.274175e-01
# 2      2.110 3.551532e-262 0.765612500  1.1514707 0.11375216     6 10.1226264 4.384884e-24
# 3     -2.169  2.392526e-79 0.069442500 -0.5682690 0.16366915     5 -3.4720591 5.164826e-04
# 4      2.061 1.442919e-276 0.749767500  0.9805250 0.10264324     6  9.5527481 1.263007e-21
# 5      2.790 5.453554e-260 0.752588300  1.1524478 0.11176449     6 10.3113946 6.258645e-25
# 6     -3.677  1.387642e-11 0.003016300 -5.4498698 3.24386571     7 -1.6800541 9.294680e-02
# 7      2.100 1.767439e-310 0.695808750  1.1798415 0.09987523     6 11.8131539 3.338006e-32
# 8      2.095 1.787905e-310 0.753228600  0.9293467 0.10168015     6  9.1399028 6.250848e-20
# 9      2.022 2.115166e-310 0.689873550  1.0132119 0.09393013     6 10.7868672 3.970958e-27
# 10    -2.123  5.370921e-83 0.079237000 -0.6461429 0.16047586     5 -4.0264179 5.663300e-05
# 11    -2.669  3.235930e-27 0.009978575 -0.9999493 0.56553306     6 -1.7681535 7.703524e-02
# 12     2.041 3.194928e-262 0.769042875  1.0739673 0.11185979     6  9.6010138 7.916195e-22
# 13    -2.281  2.244712e-11 0.012840175 -1.7011047 0.71990717     6 -2.3629501 1.813011e-02
# 14    -2.006  5.695059e-99 0.087215425 -0.8001150 0.16265184     6 -4.9191878 8.690405e-07
# 15     2.683  5.610216e-68 0.000724775  1.4633015 1.18905625     7  1.2306411 2.184571e-01
# 16     2.394 8.981408e-311 0.285767075  1.2729734 0.06876607     6 18.5116491 1.663298e-76
# 17     2.316 1.074842e-310 0.106297775  1.1299433 0.07440839     7 15.1856989 4.398500e-52
# 18    -2.158  4.619196e-44 0.032463500 -0.8967315 0.27980792     6 -3.2048110 1.351512e-03
# 19     2.039 2.033969e-310 0.066977400  1.1418232 0.09200736     7 12.4101281 2.302842e-35
