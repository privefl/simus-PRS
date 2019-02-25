library(bigsnpr)
# snp_readBed("data/data_train.bed")
train <- snp_attach("data/data_train.rds")
G <- train$genotypes

sumstats <- bigreadr::fread2("data/sumstats.txt")
CHR <- sumstats$chromosome
POS <- sumstats$physical.pos
BETA <- sumstats$beta
LPVAL <- -log10(sumstats$pval)

library(ggplot2)
# qplot(y = LPVAL) + ylim(1, NA)
THR <- exp(seq(log(0.01), log(0.999 * max(LPVAL)), length.out = 50))
hist(LPVAL); abline(v = THR, lty = 3)

ind.keep <- snp_clumping(G, infos.chr = CHR, S = LPVAL, thr.r2 = 0.1, size = 250,
                         infos.pos = POS, ncores = nb_cores())
qplot(ind.keep, LPVAL[ind.keep]) + ylim(1, NA) +
  geom_hline(yintercept = -log10(5e-8), color = "red", linetype = 2)


prs.train <- snp_PRS(train$genotypes, BETA[ind.keep], ind.keep = ind.keep,
                     lpS.keep = LPVAL[ind.keep], thr.list = THR)
auc.train <- apply(prs.train, 2, AUC, train$fam$affection)
plot(auc.train, pch = 20)
# AUCBoot(prs.train[, 1], train$fam$affection - 1)  ## 68.5 [67.0-69.9]

# snp_readBed("data/data_test.bed")
test <- snp_attach("data/data_test.rds")
prs.test <- snp_PRS(test$genotypes, BETA[ind.keep], ind.keep = ind.keep,
                    lpS.keep = LPVAL[ind.keep], thr.list = THR)
dim(prs.test) # 2000 x 50
auc.test <- apply(prs.test, 2, AUC, test$fam$affection)
plot(THR, auc.test, pch = 20, log = "x")

AUCBoot(prs.test[, which.max(auc.train)], test$fam$affection)
# 62.8 [59.8-65.7] (with thr.r2 = 0.1 and size = 250)
# 75.0 [72.3-77.6]
# 77.5 [74.9-80.1] (with thr.r2 = 0.01 and size = 10000)

#### Find beta for stacking:
b <- rnorm(ncol(prs.test))
true <- prs.test %*% b

ind_last_thr <- rowSums(outer(LPVAL[ind.keep], THR, '>'))
b2 <- c(0, cumsum(b))
try <- big_prodVec(test$genotypes, BETA[ind.keep] * b2[ind_last_thr + 1L],
                   ind.col = ind.keep)
plot(try, true, pch = 20); abline(0, 1, col = "red")
