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

ind.keep <- snp_clumping(G, infos.chr = CHR, S = LPVAL, thr.r2 = 0.05, size = 8000,
                         is.size.in.bp = TRUE, infos.pos = POS,
                         ncores = nb_cores())
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
plot(auc.test, pch = 20)

AUCBoot(prs.test[, which.max(auc.train)], test$fam$affection)
# 75.0 [72.4.8-77.6] (2000) || 75.2 [72.6-77.8] (4000)

