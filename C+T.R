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
qplot(y = LPVAL) + ylim(1, NA)
THR <- exp(seq(log(0.01), log(0.999 * max(LPVAL)), length.out = 50))
hist(LPVAL); abline(v = THR, lty = 3)

ind.keep <- snp_clumping(G, infos.chr = CHR, S = LPVAL, thr.r2 = 0.05,
                         is.size.in.bp = TRUE, infos.pos = POS,
                         ncores = nb_cores())
qplot(ind.keep, LPVAL[ind.keep]) + ylim(1, NA) +
  geom_hline(yintercept = -log10(5e-8), color = "red", linetype = 2)

prs.train <- snp_PRS(train$genotypes, -BETA[ind.keep], ind.keep = ind.keep,
                     lpS.keep = LPVAL[ind.keep], thr.list = THR)
auc.train <- apply(prs.train, 2, AUC, train$fam$affection - 1)
plot(auc.train, pch = 20)
AUCBoot(prs.train[, 1], train$fam$affection - 1)
#        Mean        2.5%       97.5%          Sd
# 0.684868419 0.670152389 0.699043396 0.007377221

# snp_readBed("data/data_test.bed")
test <- snp_attach("data/data_test.rds")
prs.test <- snp_PRS(test$genotypes, -BETA[ind.keep], ind.keep = ind.keep,
                    lpS.keep = LPVAL[ind.keep], thr.list = THR)
dim(prs.test) # 2000 x 50
auc.test <- apply(prs.test, 2, AUC, test$fam$affection)
plot(auc.test, pch = 20)

AUCBoot(prs.test[, which.max(auc.train)], test$fam$affection)
#       Mean       2.5%      97.5%         Sd
# 0.69636273 0.66797865 0.72373555 0.01428186
