library(bigsnpr)
# snp_readBed("data/data_train.bed")
train <- snp_attach("data/data_train.rds")
G <- train$genotypes
y <- train$fam$affection

sumstats <- bigreadr::fread2("data/sumstats.txt")
CHR <- sumstats$chromosome
POS <- sumstats$physical.pos
LPVAL <- -log10(sumstats$pval)
pf <- 1 / pmax(LPVAL, 1)
hist(pf)

system.time(
  mod <- big_spLogReg(G, y, alphas = 10^(-(4:8/2)), ncores = 6,
                      nlam.min = 80, dfmax = Inf, n.abort = 10)
) # 13 min
plot(mod)
summary(mod)

test <- snp_attach("data/data_test.rds")
pred.test <- predict(mod, test$genotypes)
AUCBoot(pred.test, test$fam$affection)
# 71.3 [68.4-74.1]

#### Previous simu:
## alpha = 1 -> 68.6 [65.8-71.4]
## alpha = 1e-4 - > 70.8 [68.1-73.5]
## without pf -> 66.5 [63.6-69.3]
