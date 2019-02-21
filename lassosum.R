# devtools::install_github("tshmak/lassosum")
library(lassosum)
library(bigsnpr)

sumstats <- bigreadr::fread2("data/sumstats.txt")
CHR <- sumstats$chromosome
POS <- sumstats$physical.pos
A1 <- sumstats$allele1
A2 <- sumstats$allele2
cor <- p2cor(p = sumstats$pval, sign = sumstats$beta, n = 40e3)
cor <- ifelse(is.na(cor), max(abs(cor), na.rm = TRUE) * sign(sumstats$beta), cor)
library(doParallel)
registerDoParallel(cl <- makeCluster(nb_cores()))
system.time(
  out <- lassosum.pipeline(cor = cor, chr = CHR, pos = POS, A1 = A1, A2 = A2,
                           test.bfile = "data/data_train",
                           LDblocks = "EUR.hg19",
                           cluster = cl,
                           exclude.ambiguous = FALSE)
) # 23 min
stopCluster(cl)

v <- validate(out)
test <- snp_attach("data/data_test.rds")
length(ind <- which(v$best.beta != 0))
pred <- big_prodVec(test$genotypes, v$best.beta[ind], ind.col = ind)
AUCBoot(pred, test$fam$affection)

## New simu: 72.9 [70.1-75.5]
## Previous simu: 72.7 [69.9-75.4]
