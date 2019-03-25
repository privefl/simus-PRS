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

plot(sumstats$beta[ind], v$best.beta[ind], pch = 20); abline(0, 1, col = "red")

# pROC::roc(cases = pred[test$fam$affection == 1],
#           controls =  pred[test$fam$affection == 0], ci=TRUE)
#
# mod <- glm(test$fam$affection ~ pred, family = "binomial")
# cat("Nagelkerke's R2 =", DescTools::PseudoR2(mod, "Nagelkerke"), "\n")
#
# pred2 <- 0.001 * pred
# mod2 <- glm(test$fam$affection ~ pred2, family = "binomial")
# cat("Nagelkerke's R2 =", DescTools::PseudoR2(mod2, "Nagelkerke"), "\n")
#
# pred3 <- 1 / (1 + exp(-pred))
# mod3 <- glm(test$fam$affection ~ pred3, family = "binomial")
# cat("Nagelkerke's R2 =", DescTools::PseudoR2(mod3, "Nagelkerke"), "\n")
