sumstats <- bigreadr::fread2("data/sumstats.txt")

library(dplyr)
sumstats %>%
  select(hg19chrc = chromosome, snpid = rsid,
         a1 = allele1, a2 = allele2, bp = physical.pos,
         or = beta, p = pval) %>%
  mutate(hg19chrc = paste0("chr", hg19chrc), or = exp(or),
         snpid = bigreadr::fread2("data/data_train.bim", select = 2)[[1]]) %>%
  bigreadr::fwrite2("SUM_STATS_FILE.txt", sep = "\t")

readLines("SUM_STATS_FILE.txt", n = 5)
# [1] "hg19chrc\tsnpid\ta1\ta2\tbp\tor\tp"
# [2] "chr6\t6:202076_AG_A\tAG\tA\t202076\t0.960642839169572\t0.333452277652367"
# [3] "chr6\t6:202777_C_T\tC\tT\t202777\t0.960338117401514\t0.329510213501189"
# [4] "chr6\trs11757325\tC\tT\t203397\t1.00843386744862\t0.846133552351405"
# [5] "chr6\trs80014302\tT\tG\t203722\t0.960019650560543\t0.325411512307528"

# famfile <- "data/data_test.fam"
# famfile %>%
#   bigreadr::fread2() %>%
#   mutate(V6 = V6 + 1) %>%
#   bigsnpr:::write.table2(famfile)

reticulate::use_python("/home/privef/anaconda3/bin/python3")
reticulate::py_config()
stopifnot(system("python3 --version", intern = TRUE) == "Python 3.7.0")
ldpred <- "../ldpred/LDpred.py"
unlink("OUT_COORD_FILE.hdf5")
system(glue::glue(
  "python3 {ldpred} coord",
  " --gf data/data_train",
  " --ssf SUM_STATS_FILE.txt --ssf-format BASIC",
  " --N {bigreadr::nlines('data/data_gwas.fam')}",
  " --out OUT_COORD_FILE.hdf5"
))
# How to control which are excluded?
# --skip-coordination
# --maf 0


system.time(
  system(glue::glue(
    "python3 {ldpred} gibbs",
    " --cf OUT_COORD_FILE.hdf5",
    " --ldr {round(bigreadr::nlines('data/data_gwas.bim') / 3000)}",
    " --ldf LD_FILE",
    " --f 0.01",
    " --N {bigreadr::nlines('data/data_gwas.fam')}",
    " --out OUT_WEIGHTS_FILE"
  ))
) # 12 min

# Evaluate without LDpred
weights <- bigreadr::fread2("OUT_WEIGHTS_FILE_LDpred_p1.0000e-02.txt")
weights <- bigreadr::fread2("OUT_WEIGHTS_FILE_LDpred-inf.txt")
stopifnot(!any(duplicated(weights$sid)))

train <- snp_attach("data/data_train.rds")
ind <- match(weights$sid, train$map$marker.ID)
all.equal(sumstats$beta[ind], weights$raw_beta)
prs.train <- big_prodVec(train$genotypes, weights$ldpred, ind.col = ind)
AUC(prs.train, train$fam$affection - 1)  # 60.8 with 0.01 -> 65.0 with Inf

test <- snp_attach("data/data_test.rds")
prs.test <- big_prodVec(test$genotypes, weights$ldpred, ind.col = ind)
AUC(prs.test, test$fam$affection)  # 58.7 with 0.01 -> 64.9 with Inf

maf <- snp_MAF(train$genotypes, ind.col = ind)
plot(maf, weights$ldpred, pch = 20)

# Evaluate with LDpred
system(glue::glue(
  "python3 {ldpred} score",
  " --gf data/data_test",
  " --rf OUT_WEIGHTS_FILE",
  " --out OUT_SCORE_FILE"
))

writeLines(readLines("OUT_SCORE_FILE_LDpred-inf.txt", n = 5))
scores <- readLines("OUT_SCORE_FILE_LDpred-inf.txt") %>%
  gsub(", ", " ", .) %>%
  trimws() %>%
  data.table::fread(text = ., data.table = FALSE)
AUC(scores$raw_effects_prs, scores$true_phens - 1) # 63.3
AUC(scores$pval_derived_effects_prs, scores$true_phens - 1) # 35.1?
plot(prs.test, scores$pval_derived_effects_prs, pch = 20)
plot(scores$raw_effects_prs, scores$pval_derived_effects_prs, pch = 20)
