# bim <- bigreadr::fread2("data/data_train.bim")
# bigreadr::fwrite2(bim, "data/data_train.bim", sep = "\t", col.names = FALSE)
writeLines(readLines("data/data_train.bim", n = 5))

sumstats <- bigreadr::fread2("data/sumstats.txt")
# sumstats$N <- 40e3
# sumstats$rsid <- bigreadr::fread2("data/data_train.bim", select = 2)[[1]]
# bigreadr::fwrite2(sumstats, "data/sumstats.txt", sep = "\t")
writeLines(readLines("data/sumstats.txt", n = 5))

reticulate::use_python("/home/privef/anaconda3/bin/python3")
reticulate::py_config()
stopifnot(system("python3 --version", intern = TRUE) == "Python 3.7.0")
# system("~/anaconda3/bin/pip install --user --upgrade pip")
# system("~/anaconda3/bin/pip install --user plinkio")
# system("~/anaconda3/bin/pip install --user ldpred")
ldpred <- "../ldpred/LDpred.py"
unlink("OUT_COORD_FILE.hdf5")
# system(glue::glue("python3 {ldpred} coord --help"))
system(glue::glue(
  "python3 {ldpred} coord",
  " --gf data/data_train",
  " --ssf data/sumstats.txt",
  " --maf 0 --skip-coordination",
  " --rs rsid --A1 allele1 --A2 allele2 --pos physical.pos --chr chromosome",
  " --pval pval --eff beta --beta",
  " --N 40000 --case-n 4000 --control-n 36000",
  " --out OUT_COORD_FILE.hdf5"
))

# ------------------------------ Summary statistics ------------------------------
# Num SNPs parsed from sum stats file                                       685920
# --------------------------------- Coordination ---------------------------------
# Num individuals in LD Reference data:                                       8000
# SNPs in LD Reference data:                                                686066
# Num chromosomes used:                                                          2
# SNPs common across datasets:                                              685395
# SNPs retained after filtering:                                            547355
# SNPs w ambiguous nucleotides filtered:                                     97820
# SNPs w other nucleotide discrepancies filtered:                            40220
# SNPs w allele freq discrepancy > 0.100 filtered:                               0

# system(glue::glue("python3 {ldpred} gibbs --help"))
unlink("LD_FILE_ldradius*")
system.time(
  system(glue::glue(
    "python3 {ldpred} gibbs",
    " --cf OUT_COORD_FILE.hdf5",
    " --ldr 1000",
    " --ldf LD_FILE",
    # " --f 0.01",
    " --h2 0.5",
    " --N 40000",
    " --out OUT_WEIGHTS_FILE"
  ))
) # MemoryError

train <- snp_attach("data/data_train.rds")
maf <- snp_MAF(train$genotypes, ind.col = ind)
plot(maf, weights$ldpred, pch = 20)
# test <- snp_attach("data/data_test.rds")

# Evaluate with LDpred
system(glue::glue(
  "python3 {ldpred} score",
  " --gf data/data_test",
  " --rf OUT_WEIGHTS_FILE",
  " --out OUT_SCORE_FILE"
))
sapply(list.files(pattern = "^OUT_SCORE_FILE_LDpred.*\\.txt$"), function(file) {
  bigreadr::fread2(file) %>%
  { AUC(.$PRS, .$true_phens) }
})

writeLines(readLines("OUT_SCORE_FILE_LDpred-inf.txt", n = 5))
scores <- bigreadr::fread2("OUT_SCORE_FILE_LDpred-inf.txt")
AUC(scores$PRS, scores$true_phens) # 65.8

# P+T
system.time(
  system(glue::glue(
    "python3 {ldpred} p+t",
    " --cf OUT_COORD_FILE.hdf5",
    " --ldr {round(bigreadr::nlines('data/data_train.bim') / 3000)}",
    " --r2 0.05",
    " --out OUT_WEIGHTS_FILE"
  ))
) # 5 min

# Evaluate with LDpred
system(glue::glue(
  "python3 {ldpred} score",
  " --gf data/data_test",
  " --rf OUT_WEIGHTS_FILE",
  " --rf-format P+T",
  " --r2 0.05",
  " --out OUT_SCORE_FILE"
))

sapply(list.files(pattern = "^OUT_SCORE_FILE_P\\+T_p.*\\.txt$"), function(file) {
  bigreadr::fread2(file) %>%
  { AUC(.$PRS, .$true_phens) }
})


