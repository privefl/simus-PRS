reticulate::use_python("/usr/bin/python2")
reticulate::py_config()
system("python --version")

prscs <- "../PRScs/PRScs.py"
library(glue)
system(glue("{prscs} --help"))

library(dplyr)
bigreadr::fread2("data/sumstats.txt") %>%
  select(SNP = rsid, A1 = allele1, A2 = allele2, BETA = beta, P = pval) %>%
  mutate(SNP = bigreadr::fread2("data/data_train.bim", select = 2)[[1]]) %>%
  bigreadr::fwrite2("SUM_STATS_FILE.txt", sep = "\t") %>%
  readLines(n = 5)

dir.create("PRSCS", showWarnings = FALSE)
system(glue(
  "python {prscs}",
  " --ref_dir={dirname(prscs)}",
  " --bim_prefix=data/data_train",
  " --sst_file=SUM_STATS_FILE.txt",
  " --n_gwas=40000",
  " --out_dir=PRSCS",
  " --chrom=6,8"
))
