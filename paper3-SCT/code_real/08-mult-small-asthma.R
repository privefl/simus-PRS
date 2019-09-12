# download.file("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/DemenaisF_29273806_GCST006862/TAGC_meta-analyses_results_for_asthma_risk.zip",
#               destfile = "sumstats_asthma.zip")
# unzip("sumstats_asthma.zip")
library(tidyverse)
library(bigreadr)
sumstats <- fread2(
  "TAGC meta-analyses results for asthma risk/TAGC_Multiancestry_and_European-Ancestry_Meta-analyses_Results.tsv",
  select = c(1, 3:5, 18, 20),
  col.names = c("chr", "pos", "a0", "a1", "beta", "p")
)
nrow(sumstats)  # 2,001,280

sumstats <- subset(sumstats, p < 0.5)
nrow(sumstats)  # 1,013,806
hist(sumstats$beta)
hist(sumstats$p)

info_snp_UKBB <- rbind_df(lapply(1:22, function(chr) {
  file <- paste0("data/ukb_imp_mfi/ukb_mfi_chr", chr, "_v3.txt")
  df <- fread2(file, select = c(3:5, 8), col.names = c("pos", "a0", "a1", "info"))
  cbind.data.frame(chr = chr, df)
}))
info_snp <- bigsnpr::snp_match(sumstats, info_snp_UKBB)
# 1,013,806 variants in summary statistics.
# 0 ambiguous SNPs have been removed.
# 1,013,119 variants have been matched; 0 were flipped and 1 were reversed.
info_snp <- subset(na.omit(info_snp), info > 0.3)

list_snp_id <- with(info_snp, split(paste(chr, pos, a0, a1, sep = "_"),
                                    factor(chr, levels = 1:22)))
beta <- info_snp$beta
lpval <- -log10(info_snp$p)
info <- info_snp$info


# subset samples
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

ind_asthma <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x == 1111)),
  lapply(df_ICD10,   function(x) which(substr(x, 1, 3) == "J45"))
))))
ind_respi <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x %in% 1111:1125)),
  lapply(df_ICD10,   function(x) which(substr(x, 1, 1) == "J"))
))))
y <- rep(0, nrow(df0)); y[ind_respi] <- NA; y[ind_asthma] <- 1

sub <- which(!is.na(ind.indiv) & !df0$is_rel2 & df0$is_caucasian & !is.na(y))
length(sub)  # 305,772
table(y.sub <- y[sub])
#      0      1
# 261985  43787


NCORES <- 12
system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKBB_asthma",
    ind_row = ind.indiv[sub],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = NCORES
  )
) # 3.8H


library(bigsnpr)
ukbb <- snp_attach("data/UKBB_asthma.rds")
G <- ukbb$genotypes
file.size(G$backingfile) / 1024^3  # 288 GB
CHR <- as.integer(ukbb$map$chromosome)
POS <- ukbb$map$physical.pos
Nss <- 19954 + 107715

options(bigstatsr.block.sizeGB = 50)

for (ic in 1:6) {

  res_file <- paste0("res_small/asthma_", ic, ".rds")
  if (file.exists(res_file)) next
  cat(ic, "\n")

  set.seed(ic)
  ind.train <- c(sample(which(y.sub == 0), 2000), sample(which(y.sub == 1), 500))
  ind.test <- setdiff(seq_along(sub), ind.train)

  tmp <- tempfile(tmpdir = "res_small")

  #### SCT ####
  system.time(
    all_keep <- snp_grid_clumping(
      G, CHR, POS, lpS = lpval, ind.row = ind.train, infos.imp = info,
      grid.thr.imp = c(0.3, 0.6, 0.9, 0.95),
      grid.thr.r2 = c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.95),
      grid.base.size = c(50, 100, 200, 500),
      ncores = 5
    )
  ) # 23 min

  system.time(
    multi_PRS <- snp_grid_PRS(
      G, all_keep, betas = beta, lpS = lpval, ind.row = ind.train,
      n_thr_lpS = 50, backingfile = paste0(tmp, "_scores"), ncores = NCORES
    )
  ) # 13 min

  system.time(
    final_mod <- snp_grid_stacking(
      multi_PRS, y.sub[ind.train], ncores = NCORES
    )
  ) # 2 min
  mod <- final_mod$mod
  plot(mod)
  summary(mod)

  new_beta <- final_mod$beta.G
  nb_SCT <- length(ind <- which(new_beta != 0))

  pred <- final_mod$intercept +
    big_prodVec(G, new_beta[ind], ind.row = ind.test, ind.col = ind)

  #### C+T ####
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
    slice(1)

  # Eval on test set
  ind.keep <- unlist(map(all_keep, std_prs$num))
  pred_std_prs <- snp_PRS(G, beta[ind.keep],
                          ind.test = ind.test, ind.keep = ind.keep,
                          lpS.keep = lpval[ind.keep], thr.list = std_prs$thr.lp)
  nb_std_prs <- sum(lpval[ind.keep] > std_prs$thr.lp)

  max_prs <- grid2 %>% arrange(desc(auc)) %>% slice(1)

  ind.keep <- unlist(map(all_keep, max_prs$num))
  pred_max_prs <- snp_PRS(G, beta[ind.keep],
                          ind.test = ind.test, ind.keep = ind.keep,
                          lpS.keep = lpval[ind.keep], thr.list = max_prs$thr.lp)
  nb_max_prs <- sum(lpval[ind.keep] > max_prs$thr.lp)

  #### lassosum ####
  library(lassosum)
  ukbb$map <- dplyr::mutate(ukbb$map, genetic.dist = 0, rsid = NULL,
                            chromosome = as.integer(chromosome))
  ukbb$fam <- snp_fake(nrow(G), 1)$fam
  bed <- snp_writeBed(ukbb, bedfile = paste0(tmp, ".bed"), ind.row = ind.train)
  library(doParallel)
  registerDoParallel(cl <- makeCluster(NCORES))
  system.time(
    out <- lassosum.pipeline(
      cor = p2cor(p = 10^-lpval, n = Nss, sign = beta),
      snp = ukbb$map$marker.ID,
      A1 = ukbb$map$allele1,
      test.bfile = tmp,
      LDblocks = "EUR.hg19",
      cluster = cl,
      exclude.ambiguous = FALSE
    )
  ) # 9 min
  stopCluster(cl)

  v <- validate(out, pheno = y.sub[ind.train], validate.function = AUC)
  nb_lassosum <- length(ind <- which(v$best.beta != 0))  # 290,204
  pred_lassosum <- big_prodVec(G, v$best.beta[ind], ind.row = ind.test, ind.col = ind)


  #### LDpred ####
  file_sumstats <- paste0(tmp, ".txt")
  mutate(ukbb$map, beta = beta, pval = 10^-lpval) %>%
    bigreadr::fwrite2(file_sumstats, sep = "\t")
  file_hdf5 <- paste0(tmp, ".hdf5")

  reticulate::use_python("/opt/rh/rh-python36/root/usr/bin/python")
  reticulate::py_config()
  # system("python3 --version", intern = TRUE)
  ldpred <- "../ldpred/LDpred.py"
  # system(glue::glue("python3 {ldpred} coord --help"))
  system(glue::glue(
    "python3 {ldpred} coord",
    " --gf {tmp}",
    " --ssf {file_sumstats}",
    " --skip-coordination",
    " --rs marker.ID --A1 allele1 --A2 allele2 --pos physical.pos --chr chromosome",
    " --pval pval --eff beta --beta",
    " --N {Nss}",
    " --out {file_hdf5}"
  )) # 23 min


  system.time(
    system(glue::glue(
      "python3 {ldpred} gibbs",
      " --cf {file_hdf5}",
      " --ldr {round(nrow(ukbb$map) / 3000)}",
      " --ldf {tmp}",
      " --N {Nss}",
      " --h2 0.3",
      " --out {tmp}"
    ))
  ) # 438 min

  ext <- c(sprintf("_LDpred_p%.4e.txt", c(1, 0.3, 0.1, 0.03, 0.01, 0.003, 0.001)),
           "_LDpred-inf.txt")
  files_ldpred <- paste0(tmp, ext)
  beta_ldpred <- sapply(files_ldpred, function(file) {
    res_ldpred <- bigreadr::fread2(file, select = c(3, 7))
    beta_ldpred <- numeric(nrow(ukbb$map))
    beta_ldpred[match(res_ldpred$sid, ukbb$map$marker.ID)] <- res_ldpred[[2]]
    beta_ldpred
  })

  system.time(
    pred_train_ldpred <- big_prodMat(G, beta_ldpred, ind.row = ind.train,
                                     ncores = NCORES)
  ) # 23 min
  auc_train_ldpred <- apply(-pred_train_ldpred, 2, AUC, target = y.sub[ind.train])

  beta_ldpred_max <- beta_ldpred[, which.max(auc_train_ldpred)]
  nb_LDpred <- sum(beta_ldpred_max != 0)
  system.time(
    pred_ldpred <- big_prodVec(G, ind.row = ind.test, beta_ldpred_max)) # 16 min

  #### Save results ####
  saveRDS(
    data.frame(
      auc_std_prs  = AUC(pred_std_prs,  y.sub[ind.test]), nb_std_prs,
      auc_max_prs  = AUC(pred_max_prs,  y.sub[ind.test]), nb_max_prs,
      auc_SCT      = AUC(pred,          y.sub[ind.test]), nb_SCT,
      auc_lassosum = AUC(pred_lassosum, y.sub[ind.test]), nb_lassosum,
      auc_LDpred   = AUC(-pred_ldpred,  y.sub[ind.test]), nb_LDpred
    ),
    res_file
  )
  stopifnot(file.exists(res_file))

  # Cleanup
  unlink(paste0(tmp, "*"))

}
