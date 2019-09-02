library(bigsnpr)
simu <- snp_attach("data/ukbb4simu.rds")
G <- simu$genotypes
CHR <- as.integer(simu$map$chromosome)
POS <- simu$map$physical.pos
INFO <- readRDS("data/ukbb4simu_info.rds")$info
SD <- readRDS("data/ukbb4simu_stats.rds")$scale
load("data/ukbb4simu_ind.RData")
NCORES <- nb_cores()
Nss <- length(ind.gwas)

G.train <- snp_attach("data/ukbb4simu_train.rds")$genotypes
G.test  <- snp_attach("data/ukbb4simu_test.rds")$genotypes

for (ic in 1:10) {

  res_file <- paste0("res_simu/1M_", ic, ".rds")
  if (file.exists(res_file)) next
  cat(ic, "\n")
  tmp <- tempfile(tmpdir = "res_simu")

  #### Simulate phenotype ####
  set.seed(ic)
  h2 <- 0.5; K <- 0.1; M <- 1e6
  set <- sort(sample(ncol(G), size = M))
  effects <- rnorm(M, sd = sqrt(h2 / M))
  y <- big_prodVec(G, effects / SD[set], ind.col = set)
  y <- (y - mean(y)) / sd(y) * sqrt(h2)         ## make sure that var(y) = h2
  y <- y + rnorm(nrow(G), sd = sqrt(1 - h2))
  y <- as.integer(y > qnorm(1 - K))             ## LTM

  #### GWAS ####
  system.time(
    gwas <- big_apply(G, function(G, ind, ind.ca, ind.co) {

      r <- length(ind.ca) + 0
      s <- length(ind.co) + 0
      CCa <- big_counts(G, ind.row = ind.ca, ind.col = ind)
      CCo <- big_counts(G, ind.row = ind.co, ind.col = ind)
      cco1 <- CCo[2, ]; cco01 <- 2 * CCo[1, ] + cco1; cco12 <- 2 * s - cco01
      cca1 <- CCa[2, ]; cca01 <- 2 * CCa[1, ] + cca1; cca12 <- 2 * r - cca01

      num2 <- (s * cca12 - r * cco12)
      deno <- (cco1 + 4 * CCo[3, ] + cca1 + 4 * CCa[3, ]) - (cca12 + cco12)^2 / (r + s)
      deno2 <- sqrt(r * s * deno)

      data.frame(or = cco01 / cco12 * cca12 / cca01, z = num2 / deno2)
    }, a.combine = "rbind", ncores = NCORES,
    ind.ca = intersect(which(y == 1), ind.gwas),
    ind.co = intersect(which(y == 0), ind.gwas))
  ) # 40 min

  beta_gwas <- log(gwas$or)
  # hist(beta_gwas[beta_gwas > -0.1 & beta_gwas < 0.1])
  lpval <- -pchisq(gwas$z^2, df = 1, lower.tail = FALSE, log.p = TRUE) / log(10)
  # hist(10^(-lpval))

  #### SCT ####
  system.time(
    all_keep <- snp_grid_clumping(
      G.train, CHR, POS, lpS = lpval, ncores = NCORES,
      infos.imp = INFO, grid.thr.imp = c(0.3, 0.6, 0.9, 0.95),
      grid.thr.r2 = c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.95),
      grid.base.size = c(50, 100, 200, 500)
    )
  ) # 24 min

  system.time(
    multi_PRS <- snp_grid_PRS(
      G.train, all_keep, betas = beta_gwas, lpS = lpval,
      backingfile = paste0(tmp, "_scores"),
      n_thr_lpS = 50, ncores = NCORES
    )
  ) # 1 min

  system.time(
    final_mod <- snp_grid_stacking(
      multi_PRS, y[ind.train], ncores = NCORES
    )
  ) # 13 min
  mod <- final_mod$mod
  # summary(mod)
  # plot(mod)

  new_beta <- final_mod$beta.G
  nb_SCT <- length(ind <- which(new_beta != 0))  # 560,374

  pred <- final_mod$intercept +
    big_prodVec(G.test, new_beta[ind], ind.col = ind)

  # AUCBoot(pred, y[ind.test])  # 69.3 [67.6-70.9]


  library(tidyverse)

  # ind2 <- sort(sample(ind, 10e3))
  # ggplot(data.frame(y = new_beta, x = beta_gwas)[ind, ]) +
  #   geom_abline(slope = 1, intercept = 0, color = "red") +
  #   geom_abline(slope = 0, intercept = 0, color = "blue") +
  #   geom_point(aes(x, y), size = 0.6) +
  #   theme_bigstatsr() +
  #   labs(x = "Effect sizes from GWAS", y = "Non-zero effect sizes from SCT")

  #### C+T ####
  grid2 <- attr(all_keep, "grid") %>%
    mutate(thr.lp = list(attr(multi_PRS, "grid.lpS.thr")), num = row_number()) %>%
    unnest()
  s <- nrow(grid2)
  grid2$auc <- big_apply(multi_PRS, a.FUN = function(X, ind, s, y.train) {
    # Sum over all chromosomes, for the same C+T parameters
    single_PRS <- rowSums(X[, ind + s * (0:21)])
    bigstatsr::AUC(single_PRS, y.train)
  }, ind = 1:s, s = s, y.train = y[ind.train],
  a.combine = 'c', block.size = 1, ncores = NCORES)

  std_prs <- grid2 %>%
    filter(thr.imp == 0.3, thr.r2 == 0.2, size == 500, thr.lp <= 8) %>%
    arrange(desc(auc)) %>%
    slice(1) %>%
    print()
  #   size thr.r2 thr.imp num    thr.lp       auc
  # 1  500    0.2     0.3  14 0.2399897 0.6848973

  ind.keep <- unlist(map(all_keep, std_prs$num))
  # # Verif on training set
  # AUC(
  #   snp_PRS(G.train, beta_gwas[ind.keep], ind.keep = ind.keep,
  #           lpS.keep = lpval[ind.keep], thr.list = std_prs$thr.lp),
  #   y[ind.train]
  # ) # 0.6848973
  # Eval on test set
  pred_std_prs <- snp_PRS(G.test, beta_gwas[ind.keep], ind.keep = ind.keep,
                          lpS.keep = lpval[ind.keep], thr.list = std_prs$thr.lp)
  # AUCBoot(pred_std_prs, y[ind.test])  # 69.9 [68.2-71.5]
  (nb_std_prs <- sum(lpval[ind.keep] > std_prs$thr.lp))  # 182,385

  max_prs <- grid2 %>% arrange(desc(auc)) %>% slice(1:10) %>% print() %>% slice(1)
  #    size thr.r2 thr.imp num    thr.lp       auc
  # 1  2500    0.2    0.60  44 0.2399897 0.6873083
  # 2  2500    0.2    0.60  44 0.2151112 0.6871044
  # 3  2500    0.2    0.60  44 0.2677454 0.6868253
  # 4  2500    0.2    0.30  16 0.2399897 0.6867206
  # 5  2500    0.2    0.60  44 0.1928118 0.6866413
  # 6  2500    0.2    0.95 100 0.7169479 0.6865654
  # 7  2500    0.2    0.30  16 0.2151112 0.6865352
  # 8  1000    0.2    0.60  43 0.2399897 0.6865244
  # 9  2500    0.2    0.60  44 0.1728240 0.6865223
  # 10 2500    0.2    0.60  44 0.1549083 0.6863968

  ind.keep2 <- unlist(map(all_keep, max_prs$num))
  pred_max_prs <- snp_PRS(G.test, beta_gwas[ind.keep2], ind.keep = ind.keep2,
                          lpS.keep = lpval[ind.keep2], thr.list = max_prs$thr.lp)
  # AUCBoot(pred_max_prs, y[ind.test])  # 70.0 [68.4-71.7]
  (nb_max_prs <- sum(lpval[ind.keep2] > max_prs$thr.lp))  # 175,433

  # ggplot(grid2) +
  #   geom_point(aes(thr.lp, auc)) +
  #   facet_grid(thr.imp ~ thr.r2 + size) +
  #   scale_x_log10(limits = c(0.1, 3), breaks = c(0.1, 1)) +
  #   ylim(0.65, NA) +
  #   # ylim(0.55, NA) +
  #   theme_bigstatsr(size.rel = 0.7) +
  #   labs(x = "-log10(p-value) threshold (log scale)", y = "AUC")


  #### lassosum ####
  # devtools::install_github("tshmak/lassosum")
  library(lassosum)
  library(doParallel)
  registerDoParallel(cl <- makeCluster(NCORES))
  system.time(
    out <- lassosum.pipeline(
      cor = gwas$z / sqrt(Nss),
      snp = simu$map$marker.ID,
      A1 = simu$map$allele1,
      test.bfile = "data/ukbb4simu_train",
      LDblocks = "EUR.hg19",
      cluster = cl,
      exclude.ambiguous = FALSE
    )
  ) # 19 min
  stopCluster(cl)

  v <- validate(out, pheno = y[ind.train], validate.function = AUC)
  nb_lassosum <- length(ind <- which(v$best.beta != 0))
  pred_lassosum <- big_prodVec(G.test, v$best.beta[ind], ind.col = ind)
  # AUCBoot(pred_lassosum, y[ind.test])


  #### LDpred ####
  file_sumstats <- paste0(tmp, ".txt")
  mutate(simu$map, chromosome = as.integer(chromosome),
         beta = beta_gwas, pval = 10^-lpval) %>%
    bigreadr::fwrite2(file_sumstats, sep = "\t")
  file_hdf5 <- paste0(tmp, ".hdf5")

  reticulate::use_python("/opt/rh/rh-python36/root/usr/bin/python")
  reticulate::py_config()
  # system("python3 --version", intern = TRUE)
  ldpred <- "../ldpred/LDpred.py"
  # system(glue::glue("python3 {ldpred} coord --help"))
  system(glue::glue(
    "python3 {ldpred} coord",
    " --gf data/ukbb4simu_train",
    " --ssf {file_sumstats}",
    " --skip-coordination",
    " --rs marker.ID --A1 allele1 --A2 allele2 --pos physical.pos --chr chromosome",
    " --pval pval --eff beta --beta",
    " --N {Nss}",
    " --out {file_hdf5}"
  )) # 21 min


  system.time(
    system(glue::glue(
      "python3 {ldpred} gibbs",
      " --cf {file_hdf5}",
      " --ldr {round(nrow(simu$map) / 3000)}",
      " --ldf {tmp}",
      " --N {Nss}",
      " --out {tmp}"
    ))
  ) # 127 min

  ext <- c(sprintf("_LDpred_p%.4e.txt", c(1, 0.3, 0.1, 0.03, 0.01, 0.003, 0.001)),
           "_LDpred-inf.txt")
  files_ldpred <- paste0(tmp, ext)
  beta_ldpred <- sapply(files_ldpred, function(file) {
    res_ldpred <- bigreadr::fread2(file, select = c(3, 7))
    beta_ldpred <- numeric(nrow(simu$map))
    beta_ldpred[match(res_ldpred$sid, simu$map$marker.ID)] <- res_ldpred[[2]]
    beta_ldpred
  })

  pred_train_ldpred <- big_prodMat(G.train, beta_ldpred, ncores = NCORES)
  auc_train_ldpred <- apply(-pred_train_ldpred, 2, AUC, target = y[ind.train])

  beta_ldpred_max <- beta_ldpred[, which.max(auc_train_ldpred)]
  nb_ldpred <- sum(beta_ldpred_max != 0)
  pred_ldpred <- big_prodVec(G.test, beta_ldpred_max)

  saveRDS(
    data.frame(
      auc_std_prs  = AUC(pred_std_prs,  y[ind.test]), nb_std_prs,
      auc_max_prs  = AUC(pred_max_prs,  y[ind.test]), nb_max_prs,
      auc_SCT      = AUC(pred,          y[ind.test]), nb_SCT,
      auc_lassosum = AUC(pred_lassosum, y[ind.test]), nb_lassosum,
      auc_ldpred   = AUC(-pred_ldpred,  y[ind.test]), nb_ldpred
    ),
    res_file
  )
  stopifnot(file.exists(res_file))

  unlink(paste0(tmp, "*"))

}

library(tidyverse)
list.files("res_simu", "1M_.+\\.rds$", full.names = TRUE) %>%
  map_dfr(readRDS) %>%
  print() %>%
  map_chr(~ {
    boot <- replicate(1e5, mean(sample(.x, replace = TRUE)))
    res <- signif(c(mean(boot), quantile(boot, c(0.025, 0.975))), 3)
    sprintf("%s [%s-%s]", res[1], res[2], res[3])
  }) %>%
  matrix(nrow = 2) %>%
  as.data.frame()

#    auc_std_prs nb_std_prs auc_max_prs nb_max_prs   auc_SCT nb_SCT auc_lassosum nb_lassosum auc_ldpred nb_ldpred
# 1    0.6750174     221452   0.6701758      97162 0.6832684 539697    0.6883596      546407  0.6951987    973801
# 2    0.6824284     148603   0.6834074     290575 0.6852392 558770    0.6982023      548701  0.7055711    973801
# 3    0.6972392     164720   0.7034971     317783 0.7039489 555401    0.7129182      547374  0.7215522    973801
# 4    0.6863934     154821   0.6873849     149552 0.6930221 560894    0.7019171      546225  0.7110845    973801
# 5    0.6986408     171588   0.7050889     328098 0.6994093 567386    0.7097006      546391  0.7224115    973801
# 6    0.6717313     127127   0.6748600      97052 0.6767456 559456    0.6816911      547113  0.6888234    973801
# 7    0.6825925     108826   0.6862912      90597 0.6824113 560124    0.6949600      547897  0.7043628    973801
# 8    0.6912262     221330   0.6948415     283566 0.6933723 566226    0.7016621      548215  0.7132446    973801
# 9    0.6921828     172956   0.6915479     170863 0.6962412 543071    0.7078385      544379  0.7189650    973801
# 10   0.6911168     155977   0.6985737     331550 0.6924676 557465    0.7089522      547369  0.7147439    973801

# 1    0.687 [0.681-0.692]     0.69 [0.683-0.696]    0.691 [0.686-0.696]    0.701 [0.694-0.706]     0.71 [0.703-0.716]
# 2 165000 [144000-186000] 216000 [154000-276000] 557000 [551000-562000] 547000 [546000-548000] 974000 [974000-974000]
