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

indHLA <- snp_indLRLDR(CHR, POS, subset(LD.wiki34, ID == "hild12"))
(M <- length(indHLA))  # 7105

for (ic in 1:10) {

  res_file <- paste0("res_simu/HLA_", ic, ".rds")
  if (file.exists(res_file)) next
  cat(ic, "\n")
  tmp <- tempfile(tmpdir = "res_simu")

  #### Simulate phenotype ####
  set.seed(ic)
  h2 <- 0.5; K <- 0.1
  set <- indHLA
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
  ) # 6 min
  mod <- final_mod$mod
  # summary(mod)
  # plot(mod)

  new_beta <- final_mod$beta.G
  nb_SCT <- length(ind <- which(new_beta != 0))  # 242,895

  pred <- final_mod$intercept +
    big_prodVec(G.test, new_beta[ind], ind.col = ind)

  # AUCBoot(pred, y[ind.test])  # 81.0 [79.7-82.3]


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
  #   size thr.r2 thr.imp num  thr.lp       auc
  # 1  500    0.2     0.3  14 7.89974 0.7895099

  ind.keep <- unlist(map(all_keep, std_prs$num))
  # # Verif on training set
  # AUC(
  #   snp_PRS(G.train, beta_gwas[ind.keep], ind.keep = ind.keep,
  #           lpS.keep = lpval[ind.keep], thr.list = std_prs$thr.lp),
  #   y[ind.train]
  # ) # 0.7895099
  # Eval on test set
  pred_std_prs <- snp_PRS(G.test, beta_gwas[ind.keep], ind.keep = ind.keep,
                          lpS.keep = lpval[ind.keep], thr.list = std_prs$thr.lp)
  # AUCBoot(pred_std_prs, y[ind.test])  # 79.0 [77.6-80.3]
  (nb_std_prs <- sum(lpval[ind.keep] > std_prs$thr.lp))  # 536

  max_prs <- grid2 %>% arrange(desc(auc)) %>% slice(1:10) %>% print() %>% slice(1)
  #    size thr.r2 thr.imp num   thr.lp       auc
  # 1  5000    0.1     0.6  40 51.39182 0.7990459
  # 2  2000    0.1     0.6  39 51.39182 0.7990343
  # 3  5000    0.1     0.6  40 41.73779 0.7988913
  # 4  2000    0.1     0.6  39 41.73779 0.7988509
  # 5  2000    0.1     0.6  39 33.89728 0.7987841
  # 6  2000    0.1     0.3  11 33.89728 0.7986184
  # 7  5000    0.1     0.6  40 33.89728 0.7984237
  # 8  2000    0.1     0.3  11 51.39182 0.7983256
  # 9  5000    0.1     0.3  12 41.73779 0.7983091
  # 10 5000    0.1     0.3  12 51.39182 0.7982689

  ind.keep2 <- unlist(map(all_keep, max_prs$num))
  pred_max_prs <- snp_PRS(G.test, beta_gwas[ind.keep2], ind.keep = ind.keep2,
                          lpS.keep = lpval[ind.keep2], thr.list = max_prs$thr.lp)
  # AUCBoot(pred_max_prs, y[ind.test])  # 79.9 [78.6-81.2]
  (nb_max_prs <- sum(lpval[ind.keep2] > max_prs$thr.lp))  # 139

  # ggplot(grid2) +
  #   geom_point(aes(thr.lp, auc)) +
  #   facet_grid(thr.imp ~ thr.r2 + size) +
  #   scale_x_log10(limits = c(1, 100), breaks = c(1, 10, 80), minor_breaks = 10 * 1:10) +
  #   ylim(0.76, NA) +
  #   # ylim(0.74, NA) +
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
list.files("res_simu", "HLA_.+\\.rds$", full.names = TRUE) %>%
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
# 1    0.7679828       2792   0.7866507        512 0.8222063 111960    0.7787193       14741  0.7484684    973796
# 2    0.7889066       4554   0.8050657        859 0.8359674 226256    0.8069494       14942  0.7513513    973799
# 3    0.7733389       1472   0.7938008        570 0.8254013 174081    0.8013143       15233  0.7472437    973791
# 4    0.7874589       1418   0.8020358        547 0.8345221 139853    0.8084008       14857  0.7468325    973796
# 5    0.7668200       3499   0.7797532        675 0.8095190 248518    0.7830479       14998  0.7489687    973729
# 6    0.7673479       3437   0.7870899        631 0.8165693 130741    0.7907637       14722  0.7913390    973772
# 7    0.7615027       6988   0.7766651        782 0.8145620  96042    0.7838149       15198  0.7319158    973800
# 8    0.7724445       3199   0.7857838        584 0.8161039 176073    0.7845332       14531  0.7274788    973791
# 9    0.7820110       2639   0.7893442        416 0.8250606 202433    0.7929387       15495  0.7340333    973791
# 10   0.7709909       3791   0.7790801        489 0.8120533 172338    0.7697474       14832  0.7277153    973800   972395

# 1 0.774 [0.769-0.78] 0.789 [0.783-0.794]    0.821 [0.816-0.827]  0.79 [0.783-0.797]    0.746 [0.736-0.758]
# 2   3380 [2510-4390]       607 [532-689] 168000 [139000-197000] 15000 [14800-15100] 974000 [974000-974000]
