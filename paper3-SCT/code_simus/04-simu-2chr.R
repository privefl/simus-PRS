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

ind1 <- which(CHR == 1)
ind2 <- which(CHR == 2)

for (ic in 1:10) {

  res_file <- paste0("res_simu/2chr_", ic, ".rds")
  if (file.exists(res_file)) next
  cat(ic, "\n")
  tmp <- tempfile(tmpdir = "res_simu")

  #### Simulate phenotype ####
  set.seed(ic)
  h2 <- 0.5; K <- 0.1; M1 <- 100; M2 <- length(ind2)
  set1 <- sort(sample(ind1, size = M1))
  set2 <- ind2
  effects1 <- rnorm(M1, sd = sqrt(h2 / 2 / M1))
  effects2 <- rnorm(M2, sd = sqrt(h2 / 2 / M2))
  y <- big_prodVec(G, effects1 / SD[set1], ind.col = set1) +
    big_prodVec(G, effects2 / SD[set2], ind.col = set2)
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
  ) # 5 min
  mod <- final_mod$mod
  # summary(mod)
  # plot(mod)

  new_beta <- final_mod$beta.G
  nb_SCT <- length(ind <- which(new_beta != 0))  # 179,914
  table(CHR[ind])
  # 1     2     3     4     5     6     7     8     9    10    11
  # 7537 50331 10699  9270  2383  8414  7924 29429  5758  8389   312
  # 12    13    14    15    16    17    18    19    20    21    22
  # 2452 17150  1595  1620  3600  2055  1095  1549  1408  5880  1064

  pred <- final_mod$intercept +
    big_prodVec(G.test, new_beta[ind], ind.col = ind)

  # AUCBoot(pred, y[ind.test])  # 83.1 [81.9-84.3]


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
  #   size thr.r2 thr.imp num   thr.lp       auc
  # 1  500    0.2     0.3  14 4.091176 0.7731041

  ind.keep <- unlist(map(all_keep, std_prs$num))
  # # Verif on training set
  # AUC(
  #   snp_PRS(G.train, beta_gwas[ind.keep], ind.keep = ind.keep,
  #           lpS.keep = lpval[ind.keep], thr.list = std_prs$thr.lp),
  #   y[ind.train]
  # ) # 0.7731041
  # Eval on test set
  pred_std_prs <- snp_PRS(G.test, beta_gwas[ind.keep], ind.keep = ind.keep,
                          lpS.keep = lpval[ind.keep], thr.list = std_prs$thr.lp)
  # AUCBoot(pred_std_prs, y[ind.test])  # 78.5 [77.1-79.9]
  (nb_std_prs <- sum(lpval[ind.keep] > std_prs$thr.lp))  # 1659

  max_prs <- grid2 %>% arrange(desc(auc)) %>% slice(1:10) %>% print() %>% slice(1)
  #     size thr.r2 thr.imp num   thr.lp       auc
  # 1   5000   0.01     0.6  29 4.091176 0.7908652
  # 2   5000   0.01     0.6  29 4.882081 0.7908485
  # 3  10000   0.01     0.6  30 4.882081 0.7907743
  # 4  20000   0.01     0.6  31 4.882081 0.7907743
  # 5  50000   0.01     0.6  32 4.882081 0.7907743
  # 6   5000   0.01     0.3   1 4.882081 0.7906991
  # 7  10000   0.01     0.6  30 4.091176 0.7906850
  # 8  20000   0.01     0.6  31 4.091176 0.7906850
  # 9  50000   0.01     0.6  32 4.091176 0.7906850
  # 10  5000   0.01     0.3   1 4.091176 0.7906720

  ind.keep2 <- unlist(map(all_keep, max_prs$num))
  pred_max_prs <- snp_PRS(G.test, beta_gwas[ind.keep2], ind.keep = ind.keep2,
                          lpS.keep = lpval[ind.keep2], thr.list = max_prs$thr.lp)
  # AUCBoot(pred_max_prs, y[ind.test])  # 79.8 [78.4-81.1]
  (nb_max_prs <- sum(lpval[ind.keep2] > max_prs$thr.lp))  # 675

  # ggplot(grid2) +
  #   geom_point(aes(thr.lp, auc)) +
  #   facet_grid(thr.imp ~ thr.r2 + size) +
  #   scale_x_log10(limits = c(1, 10), breaks = c(1, 2, 5, 10), minor_breaks = 1:10) +
  #   ylim(0.73, NA) +
  #   # ylim(0.65, NA) +
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
list.files("res_simu", "2chr_.+\\.rds$", full.names = TRUE) %>%
  map_dfr(readRDS) %>%
  print() %>%
  map_chr(~ {
    boot <- replicate(1e5, mean(sample(.x, replace = TRUE)))
    res <- signif(c(mean(boot), quantile(boot, c(0.025, 0.975))), 3)
    sprintf("%s [%s-%s]", res[1], res[2], res[3])
  }) %>%
  matrix(nrow = 2) %>%
  as.data.frame()

#    auc_std_prs nb_std_prs auc_max_prs nb_max_prs   auc_SCT nb_SCT auc_lassosum nb_lassosum
# 1    0.7856381       2420   0.7970973       1034 0.8358688 202518    0.7964701       14920
# 2    0.7746098       3699   0.7811054       1036 0.8197938 214078    0.7804348       14914
# 3    0.7766176       2578   0.7935147        388 0.8293244 157075    0.7965627       14526
# 4    0.7561126       2924   0.7632801        441 0.8077963 183277    0.7704855       15395
# 5    0.7584487       3550   0.7782242        483 0.8160370 154643    0.7846720       15121
# 6    0.7706646       3330   0.7828623        617 0.8176642 190917    0.7907374       15052
# 7    0.7710786       4001   0.7899965        703 0.8207598 139393    0.7933614       14982
# 8    0.7749393       2401   0.7876626        644 0.8221962 194422    0.7951758       15378
# 9    0.7779502       1854   0.7983809        358 0.8301890 174778    0.7970426       14188
# 10   0.7780627       2691   0.7886458        729 0.8241376 188311    0.7881987       15075

# 1 0.772 [0.767-0.777] 0.786 [0.78-0.792]    0.822 [0.818-0.827] 0.789 [0.784-0.794]
# 2    2950 [2550-3350]      643 [508-792] 180000 [166000-193000] 15000 [14700-15200]
