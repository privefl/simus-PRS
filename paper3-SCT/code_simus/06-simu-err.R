library(bigsnpr)
simu <- snp_attach("data/ukbb4simu.rds")
G <- simu$genotypes
CHR <- as.integer(simu$map$chromosome)
POS <- simu$map$physical.pos
INFO <- readRDS("data/ukbb4simu_info.rds")$info
SD <- readRDS("data/ukbb4simu_stats.rds")$scale
load("data/ukbb4simu_ind.RData")
NCORES <- nb_cores()

G.train <- snp_attach("data/ukbb4simu_train.rds")$genotypes
G.test  <- snp_attach("data/ukbb4simu_test.rds")$genotypes

for (ic in 1:10) {

  res_file <- paste0("res_simu/err_", ic, ".rds")
  if (file.exists(res_file)) next
  cat(ic, "\n")

  set.seed(ic)
  h2 <- 0.5; K <- 0.1; M <- 10e3
  set <- sort(sample(ncol(G), size = M))
  effects <- rnorm(M, sd = sqrt(h2 / M))
  y <- big_prodVec(G, effects / SD[set], ind.col = set)
  y <- (y - mean(y)) / sd(y) * sqrt(h2)         ## make sure that var(y) = h2
  y <- y + rnorm(nrow(G), sd = sqrt(1 - h2))
  y <- as.integer(y > qnorm(1 - K))             ## LTM

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

  beta_gwas0 <- log(gwas$or)
  err <- as.logical(rbinom(length(beta_gwas0), 1, 0.1))
  beta_gwas <- ifelse(err, -beta_gwas0, beta_gwas0)
  # hist(beta_gwas[beta_gwas > -0.1 & beta_gwas < 0.1])
  lpval <- -pchisq(gwas$z^2, df = 1, lower.tail = FALSE, log.p = TRUE) / log(10)
  # hist(10^(-lpval))

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
      backingfile = sub("\\.rds$", "_scores", res_file),
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
  nb_SCT <- length(ind <- which(new_beta != 0))  # 556,093

  pred <- final_mod$intercept +
    big_prodVec(G.test, new_beta[ind], ind.col = ind)

  # AUCBoot(pred, y[ind.test])  # 74.5 [73.0-76.0]


  library(tidyverse)

  ind2 <- sort(sample(ind, 10e3))
  ggplot(data.frame(y = new_beta, x = beta_gwas)[ind, ]) +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    geom_abline(slope = 0, intercept = 0, color = "blue") +
    geom_point(aes(x, y), size = 0.6) +
    theme_bigstatsr() +
    labs(x = "Effect sizes from GWAS", y = "Non-zero effect sizes from SCT")

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
  # 1  500    0.2     0.3  14 3.340005 0.7210221

  ind.keep <- unlist(map(all_keep, std_prs$num))
  # # Verif on training set
  # AUC(
  #   snp_PRS(G.train, beta_gwas[ind.keep], ind.keep = ind.keep,
  #           lpS.keep = lpval[ind.keep], thr.list = std_prs$thr.lp),
  #   y[ind.train]
  # ) # 0.7210221
  # Eval on test set
  pred_std_prs <- snp_PRS(G.test, beta_gwas[ind.keep], ind.keep = ind.keep,
                          lpS.keep = lpval[ind.keep], thr.list = std_prs$thr.lp)
  # AUCBoot(pred_std_prs, y[ind.test])  # 71.6 [69.9-73.1]
  (nb_std_prs <- sum(lpval[ind.keep] > std_prs$thr.lp))  # 3256

  max_prs <- grid2 %>% arrange(desc(auc)) %>% slice(1:10) %>% print() %>% slice(1)
  #     size thr.r2 thr.imp num   thr.lp       auc
  # 1  1000    0.1     0.6  38 3.340005 0.7249191
  # 2  1000    0.1     0.3  10 3.340005 0.7248890
  # 3  2000    0.1     0.6  39 3.340005 0.7239040
  # 4  5000    0.1     0.6  40 3.340005 0.7238524
  # 5  5000    0.1     0.3  12 3.340005 0.7238349
  # 6  2000    0.1     0.3  11 3.340005 0.7237971
  # 7  1000    0.2     0.3  15 3.340005 0.7233236
  # 8  1000    0.2     0.6  43 3.340005 0.7232818
  # 9   500    0.1     0.3   9 3.340005 0.7229422
  # 10  500    0.1     0.6  37 3.340005 0.7229091

  ind.keep2 <- unlist(map(all_keep, max_prs$num))
  pred_max_prs <- snp_PRS(G.test, beta_gwas[ind.keep2], ind.keep = ind.keep2,
                          lpS.keep = lpval[ind.keep2], thr.list = max_prs$thr.lp)
  # AUCBoot(pred_max_prs, y[ind.test])  # 71.5 [69.9-73.1]
  (nb_max_prs <- sum(lpval[ind.keep2] > max_prs$thr.lp))  # 2716

  # ggplot(grid2) +
  #   geom_point(aes(thr.lp, auc)) +
  #   facet_grid(thr.imp ~ thr.r2 + size) +
  #   scale_x_log10(limits = c(1, 6), breaks = c(1, 3, 10), minor_breaks = 1:10) +
  #   ylim(0.69, NA) +
  #   theme_bigstatsr(size.rel = 0.7) +
  #   labs(x = "-log10(p-value) threshold (log scale)", y = "AUC")


  # devtools::install_github("tshmak/lassosum")
  library(lassosum)
  library(doParallel)
  registerDoParallel(cl <- makeCluster(NCORES))
  system.time(
    out <- lassosum.pipeline(
      cor = gwas$z / sqrt(length(ind.gwas)),
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

  saveRDS(
    data.frame(
      auc_std_prs  = AUC(pred_std_prs,  y[ind.test]), nb_std_prs,
      auc_max_prs  = AUC(pred_max_prs,  y[ind.test]), nb_max_prs,
      auc_SCT      = AUC(pred,          y[ind.test]), nb_SCT,
      auc_lassosum = AUC(pred_lassosum, y[ind.test]), nb_lassosum
    ),
    res_file
  )
  stopifnot(file.exists(res_file))

  unlink(c(multi_PRS$rds, multi_PRS$bk))

}

library(tidyverse)
list.files("res_simu", "err_.+\\.rds$", full.names = TRUE) %>%
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
# 1    0.7208358       5990   0.7195996       5909 0.7427696 556873    0.7620051       20662
# 2    0.6740187      21344   0.6832889      28531 0.7173350 551776    0.7478481       58173
# 3    0.6889644      13998   0.7015143      18981 0.7187687 430862    0.7469321       57925
# 4    0.7014369       6862   0.7047542       5072 0.7312275 551608    0.7508955       42376
# 5    0.6937752       6455   0.7075190       2597 0.7314310 540406    0.7432223       42134
# 6    0.7084239      10220   0.7151113       4624 0.7419574 528332    0.7591377       41961
# 7    0.7190179       5516   0.7137930      23163 0.7469418 549308    0.7690968       42215
# 8    0.6866592      12846   0.7027708       2404 0.7163345 551635    0.7420381       42824
# 9    0.6948183       6692   0.7086798       2906 0.7446229 534640    0.7401974       42584
# 10   0.6871929      20441   0.7114027       2426 0.7279022 546629    0.7308000       42651
#
# 1 0.698 [0.689-0.706] 0.707 [0.7-0.712]    0.732 [0.725-0.739] 0.749 [0.743-0.756]
# 2  11000 [7740-14700] 9660 [4280-15900] 534000 [509000-550000] 43400 [37400-48700]
