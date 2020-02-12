# merged <- readRDS("tmp-data/sumstats_merged.rds")
# str(merged)
#
# merged$qual <- rowMeans(sapply(merged[c(5, 9:16)], function(x) {
#   ifelse(is.na(x), 0, pmin(pmax(0, x), 0.99))
# }))
#
# dist_to_pena <- 0
# i <- 1
# k <- 1
# while (i < 3e7) {
#   to_add <- k^2
#   dist_to_pena[i + 1] <- dist_to_pena[i]     + to_add
#   dist_to_pena[i + 2] <- dist_to_pena[i + 1] + to_add
#   k <- k + 1
#   i <- i + 2
# }
#
# LD_blocks <- bigreadr::fread2("https://bitbucket.org/nygcresearch/ldetect-data/raw/ac125e47bf7ff3e90be31f278a7b6a61daaba0dc/EUR/fourier_ls-all.bed")


library(bigparallelr)

CHR <- 15
ind_chr <- which(merged$chr == CHR)

merged_chr1 <- merged[ind_chr, ]

registerDoParallel(cl <- makeCluster(nb_cores()))

ind_keep <- foreach(BLOCK = unique(merged_chr1$block)) %dopar% {

  ind <- with(merged_chr1, which(block == BLOCK & qual > 0 & !is.na(maf)))

  block1 <- merged_chr1[ind, ]
  qual <- block1$qual
  pos <- block1$pos
  dist_to_next <- c(diff(pos), 1e5)
  dist_to_next_prev <- c(1e5, head(dist_to_next, -1))
  pena_to_rm <- dist_to_pena[dist_to_next + dist_to_next_prev] -
    (dist_to_pena[dist_to_next] + dist_to_pena[dist_to_next_prev])
  prev <- 0:(length(qual) - 1)
  prob_to_rm <- 1 / (pena_to_rm * -log(1 - qual))
  pena_tot <- 0
  pena_max <- 250e3 * (LD_blocks[BLOCK, "stop"] - LD_blocks[BLOCK, "start"])

  for (REP in seq_along(pos)) {

    # cat(REP, "")

    i <- sample.int(length(prob_to_rm), 1, prob = prob_to_rm, useHash = FALSE)
    pena_tot <- pena_tot + pena_to_rm[i]
    if (pena_tot > pena_max) break
    prob_to_rm[i] <- 0

    # need to update prev and next
    i_next <- i + 1
    while (prev[i_next] == 0) i_next <- i_next + 1
    stopifnot(prev[i_next] == i)
    prev[i_next] <- i_prev <- prev[i]
    prev[i] <- 0

    dist_to_next_prev[i_next] <- dist_to_next[i_prev] <- pos[i_next] - pos[i_prev]

    pena_to_rm[i_next] <- dist_to_pena[dist_to_next[i_next] + dist_to_next_prev[i_next]] -
      (dist_to_pena[dist_to_next[i_next]] + dist_to_pena[dist_to_next_prev[i_next]])
    pena_to_rm[i_prev] <- dist_to_pena[dist_to_next[i_prev] + dist_to_next_prev[i_prev]] -
      (dist_to_pena[dist_to_next[i_prev]] + dist_to_pena[dist_to_next_prev[i_prev]])
    prob_to_rm[i_next] <- 1 / (pena_to_rm[i_next] * -log(1 - qual[i_next]))
    prob_to_rm[i_prev] <- 1 / (pena_to_rm[i_prev] * -log(1 - qual[i_prev]))
  }

  ind[prob_to_rm != 0]
}

stopCluster(cl)

ind_chr1 <- ind_chr[unlist(ind_keep)]
length(ind_chr1) / nrow(merged_chr1) * nrow(merged)


closest_dist <- function(from, to) {
  drop(nabor::knn(as.matrix(to), as.matrix(from), k = 1)$nn.dist)
}

intersect_all <- function(l) Reduce(intersect, l)

roll_mean_sq <- function(x, size) sqrt(bigutilsr::rollmean(x ** 2, size))


subset_ukbb <- with(merged, which(info_ukbb >= 0.8 & chr == CHR))

map_hapmap3 <- bigreadr::fread2(
  "ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/hapmap3/plink_format/draft_2/hapmap3_r2_b36_fwd.consensus.qc.poly.map.bz2")
inter2 <- with(merged, which(rsID %in% map_hapmap3[[2]] & chr == CHR))


pos_chr <- merged$pos[ind_chr]
pos_seq_chr <- seq(min(pos_chr), max(pos_chr), by = 1000L)
dist_ukbb  <- closest_dist(
  from = pos_seq_chr, to = merged$pos[subset_ukbb])
dist_keep <- closest_dist(
  from = pos_seq_chr, to = merged$pos[ind_chr1])
dist_hm3   <- closest_dist(
  from = pos_seq_chr, to = merged$pos[inter2])

(list_var <- tibble::tibble(chr = 1, pos_seq_chr, dist_ukbb, dist_keep, dist_hm3))

library(ggplot2)
ggplot(list_var) +
  bigstatsr::theme_bigstatsr(0.7) +
  geom_line(aes(pos_seq_chr, roll_mean_sq(dist_ukbb, 100))) +
  geom_line(aes(pos_seq_chr, roll_mean_sq(dist_keep, 100)), color = "green") +
  geom_line(aes(pos_seq_chr, roll_mean_sq(dist_hm3,  100)), color = "red") +
  scale_y_log10() +
  labs(x = "BP position", y = "Distance to closest (smoothed)")

qual <- merged$qual

print(log(sum(dist_ukbb ** 2))) +
  print(log(1 + sum(1 - qual[subset_ukbb], na.rm = TRUE)))  # 41.2 + 13.2 -> 54.4
print(log(sum(dist_keep ** 2))) +
  print(log(1 + sum(1 - qual[ind_chr1], na.rm = TRUE)))     # 41.2 + 10.2 -> 51.4
print(log(sum(dist_hm3  ** 2))) +
  print(log(1 + sum(1 - qual[inter2], na.rm = TRUE)))       # 41.4 + 10.0 -> 51.4

mean(qual[subset_ukbb], na.rm = TRUE)  # 0.553
mean(qual[ind_chr1],    na.rm = TRUE)  # 0.842
mean(qual[inter2],      na.rm = TRUE)  # 0.798 -> but higher maf?


merged$grp_maf <- cut(merged$maf, c(0, 0.02, 0.1, 0.25, 0.5))

library(ggplot2)
ggplot() +
  bigstatsr::theme_bigstatsr() +
  geom_histogram(aes(qual, y = ..density..), data = merged[subset_ukbb, ],
                 breaks = seq(0, 1, by = 0.01),
                 color = "black", fill = "black", size = 1.2, alpha = 0.2) +
  geom_histogram(aes(qual, y = ..density..), data = merged[inter2, ],
                 breaks = seq(0, 1, by = 0.01),
                 color = "red", fill = "red", size = 1.2, alpha = 0.2) +
  geom_histogram(aes(qual, y = ..density..), data = merged[ind_chr1, ],
                 breaks = seq(0, 1, by = 0.01),
                 color = "green", fill = "green", size = 1.2, alpha = 0.2) +
  facet_wrap(~ grp_maf) +
  xlim(0.6, 1) +
  labs(x = "Mean INFO score (0 if missing)")
