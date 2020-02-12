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

which.max(LD_blocks[, 3] - LD_blocks[, 2])

BLOCK <- 3

block1 <- dplyr::filter(merged, block == BLOCK, qual > 0)
str(block1) # 14708

# ind_seq <- seq(LD_blocks[BLOCK, "start"], LD_blocks[BLOCK, "stop"])
#
# closest_dist <- function(from, to) {
#   drop(nabor::knn(as.matrix(to), as.matrix(from), k = 1)$nn.dist)
# }
#
# coverage_ind <- function(pos) {
#   log(sum(closest_dist(ind_seq, pos) ** 2))
# }
#
# coverage_ind(block1$pos)                # 31.9 -> best we can do
# coverage_ind(na.omit(block1)$pos)       # 40.1
#
# log(1 + sum(1 - block1$qual))           #  9.3 -> worst we can do
# log(1 + sum(1 - na.omit(block1)$qual))  #  4.6


qual <- block1$qual
pos <- block1$pos
dist_to_next <- c(diff(pos), 1e5)
dist_to_next_prev <- c(1e5, head(dist_to_next, -1))
pena_to_rm <- dist_to_pena[dist_to_next + dist_to_next_prev] -
  (dist_to_pena[dist_to_next] + dist_to_pena[dist_to_next_prev])
prev <- 0:(length(qual) - 1)
prob_to_rm <- 1 / (pena_to_rm * -log(1 - qual))
pena_tot <- 0
pena_max <- 100000 * (LD_blocks[BLOCK, "stop"] - LD_blocks[BLOCK, "start"])

for (REP in seq_along(pos)) {

  cat(REP, "")

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

# plot(cumsum(unlist(pena_tot)), log = "xy")
# abline(h = 10000 * pena_max, col = "red")

hist(qual, "FD")
hist(qual[prob_to_rm != 0], breaks = .Last.value$breaks,
     add = TRUE, col = scales::alpha("red", 0.2))
mean(qual)
mean(qual[prob_to_rm != 0])
mean(prob_to_rm != 0)
