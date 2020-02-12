closest_dist <- function(from, to) {
  drop(nabor::knn(as.matrix(to), as.matrix(from), k = 1)$nn.dist)
}

intersect_all <- function(l) Reduce(intersect, l)

roll_mean_sq <- function(x, size) sqrt(bigutilsr::rollmean(x ** 2, size))


merged <- readRDS("tmp-data/sumstats_merged.rds")
str(merged)

merged$qual <- qual <- rowMeans(sapply(merged[c(5, 9:16)], function(x) {
  ifelse(is.na(x), 0, pmin(pmax(0, x), 1))
}))
hist(qual)

subset_ukbb <- which(merged$info_ukbb >= 0.8)

s <- lapply(merged[c(5, 9:16)], function(x) which(x >= 0.8))
str(inter <- intersect_all(s)) # 2,545,046
table(merged$chr[inter])
#      1      2      3      4      5      6      7      8      9     10     11
# 250220 297503 249725 238833   1343 216944 184877 199240 141400  61061 171654
#     12     13     14     15     16     17     18     19     20     21     22
#    674 132535   2693  87587  23372  63422  88688  34740  49693  40150   8692

map_hapmap3 <- bigreadr::fread2(
  "ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/hapmap3/plink_format/draft_2/hapmap3_r2_b36_fwd.consensus.qc.poly.map.bz2")
inter2 <- which(merged$rsID %in% map_hapmap3[[2]]) # 1,336,684
table(merged$chr[inter2])
#      1      2      3      4      5      6      7      8      9     10     11
# 111446 112303  93286  83069  85011  88141  72569  72867  61097  71460  68313
#     12     13     14     15     16     17     18     19     20     21     22
#  65954  50151  43916  40762  42695  36379  39690  24620  34978  18779  19198


pos_seq <- dist_ukbb <- dist_inter <- dist_hm3 <- vector("list", 22)

for (chr in 1:22) {
  cat(chr, "\n")
  ind_chr <- which(merged$chr == chr)
  pos_chr <- merged$pos[ind_chr]
  pos_seq[[chr]] <- pos_seq_chr <- seq(min(pos_chr), max(pos_chr), by = 1000L)
  dist_ukbb[[chr]]  <- closest_dist(
    from = pos_seq_chr, to = merged$pos[intersect(ind_chr, subset_ukbb)])
  dist_inter[[chr]] <- closest_dist(
    from = pos_seq_chr, to = merged$pos[intersect(ind_chr, inter)])
  dist_hm3[[chr]]   <- closest_dist(
    from = pos_seq_chr, to = merged$pos[intersect(ind_chr, inter2)])
}

(list_var <- tibble::tibble(chr = 1:22, pos_seq, dist_ukbb, dist_inter, dist_hm3))

library(ggplot2)
ggplot(subset(tidyr::unnest(list_var))) +
  bigstatsr::theme_bigstatsr(0.7) +
  geom_line(aes(pos_seq, roll_mean_sq(dist_ukbb,  50))) +
  geom_line(aes(pos_seq, roll_mean_sq(dist_inter, 50)), color = "green") +
  geom_line(aes(pos_seq, roll_mean_sq(dist_hm3,   50)), color = "red") +
  scale_y_log10() +
  facet_wrap(~chr, scales = "free") +
  labs(x = "BP position", y = "Distance to closest (smoothed)")


print(log(sum(unlist(dist_ukbb)  ** 2))) +
  print(log(1 + sum(1 - qual[subset_ukbb])))  # 41.8 + 15.8 -> 57.6
print(log(sum(unlist(dist_inter) ** 2))) +
  print(log(1 + sum(1 - qual[inter])))        # 49.5 + 11.7 -> 61.2
print(log(sum(unlist(dist_hm3)   ** 2))) +
  print(log(1 + sum(1 - qual[inter2])))       # 42.0 + 12.6 -> 54.5

mean(qual[subset_ukbb])
mean(qual[inter])
mean(qual[inter2])
