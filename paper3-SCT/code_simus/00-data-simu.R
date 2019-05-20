# Get variant info
library(doParallel)
registerDoParallel(cl <- makeCluster(22))
info_all_chr <- foreach(chr = 1:22) %dopar% {
  file <- paste0("data/ukb_imp_mfi/ukb_mfi_chr", chr, "_v3.txt")
  info_chr <- bigreadr::fread2(file, nThread = 1)
  info_chr_sub <- subset(info_chr, V6 > 0.01 & V8 > 0.3, -c(V1, V2, V7))
  cbind.data.frame(chr = chr, setNames(info_chr_sub, c("pos", "a1", "a2", "maf", "info")))
}
stopCluster(cl)

# Keep 1M variants
nb_snp <- sum(sapply(info_all_chr, nrow)) # 9,926,099
set.seed(1); ind_snp <- sort(sample(nb_snp, 1e6))
info_final <- do.call(rbind, info_all_chr)[ind_snp, ]
stopifnot(nrow(info_final) == 1e6)
head(info_final)
saveRDS(info_final, "data/ukbb4simu_info.rds")


# system("./ukbgene imp -c1 -m")
library(bigreadr)
sample <- fread2("ukb25589_imp_chr1_v3_s487327.sample")
str(sample)
(N <- readBin("data/ukb_imp_chr1_v3.bgen", what = 1L, size = 4, n = 4)[4])
sample <- sample[-1, ]
stopifnot(nrow(sample) == N)

csv <- "ukb22544.csv"
df0 <- fread2(csv, select = c("eid", "22006-0.0"),
              col.names = c("eid", "is_caucasian"))
# sample still in data
ind.indiv <- match(df0$eid, sample$ID_2)
# related samples
rel <- fread2("ukb25589_rel_s488346.dat")
df0$is_rel2 <- df0$eid %in% rel$ID2
# + keep caucasian only
sub <- which(!is.na(ind.indiv) & !df0$is_rel2 & df0$is_caucasian)
length(sub)  # 335,609

nPC <- 10
PC <- fread2(csv, select = paste0("22009-0.", 1:nPC),
             col.names = paste0("PC", 1:nPC))
saveRDS(as.matrix(PC[sub, ]), "PC_sub.rds")

list_snp_id <- with(info_final, split(paste(chr, pos, a1, a2, sep = "_"), chr))
stopifnot(sum(lengths(list_snp_id)) == 1e6)

library(bigsnpr)
# Read all data as random hard calls
system.time(
  snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/ukbb4simu",
    ind_row = ind.indiv[sub],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = 10,
    read_as = "random"
  )
) # 4H

file.size("data/ukbb4simu.bk") / 1024^3  # 313 GB
simu <- snp_attach("data/ukbb4simu.rds")
G <- simu$genotypes
G[, 1]
G[, ncol(G)]
stats <- big_scale()(G)
saveRDS(stats, "data/ukbb4simu_stats.rds")

set.seed(1)
ind.train <- sort(sample(nrow(G), 10e3))
ind.test <- sort(sample(setdiff(rows_along(G), ind.train), 10e3))
ind.gwas <- setdiff(rows_along(G), c(ind.train, ind.test))
save(ind.gwas, ind.train, ind.test, file = "data/ukbb4simu_ind.RData")

# Read train data as dosages
system.time(
  snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/ukbb4simu_train",
    ind_row = ind.indiv[sub][ind.train],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = 10
  )
) # 3.5H

# Read test data as dosages
system.time(
  snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/ukbb4simu_test",
    ind_row = ind.indiv[sub][ind.test],
    bgi_dir = "data/ukb_imp_bgi",
    ncores = 10
  )
) # 3.3H

# Write bed/bim/fam for LDpred and lassosum
library(dplyr)
train <- snp_attach("data/ukbb4simu_train.rds")
train$map <- mutate(train$map, genetic.dist = 0, rsid = NULL,
                    chromosome = as.integer(chromosome))
train$fam <- snp_fake(nrow(G), 1)$fam
snp_writeBed(train, bedfile = "data/ukbb4simu_train.bed")
test <- snp_attach("data/ukbb4simu_test.rds")
test$map <- mutate(test$map, genetic.dist = 0, rsid = NULL,
                   chromosome = as.integer(chromosome))
test$fam <- snp_fake(nrow(G), 1)$fam
snp_writeBed(test, bedfile = "data/ukbb4simu_test.bed")

# Verif GWAS
ind_first <- seq_len(20e3)
# devtools::install_github("privefl/paper2-PRS/pkg.paper.PRS")
pheno <- pkg.paper.PRS::get_pheno(G, h2 = 0.5, M = 100, K = 0.1,
                                  ind.possible = ind_first)
y <- pheno$pheno
PC <- readRDS("PC_sub.rds")
system.time(
  gwas <- big_univLogReg(G, y[ind.gwas], ind.gwas, covar.train = PC[ind.gwas, ],
                         ind.col = ind_first, ncores = nb_cores())
) # 27 min

system.time(
  gwas2 <- big_apply(G, function(G, ind, ind.ca, ind.co) {

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
  }, a.combine = "rbind", ind = ind_first, ncores = nb_cores(),
  ind.ca = intersect(which(y == 1), ind.gwas),
  ind.co = intersect(which(y == 0), ind.gwas))
) # 18 sec

plot(gwas$estim, log(gwas2$or), pch = 20); abline(0, 1, col = "red")
all.equal(gwas$estim, log(gwas2$or))  # Mean relative difference: 0.01833519
summary(gwas2$z)
plot(gwas$score, gwas2$z, pch = 20); abline(0, 1, col = "red")
all.equal(gwas$score, gwas2$z)  # Mean relative difference: 0.02131113

library(ggplot2)
p1 <- ggplot(cbind(gwas, gwas2)) +
  theme_bigstatsr() +
  geom_point(aes(estim, log(or))) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(x = "Effect sizes from logistic regression with 10 PCs",
       y = "Effect sizes from Cochran-Armitage additive test") +
  ggtitle(all.equal(gwas$estim, log(gwas2$or)))

p2 <- ggplot(cbind(gwas, gwas2)) +
  theme_bigstatsr() +
  geom_point(aes(score, z)) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(x = "Z-scores from logistic regression with 10 PCs",
       y = "Z-scores from Cochran-Armitage additive test") +
  ggtitle(all.equal(gwas$score, gwas2$z))

cowplot::plot_grid(p1, p2, scale = 0.95, labels = LETTERS[1:2], label_size = 25)
