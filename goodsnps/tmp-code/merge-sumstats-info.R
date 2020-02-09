wd <- setwd("tmp-data")

library(bigreadr)
library(data.table)

#### UKBB ####
## Download and filter data
# sumstats_ukbb <- do.call("rbind", lapply(1:22, function(chr) {
#   url <- "http://biobank.ctsu.ox.ac.uk/crystal/crystal/auxdata/"
#   mfi <- fread2(paste0(url, "ukb_mfi_chr", chr, "_v3.txt"),
#                 select = c(3, 2, 6, 8),
#                 col.names = c("pos", "rsID", "maf", "info_ukbb"))
#   cbind(chr = chr, subset(mfi, maf > 0.001))
# }))
## Provide positions in different builds
# sumstats_ukbb$pos37 <- sumstats_ukbb$pos
# sumstats_ukbb <- bigsnpr::snp_modifyBuild(sumstats_ukbb, "./liftOver",
#                                           from = "hg19", to = "hg18")
# sumstats_ukbb$pos36 <- sumstats_ukbb$pos
# sumstats_ukbb$pos <- sumstats_ukbb$pos37
# sumstats_ukbb <- bigsnpr::snp_modifyBuild(sumstats_ukbb, "./liftOver",
#                                           from = "hg19", to = "hg38")
# sumstats_ukbb$pos38 <- sumstats_ukbb$pos
# sumstats_ukbb$pos <- sumstats_ukbb$pos37
## Save to compute only once
# saveRDS(sumstats_ukbb, "sumstats_UKBB.rds")
sumstats_ukbb <- readRDS("sumstats_UKBB.rds")
str(sumstats_ukbb)

merged <- as.data.table(sumstats_ukbb)
setkey(merged, chr, pos)

#### BRCA ####
# download.file("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/MichailidouK_29059683_GCST004988/oncoarray_bcac_public_release_oct17%20(1).txt.gz",
#               destfile = "sumstats_BRCA.txt.gz")
sumstats_brca <- fread2(
  "sumstats_BRCA.txt.gz", na.strings = "NULL",
  select = c("chr", "position_b37", "bcac_onco2_r2", "bcac_icogs2_r2"),
  col.names = c("chr", "pos", "info_onco_brca", "info_icogs_brca"))
str(sumstats_brca)

merged <- merge(merged, sumstats_brca, all = TRUE)

#### PRCA ####
# download.file("http://practical.icr.ac.uk/blog/wp-content/uploads/uploadedfiles/oncoarray/MetaSummaryData/meta_v3_onco_euro_overall_ChrAll_1_release.zip",
#               destfile = "sumstats_PRCA.zip")
# unzip("sumstats_PRCA.zip"); file.remove("sumstats_PRCA.zip")
sumstats_prca <- fread2(
  "meta_v3_onco_euro_overall_ChrAll_1_release.txt",
  select = c("Chr", "position", "OncoArray_imputation_r2"),
  col.names = c("chr", "pos", "info_onco_prca")
)
str(sumstats_prca)

merged <- merge(merged, sumstats_prca, all = TRUE)

#### CAD ####
# download.file("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/NikpayM_26343387_GCST003116/cad.add.160614.website.txt",
#               destfile = "sumstats_CAD.txt")
sumstats_cad <- fread2("sumstats_CAD.txt",
                       select = c("chr", "bp_hg19", "median_info"),
                       col.names = c("chr", "pos", "info_cad"))

merged <- merge(merged, sumstats_cad, all = TRUE)

#### T1D ####
# https://datadryad.org//resource/doi:10.5061/dryad.ns8q3
# untar("65009084.tar.gz", exdir = "T1D"); file.remove("65009084.tar.gz")
sumstats_t1d <- fread2(paste0("T1D/meta_chr_", 1:22),
                       select = c("chromosome", "position", "info_score.I", "info_score.A"),
                       col.names = c("chr", "pos", "info_illu_t1d", "info_affy_t1d"))

merged <- merge(merged, sumstats_t1d, all = TRUE)

#### ADHD ####
# https://ipsych.dk/en/research/downloads/data-download-agreement-adhd-european-ancestry-gwas-june-2017/
sumstats_adhd <- fread2("adhd_eur_jun2017.gz",
                        select = c("CHR", "BP", "INFO"),
                        col.names = c("chr", "pos", "info_adhd"))

merged <- merge(merged, sumstats_adhd, all = TRUE)

#### IHPS ####
# download.file("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/FadistaJ_30281099_GCST006705/MEGA_CIDR_IHPS_summaryStats.txt.gz",
#               destfile = "sumstats_IHPS.txt.gz")
sumstats_ihps <- fread2("sumstats_IHPS.txt.gz",
                        select = c("CHR", "BP"),
                        col.names = c("chr", "pos"))
sumstats_ihps$info_ihps <- 0.8

merged <- merge(merged, sumstats_ihps, all = TRUE)

## Remove multiallelic and sort by chr/pos
merged <- dplyr::arrange(merged, chr, pos)
merged <- merged[!vctrs::vec_duplicate_detect(merged[, c("chr", "pos")]), ]
str(merged)

## Add Berisa LD blocks
library(bigsnpr)
LD_blocks <- fread2("https://bitbucket.org/nygcresearch/ldetect-data/raw/ac125e47bf7ff3e90be31f278a7b6a61daaba0dc/EUR/fourier_ls-all.bed")
LD_blocks$chr <- as.integer(sub("^chr([0-9]+)$", "\\1", LD_blocks$chr))

merged$block <- NA_integer_
for (i in 1:nrow(LD_blocks)) {
  cat(i, "")
  merged[chr == LD_blocks[i, "chr"] & pos %between% LD_blocks[i, 2:3], block := i]
}
sort(table(merged$block, exclude = NULL))
sort(table(merged$block, exclude = NULL), decreasing = TRUE)


# saveRDS(as.data.frame(merged[!is.na(block)]), "sumstats_merged.rds")

setwd(wd)
