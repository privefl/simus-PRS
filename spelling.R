WORDS_TO_IGNORE <- c(
  "Aarhus", "almsson", "Alpes", "Aschard", "BGEN", "bigsnpr", "bigstatsr",
  "Biobank", "Bioinformatique", "Biologie", "Biostatistique", "Bjarni", "Blum",
  "BRCA", "Centre", "chr", "CNRS", "de", "egrative", "et", "Florian", "GWAS",
  "Hugues", "IMAG", "Institut", "kb", "Laboratoire", "LD", "MDD", "phenotypes",
  "PLINK", "polygenic", "polygenicity", "PRCA", "Priv", "PRS", "PRSice",
  "readBGEN", "SCT", "snp", "th", "thr", "TIMC", "Tronche", "UKBB", "UMR", "Vilhj")

spelling::spell_check_files(
  list.files("paper3-SCT/", pattern = "\\.tex$", full.names = TRUE),
  ignore = WORDS_TO_IGNORE
)
