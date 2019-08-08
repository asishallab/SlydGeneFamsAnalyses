require(GeneFamilies)
options(mc.cores=detectCores())

message("USAGE: Rscript path/2/GeneFamilies/exec/parseCafeResult.R path/2/summary_run_2_cafe_fams.txt path/2/SlydGeneFamsAnalyses/data")

input.args <- commandArgs(trailingOnly = TRUE)

spec.names <- c("Pi", "Cag", "Dati", "Cam", "Nat", "Ntab", "Date", "Stu", "Spe", "Sly", "Nto", "Nsy", "Spi")

cafe.result.df <- parseCafeSignExpContrFamFile(input.args[[1]], spec.names)

#' Save results:
save(cafe.result.df, file = file.path(input.args[[2]], "cafe_result.RData")) 
