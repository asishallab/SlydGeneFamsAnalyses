require(SlydGeneFamsAnalyses)

option_list = list(
  make_option(
    c("-w", "--workDir"),
    type = "character",
    default = NULL,
    help = "The directory in which to find each family's work-dir and the respective MEME results in it.",
    metavar = "character"
  ),
  make_option(
    c("-o", "--outDir"),
    type = "character",
    default = detectCores(),
    help = "The path to the directory into which to write the binary results.",
    metavar = "character"
  ),
  make_option(
    c("-c", "--cores"),
    type = "integer",
    default = detectCores(),
    help = "The number of cores to use in parallel. Default is the result of 'detectCores()', i.e. ALL.",
    metavar = "integer"
  )
)

opt_parser = OptionParser(option_list = option_list)
script.args = parse_args(opt_parser)

# Validate input:
if (is.null(script.args$workDir)) {
    stop("Please specify the required argument '--workDir'.")
}

# Prepare multi core analysis:
options(mc.cores = script.args$cores)
message("Set mc.cores to ", script.args$cores)

# For all alignments compute the valdar score:
msa.files <- system(paste("find", script.args$workDir, "-type f", 
    "-name '*_AA_MSA_orig_gene_ids.fa'"), intern = TRUE)
families.msa.scores.df <- do.call(rbind, mclapply(msa.files, function(msa.file) {
    fam.name <- sub("^.*/", "", sub("_AA_MSA_orig_gene_ids.fa$", "", msa.file))
    data.frame(Family = fam.name, Valdar.Score = valdarMultipleAlignmentScore(as.matrix(read.alignment(msa.file, 
        format = "fasta"))), stringsAsFactors = FALSE)
}))

# Save results:
save(families.msa.scores.df, file = file.path(script.args$outDir, 
    "valdarScores.RData"))

message("DONE")
