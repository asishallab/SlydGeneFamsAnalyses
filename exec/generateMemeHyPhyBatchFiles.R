require(SlydGeneFamsAnalyses)

option_list = list(
  make_option(
    c("-w", "--workDir"),
    type = "character",
    default = NULL,
    help = "The directory into which to write each family's parsed phylogenetic tree and multiple amino acid sequence alignment.",
    metavar = "character"
  ),
  make_option(
    c("-b", "--HyPhyBatchFilesDir"),
    type = "character",
    default = NULL,
    help = "The sub-directory of the HyPhy installation in which to find HyPhy's Batch-Files.",
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
args = parse_args(opt_parser)

if (is.null(args$workDir) || is.null(args$HyPhyBatchFilesDir)) {
    stop("Please specify the required argument '--workDir' and '--HyPhyBatchFilesDir'.")
}

options(mc.cores = args$cores)
message("Set mc.cores to ", args$cores)

# Generate each family's multiple coding sequence alignment:
generateMemeBatchFiles(args$workDir, args$HyPhyBatchFilesDir)

# DONE
message("DONE")
