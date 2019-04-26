require(SlydGeneFamsAnalyses)

option_list = list(
  make_option(
    c("-o", "--orthogroupsTxt"),
    type = "character",
    default = NULL,
    help = "The valid file path to Orthofinder's output 'Orthogroups.txt'. Use the txt and not the tsv output.",
    metavar = "character"
  ),
  make_option(
    c("-m", "--orthogroupsMSADir"),
    type = "character",
    default = NULL,
    help = "The valid file path to Orthofinder's output in which the respective multiple amino acid sequence alignments are stored.",
    metavar = "character"
  ),
  make_option(
    c("-p", "--orthogroupsTreesDir"),
    type = "character",
    default = NULL,
    help = "The valid file path to Orthofinder's output directory in which the respective phylogenetic trees are stored.",
    metavar = "character"
  ),
  make_option(
    c("-d", "--dataDir"),
    type = "character",
    default = './data',
    help = "The directory into which to write the binary RData holding the parsed gene families objects. Default is './data'",
    metavar = "character"
  ),
  make_option(
    c("-w", "--workDir"),
    type = "character",
    default = NULL,
    help = "The directory into which to write each family's parsed phylogenetic tree and multiple amino acid sequence alignment.",
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

if (is.null(args$orthogroupsTxt) || is.null(args$orthogroupsMSADir) || 
    is.null(args$orthogroupsTreesDir) || is.null(args$workDir)) {
    stop("Please specify the required arguments '--orthogroupsTxt', '--orthogroupsMSADir', '--orthogroupsTreesDir', and '--workDir'.")
}

options(mc.cores = args$cores)
message("Set mc.cores to ", args$cores)

# Parse and load Gene Families:
gene.families <- parseOrthogroupsTxt(args$orthogroupsTxt)
gene.families.sizes <- lapply(gene.families, length)
gene.families.non.singletons.names <- names(gene.families.sizes[which(unlist(gene.families.sizes) > 
    1)])

# Parse and load the Gene Families' multiple amino acid
# alignments, write them into the families' respective work dir:
loadAndSanitizeAAMsas(args$orthogroupsMSADir, args$workDir)

# Parse and load the Gene Families' phylogenetic trees, write them
# into the families' respective work dir:
loadAndSanitizeTrees(args$orthogroupsTreesDir, args$workDir)

# Save results
save(gene.families, gene.families.sizes, gene.families.non.singletons.names, 
    file = file.path(args$dataDir, "geneFamilies.RData"))
