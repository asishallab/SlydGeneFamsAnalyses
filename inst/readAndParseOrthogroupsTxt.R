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
    c("-d", "--dataDir"),
    type = "character",
    default = './data',
    help = "The directory into which to write the binary RData holding the parsed gene families objects. Default is './data'",
    metavar = "character"
  )
)

opt_parser = OptionParser(option_list = option_list)
args = parse_args(opt_parser)

if (is.null(args$orthogroupsTxt)) {
    stop("Please specify the required argument '--orthogroupsTxt'")
}

gene.families <- parseOrthogroupsTxt(args$orthogroupsTxt)
gene.families.sizes <- lapply(gene.families, 
    length)
gene.families.non.singletons.names <- names(gene.families.sizes[which(unlist(gene.families.sizes) > 
    1)])

save(gene.families, gene.families.sizes, 
    gene.families.non.singletons.names, 
    file = file.path(args$dataDir, 
        "geneFamilies.RData"))
