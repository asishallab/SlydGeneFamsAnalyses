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
    c("-m", "--mercatorMapManAnnotations"),
    type = "character",
    default = detectCores(),
    help = "The tabular output of running Mercator (MapMan 4) on the proteins of the gene families.",
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

# Read Mercator (MapMan 4) gene function annotations:
mm4.tbl <- if (!is.null(script.args$mercatorMapManAnnotations)) {
    read.table(script.args$mercatorMapManAnnotations, sep = "\t", 
        header = TRUE, quote = "'", stringsAsFactors = FALSE)
} else NULL

# Annotate those Families for which we found evidence of positive
# selection:
fams.meme.anno.df <- Reduce(rbind, mclapply(names(fams.meme[fams.meme.i]), 
    function(fam.name) {
        fam.dir <- normalizePath(file.path(script.args$workDir, fam.name))
        fam.gene.ids <- gene.families[[fam.name]]
        fam.pfam.hmmer3.domtbl <- parseHmmer3DomTableOut(file.path(fam.dir, 
            paste0(fam.name, "_HMMER3_PfamA_domtblout.txt")))
        annotateGeneFamily(fam.gene.ids, fam.name, fam.pfam.hmmer3.domtbl, 
            mm4.tbl)
    }))

# Save results
write.table(fams.meme.anno.df, file.path(script.args$outDir, "geneFamilyAnnotations.txt"), 
    sep = "\t", row.names = FALSE, quote = FALSE)

message("DONE")
