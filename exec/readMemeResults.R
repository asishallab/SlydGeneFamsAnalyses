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

# MEME results as a list of tables. Two per family:
fams.meme <- readMemeResults(script.args$workDir)
fams.meme.i <-
  which(unlist(lapply(fams.meme, function(x)
    ! is.null(x))))

# Table of Slyd genes that are found to be subject of positive
# selection:
fams.meme.slyd.genes.df <-
  do.call(rbind, lapply(names(fams.meme[fams.meme.i]),
                       function(g.f) {
                         x <- fams.meme[[g.f]]$MEME.branches
                         d.f <- x[which(x$Protein %in% names(slyd.cds)),]
                         if (!is.null(d.f) && nrow(d.f) > 0) {
                           d.f
                         } else {
                           NULL
                         }
                       }))

# Assign the Mercator MapMan4 gene function annotations to the
# MEME results:
if (!is.null(all.prots.mm4)) {
  fams.meme.slyd.genes.df$BINCODE <- c()
  fams.meme.slyd.genes.df$NAME <- c()
  for (i in 1:nrow(fams.meme.slyd.genes.df)) {
    prot.id <- tolower(fams.meme.slyd.genes.df[[i, "Protein"]])
    if (prot.id %in% all.prots.mm4$IDENTIFIER) {
      all.prots.mm4.prot <- all.prots.mm4[which(all.prots.mm4$IDENTIFIER ==
                                      prot.id),]
      fams.meme.slyd.genes.df$BINCODE[[i]] <-
        paste(all.prots.mm4.prot$BINCODE,
              collapse = ",")
      fams.meme.slyd.genes.df$NAME[[i]] <-
        paste(all.prots.mm4.prot$NAME,
              collapse = ",")
    }
  }
}

# Create the plots with the MEME results:
for (fam.nm in names(fams.meme[fams.meme.i])) {
  fam.meme.branches.df <- fams.meme[[fam.nm]]$MEME.branches
  if (!is.null(fam.meme.branches.df)) {
    tryCatch({
      fam.dir <- normalizePath(file.path(script.args$workDir,
                                         fam.nm))
      fam.aa.msa.fasta <- file.path(fam.dir, paste0(fam.nm,
                                                    "_AA_MSA_orig_gene_ids.fa"))
      fam.pfam.results <-
        parseHmmer3DomTableOut(file.path(fam.dir,
                                         paste0(
                                           fam.nm, "_HMMER3_PfamA_domtblout.txt"
                                         )))
      generateInteractiveMsaPlot(fam.aa.msa.fasta,
                                 fam.meme.branches.df,
                                 fam.pfam.results,
                                 fam.dir)
    }, error = function(e) {
      message(
        "An error occurred when plotting MEME results for family '",
        fam.nm,
        "'. Will continue with next family.\n",
        e
      )
    })
  }
}

# Save results
save(
  fams.meme,
  fams.meme.i,
  fams.meme.slyd.genes.df,
  file = file.path(script.args$outDir,
                   "memeResults.RData")
)

message("DONE")
