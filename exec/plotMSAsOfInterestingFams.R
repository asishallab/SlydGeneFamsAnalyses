require(SlydGeneFamsAnalyses)

option_list <- list(
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
    help = "The path to the directory into which to write the plots.",
    metavar = "character"
  )
)

opt_parser = OptionParser(option_list = option_list)
script.args = parse_args(opt_parser)

# Validate input:
if (is.null(script.args$outDir) || is.null(script.args$workDir)) {
    stop("Please specify the required arguments '--outDir' and '--workDir'.")
}

# Plot the MSAs of the families with intersting results for Slyd:
pos.sel.fams <- unique(fams.meme.slyd.genes.sel.pfam.doms.df$Family)
for (fam.name in pos.sel.fams) {
    tryCatch({
        printAaMsaWithSelection(file.path(script.args$workDir, fam.name), 
            paste0(fam.name, "_AA_MSA_orig_gene_ids.fa"), fams.meme[[fam.name]]$MEME.branches, 
            NULL)
    }, error = function(e) {
        warning("Could not plot the MSA for gene family '", fam.name, 
            "'.\n", e)
    })
}

fams.with.slyd.divergent.aa <- unique(fams.mapp.slyd.pfam.doms.df$Family)
for (fam.name in fams.with.slyd.divergent.aa) {
    tryCatch({
        fam.mock.meme.branches <- fams.mapp.slyd.pfam.doms.df[which(fams.mapp.slyd.pfam.doms.df$Family == 
            fam.name), ]
        fam.mock.meme.branches$aligned.pos.sel.codon <- fam.mock.meme.branches$aligned.divergent.site
        fam.mock.meme.branches$unaligned.pos.sel.codon <- unlist(lapply(1:nrow(fam.mock.meme.branches), 
            function(i) {
                unalignedAAforAlignedAAPos(fam.mock.meme.branches[i,]$Protein, fam.mock.meme.branches[i,]$aligned.divergent.site, 
                  read.fasta(file.path(script.args$workDir, fam.name, 
                    paste0(fam.name, "_AA_MSA_orig_gene_ids.fa")), 
                    seqtype = "AA", strip.desc = TRUE, as.string = TRUE))
            }))
        printAaMsaWithSelection(file.path(script.args$workDir, fam.name), 
            paste0(fam.name, "_AA_MSA_orig_gene_ids.fa"), fam.mock.meme.branches, 
            NULL)
    }, error = function(e) {
        warning("Could not plot the MSA for gene family '", fam.name, 
            "'.\n", e)
    })
}


# Use absolute path:
script.args$outDir <- normalizePath(script.args$outDir)

message("DONE")
