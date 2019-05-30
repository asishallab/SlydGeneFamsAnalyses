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

# Read in MAPP results for all families with positive selection:
fams.meme.mapp.df <- Reduce(rbind, mclapply(names(fams.meme.i), function(fam.name) {
    fam.mapp.tbl <- file.path(script.args$workDir, fam.name, paste0(fam.name, 
        "_MAPP_out.txt"))
    readMappResult(fam.mapp.tbl, fam.name)
}))
fams.meme.mapp.df$Column.p.adjusted <- p.adjust(fams.meme.mapp.df$Column.p.value, 
    method = "fdr")

# Find Slyd genes of families with good MSAs and with residues at
# sites with significant divergency:
fams.good.msa.slyd.genes <- intersect(names(slyd.cds), unlist(gene.families[families.msa.scores.df[which(families.msa.scores.df$Valdar.Score >= 
    0.6), "Family"]]))
fams.meme.mapp.slyd.df <- Reduce(rbind, mclapply(unique(fams.meme.mapp.df$Family), 
    function(fam.name) {
        fam.aa.msa.mtrx <- readMultipleSequenceAlignmentAsMatrix(file.path(script.args$workDir, 
            fam.name, paste0(fam.name, "_AA_MSA_orig_gene_ids.fa")))
        fam.mapp.df <- fams.meme.mapp.df[which(fams.meme.mapp.df$Family == 
            fam.name), ]
        findGenesWithPhysicoChemicalDivergentAA(fam.mapp.df, fam.aa.msa.mtrx, 
            genes.of.interest = fams.good.msa.slyd.genes)
    }))
fams.meme.mapp.slyd.df$AA.p.adjusted <- p.adjust(fams.meme.mapp.slyd.df$AA.p.value, 
    method = "BH")

# Find Pfam domains overlapping with significantly divergent
# residues in genes of families with good alignments:
fams.meme.mapp.slyd.sel.pfam.doms.df <- Reduce(rbind, mclapply(1:nrow(fams.meme.mapp.slyd.df), 
    function(row.i) {
        fmms.row <- fams.meme.mapp.slyd.df[row.i, ]
        d.f.p <- domainsForPos(fmms.row$Protein, fmms.row$Site, fams.meme.hmmer3.pfam.df, 
            gene.col = "query.name", start.col = "ali.coord.from", 
            end.col = "ali.coord.to", ipr.col = "target.accession")
        if (length(d.f.p) > 0) {
            data.frame(Protein = fmms.row$Protein, aligned.divergent.site = fmms.row$Site, 
                Divergent.AA = fmms.row$Divergent.AA, AA.p.value = fmms.row$AA.p.value, 
                AA.p.adjusted = fmms.row$AA.p.adjusted, Pfam.accession = d.f.p, 
                Pfam.description = unlist(lapply(d.f.p, function(x) {
                  fams.meme.hmmer3.pfam.df[which(fams.meme.hmmer3.pfam.df$target.accession == 
                    x), "description.of.target"][[1]]
                })), stringsAsFactors = FALSE)
        } else NULL
    }))

# Save results:
save(fams.meme.mapp.df, fams.good.msa.slyd.genes, fams.meme.mapp.slyd.df, 
    fams.meme.mapp.slyd.sel.pfam.doms.df, file = file.path(script.args$outDir, 
        "mappResults.RData"))

message("DONE")
