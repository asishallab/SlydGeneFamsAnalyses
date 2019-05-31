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

# Find the families for the pos sel Slyd genes:
fams.meme.hmmer3.pfam.df <- do.call(rbind, mclapply(names(fams.meme.i), 
    function(fam.name) {
        parseHmmer3DomTableOut(file.path(script.args$workDir, fam.name, 
            paste0(fam.name, "_HMMER3_PfamA_domtblout.txt")))
    }))
fams.meme.slyd.genes.sel.pfam.doms.df <- do.call(rbind, mclapply(1:nrow(fams.meme.slyd.genes.df), 
    function(row.i) {
        fmsg.row <- fams.meme.slyd.genes.df[row.i, ]
        d.f.p <- domainsForPos(fmsg.row$Protein, fmsg.row$aligned.pos.sel.codon, 
            fams.meme.hmmer3.pfam.df, gene.col = "query.name", start.col = "ali.coord.from", 
            end.col = "ali.coord.to", ipr.col = "target.accession")
        if (length(d.f.p) > 0) {
            prot.mm4 <- all.prots.mm4[which(all.prots.mm4$IDENTIFIER == 
                tolower(fmsg.row$Protein)), ]
            prot.mm4.BINCODES <- paste(prot.mm4$BINCODE, collapse = ",")
            prot.mm4.NAMES <- paste(prot.mm4$NAME, collapse = ",")
            data.frame(Protein = fmsg.row$Protein, aligned.pos.sel.codon = fmsg.row$aligned.pos.sel.codon, 
                unaligned.pos.sel.codon = fmsg.row$unaligned.pos.sel.codon, 
                MEME.site.p.value = fmsg.row$MEME.site.p.value, MEME.site.p.adj = fmsg.row$MEME.site.p.adj, 
                Pfam.accession = d.f.p, Pfam.description = unlist(lapply(d.f.p, 
                  function(x) {
                    fams.meme.hmmer3.pfam.df[which(fams.meme.hmmer3.pfam.df$target.accession == 
                      x), "description.of.target"][[1]]
                  })), BINCODE = prot.mm4.BINCODES, NAME = prot.mm4.NAMES, 
                stringsAsFactors = FALSE)
        } else NULL
    }))

# Find Slyd genes in pos sel families with good AA MSAs:
fams.meme.good.msas <- families.msa.scores.df[which(families.msa.scores.df$Family %in% 
    names(fams.meme.i) & families.msa.scores.df$Valdar.Score >= 0.6), 
    "Family"]
fams.meme.good.msas.slyd.genes <- intersect(fams.meme.slyd.genes.df$Protein, 
    intersect(names(slyd.cds), unlist(gene.families[fams.meme.good.msas])))


# Save results:
save(fams.meme.hmmer3.pfam.df, fams.meme.slyd.genes.sel.pfam.doms.df, 
    fams.meme.good.msas, fams.meme.good.msas.slyd.genes, file = file.path(script.args$outDir, 
        "domainsAtSelectedSites.RData"))

message("DONE")
