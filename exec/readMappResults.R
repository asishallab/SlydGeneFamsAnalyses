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

# Use absolute path:
script.args$workDir <- normalizePath(script.args$workDir)

# Read in MAPP results for all families:
mapp.result.files <- system(paste("find", script.args$workDir, "-type f", 
    "! -size 0", "-name 'OG*_MAPP_out.txt'"), intern = TRUE)
names(mapp.result.files) <- sub("^.*/", "", sub("_MAPP_out.txt$", 
    "", mapp.result.files))
fams.mapp.df <- do.call(rbind, mclapply(names(mapp.result.files), 
    function(fam.name) {
        readMappResult(mapp.result.files[[fam.name]], fam.name)
    }))

# Adjust P values for multiple hypothesis testing:
p.value.cols <- c("Column.p.value", "A.1", "C.1", "D.1", "E.1", "F.1", 
    "G.1", "H.1", "I.1", "K.1", "L.1", "M.1", "N.1", "P.1", "Q.1", 
    "R.1", "S.1", "T.1", "V.1", "W.1", "Y.1")
fams.mapp.df <- cbind(fams.mapp.df, matrix(p.adjust(unlist(fams.mapp.df[, 
    p.value.cols]), method = "BH"), ncol = length(p.value.cols), dimnames = list(c(), 
    paste0(p.value.cols, ".adj"))))


# Find Slyd genes of families with good MSAs and with residues at
# sites with significant divergency:
fams.good.msa.slyd.genes <- intersect(names(slyd.cds), unlist(gene.families[families.msa.scores.df[which(families.msa.scores.df$Valdar.Score >= 
    0.6), "Family"]]))
fams.mapp.slyd.df <- do.call(rbind, mclapply(unique(fams.mapp.df$Family), 
    function(fam.name) {
        fam.aa.msa.mtrx <- readMultipleSequenceAlignmentAsMatrix(file.path(script.args$workDir, 
            fam.name, paste0(fam.name, "_AA_MSA_orig_gene_ids.fa")))
        fam.mapp.df <- fams.mapp.df[which(fams.mapp.df$Family == fam.name), 
            ]
        findGenesWithPhysicoChemicalDivergentAA(fam.mapp.df, fam.aa.msa.mtrx, 
            fam.name, genes.of.interest = fams.good.msa.slyd.genes, p.adjusted.cutoff = 0.01)
    }))

# Find Pfam domains overlapping with significantly divergent
# residues in genes of families with good alignments:
fams.mapp.slyd.pfam.doms.df <- do.call(rbind, mclapply(1:nrow(fams.mapp.slyd.df), 
    function(row.i) {
        fmms.row <- fams.mapp.slyd.df[row.i, ]
        fams.hmmer3.pfam.df <- parseHmmer3DomTableOut(file.path(script.args$workDir, 
            fmms.row$Family, paste0(fmms.row$Family, "_HMMER3_PfamA_domtblout.txt")))
        d.f.p <- domainsForPos(fmms.row$Protein, fmms.row$Site, fams.hmmer3.pfam.df, 
            gene.col = "query.name", start.col = "ali.coord.from", 
            end.col = "ali.coord.to", ipr.col = "target.accession")
        if (length(d.f.p) > 0) {
            pos.sel.site <- fmms.row$Site %in% fams.meme.slyd.genes.df[which(fams.meme.slyd.genes.df$Protein == 
                fmms.row$Protein), "aligned.pos.sel.codon"]
            prot.mm4 <- all.prots.mm4[which(all.prots.mm4$IDENTIFIER == 
                tolower(fmms.row$Protein)), ]
            prot.mm4.BINCODES <- paste(prot.mm4$BINCODE, collapse = ",")
            prot.mm4.NAMES <- paste(prot.mm4$NAME, collapse = ",")
            data.frame(Protein = fmms.row$Protein, Family=fmms.row$Family, aligned.divergent.site = fmms.row$Site, 
                Divergent.AA = fmms.row$Divergent.AA, AA.p.value.adj = fmms.row$AA.p.value.adj, 
                is.pos.sel.site = pos.sel.site, Pfam.accession = d.f.p, 
                Pfam.description = unlist(lapply(d.f.p, function(x) {
                  fams.hmmer3.pfam.df[which(fams.hmmer3.pfam.df$target.accession == 
                    x), "description.of.target"][[1]]
                })), BINCODE = prot.mm4.BINCODES, NAME = prot.mm4.NAMES, 
                stringsAsFactors = FALSE)
        } else NULL
    }))

# Save results:
save(fams.mapp.df, fams.good.msa.slyd.genes, fams.mapp.slyd.df, fams.mapp.slyd.pfam.doms.df, 
    file = file.path(script.args$outDir, "mappResults.RData"))

message("DONE")
