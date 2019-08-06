require(SlydGeneFamsAnalyses)

option_list = list(make_option(c("-w", "--workDir"), type = "character", 
    default = NULL, help = "The directory in which to find each family's work-dir and the respective MEME results in it.", 
    metavar = "character"), make_option(c("-o", "--outDir"), 
    type = "character", default = NULL, help = "The path to the directory into which to write the binary results.", 
    metavar = "character"), make_option(c("-p", "--allInterProAnnosTbl"), 
    type = "character", default = NULL, help = "The path to the tabular InterProScan results for all proteomes analyzed in this package.", 
    metavar = "character"), make_option(c("-c", "--cores"), type = "integer", 
    default = detectCores(), help = "The number of cores to use in parallel. Default is the result of 'detectCores()', i.e. ALL.", 
    metavar = "integer"))

opt_parser = OptionParser(option_list = option_list)
script.args = parse_args(opt_parser)

# Validate input:
if (is.null(script.args$workDir) || is.null(script.args$allInterProAnnosTbl) || 
    is.null(script.args$outDir)) {
    stop("Please specify the required arguments '--workDira, '--allInterProAnnosTbl', and '--outDir'.")
}

# Prepare multi core analysis:
options(mc.cores = script.args$cores)
message("Set mc.cores to ", script.args$cores)

# Use absolute path:
script.args$workDir <- normalizePath(script.args$workDir)
script.args$outDir <- normalizePath(script.args$outDir)
script.args$allInterProAnnosTbl <- normalizePath(script.args$allInterProAnnosTbl)

# Read InterProScan results table for all analyzed Proteomes:
all.ipr <- readInterProScanResultTable(script.args$allInterProAnnosTbl)

# Read in MAPP results for all families:
mapp.result.files <- system(paste("find", script.args$workDir, 
    "-type f", "! -size 0", "-name 'OG*_mapp.out'"), intern = TRUE)
names(mapp.result.files) <- sub("^.*/", "", sub("_mapp.out$", 
    "", mapp.result.files))
fams.mapp.df <- do.call(rbind, mclapply(names(mapp.result.files), 
    function(fam.name) {
        readMappResult(mapp.result.files[[fam.name]], fam.name)
    }))

# Adjust P values for multiple hypothesis testing:
p.value.cols <- c("Column.p.value", "A.1", "C.1", "D.1", "E.1", 
    "F.1", "G.1", "H.1", "I.1", "K.1", "L.1", "M.1", "N.1", "P.1", 
    "Q.1", "R.1", "S.1", "T.1", "V.1", "W.1", "Y.1")
fams.mapp.df <- cbind(fams.mapp.df, matrix(p.adjust(unlist(fams.mapp.df[, 
    p.value.cols]), method = "BH"), ncol = length(p.value.cols), 
    dimnames = list(c(), paste0(p.value.cols, ".adj"))))


# Find Datura genes of families with good MSAs and with
# residues at sites with significant divergency:
fams.good.msa.slyd.genes <- intersect(c(names(date.cds), names(dati.cds)), 
    unlist(gene.families[families.msa.scores.df[which(families.msa.scores.df$Valdar.Score >= 
        0.6), "Family"]]))
fams.mapp.slyd.df <- do.call(rbind, mclapply(unique(fams.mapp.df$Family), 
    function(fam.name) {
        fam.aa.msa.mtrx <- readMultipleSequenceAlignmentAsMatrix(file.path(script.args$workDir, 
            fam.name, paste0(fam.name, "_AA_MSA_orig_gene_ids.fa")))
        fam.mapp.df <- fams.mapp.df[which(fams.mapp.df$Family == 
            fam.name), ]
        findGenesWithPhysicoChemicalDivergentAA(fam.mapp.df, 
            fam.aa.msa.mtrx, fam.name, genes.of.interest = fams.good.msa.slyd.genes, 
            p.adjusted.cutoff = 0.01)
    }))

# Find ipr domains overlapping with significantly divergent
# residues in genes of families with good alignments:
fams.mapp.slyd.ipr.doms.df <- do.call(rbind, mclapply(1:nrow(fams.mapp.slyd.df), 
    function(row.i) {
        fmms.row <- fams.mapp.slyd.df[row.i, ]
        d.f.p <- domainsForPos(fmms.row$Protein, fmms.row$Site, 
            all.ipr)
        if (length(d.f.p) > 0) {
            x <- data.frame(Protein = fmms.row$Protein, Family = fmms.row$Family, 
                aligned.divergent.site = fmms.row$Site, Divergent.AA = fmms.row$Divergent.AA, 
                AA.p.value.adj = fmms.row$AA.p.value.adj, ipr.accession = d.f.p, 
                ipr.name = unlist(lapply(d.f.p, function(x) {
                  if (x %in% names(ipr.db)) {
                    ipr.entry <- ipr.db[[x]]
                    ipr.entry$NAME
                  } else NA
                })), stringsAsFactors = FALSE)
        } else NULL
    }))


# Enrichment test for Datura proteins with InterPro Domains
# overlapping MAPP significant amino acids:
date.dati.mapp.ipr.prots <- unique(fams.mapp.slyd.ipr.doms.df$Protein)
all.ipr.4.enrich <- unique(all.ipr[which(grepl("^(date)|(dati)", 
    all.ipr$V1, perl = TRUE) & !is.na(all.ipr$V12) & !is.na(all.ipr$V1)), 
    c("V1", "V12")])
colnames(all.ipr.4.enrich) <- paste0("V", 1:2)
fams.mapp.slyd.ipr.doms.enriched <- enrichedAnnotations(date.dati.mapp.ipr.prots, 
    all.ipr.4.enrich, annos.2.test = unique(fams.mapp.slyd.ipr.doms.df$ipr.accession))
fams.mapp.slyd.ipr.doms.enriched.sign <- selectSignificantlyEnrichedAnnotations(fams.mapp.slyd.ipr.doms.enriched)

# Save results:
save(fams.mapp.df, fams.good.msa.slyd.genes, fams.mapp.slyd.df, 
    fams.mapp.slyd.ipr.doms.df, fams.mapp.slyd.ipr.doms.enriched.sign, 
    file = file.path(script.args$outDir, "mappResults.RData"))

message("DONE")
