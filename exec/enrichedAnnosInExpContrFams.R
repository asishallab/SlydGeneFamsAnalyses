require(SlydGeneFamsAnalyses)

option_list = list(
  make_option(
    c("-w", "--mappResults"),
    type = "character",
    default = NULL,
    help = "The binary (RData) file in which the parsed and loaded MAPP results have been stored. See script ./exec/readMappResults.R for details.",
    metavar = "character"
  ),
  make_option(
    c("-o", "--outDir"),
    type = "character",
    default = NULL,
    help = "The path to the directory into which to write the table of enriched MapMan4 Annos found in Datura genes with significant divergent sites in conserved protein domains.",
    metavar = "character"
  ),
  make_option(
    c("-p", "--allMapMan4Annos"),
    type = "character",
    default = NULL,
    help = "The path to the tabular InterProScan results for all proteomes analyzed in this package.",
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
if (is.null(script.args$mappResults) ||
    is.null(script.args$allMapMan4Annos) ||
    is.null(script.args$outDir)) {
  stop(
    "Please specify the required arguments '--mappResultsa, '--allMapMan4Annos', and '--outDir'."
  )
}

# Prepare multi core analysis:
options(mc.cores = script.args$cores)
message("Set mc.cores to ", script.args$cores)

# Use absolute path:
script.args$mappResults <- normalizePath(script.args$mappResults)
script.args$outDir <- normalizePath(script.args$outDir)
script.args$allMapMan4Annos <-
  normalizePath(script.args$allMapMan4Annos)

# Read Mapp results
load(script.args$mappResults)

# Read all MapMan4 Annotations
all.prots.mm4 <-
  readMercatorResultTable(script.args$allMapMan4Annos,
                          FALSE)

# Enrichment test for Datura proteins with MapMan4
# Annotations overlapping MAPP significant amino acids:
date.dati.mapp.ipr.prots <-
  unique(fams.mapp.slyd.ipr.doms.df$Protein)
all.prots.mm4.4.enrich <-
  as.data.frame(unique(all.prots.mm4[which(grepl("^(date)|(dati)",
                                                 all.prots.mm4$IDENTIFIER, perl = TRUE)),]))
mm4.to.test <-
  sort(unique(all.prots.mm4.4.enrich[which(all.prots.mm4.4.enrich$IDENTIFIER %in%
                                             date.dati.mapp.ipr.prots),]$BINCODE))
fams.mapp.slyd.mm4.doms.enriched <-
  enrichedAnnotations(
    date.dati.mapp.ipr.prots,
    all.prots.mm4.4.enrich,
    univ.gene.col = "IDENTIFIER",
    univ.anno.col = "BINCODE",
    annos.2.test = mm4.to.test
  )
fams.mapp.slyd.mm4.doms.enriched.sign.df <-
  unique(all.prots.mm4[which(all.prots.mm4$BINCODE %in%
                               names(fams.mapp.slyd.mm4.doms.enriched[which(fams.mapp.slyd.mm4.doms.enriched <
                                                                              0.05)])), c("BINCODE", "NAME")])

# Save results:
write.table(
  fams.mapp.slyd.mm4.doms.enriched.sign.df,
  file.path(
    script.args$outDir,
    "enrichedMapMan4BinsInDaturaGenesWithSignDivergentSitesInConservedDomains.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("DONE")
