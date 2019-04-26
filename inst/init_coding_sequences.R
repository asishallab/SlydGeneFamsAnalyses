require(SlydGeneFamsAnalyses)

option_list = list(
  make_option(
    c("-o", "--outdir"),
    type = "character",
    default = NULL,
    help = "The working directory in which the resulting binary database of coding sequences is to be stored.",
    metavar = "character"
  ),
  make_option(
    c("-i", "--inputdir"),
    type = "character",
    default = NULL,
    help = "The 'material' directory in which to lookup the respective coding sequence FASTA files. See script './inst/init_coding_sequences.R' for details on which files are expected to be present there.",
    metavar = "character"
  )
)

opt_parser = OptionParser(option_list = option_list)
args = parse_args(opt_parser)

if (is.null(args$inputdir) || is.null(args$outdir)) stop("Please specify the required arguments")

ath.cds <- read.fasta(file.path(args$inputdir, "TAIR10_cds_20110103_representative_gene_model_updated"), 
    as.string = TRUE, strip.desc = TRUE, seqtype = "DNA")
vvi.cds <- read.fasta(file.path(args$inputdir, "Vitis_vinifera_mRNA.fa"), 
    strip.desc = TRUE, as.string = TRUE, seqtype = "DNA")
slyd.cds <- read.fasta(file.path(args$inputdir, "SlydLA2951_v0.6_cds_all.fasta"), 
    strip.desc = TRUE, as.string = TRUE, seqtype = "DNA")
stu.cds <- read.fasta(file.path(args$inputdir, "Solanum_tuberosum_v3.4_cds_matching_peptides.fa"), 
    strip.desc = TRUE, as.string = TRUE, seqtype = "DNA")
slyc.cds <- read.fasta(file.path(args$inputdir, "ITAG3.2_CDS.fasta"), 
    strip.desc = TRUE, as.string = TRUE, seqtype = "DNA")
spe.cds <- read.fasta(file.path(args$inputdir, "Spenn-v2-cds-annot.fa"), 
    strip.desc = TRUE, as.string = TRUE, seqtype = "DNA")

save(ath.cds, slyc.cds, slyd.cds, spe.cds, stu.cds, 
    vvi.cds, file = file.path(args$outdir, "all_CDS.RData"))
