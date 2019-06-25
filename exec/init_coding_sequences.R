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

cag.cds <- read.fasta(file.path(args$inputdir, "Capsicumannuumglabriusculum.fasta"),
    as.string = TRUE, strip.desc = TRUE, seqtype = "DNA")
cam.cds <- read.fasta(file.path(args$inputdir, "Capsicumannuummorelia.fasta"),
    strip.desc = TRUE, as.string = TRUE, seqtype = "DNA")
date.cds <- read.fasta(file.path(args$inputdir, "DaturastramoniumTeo1.fasta"),
    strip.desc = TRUE, as.string = TRUE, seqtype = "DNA")
dati.cds <- read.fasta(file.path(args$inputdir, "DaturastramoniumTic23.fasta"),
    strip.desc = TRUE, as.string = TRUE, seqtype = "DNA")
nat.cds <- read.fasta(file.path(args$inputdir, "Nicotianaattenuata.fasta"),
    strip.desc = TRUE, as.string = TRUE, seqtype = "DNA")
nsy.cds <- read.fasta(file.path(args$inputdir, "Nicotianasylvestris.fasta"),
    strip.desc = TRUE, as.string = TRUE, seqtype = "DNA")
nta.cds <- read.fasta(file.path(args$inputdir, "Nicotianatabacum.fasta"),
    strip.desc = TRUE, as.string = TRUE, seqtype = "DNA")
nto.cds <- read.fasta(file.path(args$inputdir, "Nicotianatomentosiformis.fasta"),
    strip.desc = TRUE, as.string = TRUE, seqtype = "DNA")
pi.cds <- read.fasta(file.path(args$inputdir, "Petuniainflata.fasta"),
    strip.desc = TRUE, as.string = TRUE, seqtype = "DNA")
sly.cds <- read.fasta(file.path(args$inputdir, "Solanumlycopersicum.fasta"),
    strip.desc = TRUE, as.string = TRUE, seqtype = "DNA")
spe.cds <- read.fasta(file.path(args$inputdir, "Solanumpennellii.fasta"),
    strip.desc = TRUE, as.string = TRUE, seqtype = "DNA")
spi.cds <- read.fasta(file.path(args$inputdir, "Solanumpimpinellifolium.fasta"),
    strip.desc = TRUE, as.string = TRUE, seqtype = "DNA")
stu.cds <- read.fasta(file.path(args$inputdir, "Solanumtuberosum.fasta"),
    strip.desc = TRUE, as.string = TRUE, seqtype = "DNA")

save(cag.cds, cam.cds, date.cds, dati.cds, nat.cds, nsy.cds, nta.cds, nto.cds, pi.cds, sly.cds, spe.cds, spi.cds, stu.cds,
         file = file.path(args$outdir, "all_CDS.RData"))
