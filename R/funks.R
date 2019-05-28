#' Reads and parses the Orthofinder gene families output file.
#'
#' @param path.2.orthogroups.txt - The valid file path to the Orthofinder
#' output 'orthogroups.txt'. Use the txt version, not the tsv, please.
#'
#' @return A list of gene families, names are the family names, and values are
#' the character vectors of each respective family's gene members.
#' @export
parseOrthogroupsTxt <- function(path.2.orthogroups.txt) {
    fams.raw <- readLines(path.2.orthogroups.txt)
    fams.splt <- strsplit(fams.raw, ":")
    fam.names <- unlist(lapply(fams.splt, function(x) x[[1]]))
    setNames(lapply(fams.splt, function(fam.str) {
        setdiff(strsplit(fam.str[[2]], " ")[[1]], "")
    }), fam.names)
}

#' Orthofinder puts the species name as a prefix before all gene identifiers.
#' This functions removes them by deleting everything and including to the last
#' '_' underscore character. Note, make sure your gene identifiers do NOT use
#' underscores.
#'
#' @param gene.ids - A character vector of gene identifiers modified by
#' Orthofinder.
#'
#' @return A modified version of the argument 'gene.ids'.
#' @export
sanitizeOrthofinderGeneIds <- function(gene.ids) {
    sub("^.*_", "", gene.ids)
}

#' Creates the argument gene family 'gene.family.name' working directory in the
#' argument 'base.dir' if it does not already exists.
#'
#' @param base.dir - The valid file path to the base working directory
#' @param gene.family.name - The name of the gene family which defines the name
#' of the directory to be created within the argument 'base.dir'.
#'
#' @return The full path to the gene family's working directy, created if not
#' already existed.
#' @export
createGeneFamilyWorkingDirIfNotExists <- function(base.dir, gene.family.name) {
    gene.fam.wd <- normalizePath(file.path(base.dir, gene.family.name))
    system(paste("mkdir -p", gene.fam.wd))
    gene.fam.wd
}

#' Creates the argument gene family 'gene.family.name' mapping of original gene
#' identifiers to sanitized 'PROT1' ... 'PROTN' identifiers.
#'
#' @param gene.family.dir - The valid file path to the gene family's working
#' directory.
#' @param gene.family.name - The name of the gene family which defines the name
#' of the directory to be created within the argument 'base.dir'.
#' @param gene.fams.lst - The list of gene families in which to lookup the
#' argument 'gene.family.name'. Default is SlydGeneFamsAnalyses's
#' \code{gene.families}.
#'
#' @return A data.frame holding two columns: original and sanitized. The
#' dataframe is either loaded from a file 'name_mappings.txt' or created,
#' stored in that file, and then returned.
#' @export
geneNameMappings <- function(gene.family.dir, gene.family.name, gene.fams.lst = gene.families) {
    gene.fam.name.maps.path <- file.path(gene.family.dir, paste0(gene.family.name, 
        "_name_mappings.txt"))
    if (file.exists(gene.fam.name.maps.path)) {
        read.table(gene.fam.name.maps.path, header = TRUE, sep = "\t", 
            stringsAsFactors = FALSE)
    } else {
        gf.name.maps.df <- data.frame(original = gene.fams.lst[[gene.family.name]], 
            stringsAsFactors = FALSE)
        gf.name.maps.df$sanitized <- paste0("PROT", 1:length(gene.fams.lst[[gene.family.name]]))
        write.table(gf.name.maps.df, gene.fam.name.maps.path, row.names = FALSE, 
            quote = FALSE, sep = "\t")
        gf.name.maps.df
    }
}

#' Returns the sanitized gene names mapped to the argument 'originals' as
#' assigned in argument 'mappings' table.
#'
#' @param originals - Character vector of gene identifiers.
#' @param mappings - Instance of base::data.frame with two columns: 'original'
#' holding the original gene identifiers, and 'sanitized' holding the mapped
#' simplified gene identifiers. See function \code{geneNameMappings} for more
#' details.
#'
#' @return A character vector of same length as argument 'originals' holding
#' the respective sanitized gene identifiers.
#' @export
mapOriginalToSanitizedNames <- function(originals, mappings) {
    rownames(mappings) <- mappings$original
    unlist(lapply(originals, function(x) mappings[[x, "sanitized"]]))
}

#' Returns the original gene names mapped to the argument 'sanitized' as
#' assigned in argument 'mappings' table.
#'
#' @param sanitized - Character vector of gene identifiers.
#' @param mappings - Instance of base::data.frame with two columns: 'original'
#' holding the original gene identifiers, and 'sanitized' holding the mapped
#' simplified gene identifiers. See function \code{geneNameMappings} for more
#' details.
#'
#' @return A character vector of same length as argument 'sanitized' holding
#' the respective original gene identifiers.
#' @export
mapSanitizedToOriginalNames <- function(sanitized, mappings) {
    rownames(mappings) <- mappings$sanitized
    unlist(lapply(sanitized, function(x) mappings[[x, "original"]]))
}

#' Loads all multiple amino acid sequence alignments (MSAs) for the respective
#' gene families reconstructed by Orthofinder. They are read in using
#' \code{seqinr::read.fasta} and gene identifiers are sanitized, removing the
#' species name prefixes added by Orthofinder. See function
#' \code{sanitizeOrthofinderGeneIds} for more details. Note, that this function
#' uses \code{mclapply}.
#'
#' @param path.2.msa.dir - The valid file path to the directory the MSAs are
#' stored in. Make sure their file-extensions are '.fa'.
#' @param working.dir - The valid file path to the directory in which to store
#' the sanitized and parsed MSAs. Each will be stored in its respective gene
#' family's sub-directory, identified by the family's name itself.
#' @param gene.families.names - A character vector of gene family names that
#' are to be expected to contain the file names in argument 'path.2.msa.dir'
#' excluding the file-extensions. Default is this package's
#' \code{names(gene.families)}.
#'
#' @return TRUE if and only if no error has occurred.
#' @export
loadAndSanitizeAAMsas <- function(path.2.msa.dir, working.dir, gene.families.names = names(gene.families)) {
    msa.fls <- system(paste("ls", path.2.msa.dir, "| grep -P \"^OG.*\\.fa$\""), 
        intern = TRUE)
    names(msa.fls) <- sub("\\.fa$", "", msa.fls)
    if (!all(names(msa.fls) %in% gene.families.names)) {
        warning("Some of the multiple amino acid sequence alignments' file names do not appear in the argument 'gene.families.names'.")
    }
    mclapply(names(msa.fls), function(gene.fam) {
        msa.file <- msa.fls[[gene.fam]]
        msa <- read.fasta(file.path(path.2.msa.dir, msa.file), seqtype = "AA", 
            as.string = TRUE, strip.desc = TRUE)
        names(msa) <- sanitizeOrthofinderGeneIds(names(msa))
        if (!all(names(msa) %in% names(all.cds))) {
            warning("Some of the amino acid sequences in MSA '", msa.file, 
                "' are not in the database of coding sequences ('all.cds')!")
        }
        gene.fam.wd <- createGeneFamilyWorkingDirIfNotExists(working.dir, 
            gene.fam)
        write.fasta(msa, names(msa), file.path(gene.fam.wd, paste0(gene.fam, 
            "_AA_MSA_orig_gene_ids.fa")))
        fam.gene.id.maps <- geneNameMappings(gene.fam.wd, gene.fam)
        names(msa) <- mapOriginalToSanitizedNames(names(msa), fam.gene.id.maps)
        write.fasta(msa, names(msa), file.path(gene.fam.wd, paste0(gene.fam, 
            "_AA_MSA.fa")))
        NULL
    })
    TRUE
}

#' Loads all phylogenetic trees for the respective gene families reconstructed
#' by Orthofinder. They are read in using \code{ape::read.tree} and gene
#' identifiers are sanitized, removing the species name prefixes added by
#' Orthofinder. Also node support values are removed to prepare the tree for
#' usage with HyPhy. See function \code{sanitizeOrthofinderGeneIds} for more
#' details. Note, that this function uses \code{mclapply}.
#'
#' @param path.2.trees.dir - The valid file path to the directory the MSAs are
#' stored in. Make sure their file-extensions are '.fa'.
#' @param working.dir - The valid file path to the directory in which to store
#' the sanitized and parsed MSAs. Each will be stored in its respective gene
#' family's sub-directory, identified by the family's name itself.
#' @param gene.families.names - A character vector of gene family names that
#' are to be expected to contain the file names in argument 'path.2.trees.dir'
#' excluding the file-extensions. Default is this package's
#' \code{names(gene.families)}.
#'
#' @return TRUE if and only if no error has occurred.
#' @export
loadAndSanitizeTrees <- function(path.2.trees.dir, working.dir, gene.families.names = names(gene.families)) {
    tree.fls <- system(paste("ls", path.2.trees.dir, "| grep -P \"^OG.*\\_tree.txt$\""), 
        intern = TRUE)
    names(tree.fls) <- sub("\\_tree.txt$", "", tree.fls)
    if (!all(names(tree.fls) %in% gene.families.names)) {
        warning("Some of the phylogenetic tree file names do not appear in the argument 'gene.families.names'.")
    }
    mclapply(names(tree.fls), function(gene.fam) {
        tree.file <- tree.fls[[gene.fam]]
        fam.tree <- read.tree(file.path(path.2.trees.dir, tree.file))
        fam.tree$tip.label <- sanitizeOrthofinderGeneIds(fam.tree$tip.label)
        if (!all(fam.tree$tip.label %in% names(all.cds))) {
            warning("Some of the gene identifier in tree '", tree.file, 
                "' are not in the database of coding sequences ('all.cds')!")
        }
        fam.tree$node.label <- NULL
        gene.fam.wd <- createGeneFamilyWorkingDirIfNotExists(working.dir, 
            gene.fam)
        write.tree(fam.tree, file.path(gene.fam.wd, paste0(gene.fam, 
            "_phyl_tree_orig_gene_ids.newick")))
        fam.gene.id.maps <- geneNameMappings(gene.fam.wd, gene.fam)
        fam.tree$tip.label <- mapOriginalToSanitizedNames(fam.tree$tip.label, 
            fam.gene.id.maps)
        write.tree(fam.tree, file.path(gene.fam.wd, paste0(gene.fam, 
            "_phyl_tree.newick")))
        NULL
    })
    TRUE
}

#' Generates the coding sequence alignments for all families found in
#' sub-directories of argument 'working.dir'.
#'
#' @param working.dir - The valid file path to the directory in which to store
#' the sanitized and parsed MSAs. Each will be stored in its respective gene
#' family's sub-directory, identified by the family's name itself.
#' @param gene.families.names - A character vector of gene family names that
#' are to be expected to contain the sub-directory names in argument
#' 'working.dir' excluding the file-extensions. Default is this package's
#' \code{names(gene.families)}.
#'
#' @return TRUE if and only if no error has occurred.
#' @export
generateCodingSequenceMSAs <- function(working.dir, gene.families.names = names(gene.families)) {
    gene.fams.dirs <- system(paste("ls", working.dir, "| grep -P \"^OG\""), 
        intern = TRUE)
    if (!all(gene.fams.dirs %in% gene.families.names)) {
        warning("Some of the gene families stored in working directory '", 
            working.dir, "' do not appear in the argument gene.families.names!")
    }
    mclapply(gene.fams.dirs, function(gene.fam) {
        gene.fam.wd <- normalizePath(file.path(working.dir, gene.fam))
        fam.cds <- all.cds[gene.families[[gene.fam]]]
        fam.aa.msa <- seqinr::read.fasta(file.path(gene.fam.wd, paste0(gene.fam, 
            "_AA_MSA_orig_gene_ids.fa")), seqtype = "AA", as.string = TRUE, 
            strip.desc = TRUE)
        fam.cds.msa <- GeneFamilies::alignCDSSetWithAlignedAAsAsGuide(fam.cds, 
            fam.aa.msa)
        fam.gene.id.maps <- geneNameMappings(gene.fam.wd, gene.fam)
        names(fam.cds.msa) <- mapOriginalToSanitizedNames(names(fam.cds.msa), 
            fam.gene.id.maps)
        seqinr::write.fasta(fam.cds.msa, names(fam.cds.msa), file.path(gene.fam.wd, 
            paste0(gene.fam, "_CDS_MSA.fa")))
        NULL
    })
    TRUE
}

#' Generates the MEME HyPhy batch-files for all families found in
#' sub-directories of argument 'working.dir'.
#'
#' @param working.dir - The valid file path to the directory in which to store
#' the sanitized and parsed MSAs. Each will be stored in its respective gene
#' family's sub-directory, identified by the family's name itself.
#' @param hyphy.batch.files.dir - The valid file path to the directory of the
#' HyPhy installation in which to find the HyPhy batch files.
#' @param gene.families.names - A character vector of gene family names that
#' are to be expected to contain the sub-directory names in argument
#' 'working.dir' excluding the file-extensions. Default is this package's
#' \code{names(gene.families)}.
#'
#' @return TRUE if and only if no error has occurred.
#' @export
generateMemeBatchFiles <- function(working.dir, hyphy.batch.files.dir, 
    gene.families.names = names(gene.families)) {
    hyphy.batch.files.dir <- normalizePath(hyphy.batch.files.dir)
    gene.fams.dirs <- system(paste("ls", working.dir, "| grep -P \"^OG\""), 
        intern = TRUE)
    if (!all(gene.fams.dirs %in% gene.families.names)) {
        warning("Some of the gene families stored in working directory '", 
            working.dir, "' do not appear in the argument gene.families.names!")
    }
    mclapply(gene.fams.dirs, function(gene.fam) {
        gene.fam.wd <- normalizePath(file.path(working.dir, gene.fam))
        fam.cds.msa.path <- file.path(gene.fam.wd, paste0(gene.fam, 
            "_CDS_MSA.fa"))
        fam.tree.no.node.labels.path <- file.path(gene.fam.wd, paste0(gene.fam, 
            "_phyl_tree.newick"))
        fam.hyphy.meme.log.path <- file.path(gene.fam.wd, paste0(gene.fam, 
            "_HyPhy_MEME_log.txt"))
        fam.hyphy.meme.output.path <- file.path(gene.fam.wd, paste0(gene.fam, 
            "_HyPhy_MEME_output.txt"))
        fam.hyphy.meme.batch.file.path <- file.path(gene.fam.wd, paste0(gene.fam, 
            "_HyPhy_MEME_input.bf"))
        brew(text = GeneFamilies::hyphy.meme.bf, output = fam.hyphy.meme.batch.file.path)
        NULL
    })
    TRUE
}

#' Within the argument 'work.dir' all gene family work dirs are looked up that
#' have HyPhy MEME result files. These are parsed and sites showing significant
#' evidence for positive selection are identified. Subsequently leaves of the
#' phylogenetic tree (proteins) that show significant evidence for positive
#' selection at the identified sites are also identified.
#'
#' @param work.dir - The valid file path to the directory in which each
#' family's MEME results are stored.
#' @param p.adjust.method - The argument passed to \code{base::p.adjust} as the
#' method to apply for correcting p-values for multiple hypotheses testing.
#' @param pval.cutoff - The significance level for corrected p-values. Default
#' is \code{0.05}.
#' @param empirical.bayes.factor.cutoff - The significance level to be apliad
#' to empirical Bayes Factor values. Default is \code{100}.
#' @param phyl.tree.leaf.regex - The regular expression to be applied to
#' distinguish inner tree nodes from leaf nodes. Default is
#' \code{'^PROT\\d+$'}.
#' @param report.original.gene.names - Boolean switch to control whether the
#' genes found to be subject to positive selection whall be reported with their
#' original and not their sanitized names. See
#' \code{SolycGeneFams::mapOriginalToSanitizedNames} and
#' \code{SolycGeneFams::geneNameMappings} for more details. Default of this
#' argument is \code{TRUE}.
#'
#' @return An instance of \code{base::list} with names the families for which
#' MEME results were found and values lists with two entries 'MEME.sites' and
#' 'MEME.branches', the respectively found subjects to positive selection
#' passing the significance levels.
#' @export
readMemeResults <- function(work.dir, p.adjust.method = "fdr", pval.cutoff = 0.05, 
    empirical.bayes.factor.cutoff = 100, phyl.tree.leaf.regex = "^PROT\\d+$", 
    report.original.gene.names = TRUE, gene.family.sizes = gene.families.sizes) {
    norm.wd <- normalizePath(work.dir)
    fams.meme.outs <- system(paste("find", norm.wd, "-type f", "-name '*_HyPhy_MEME_output.txt'"), 
        intern = TRUE)
    fams.meme.branches <- paste0(fams.meme.outs, ".branches")
    fam.nms <- sub("^.*/", "", sub("_HyPhy_MEME_output.txt$", "", 
        fams.meme.outs))
    names(fams.meme.outs) <- fam.nms
    names(fams.meme.branches) <- fam.nms
    setNames(mclapply(fam.nms, function(fam) {
        m.o <- read.table(fams.meme.outs[[fam]], sep = ",", header = TRUE, 
            stringsAsFactors = FALSE)
        m.o$p.adj <- p.adjust(m.o$p.value, method = p.adjust.method)
        m.sites <- which(m.o$p.adj <= pval.cutoff)
        if (length(m.sites) > 0) {
            m.sites.tbl <- data.frame(MEME.site = m.sites, P.val = m.o[m.sites, 
                "p.value"], P.adj = m.o[m.sites, "p.adj"], stringsAsFactors = FALSE)
            m.b <- read.table(fams.meme.branches[[fam]], sep = ",", 
                header = TRUE, stringsAsFactors = FALSE, skip = 2)
            m.branches <- unique(m.b[which(m.b$Site %in% m.sites & 
                m.b$EmpiricalBayesFactor >= empirical.bayes.factor.cutoff & 
                grepl(phyl.tree.leaf.regex, m.b$Branch)), c("Site", 
                "Branch")])
            meme.branches.df <- if (nrow(m.branches) > 0) {
                fam.wd <- sub("/OG\\d+_HyPhy_MEME_output.txt$", "", 
                  fams.meme.outs[[fam]])
                if (report.original.gene.names) {
                  m.branches$Branch <- mapSanitizedToOriginalNames(m.branches$Branch, 
                    geneNameMappings(fam.wd, fam))
                  fam.aa.msa.file <- normalizePath(file.path(fam.wd, 
                    paste0(fam, "_AA_MSA_orig_gene_ids.fa")))
                } else {
                  fam.aa.msa.file <- normalizePath(file.path(fam.wd, 
                    paste0(fam, "_AA_MSA.fa")))
                }
                m.b.unq <- unique(m.branches$Branch)
                fam.aa.msa <- read.fasta(fam.aa.msa.file, seqtype = "AA", 
                  as.string = TRUE, strip.desc = TRUE)
                Reduce(rbind, lapply(m.b.unq, function(meme.branch) {
                  align.sites <- sort(unique(m.branches[which(m.branches$Branch == 
                    meme.branch), "Site"]))
                  Reduce(rbind, lapply(align.sites, function(align.site) {
                    unalign.site <- unalignedAAforAlignedAAPos(meme.branch, 
                      align.site, fam.aa.msa)
                    mst.i <- which(m.sites.tbl$MEME.site == align.site)
                    meme.site.p.val <- m.sites.tbl[[mst.i, "P.val"]]
                    meme.site.p.adj <- m.sites.tbl[[mst.i, "P.adj"]]
                    data.frame(Protein = meme.branch, Gene.Family = fam, 
                      Gene.Family.Size = gene.family.sizes[[fam]], 
                      unaligned.pos.sel.codon = unalign.site, aligned.pos.sel.codon = align.site, 
                      MEME.site.p.value = meme.site.p.val, MEME.site.p.adj = meme.site.p.adj, 
                      stringsAsFactors = FALSE)
                  }))
                }))
            } else NULL
            list(MEME.sites = m.sites.tbl, MEME.branches = meme.branches.df)
        } else NULL
    }), fam.nms)
}

#' Parses the output table of a HMMER3 run with command line switch
#' \code{--domtblout}.
#'
#' @param path.2.hmmer3.domtbl.out - Valid file path to the output of hmmer3
#' run with the command line switch \code{--domtblout}.
#'
#' @return An instance of \code{data.frame} with 23 columns containing the
#' output of HMMER3.
#' @export
parseHmmer3DomTableOut <- function(path.2.hmmer3.domtbl.out) {
    dom.tbl <- fread(cmd = paste("sed", "-e '/^#/d'", paste(rep("-e 's/ \\+/\\t/'", 
        22), collapse = " "), path.2.hmmer3.domtbl.out), data.table = FALSE, 
        sep = "\t", header = FALSE)
    colnames(dom.tbl) <- c("target.name", "target.accession", "tlen", 
        "query.name", "query.accession", "qlen", "full.sequence.E-value", 
        "full.sequence.score", "full.sequence.bias", "this.domain.#", 
        "this.domain.of", "this.domain.c-Evalue", "this.domain.i-Evalue", 
        "this.domain.score", "this.domain.bias", "hmm.coord.from", 
        "hmm.coord.to", "ali.coord.from", "ali.coord.to", "env.coord.from", 
        "env.coord.to", "acc", "description.of.target")
    dom.tbl
}


#' Generates an interactive plot of the family's multiple amino acid sequence
#' alignment. The output is an HTML file which requires this packages's
#' \code{./inst/assets} folder to be at the same location as the output file. A
#' modified version of the Javascript library
#' \code{github.com/AndrewCRMartin/JSAV} is used to display the alignment.
#' Warning: There are more beatifull ways to interactively show an alignment.
#'
#' @param path.2.msa.fasta - The valid file path to the family's multiple amino
#' acid alignment in fasta format. The filename is also used to extract the
#' family's name and the name of the output file.
#' @param meme.branches.tbl - An instance of \code{base::data.frame} holding
#' the results of HyPhy's MEME analysis (branches). See
#' \code{SlydGeneFamsAnalyses::readMemeResults} for more details.
#' @param pfam.hmmer3.domtbl - An instance of \code{base::data.frame} holding
#' the results of running HMMER3 with the \code{--domtblout} command line
#' switch. Use \code{SlydGeneFamsAnalyses::parseHmmer3DomTableOut} to read in
#' the resulting table.
#' @param output.dir - The directory in which to write the resulting HTML file.
#' @param brew.template - The brew template to be used to render the HTML file.
#' Default is
#' \code{normalizePath(file.path(path.package("SlydGeneFamsAnalyses"),
#' "msa_template_brew.html"))}.
#'
#' @return The result of invoking \code{brew::brew(...)}.
#' @export
generateInteractiveMsaPlot <- function(path.2.msa.fasta, meme.branches.tbl, 
    pfam.hmmer3.domtbl, output.dir, brew.template = normalizePath(file.path(path.package("SlydGeneFamsAnalyses"), 
        "msa_template_brew.html"))) {
    aa.msa <- read.fasta(path.2.msa.fasta, seqtype = "AA", strip.desc = TRUE, 
        as.string = TRUE)
    # Javascript starts counting at 0, not 1:
    meme.branches.tbl$aligned.pos.sel.codon <- meme.branches.tbl$aligned.pos.sel.codon - 
        1
    pfam.hmmer3.domtbl$aligned.ali.coord.from <- c()
    pfam.hmmer3.domtbl$aligned.ali.coord.to <- c()
    for (i in 1:nrow(pfam.hmmer3.domtbl)) {
        gene.id <- pfam.hmmer3.domtbl$query.name[[i]]
        # Javascript starts counting at 0, not 1:
        pfam.hmmer3.domtbl$aligned.ali.coord.from[[i]] <- alignedForUnalignedAAPos(gene.id, 
            pfam.hmmer3.domtbl$ali.coord.from[[i]], aa.msa) - 1
        pfam.hmmer3.domtbl$aligned.ali.coord.to[[i]] <- alignedForUnalignedAAPos(gene.id, 
            pfam.hmmer3.domtbl$ali.coord.to[[i]], aa.msa) - 1
    }
    fam.name <- sub("\\.[^.]+$", "", sub("^.*/", "", path.2.msa.fasta))
    out.html <- sub("\\.[^.]+$", ".html", sub("^.*/", "", path.2.msa.fasta))
    brew(file = brew.template, output = file.path(output.dir, out.html))
}
