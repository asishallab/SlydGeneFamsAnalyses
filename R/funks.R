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

#' Uses Bioconductor's 'msa' package to pretty print the alignment of the
#' argument gene family highlighting the genes and amino acids found to be
#' subject to positive selection.
#'
#' @param fam.fir - The path to the working directory of the argument gene
#' family. This dir should contain the multiple sequence alignment file and
#' will be used to write the alignment plots into.
#' @param aa.msa.file - The name, not path, of the multiple sequence alignment
#' file to be used. It is of key importance that the protein identifiers in the
#' alignment file are identical to those in the argument
#' 'selected.pos.and.genes.tbl' column 'Protein'.
#' @param selected.pos.and.genes.tbl - An instance of \code{base::data.frame}
#' holding information about positively selected proteins, aligned amino acid
#' positions and their corresponding unaligned positions. See function
#' \code{SlydGeneFamsAnalyses::readMemeResults} for details on how this table
#' can be generated.
#' @param pfam.tbl - An instance of \code{base::data.frame} the result from
#' invoking \code{SlydGeneFamsAnalyses::parseHmmer3DomTableOut}. If NULL this
#' will be ignored. If set and some protein domains fall within the displayed
#' region of the MSA the domains will be tinted. See CTAN texshade
#' documentation for more details.
#' @param plot.pdf.prefix - The file name with file extension ('.pdf') of the
#' plot that will be generated and saved in the argument 'fam.dir'. Default is
#' \code{sub('.fa', '', aa.msa.file, fixed=TRUE)}.
#' @param before.sel.pos.offset - The offset before the first selected position
#' to be included in the plot. Default is \code{5}.
#' @param after.sel.pos.offset - The offset after the first selected position
#' to be included in the plot. Default is \code{5}.
#'
#' @return The result of invoking \code{msa::msaPrettyPrint(...)}.
#' @export
printAaMsaWithSelection <- function(fam.dir, aa.msa.file, selected.pos.and.genes.tbl, 
    pfam.tbl, plot.pdf.prefix = sub(".fa", "", aa.msa.file, fixed = TRUE), 
    before.sel.pos.offset = 5, after.sel.pos.offset = 5) {
    begin.codon <- min(selected.pos.and.genes.tbl$aligned.pos.sel.codon) - 
        before.sel.pos.offset
    end.codon <- max(selected.pos.and.genes.tbl$aligned.pos.sel.codon) + 
        after.sel.pos.offset
    msa.fasta.path <- normalizePath(file.path(fam.dir, aa.msa.file))
    m <- Biostrings::readAAMultipleAlignment(msa.fasta.path)
    aa.msa.prot.ids <- names(unmasked(m))
    gene.name.color.tex <- paste(unlist(lapply(selected.pos.and.genes.tbl$Protein, 
        function(sel.gene) {
            i.gene <- which(aa.msa.prot.ids == sel.gene)
            paste0("\\namecolor{", i.gene, "}{Magenta}")
        })), collapse = "\n")
    
    frameblocks.emphregions.tex <- paste(unlist(lapply(selected.pos.and.genes.tbl$Protein, 
        function(sel.gene) {
            gene.no <- which(aa.msa.prot.ids == sel.gene)
            unalgnd.pos <- sort(unique(selected.pos.and.genes.tbl[which(selected.pos.and.genes.tbl$Protein == 
                sel.gene), "unaligned.pos.sel.codon"]))
            framebl.tex <- paste0("\\frameblock{", gene.no, "}{", 
                paste(sapply(unalgnd.pos, function(u.p) {
                  paste(u.p, u.p, sep = "..")
                }), collapse = ","), "}{Cyan[1pt]}")
            emphreg.tex <- paste0("\\emphregion{", gene.no, "}{", 
                paste(sapply(unalgnd.pos, function(u.p) {
                  paste(u.p, u.p, sep = "..")
                }), collapse = ","), "}")
            paste(framebl.tex, emphreg.tex, sep = "\n")
        })), collapse = "\n")
    
    tinted.tex <- ""
    feature.tex <- ""
    caption.tex <- ""
    if (!is.null(pfam.tbl)) {
        msa.fa <- read.fasta(msa.fasta.path, seqtype = "AA", strip.desc = TRUE, 
            as.string = TRUE)
        prot.dom.lst <- list()
        tinted.tex <- paste(unlist(lapply(1:nrow(pfam.tbl), function(i) {
            prot <- pfam.tbl[[i, "query.name"]]
            gene.no <- which(aa.msa.prot.ids == prot)
            unalign.start <- pfam.tbl[[i, "ali.coord.from"]]
            unalign.end <- pfam.tbl[[i, "ali.coord.to"]]
            aligned.start <- alignedForUnalignedAAPos(prot, unalign.start, 
                msa.fa)
            aligned.end <- alignedForUnalignedAAPos(prot, unalign.end, 
                msa.fa)
            if (prot %in% aa.msa.prot.ids && aligned.start >= begin.codon && 
                aligned.end <= end.codon) {
                pdl.entry <- paste(pfam.tbl[[i, "target.accession"]], 
                  pfam.tbl[[i, "this.domain.#"]], sep = " No.")
                prot.dom.lst[[pdl.entry]] <<- c(min(prot.dom.lst[[pdl.entry]], 
                  aligned.start), max(prot.dom.lst[[pdl.entry]], aligned.end))
                paste0("\\tintregion{", gene.no, "}{", unalign.start, 
                  "..", unalign.end, "}")
            } else NULL
        })), collapse = "\n")
        first.prot <- aa.msa.prot.ids[[1]]
        feat.pos <- c("top", "ttop", "tttop", "ttttop", "bottom", 
            "bbottom", "bbbottom", "bbbbottom")
        feat.pos.i <- 0
        feature.tex <- paste(unlist(lapply(names(prot.dom.lst), function(prot.dom) {
            first.prot.unalign.start <- unalignedAAforAlignedAAPos(first.prot, 
                prot.dom.lst[[prot.dom]][[1]], msa.fa)
            first.prot.unalign.end <- unalignedAAforAlignedAAPos(first.prot, 
                prot.dom.lst[[prot.dom]][[2]], msa.fa)
            feat.pos.i <<- if (feat.pos.i + 1 > length(feat.pos)) {
                1
            } else {
                feat.pos.i + 1
            }
            paste0("\\feature{", feat.pos[[feat.pos.i]], "}{1}{", 
                first.prot.unalign.start, "..", first.prot.unalign.end, 
                "}{brace}{", prot.dom, "}")
        })), collapse = "\n")
        caption.tex <- paste0("\\showcaption[bottom]{Protein Domains:\\\\", 
            paste(unlist(lapply(unique(names(prot.dom.lst)), function(prot.dom.nm) {
                prot.dom.acc <- sub(" No\\.\\d+$", "", prot.dom.nm)
                prot.dom.coords <- prot.dom.lst[[prot.dom.nm]]
                paste0(prot.dom.nm, " (", prot.dom.coords[[1]], "-", 
                  prot.dom.coords[[2]], ") - ", pfam.tbl[which(pfam.tbl$target.accession == 
                    prot.dom.acc), "description.of.target"][[1]])
            })), collapse = "\\\\"), "}")
    }
    msa::msaPrettyPrint(m, c(begin.codon, end.codon), file = file.path(fam.dir, 
        paste0(plot.pdf.prefix, "_", begin.codon, "-", end.codon, 
            ".pdf")), shadingMode = "functional", askForOverwrite = FALSE, 
        shadingModeArg = "chemical", output = "pdf", furtherCode = paste("\\tintdefault{weak}", 
            gene.name.color.tex, frameblocks.emphregions.tex, tinted.tex, 
            caption.tex, caption.tex, feature.tex, sep = "\n"), showNumbering = "none")
}
