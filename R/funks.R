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
        fam.tree$tip.label <- mapOriginalToSanitizedNames(fam.tree$tip.label, fam.gene.id.maps)
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
            "_AA_MSA_orig_gene_ids.fa")), seqtype = "AA", as.string = TRUE, strip.desc = TRUE)
        fam.cds.msa <- GeneFamilies::alignCDSSetWithAlignedAAsAsGuide(fam.cds, 
            fam.aa.msa)
        fam.gene.id.maps <- geneNameMappings(gene.fam.wd, gene.fam)
        names(fam.cds.msa) <- mapOriginalToSanitizedNames(names(fam.cds.msa), fam.gene.id.maps)
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
