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
            "_phyl_tree.newick")))
        NULL
    })
    TRUE
}

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
            "_AA_MSA.fa")), seqtype = "AA", as.string = TRUE, strip.desc = TRUE)
        fam.cds.msa <- GeneFamilies::alignCDSSetWithAlignedAAsAsGuide(fam.cds, 
            fam.aa.msa)
        seqinr::write.fasta(fam.cds.msa, names(fam.cds.msa), file.path(gene.fam.wd, 
            paste0(gene.fam, "_CDS_MSA.fa")))
        NULL
    })
    TRUE
}
