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
    fams.splt <- strsplit(fams.raw, 
        ":")
    fam.names <- unlist(lapply(fams.splt, 
        function(x) x[[1]]))
    setNames(lapply(fams.splt, 
        function(fam.str) {
            setdiff(strsplit(fam.str[[2]], 
                " ")[[1]], "")
        }), fam.names)
}
