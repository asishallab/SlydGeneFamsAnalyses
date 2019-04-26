.onLoad <- function(libname = find.package("SlydGeneFamsAnalyses"), 
    pkgname = "SlydGeneFamsAnalyses") {
    data("all_CDS", package = "SlydGeneFamsAnalyses")
    # Define 'all.cds' constant:
    if (exists("ath.cds") && exists("slyc.cds") && 
        exists("slyd.cds") && exists("spe.cds") && 
        exists("stu.cds") && exists("vvi.cds")) {
        assign("all.cds", c(ath.cds, 
            slyc.cds, slyd.cds, 
            spe.cds, stu.cds, vvi.cds), 
            envir = globalenv())
    }
}
