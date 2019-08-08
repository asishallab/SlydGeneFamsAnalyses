.onLoad <- function(libname = find.package("SlydGeneFamsAnalyses"), 
    pkgname = "SlydGeneFamsAnalyses") {
    data("all_CDS", package = "SlydGeneFamsAnalyses")
    data("geneFamilies", package = "SlydGeneFamsAnalyses")
    data("valdarScores", package = "SlydGeneFamsAnalyses")
    data("interProDb", package = "SlydGeneFamsAnalyses")
    data("valdarScores", package = "SlydGeneFamsAnalyses")
    data("H6H", package = "SlydGeneFamsAnalyses")
    data("allInterProAnnotations", package = "SlydGeneFamsAnalyses")
#
H6H.RData
allInterProAnnotations.RData
allProteomesInterProAnnos.RData

    # data("mappResults", package = "SlydGeneFamsAnalyses")
    # Define 'all.cds' constant:
    if (exists("cag.cds") && exists("cam.cds") &&
        exists("date.cds") && exists("dati.cds") &&
        exists("nat.cds") && exists("nsy.cds") &&
        exists("nta.cds") && exists("nto.cds") &&
        exists("pi.cds") && exists("sly.cds") &&
        exists("spe.cds") && exists("spi.cds")&&
        exists("stu.cds")) {
        assign("all.cds", c(cag.cds, cam.cds, date.cds,
            dati.cds, nat.cds, nsy.cds, nta.cds, nto.cds,
            pi.cds, sly.cds, spe.cds, spi.cds, stu.cds),
            envir = globalenv())
    }
}
