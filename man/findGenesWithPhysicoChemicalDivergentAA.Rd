% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funks.R
\name{findGenesWithPhysicoChemicalDivergentAA}
\alias{findGenesWithPhysicoChemicalDivergentAA}
\title{For a given MAPP result table find those argument genes of interest that
have a significantly divergent amino acid at the positions that MAPP
identified as divergent.}
\usage{
findGenesWithPhysicoChemicalDivergentAA(mapp.tbl, fam.aa.msa, fam.name,
  genes.of.interest = names(slyd.cds), p.adjusted.cutoff = 0.05,
  mapp.aa.p.val.cols = setNames(c("A.1.adj", "C.1.adj", "D.1.adj",
  "E.1.adj", "F.1.adj", "G.1.adj", "H.1.adj", "I.1.adj", "K.1.adj",
  "L.1.adj", "M.1.adj", "N.1.adj", "P.1.adj", "Q.1.adj", "R.1.adj",
  "S.1.adj", "T.1.adj", "V.1.adj", "W.1.adj", "Y.1.adj"), c("A", "C", "D",
  "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T",
  "V", "W", "Y")))
}
\arguments{
\item{mapp.tbl}{- The result of invoking
\code{SlydGeneFamsAnalyses::readMappResult}.}

\item{fam.aa.msa}{- The result of invoking
\code{SlydGeneFamsAnalyses::readMultipleSequenceAlignmentAsMatrix}.}

\item{fam.name}{- The name of the gene family}

\item{genes.of.interest}{- A character vector of gene identifier. Set to
rownames(fam.aa.msa) if you want to get results for all genes in the
alignment. Default is \code{names(slyd.cds)}.}

\item{p.adjusted.cutoff}{- The cutoff for the adjusted P values to infer
significancy. Default is \code{.05}.}

\item{mapp.aa.p.val.cols}{- A named character or integer vector identifying
those columns of argument 'mapp.tbl' in which to lookup the respective amino
acid P values. Default is \code{setNames(c('A.1.adj', 'C.1.adj', 'D.1.adj',
'E.1.adj', 'F.1.adj', 'G.1.adj', 'H.1.adj', 'I.1.adj', 'K.1.adj', 'L.1.adj',
'M.1.adj', 'N.1.adj', 'P.1.adj', 'Q.1.adj', 'R.1.adj', 'S.1.adj', 'T.1.adj',
'V.1.adj', 'W.1.adj', 'Y.1.adj'), c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I',
'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'))}.}
}
\value{
An instance of \code{base::data.frame} with the following columns:
'Protein', 'Site', 'Divergent.AA', and 'AA.p.value.adj'. Returns NULL if no
matches were found.
}
\description{
For a given MAPP result table find those argument genes of interest that
have a significantly divergent amino acid at the positions that MAPP
identified as divergent.
}
