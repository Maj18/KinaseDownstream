\name{runPTMSEA}
\alias{runPTMSEA}
\title{
Run PTMSEA analysis from limma results
}
\description{
Prepares input data from limma results, downloads the PTMsigDB gene set database, and runs PTMSEA analysis using ssGSEA2. Writes GCT files and PTMSEA output to the specified directory.
}
\usage{
runPTMSEA(limma_rslt, PTMSEA_OUTDIR)
}
\arguments{
  \item{limma_rslt}{A list of data frames, each containing phosphosite-level statistics from limma analysis.}
  \item{PTMSEA_OUTDIR}{Directory path for saving PTMSEA input and output files.}
}
\details{
This function formats limma results for PTMSEA, writes a GCT file, downloads the PTMsigDB gene set file if needed, and runs ssGSEA2 for PTMSEA analysis. The function is designed to automate the PTMSEA workflow for phosphoproteomics data.
}
\value{
A data frame containing the formatted input for PTMSEA (GCT format).
}
\references{
}
\author{
Yuan Li \email{yuan.li@nbis.se}
}
\note{
Requires the \code{dplyr}, \code{tibble}, \code{cmapR}, and \code{ssGSEA2} packages. Downloads PTMsigDB from the official repository if not present.
}

\examples{
## Not run:
## Example usage (requires appropriate input data):
# runPTMSEA(limma_rslt, PTMSEA_OUTDIR = "results/PTMSEA/")

## The function is currently defined as
function (limma_rslt, PTMSEA_OUTDIR) 
{
    limma_rslt_4gct = lapply(seq_along(limma_rslt), function(i) {
        rslt = limma_rslt[[i]] \%>\% arrange(-logFC)
        rownames(rslt) = NULL
        temp = rslt[, c("PTM.FlankingRegion", "logFC")]
        colnames(temp) = c("PTM.FlankingRegion", names(limma_rslt)[i])
        temp$PTM.FlankingRegion = paste0(temp$PTM.FlankingRegion, 
            "_p")
        temp
    }) \%>\% Reduce(full_join, .) \%>\% .[!duplicated(.$PTM.FlankingRegion), 
        ]
    rownames(limma_rslt_4gct) = NULL
    limma_rslt_4gct = tibble::column_to_rownames(limma_rslt_4gct, 
        var = "PTM.FlankingRegion")
    write_gct(as.matrix(limma_rslt_4gct) \%>\% cmapR::new("GCT", 
        mat = .), paste0(PTMSEA_OUTDIR, "/PTM-SEA"), precision = 2)
    input_gct_file = list.files(path = PTMSEA_OUTDIR, pattern = "PTM-SEA", 
        full.names = TRUE)
    download.file(url = "https://raw.githubusercontent.com/nicolerg/ssGSEA2/refs/heads/master/db/ptmsigdb/ptm.sig.db.all.flanking.human.v2.0.0.gmt", 
        destfile = paste0(PTMSEA_OUTDIR, "/ptm.sig.db.all.flanking.human.v2.0.0.gmt"))
    set.seed(123)
    invisible(capture.output(ssGSEA2::run_ssGSEA2(input_gct_file, 
        output.prefix = "PTMSEA_OUTPUT", gene.set.databases = paste0(PTMSEA_OUTDIR, 
            "/ptm.sig.db.all.flanking.human.v2.0.0.gmt"), output.directory = PTMSEA_OUTDIR, 
        sample.norm.type = "none", weight = 1, correl.type = "rank", 
        statistic = "area.under.RES", spare.cores = 4, output.score.type = "NES", 
        nperm = 1000, min.overlap = 5, extended.output = TRUE, 
        global.fdr = FALSE, export.signat.gct = T, param.file = T, 
        log.file = paste0(PTMSEA_OUTDIR, "/run.log"))))
    return(limma_rslt_4gct)
  }
}
\keyword{PTMSEA}
\keyword{phosphoproteomics}
\keyword{enrichment}
\keyword{ssGSEA2}
\keyword{data processing}