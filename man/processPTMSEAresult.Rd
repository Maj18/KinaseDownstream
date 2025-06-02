\name{processPTMSEAresult}
\alias{processPTMSEAresult}
\title{
Process and visualize PTMSEA results
}
\description{
Reads PTMSEA result files, prepares it for network analyses, and generates dot plots for significant kinase sets.
}
\usage{
processPTMSEAresult(PTMSEA_FILE_PATH, output.score.type = "NES", sig.thresh = 0.05)
}
\arguments{
  \item{PTMSEA_FILE_PATH}{Path to the PTMSEA result file (GCT format) to be processed.}
  \item{output.score.type}{Score type to use for output columns (default is \code{"NES"}).}
  \item{sig.thresh}{Significance threshold for filtering results (default is 0.05).}
}
\details{
This function reads a PTMSEA result file, prepares it for network analyses, and generates a dot plot of significant kinase sets using ggplot2. The function returns a list containing the processed PTMSEA results and the data used for plotting.
}
\value{
A list with components:
  \item{ptmsea_rslt}{A named list of data frames, each containing PTMSEA results for a comparison.}
  \item{gct4plot}{A data frame containing the filtered and summarized results used for plotting.}
}
\references{
}
\author{
Yuan Li \email{yuan.li@nbis.se}
}
\note{
Requires the \code{dplyr}, \code{tidyr}, \code{ggplot2}, and \code{xlsx} packages.
}

\examples{
## Not run:
## Example usage (requires appropriate input data):
# result <- processPTMSEAresult("PTMSEA_results.gct")
# head(result$gct4plot)

## The function is currently defined as
function (PTMSEA_FILE_PATH, output.score.type = "NES", sig.thresh = 0.05) 
{
    gct_data = read.delim(PTMSEA_FILE_PATH, skip = 2, header = TRUE, 
        sep = "\t", check.names = FALSE)
    pairs = grep("fdr.pvalue", colnames(gct_data), value = T) \%>\% 
        gsub("fdr.pvalue.", "", .)
    colnames(gct_data)[colnames(gct_data) \%in\% pairs] = paste0(output.score.type, 
        "_", pairs)
    xlsx::write.xlsx(gct_data, gsub("gct", "xlsx", PTMSEA_FILE_PATH), 
        append = F, sheetName = "PTMSEAonDiff", row.names = F)
    ptmsea_rslt = lapply(pairs, function(Pair) {
        temp = gct_data[, c("id", "Signature.set.description", 
            "Signature.set.size", grep(Pair, colnames(gct_data), 
                value = T))]
        colnames(temp) = gsub(paste0("[.]", Pair), "", colnames(temp)) \%>\% 
            gsub(paste0("_", Pair), "", .)
        temp
    }) \%>\% setNames(pairs)
    gct4plot_NES = gct_data[, c("id", grep("NES_", colnames(gct_data), 
        value = T))]
    colnames(gct4plot_NES) = gsub("NES_", "", colnames(gct4plot_NES))
    gct4plot_fdr.pvalue = gct_data[, c("id", grep("fdr.pvalue.", 
        colnames(gct_data), value = T))]
    colnames(gct4plot_fdr.pvalue) = gsub("fdr.pvalue.", "", colnames(gct4plot_fdr.pvalue))
    gct4plot = full_join(gather(gct4plot_NES, key = "Pair", value = "NES", 
        2:4), gather(gct4plot_fdr.pvalue, key = "Pair", value = "fdr.pvalue", 
        2:4)) \%>\% filter(fdr.pvalue < 0.05) \%>\% arrange(-NES)
    gct4plot$Regulation = ifelse(gct4plot$NES > 0, "Up", "Down")
    gct4plot$Pair = factor(gct4plot$Pair)
    gct4plot$id = factor(gct4plot$id, levels = unique(rev(gct4plot$id)))
    ggplot(gct4plot, aes(x = NES, y = id)) + geom_point(shape = 21, 
        aes(size = fdr.pvalue, fill = Pair)) + facet_wrap(~Pair, 
        ncol = 3) + ggtitle("PTM-SEA_PTMsigDB2") + xlab("NES") + 
        ylab("") + labs(size = "fdr.pvalue") + theme_bw(base_size = 15) + 
        scale_size(trans = "reverse") + guides(alpha = "none", 
        size = guide_legend(override.aes = list(shape = 21))) + 
        theme(strip.text.x = element_text(size = 14)) + labs(caption = paste0("fdr.pvalueCutoff = ", 
        sig.thresh, ", pAdjustMethod = fdr"))
    return(list(ptmsea_rslt = ptmsea_rslt, gct4plot = gct4plot))
  }
}
\keyword{PTMSEA}
\keyword{visualization}
\keyword{kinase}
\keyword{phosphoproteomics}
\keyword{data processing}
