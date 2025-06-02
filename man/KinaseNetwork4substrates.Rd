\name{KinaseNetwork4substrates}
\alias{KinaseNetwork4substrates}
\title{
Visualize kinase-substrate networks from PTMSEA results}
\description{
Creates and visualizes kinase-substrate networks based on significant PTMSEA results and kinase-substrate mappings. The function generates network plots for each significant kinase, showing its downstream substrates and their regulation.
}
\usage{
KinaseNetwork4substrates(pair, limma_output, PTMSEA_output, significance_cutoff = 0.05, significance_statistic = "fdr.pvalue", mapping_ID, PTMSEA_OUTDIR, mapping_regulation, output_file_suffix = "", PTMsigDB_collection_file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pair}{A character string specifying the comparison or experimental condition.}
  \item{limma_output}{A data frame containing phosphosite-level statistics, typically from limma analysis.}
  \item{PTMSEA_output}{A data frame with PTMSEA results, including significance statistics for kinase sets.}
  \item{significance_cutoff}{Numeric threshold for significance (default 0.05).}
  \item{significance_statistic}{Column name in \code{PTMSEA_output} to use for significance (e.g., "fdr.pvalue" or "pvalue").}
  \item{mapping_ID}{A named character vector mapping phosphosite flanking regions to UniProt IDs or other identifiers.}
  \item{PTMSEA_OUTDIR}{Directory path for saving output network plots.}
  \item{mapping_regulation}{A named numeric vector mapping phosphosites to regulation values (e.g., logFC or t-statistics).}
  \item{output_file_suffix}{Optional suffix for output file names (default "").}
  \item{PTMsigDB_collection_file}{Path to the PTMsigDB collection file in GMT format.}
}
\details{
This function identifies significant kinase sets from PTMSEA results, extracts their downstream substrates, and visualizes the kinase-substrate network using \code{igraph}. Each network is saved as a PDF in the specified output directory.
}
\value{
No return value. Generates and saves network plots as PDF files.
}
\references{

}
\author{
Yuan Li \email{yuan.li@nbis.se}
}
\note{
Requires the \code{igraph}, \code{dplyr}, and \code{GSEABase} packages.
}

\examples{
## Not run:
## Example usage (requires appropriate input data):
# KinaseNetwork4substrates(
#   pair = "ConditionA_vs_ConditionB",
#   limma_output = limma_df,
#   PTMSEA_output = ptmsea_df,
#   significance_cutoff = 0.05,
#   significance_statistic = "fdr.pvalue",
#   mapping_ID = mapping_vec,
#   PTMSEA_OUTDIR = "results/",
#   mapping_regulation = reg_vec,
#   output_file_suffix = "",
#   PTMsigDB_collection_file = "PTMsigDB.gmt"
# )

## The function is currently defined as
function (pair, limma_output, PTMSEA_output, significance_cutoff = 0.05, 
    significance_statistic = "fdr.pvalue", mapping_ID, PTMSEA_OUTDIR, 
    mapping_regulation, output_file_suffix = "", PTMsigDB_collection_file) 
{
    limma_output$PTM.FlankingRegion2 = paste0(limma_output$PTM.FlankingRegion, 
        "-p")
    limma_output$uniprotID = limma_output$Phosphosite \%>\% sub("_.*", 
        "", .)
    limma_output$PhosLocation = limma_output$Phosphosite \%>\% 
        sub(".*_", "", .)
    limma_output = limma_output \%>\% dplyr::select(uniprotID, 
        PhosLocation, PTM.FlankingRegion2, everything())
    sigIDs = PTMSEA_output[PTMSEA_output[, significance_statistic, 
        drop = T] < significance_cutoff, ] \%>\% pull(id) \%>\% .[!is.na(.)]
    ptmsigdb = GSEABase::getGmt(PTMsigDB_collection_file)
    ptmsigdb0 = lapply(ptmsigdb, function(K) {
        data.frame(Term = K@setName, Phosphosite = K@geneIds)
    }) \%>\% Reduce(rbind, .)
    ptmsigdb2 = as.data.frame(stringr::str_split_fixed(ptmsigdb0$Phosphosite, 
        ";", 2)) \%>\% setNames(c("Phosphosite", "Regulation"))
    ptmsigdb3 = data.frame(Term = ptmsigdb0$Term, Phosphosite = ptmsigdb2$Phosphosite, 
        Regulation = ptmsigdb2$Regulation)
    lapply(sigIDs, function(ID) {
        substrates = ptmsigdb3 \%>\% filter(Term == ID) \%>\% pull(Phosphosite)
        shared = intersect(substrates, limma_output$PTM.FlankingRegion2)
        regulation = ptmsigdb3 \%>\% filter(Term == ID) \%>\% filter(Phosphosite \%in\% 
            shared) \%>\% pull(Regulation)
        edge = data.frame(From = rep(ID, length(shared)), To = shared)
        g = graph_from_data_frame(d = edge, directed = FALSE)
        V(g)$label = c(V(g)[[1]]$name, mapping_ID[V(g)[2:length(V(g))]$name] \%>\% 
            as.character())
        mapping_regulation = mapping_regulation[shared]
        E(g)$color = ifelse(regulation == "u", alpha("orange", 
            0.25), alpha("darkturquoise", 0.25))
        V(g)$color = c(alpha("black", 0.25), ifelse(mapping_regulation < 
            0, alpha("darkturquoise", 0.5), alpha("orange", 0.5)))
        layout = layout_with_fr(g)
        dir.create(paste0(PTMSEA_OUTDIR, "/Network/kinase_substrates/", 
            pair, "/"), showWarnings = FALSE, recursive = TRUE)
        L = length(V(g))
        ID2 = gsub("/", "_", ID)
        if (output_file_suffix != "") {
            suffix = paste0("_", output_file_suffix)
            pdf(paste0(PTMSEA_OUTDIR, "/Network/kinase_substrates/", 
                pair, "/", ID2, suffix, ".pdf"), h = (L/2.5) + 
                3, w = (L/2.5) + 3)
        }
        else {
            suffix = ""
            pdf(paste0(PTMSEA_OUTDIR, "/Network/kinase_substrates/", 
                pair, "/", ID2, suffix, ".pdf"), h = (L/4) + 
                2.5, w = (L/4) + 2.8)
        }
        print(plot(g, vertex.label.dist = 0, layout = layout, 
            vertex.label.cex = 0.8, vertex.label = V(g)$label, 
            vertex.color = V(g)$color, vertex.frame.color = V(g)$color, 
            vertex.size = abs(round(c(6, as.numeric(mapping_regulation * 
                4)))), vertex.label.color = "black", edge.color = E(g)$color))
        dev.off()
        print(plot(g, vertex.label.dist = 0, layout = layout, 
            vertex.label.cex = 0.8, vertex.label = V(g)$label, 
            vertex.color = V(g)$color, vertex.frame.color = V(g)$color, 
            vertex.size = abs(round(c(6, as.numeric(mapping_regulation * 
                4)))), vertex.label.color = "black", edge.color = E(g)$color))
    })
  }
}
\keyword{kinase}
\keyword{network}
\keyword{visualization}
\keyword{substrate}
\keyword{PTMSEA}
\keyword{PTMsigDB}
\keyword{phosphoproteomics}
\keyword{proteomics}
