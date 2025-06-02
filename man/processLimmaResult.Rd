\name{processLimmaResult}
\alias{processLimmaResult}
\title{
Process limma results for PTMSEA analysis
}
\description{
Processes a list of limma result data frames, filters for phosphosites present in a given set, and generates mapping vectors for downstream PTMSEA analysis.
}
\usage{
processLimmaResult(limma_rslt, PTM.FlankingRegion4PTMSEAanalysis)
}
\arguments{
  \item{limma_rslt}{A list of data frames, each containing phosphosite-level statistics from limma analysis.}
  \item{PTM.FlankingRegion4PTMSEAanalysis}{A character vector of phosphosite flanking regions to retain for PTMSEA analysis.}
}
\details{
This function filters each limma result data frame to retain only phosphosites present in the specified set. It then creates mapping vectors for UniProt ID and phosphosite location, as well as log fold change, for use in downstream PTMSEA and network analysis.
}
\value{
A list of lists, each containing:
  \item{mapUniprotPhosLocation2FlankingRegion2}{Named character vector mapping flanking region to UniProt ID and phosphosite location.}
  \item{mapUniprotPhosLocationFlankingRegion2FlankingRegion2}{Named character vector mapping flanking region to UniProt ID, phosphosite location, and original flanking region.}
  \item{mapLogFC2FlankingRegion2}{Named numeric vector mapping flanking region to log fold change.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Yuan Li \email{yuan.li@nbis.se}
}
\note{
Requires the \code{dplyr} and \code{stringr} packages.
}
\examples{
## Example usage (requires appropriate input data):
# result <- processLimmaResult(limma_rslt, PTM.FlankingRegion4PTMSEAanalysis)

## The function is currently defined as
function (limma_rslt, PTM.FlankingRegion4PTMSEAanalysis) 
{
    maps = lapply(limma_rslt, function(rslt) {
        rslt = filter(rslt, PTM.FlankingRegion \%in\% PTM.FlankingRegion4PTMSEAanalysis)
        rslt$uniprotID = rslt$Phosphosite \%>\% sub("_.*", "", 
            .)
        rslt$PTM.FlankingRegion2 = paste0(rslt$PTM.FlankingRegion, 
            "-p")
        rslt$PhosLocation = rslt$Phosphosite \%>\% sub(".*_", "", 
            .)
        mapUniprotPhosLocation2FlankingRegion2 = setNames(paste0(rslt$uniprotID, 
            "\n", rslt$PhosLocation), rslt$PTM.FlankingRegion2)
        mapUniprotPhosLocationFlankingRegion2FlankingRegion2 = setNames(paste0(rslt$uniprotID, 
            "\n", rslt$PhosLocation, "\n", rslt$PTM.FlankingRegion), 
            rslt$PTM.FlankingRegion2)
        mapLogFC2FlankingRegion2 = setNames(rslt$logFC, rslt$PTM.FlankingRegion2)
        list(mapUniprotPhosLocation2FlankingRegion2 = mapUniprotPhosLocation2FlankingRegion2, 
            mapUniprotPhosLocationFlankingRegion2FlankingRegion2 = mapUniprotPhosLocationFlankingRegion2FlankingRegion2, 
            mapLogFC2FlankingRegion2 = mapLogFC2FlankingRegion2)
    })
    return(maps)
  }
}

\keyword{PTMSEA}
\keyword{phosphoproteomics}
\keyword{limma}
\keyword{mapping}
\keyword{data processing}
