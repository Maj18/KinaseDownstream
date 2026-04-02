![Version](https://img.shields.io/badge/Version-2.0.0-green)  [![DOI](https://zenodo.org/badge/994653217.svg)](https://doi.org/10.5281/zenodo.19384525)

# KinaseDownstream

This package provides functions to visualize the significantly enriched kinases in the phosphoproteomics data, the networks between each kinase, the downstream targets and their immediate String_db interacting neighbors that are observed from the matching proteomics data, facilitating the understanding of kinase interactions and their biological implications.

![Diagram of the KinaseDownstream package](https://github.com/Maj18/KinaseDownstream/blob/main/man/figures/Diagram.png)

# How to use this package

To install the permuDE package in a R session:

```
install.packages("remotes")
remotes::install_github("Maj18/KinaseDownstream", dependencies = TRUE)
```

Alternatively, you can download the package locally, and in the directory where the package directory is loaded, open a terminal and run

```
R CMD INSTALL KinaseDownstream
```


To check whether the package has been properly installed, start a R session and type:

```
library(KinaseDownstream)
```

## Example run:
```
library(KinaseDownstream)
INDIR = "./"
limma_rslt = readRDS(paste0(INDIR, "/Limma_differentialAnalysis_result.RDS"))
```
The limma input should look like this:

$PNLvsPL
|   |Phosphosite  | PTM.FlankingRegion |logFC   |AveExpr    |    t    |   P.Value  |  adj.P.Val |    B    |
|---|-------------|--------------------|--------|-----------|---------|------------|------------|---------|
|1  |Q9GZM8_T245  |  GFGTSPLTPSARISA   |7.939980| 6.880564  | 5.174157|3.060240e-05|1.160624e-04| 2.282851|
|2  |Q09666_S2397 |  DLDLHLKSPKAKGEV   |7.768347| 5.267173  |13.931695|5.022985e-13|1.024772e-11|16.576525|

$HSvsPNL
|  |Phosphosite   |PTM.FlankingRegion |logFC    |AveExpr  |    t     |  P.Value  |  adj.P.Val |    B     |
|--|--------------|-------------------|---------|---------|----------|-----------|------------|----------|
|1 |Q9UPU9_S420   |TPIKAYSSPSTTPEA    |41.798077| 4.815310| 0.7234658| 0.51001991| 0.72947446 | -4.595267|
|2 |Q4KMP7_S678   |RAAGGAPSPPPPVRR    |37.988089| 2.888568| 3.6986481| 0.02143313| 0.09818054 | -4.564233|

#### Run PTMSEA analysis
```
PTMSEA_OUTDIR = "./"
PTMSEA_rslt = runPTMSEA(limma_rslt, PTMSEA_OUTDIR)
# Get the PTM substrates that have been used for PTMSEA analysis
PTMsubstrates4PTMSEAanalysis = getPTMsubstrates4PTMSEAanalysis(limma_rslt)
```

**Note:** For the moment, we only support the PTMSEA analysis (based on PTMsigDB2.0) with the flanking regions of the phosphosites, which are defined as 7 amino acids upstream and downstream of the phosphosite, in human. The PTMSEA analysis is performed on the phosphoproteomics data that has been processed by limma. 

#### Prepare the PTMSEA output for the kinase-substrate network analysis
```
PTMSEA_FILE_PATH = paste0(PTMSEA_OUTDIR,"/PTMSEA_OUTPUT-combined.gct")
significance_cutoff = 0.01
ptmsea_rslt_all = processPTMSEAresult(PTMSEA_FILE_PATH, output.score.type = "NES", 
    sig.thresh=significance_cutoff, PTMSEA_OUTDIR = PTMSEA_OUTDIR)
ptmsea_rslt = ptmsea_rslt_all$ptmsea_rslt
```

Check the dotplots for the significant kinases
```
SignificantKinaseDotplots = ptmsea_rslt_all$plot
SignificantKinaseDotplots
```
![Example dotplot](https://github.com/Maj18/KinaseDownstream/blob/main/man/figures/KinaseDotplot_example.png)

One can also map the significant kinases to manning kinase tree at http://kinhub.org/kinmap/index.html

The significant kinases can be found in the SigKinaseList_*.txt within PTMSEA_OUTDIR.

![Example manning kinase tree mapping](https://github.com/Maj18/KinaseDownstream/blob/main/man/figures/ExampleManningKinaseTreeMapping.png)

#### Process the limma output to build maps among FlankingRegions, uniprot IDs and Phosphosites, for the kinase-substrate network analysis
```
maps = processLimmaResult(limma_rslt, PTMsubstrates4PTMSEAanalysis)
```

#### Run the kinase-substrate network analysis for each pair of conditions
```
significance_statistic = "fdr.pvalue"
invisible(capture.output(lapply(seq_along(limma_rslt), function(i) {
  pair = names(limma_rslt)[i]
  limma_output = limma_rslt[[i]]
  PTMSEA_output = ptmsea_rslt[[i]]
  mapping_regulation = maps[[i]]$mapLogFC2FlankingRegion2
  mapping_ID = maps[[i]]$mapUniprotPhosLocation2FlankingRegion2
  KinaseNetwork4substrates(pair,
                           PTMsubstrates4PTMSEAanalysis,
                           limma_output,
                           PTMSEA_output,
                           significance_cutoff = significance_cutoff, 
                           significance_statistic = significance_statistic, 
                           mapping_ID = mapping_ID,
                           mapping_regulation = mapping_regulation,
                           PTMSEA_OUTDIR = PTMSEA_OUTDIR,
                           PTMsigDB_collection_file = 
                            paste0(PTMSEA_OUTDIR, "/ptm.sig.db.all.flanking.human.v2.0.0.gmt"))

  # Add flanking region to the kinase-substrate networks:
  mapping_ID = maps[[i]]$mapUniprotPhosLocationFlankingRegion2FlankingRegion2
  KinaseNetwork4substrates(pair, PTMsubstrates4PTMSEAanalysis,
                           limma_output,
                           PTMSEA_output,
                           significance_cutoff = significance_cutoff, 
                           significance_statistic = "fdr.pvalue", 
                           mapping_ID = mapping_ID,
                           mapping_regulation = mapping_regulation,
                           output_file_suffix = "wPhosphosites",
                           PTMSEA_OUTDIR = PTMSEA_OUTDIR,
                           PTMsigDB_collection_file = 
                            paste0(PTMSEA_OUTDIR, "/ptm.sig.db.all.flanking.human.v2.0.0.gmt")) # This is needed when flanking regions are added to the figures, otherwise, please DO NOT set outputFileSuffix.
})))

```
![Example kinase-substrate network](https://github.com/Maj18/KinaseDownstream/blob/main/man/figures/Kinase_substrates_network_example.png)

> In the example kinase-substrate network plot above, the nodes represent the kinases (grey) and their downstream targets, while the edges represent the regulation status of interactions in the PTMsigDB (upregulated: orange or downregulated: turquoise). The color of the substrate nodes indicates the regulation status (upregulated: orange or downregulated: turquoise) based on the phosphoproteomics data.

Save the kinase-substrate network as csv file for collective kinase-substrate network analysis using PhosNetVis.
```
invisible(capture.output(lapply(seq_along(limma_rslt), function(i) {
  pair = names(limma_rslt)[i]
  limma_output = limma_rslt[[i]]
  PTMSEA_output = ptmsea_rslt[[i]]
  prepInput4PhosNetVis(pair, limma_output, PTMsubstrates4PTMSEAanalysis, 
                                PTMSEA_output, significance_cutoff=0.05, 
                                significance_statistic="fdr.pvalue", 
                                PTMsigDB_collection_file = 
                                paste0(OUTDIR_PTMSEA, "/ptm.sig.db.all.flanking.human.v2.0.0.gmt"))
})))

```

##### Prepare proteomics dataset that matches the phosphoproteomics data
```
# Download the psoriasis-associated proteins from the reference
proteomics_dat = openxlsx::read.xlsx(paste0(INDIR, "/12014_2020_9293_MOESM5_ESM.xlsx"), sheet = 1)
# Let's get significantly upregulated proteins
up_proteomics = proteomics_dat %>% filter(`pso/con.Ratio`>1.3) %>% pull(Protein.accession) # This is a list protein uniprot id
# [1] "A0A075B6I9" "A0A075B6K0" "A0A075B6K4" "A0A0C4DH38" "A8K2U0"     "O00148"
# Let's get significantly downregulated proteins
down_proteomics = proteomics_dat %>% filter(`pso/con.Ratio`< (1/1.3)) %>% pull(Protein.accession)

```

#### Run the PPI network analysis for each pair of conditions

The PPI network will connect each kinase to its downstream targets and their immediate String_db interacting neighbors that are observed from the matching proteomics data.

```
significance_statistic = "fdr.pvalue"
significance_cutoff = 0.01
invisible(capture.output(lapply(seq_along(limma_rslt), function(i) {
  pair = names(limma_rslt)[i]
  outdir_ppi = paste0(PTMSEA_OUTDIR, "/Network/PPI/", pair, "/")
  dir.create(outdir_ppi, recursive=T, showWarnings=F)
  limma_output = limma_rslt[[i]]
  PTMSEA_output = ptmsea_rslt[[i]]
  mapping_regulation = maps[[i]]$mapLogFC2FlankingRegion2
  mapping_ID = maps[[i]]$mapUniprotPhosLocation2FlankingRegion2
  ppiNetwork4substrates(limma_output, PTMsubstrates4PTMSEAanalysis,
                        PTMSEA_output,
                        PTMSEA_OUTDIR = PTMSEA_OUTDIR,
                        significance_cutoff = significance_cutoff, 
                        significance_statistic = significance_statistic, 
                        mapping_ID = mapping_ID,
                        mapping_regulation = mapping_regulation,
                        proteomics = up_proteomics,
                        outdir_ppi = outdir_ppi,
                        PTMsigDB_collection_file = 
                            paste0(PTMSEA_OUTDIR, "/ptm.sig.db.all.flanking.human.v2.0.0.gmt"),
                        uniprot_provided = T)
})))  

```

![Example PPI network](https://github.com/Maj18/KinaseDownstream/blob/main/man/figures/PPI_example.png)

> In the example PPI network plot above, the nodes represent the kinases (purple), their downstream targets (orange or turquoise) and the immediate String_db interacting neighbors (gray) of the targeted phosphoproteins. The immediate neighbors are required to be significantly regulated in the proteomics data. The edges represent the interactions between them. The color of the nodes for the substrates indicates the regulation status (upregulated: orange or downregulated: turquoise) based on the phosphoproteomics data. The color of the edges of the kinase-substrate interactions indicates the regulation status of the targets (upregulated: purple or downregulated: turquoise) in PTMsigDB, while the color of the edges of the substrate-neighbor interactions are in orange.
