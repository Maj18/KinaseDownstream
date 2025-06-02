\name{ppiNetwork4substrates}
\alias{ppiNetwork4substrates}

\title{
Build and visualize protein-protein interaction networks for kinases, their substrates and the immediate String_db neighbors observed in matching proteomics data. 
}
\description{
Constructs and visualizes protein-protein interaction (PPI) networks for kinase substrates identified as significant in PTMSEA analysis. The function integrates STRINGdb interactions and proteomics data, generating network plots for each kinase.
}
\usage{
ppiNetwork4substrates(limma_output, PTMSEA_output, significance_cutoff = 0.05, significance_statistic = "fdr.pvalue", mapping_ID, mapping_regulation, proteomics, outdir_ppi, PTMsigDB_collection_file, uniprot_provided = T)
}
\arguments{
  \item{limma_output}{A data frame containing phosphosite-level statistics, typically from limma analysis.}
  \item{PTMSEA_output}{A data frame with PTMSEA results, including significance statistics for kinase sets.}
  \item{significance_cutoff}{Numeric threshold for significance (default 0.05).}
  \item{significance_statistic}{Column name in \code{PTMSEA_output} to use for significance (e.g., "fdr.pvalue" or "pvalue").}
  \item{mapping_ID}{A named character vector mapping phosphosite flanking regions to UniProt IDs or other identifiers.}
  \item{mapping_regulation}{A named numeric vector mapping phosphosites to regulation values (e.g., logFC or t-statistics).}
  \item{proteomics}{A character vector of protein identifiers observed in the proteomics dataset.}
  \item{outdir_ppi}{Directory path for saving output PPI network plots.}
  \item{PTMsigDB_collection_file}{Path to the PTMsigDB collection file in GMT format.}
  \item{uniprot_provided}{Logical; if \code{TRUE}, UniProt IDs are provided for mapping (default \code{TRUE}).}
}
\details{
This function identifies significant kinase sets from PTMSEA results, maps their substrates to STRINGdb protein-protein interaction networks, and visualizes the resulting networks. Each network is saved as a PDF in the specified output directory. Node and edge attributes reflect substrate regulation and interaction types.
}
\value{
No return value. Generates and saves PPI network plots as PDF files.
}
\references{
}
\author{
Yuan Li \email{yuan.li@nbis.se}
}
\note{
Requires the \code{STRINGdb}, \code{igraph}, \code{dplyr}, \code{biomaRt}, and \code{future.apply} packages.
}

\examples{
## Not run:
## Example usage (requires appropriate input data):
# ppiNetwork4substrates(
#   limma_output = limma_df,
#   PTMSEA_output = ptmsea_df,
#   significance_cutoff = 0.05,
#   significance_statistic = "fdr.pvalue",
#   mapping_ID = mapping_vec,
#   mapping_regulation = reg_vec,
#   proteomics = proteomics_vec,
#   outdir_ppi = "results/PPI/",
#   PTMsigDB_collection_file = "PTMsigDB.gmt",
#   uniprot_provided = TRUE
# )

## The function is currently defined as
function (limma_output, PTMSEA_output, significance_cutoff = 0.05, 
    significance_statistic = "fdr.pvalue", mapping_ID, mapping_regulation, 
    proteomics, outdir_ppi, PTMsigDB_collection_file, uniprot_provided = T) 
{
    dir.create(outdir_ppi, recursive = T, showWarnings = F)
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
    string_db = STRINGdb::STRINGdb$new(version = "11.5", species = 9606, 
        score_threshold = 700)
    mapped_proteomics = string_db$map(data.frame(protein = proteomics), 
        "protein", removeUnmappedRows = TRUE)
    substrates_all = ptmsigdb3 \%>\% filter(Term \%in\% sigIDs) \%>\% 
        pull(Phosphosite)
    shared_all = intersect(substrates_all, limma_output$PTM.FlankingRegion2)
    shared_uniprotID_all = mapping_ID[shared_all] \%>\% as.character() \%>\% 
        sub("\n.*", "", .)
    mapping_substrates = string_db$map(data.frame(protein = shared_uniprotID_all), 
        "protein", removeUnmappedRows = TRUE)
    print("When one uniprot ID mapped to multiple STRING IDs, only keep the 1st one!")
    mapping_substrates = mapping_substrates[!duplicated(mapping_substrates$protein), 
        ]
    neighbors_all = string_db$get_neighbors(mapping_substrates$STRING_id)
    proteomics_neighbors_all = intersect(neighbors_all, mapped_proteomics$STRING_id)
    all_proteins_all = c(mapping_substrates$STRING_id, proteomics_neighbors_all)
    mapped_back_all_proteins_all = string_db$add_proteins_description(data.frame(STRING_id = all_proteins_all)) \%>\% 
        mutate(hgnc_symbol = preferred_name)
    if (!file.exists(paste0(PTMSEA_OUTDIR, "/Network/PPI/ensembl.rds"))) {
        ensembl = biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", 
            mirror = "useast")
        saveRDS(ensembl, paste0(PTMSEA_OUTDIR, "/Network/PPI/ensembl.rds"))
    }
    else {
        ensembl = readRDS(paste0(PTMSEA_OUTDIR, "/Network/PPI/ensembl.rds"))
    }
    temp_all_proteins_all = biomaRt::getBM(attributes = c("hgnc_symbol", 
        "uniprotswissprot", "uniprot_gn_symbol"), filters = "hgnc_symbol", 
        values = mapped_back_all_proteins_all$preferred_name, 
        mart = ensembl) \%>\% filter(uniprotswissprot != "") \%>\% 
        as.data.frame() \%>\% dplyr::select(hgnc_symbol, uniprotswissprot) \%>\% 
        unique()
    duplicated_hgnc_symbol = temp_all_proteins_all$hgnc_symbol[duplicated(temp_all_proteins_all$hgnc_symbol)]
    if (length(duplicated_hgnc_symbol) > 0) {
        uniprot4duplicated_hgnc_symbol = temp_all_proteins_all \%>\% 
            filter(temp_all_proteins_all$hgnc_symbol \%in\% duplicated_hgnc_symbol) \%>\% 
            pull(uniprotswissprot)
        all_IDs_from_dat = sapply(mapping_ID, function(I) {
            strsplit(I, "\n")
        }) \%>\% unlist() \%>\% as.character()
        uniprot4duplicated_hgnc_symbol_in_our_dat = uniprot4duplicated_hgnc_symbol \%>\% 
            .[. \%in\% all_IDs_from_dat]
        if (length(uniprot4duplicated_hgnc_symbol_in_our_dat) > 
            0) {
            dup_hgnc_symbol_inOurData = filter(temp_all_proteins_all, 
                uniprotswissprot \%in\% uniprot4duplicated_hgnc_symbol_in_our_dat) \%>\% 
                pull(hgnc_symbol)
            rest = filter(temp_all_proteins_all, !hgnc_symbol \%in\% 
                dup_hgnc_symbol_inOurData)
            hit = filter(temp_all_proteins_all, hgnc_symbol \%in\% 
                dup_hgnc_symbol_inOurData)
            temp_all_proteins_all = rbind(rest[!duplicated(rest$hgnc_symbol), 
                ], filter(hit, uniprotswissprot \%in\% uniprot4duplicated_hgnc_symbol_in_our_dat))
        }
        else {
            temp_all_proteins_all = temp_all_proteins_all[!duplicated(temp_all_proteins_all$hgnc_symbol), 
                ]
        }
    }
    temp_all_proteins_all = full_join(temp_all_proteins_all, 
        mapped_back_all_proteins_all, by = "hgnc_symbol")
    mapping_substrates = setNames(mapping_substrates$protein, 
        mapping_substrates$STRING_id)
    which_substrates = temp_all_proteins_all$STRING_id \%in\% names(mapping_substrates)
    temp_all_proteins_all$uniprotswissprot[which_substrates] = mapping_substrates[temp_all_proteins_all$STRING_id[which_substrates]]
    temp_all_proteins_all = na.omit(temp_all_proteins_all)
    interactions_all = string_db$get_interactions(temp_all_proteins_all$STRING_id)
    future::plan("multisession", workers = 4)
    future.apply::future_lapply(sigIDs, function(ID) {
        substrates = ptmsigdb3 \%>\% filter(Term == ID) \%>\% pull(Phosphosite)
        shared = intersect(substrates, limma_output$PTM.FlankingRegion2)
        shared_uniprotID = mapping_ID[shared] \%>\% as.character() \%>\% 
            sub("\n.*", "", .)
        mapped_interest = mapping_substrates[mapping_substrates \%in\% 
            shared_uniprotID]
        neighbors = string_db$get_neighbors(names(mapped_interest))
        proteomics_neighbors = intersect(neighbors, mapped_proteomics$STRING_id)
        all_proteins = c(names(mapped_interest), proteomics_neighbors)
        mapped_back_all_proteins = mapped_back_all_proteins_all \%>\% 
            filter(STRING_id \%in\% all_proteins)
        if (uniprot_provided) {
            temp_all_proteins = filter(temp_all_proteins_all, 
                hgnc_symbol \%in\% mapped_back_all_proteins$preferred_name)
            mapping_back_all_proteins = setNames(paste0(temp_all_proteins$hgnc_symbol, 
                "\n(", temp_all_proteins$uniprotswissprot, ")"), 
                temp_all_proteins$STRING_id)
            mapping_back_all_proteins_forKinaseSubstrates = setNames(paste0(temp_all_proteins$hgnc_symbol, 
                "\n(", temp_all_proteins$uniprotswissprot, ")"), 
                temp_all_proteins$uniprotswissprot)
        }
        else {
            mapping_back_all_proteins = setNames(mapped_back_all_proteins$preferred_name, 
                mapped_back_all_proteins$STRING_id)
            mapping_back_all_proteins_forKinaseSubstrates = setNames(paste0(temp_all_proteins$hgnc_symbol), 
                temp_all_proteins$uniprotswissprot)
        }
        kinase_substrate_ineractions = data.frame(from = rep(ID, 
            length(shared_uniprotID)), to = shared_uniprotID, 
            combined_score = "NA")
        kinase_substrate_ineractions = kinase_substrate_ineractions[kinase_substrate_ineractions$to \%in\% 
            names(mapping_back_all_proteins_forKinaseSubstrates), 
            ]
        if (uniprot_provided) {
            kinase_substrate_ineractions$to = mapping_back_all_proteins_forKinaseSubstrates[kinase_substrate_ineractions$to]
        }
        kinase_substrate_ineractions = unique(kinase_substrate_ineractions)
        interactions = interactions_all[interactions_all$from \%in\% 
            names(mapping_back_all_proteins) & interactions_all$to \%in\% 
            names(mapping_back_all_proteins), ] \%>\% unique()
        interactions$from = mapping_back_all_proteins[as.character(interactions$from)]
        interactions$to = mapping_back_all_proteins[as.character(interactions$to)]
        interactions = rbind(kinase_substrate_ineractions, interactions)
        edge.colors = rep("snow3", nrow(interactions))
        edge.colors[interactions$from \%in\% as.character(mapping_back_all_proteins[names(mapped_interest)]) | 
            interactions$to \%in\% as.character(mapping_back_all_proteins[names(mapped_interest)])] = alpha("orange", 
            0.5)
        t1 = ptmsigdb3 \%>\% filter(Term == ID)
        t2 = mapping_ID[t1$Phosphosite] \%>\% na.omit()
        interactions2 = filter(t1, t1$Phosphosite \%in\% names(t2)) \%>\% 
            left_join(data.frame(Phosphosite = names(t2), Names = t2, 
                stringsAsFactors = FALSE)) \%>\% mutate(uniprotID = gsub("\n.*", 
            "", Names)) \%>\% right_join(filter(interactions, from \%in\% 
            ID) \%>\% mutate(uniprotID = gsub(".*\n", "", to)) \%>\% 
            mutate(uniprotID = gsub("[(]|[)]", "", uniprotID))) \%>\% 
            left_join(data.frame(Phosphosite = names(mapping_regulation), 
                effect = mapping_regulation, stringsAsFactors = FALSE), 
                )
        mapping_TO_regualtion = setNames(interactions2$Regulation, 
            interactions2$to)
        regulations_ID = mapping_TO_regualtion[filter(interactions, 
            from \%in\% ID) \%>\% pull(to)]
        edge.colors[interactions$from \%in\% ID] = ifelse(regulations_ID == 
            "u", alpha("purple", 0.5), alpha("darkturquoise", 
            0.5))
        edge.widths = rep(0.1, nrow(interactions))
        network_vertexIDs_4substrates = as.character(mapping_back_all_proteins[names(mapped_interest)]) \%>\% 
            na.omit()
        edge.widths[interactions$from \%in\% network_vertexIDs_4substrates \%>\% 
            na.omit() | interactions$to \%in\% network_vertexIDs_4substrates \%>\% 
            na.omit()] = 0.75
        ppi_network = graph_from_data_frame(interactions[, 1:2], 
            directed = FALSE)
        vertex.label.colors = rep(alpha("black", 0.5), length(V(ppi_network)$name))
        vertex.label.colors[V(ppi_network)$name \%in\% ID] = alpha("black", 
            1)
        vertex.label.colors[V(ppi_network)$name \%in\% network_vertexIDs_4substrates] = alpha("black", 
            0.8)
        mapping_TO_effect = setNames(interactions2$effect, interactions2$to)
        effect = mapping_TO_effect[V(ppi_network)$name \%>\% .[. \%in\% 
            network_vertexIDs_4substrates]]
        vertex.colors = rep("snow3", length(V(ppi_network)$name))
        vertex.colors[V(ppi_network)$name \%in\% ID] = alpha("purple", 
            0.5)
        vertex.colors[V(ppi_network)$name \%in\% network_vertexIDs_4substrates] = ifelse(effect < 
            0, alpha("darkturquoise", 0.5), alpha("orange", 0.5))
        vertex.sizes = rep(1, length(V(ppi_network)$name))
        vertex.sizes[V(ppi_network)$name \%in\% network_vertexIDs_4substrates] = abs(round(as.numeric(effect)) * 
            2)
        vertex.sizes[V(ppi_network)$name \%in\% ID] = 10
        coords = layout_in_circle(ppi_network)
        L = length(V(ppi_network))/60
        ID2 = sub("/", "_", ID)
        pdf(paste0(outdir_ppi, "/", ID2, "_substrate_PPI.pdf"), 
            h = 7 * L + 2, w = 10 * L + 2)
        print(plot(ppi_network, vertex.shape = "circle", vertex.label.cex = 0.7, 
            vertex.label.color = vertex.label.colors, edge.color = edge.colors, 
            vertex.label.dist = 0, edge.width = edge.widths, 
            layout = coords, main = ID, asp = 0.7, vertex.color = vertex.colors, 
            vertex.frame.color = "white", vertex.size = vertex.sizes, 
            edge.curved = 0))
        dev.off()
        print(plot(ppi_network, vertex.shape = "circle", vertex.label.cex = 0.7, 
            vertex.label.color = vertex.label.colors, edge.color = edge.colors, 
            vertex.label.dist = 0, edge.width = edge.widths, 
            layout = coords, main = ID, asp = 0.7, vertex.color = vertex.colors, 
            vertex.frame.color = "white", vertex.size = vertex.sizes, 
            edge.curved = 0))
    })
  }
}
\keyword{network}
\keyword{visualization}
\keyword{protein-protein interaction}
\keyword{kinase}
\keyword{substrate}
\keyword{PTMSEA}
\keyword{PTMsigDB}
\keyword{phosphoproteomics}
\keyword{proteomics}
\keyword{PPI network}
