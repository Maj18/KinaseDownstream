#' Run PTM-SEA on the limma results
#' @import dplyr
#' 
runPTMSEA = function(limma_rslt, PTMSEA_OUTDIR) {
  limma_rslt_4gct = lapply(seq_along(limma_rslt), function(i){
    rslt = limma_rslt[[i]] %>% arrange(-logFC)
    rownames(rslt) = NULL
    temp = rslt[, c("PTM.FlankingRegion", "logFC")]
    colnames(temp) = c("PTM.FlankingRegion", names(limma_rslt)[i])
    temp$PTM.FlankingRegion = paste0(temp$PTM.FlankingRegion, "-p")
    temp
  }) %>% Reduce(full_join, .) %>%
    .[!duplicated(.$PTM.FlankingRegion), ] 
  rownames(limma_rslt_4gct) = NULL
  limma_rslt_4gct = 
    tibble::column_to_rownames(limma_rslt_4gct, var="PTM.FlankingRegion")
  cmapR::write_gct(as.matrix(limma_rslt_4gct)%>%
  cmapR::GCT(.),
  #cmapR::new("GCT",mat=.), 
    paste0(PTMSEA_OUTDIR ,"/PTM-SEA"), precision=2)
  # run PTM-SEA
  input_gct_file = list.files(path=PTMSEA_OUTDIR, pattern = "PTM-SEA", full.names = TRUE)

  # # Download gene set database 
  download.file(url = "https://raw.githubusercontent.com/nicolerg/ssGSEA2/refs/heads/master/db/ptmsigdb/ptm.sig.db.all.flanking.human.v2.0.0.gmt",
                destfile = paste0(PTMSEA_OUTDIR, "/ptm.sig.db.all.flanking.human.v2.0.0.gmt"))
  set.seed(123)
  invisible(capture.output(ssGSEA2::run_ssGSEA2(input_gct_file,
                    output.prefix = "PTMSEA_OUTPUT",
                    gene.set.databases = paste0(PTMSEA_OUTDIR, "/ptm.sig.db.all.flanking.human.v2.0.0.gmt"),
                    output.directory = PTMSEA_OUTDIR,
                    sample.norm.type = "none", #uses actual expression values
                    weight = 1, 
                    correl.type = "rank", #genes are weighted by actual values 
                    statistic = "area.under.RES",
                    spare.cores = 4,
                    output.score.type = "NES", 
                    nperm = 1000, 
                    min.overlap = 5, 
                    extended.output = TRUE, 
                    global.fdr = FALSE,
                    export.signat.gct = T,
                    param.file=T,
                    log.file = paste0(PTMSEA_OUTDIR, "/run.log")))) 
                    
  return(limma_rslt_4gct)
}




#' Prepare the PTM-SEA results
#' @import tidyr
#' @import ggplot2
#' @import dplyr
#' @param PTMSEA_FILE_PATH the path to the PTM-SEA output file: PTMSEA_OUTPUT-combined.gct.
#' @param output.score.type the type of score to be used in the output, default is "NES".
#' @param sig.thresh the fdr.pvalue threshold for filtering the PTMSEA results, default is 0.05.
#' 
processPTMSEAresult = function(PTMSEA_FILE_PATH, output.score.type = "NES", sig.thresh = 0.05, PTMSEA_OUTDIR) {
  gct_data=read.delim(PTMSEA_FILE_PATH , skip = 2, 
                      header = TRUE, sep = "\t", check.names = FALSE) 
  pairs = grep("fdr.pvalue", colnames(gct_data), value=T) %>% gsub("fdr.pvalue.", "", .)
  colnames(gct_data)[colnames(gct_data) %in% pairs] =
    paste0(output.score.type, "_", pairs)
  # Save PTMSEA result to an excel file
  xlsx::write.xlsx(gct_data,
            gsub("gct", "xlsx", PTMSEA_FILE_PATH),
            append=F, sheetName="PTMSEAonDiff",row.names = F)
  ptmsea_rslt = lapply(pairs, function(Pair) {
    temp = gct_data[, c("id", "Signature.set.description", "Signature.set.size", 
                        grep(Pair, colnames(gct_data), value=T))]
    colnames(temp) =
      gsub(paste0("[.]",Pair),"", colnames(temp)) %>% gsub(paste0("_",Pair),"",.)
    temp
  }) %>% setNames(pairs)

  # Make a dot plot
  # library(tidyr)
  gct4plot_NES =  gct_data[, c("id", grep("NES_", colnames(gct_data), value=T))]
  colnames(gct4plot_NES) = gsub("NES_", "", colnames(gct4plot_NES))
  gct4plot_fdr.pvalue =  gct_data[, c("id", grep("fdr.pvalue.", colnames(gct_data), value=T))]
  colnames(gct4plot_fdr.pvalue) = gsub("fdr.pvalue.", "", colnames(gct4plot_fdr.pvalue))
  gct4plot = full_join(tidyr::gather(gct4plot_NES,key="Pair",value="NES",2:ncol(gct4plot_NES)), 
    tidyr::gather(gct4plot_fdr.pvalue,key="Pair",value="fdr.pvalue",2:ncol(gct4plot_fdr.pvalue))) %>% 
    filter(fdr.pvalue<sig.thresh) %>% arrange(-NES) 
  gct4plot$Regulation = ifelse(gct4plot$NES>0, "Up", "Down")
  gct4plot$Pair = factor(gct4plot$Pair)
  gct4plot$id = 
    factor(gct4plot$id, levels=unique(rev(gct4plot$id)))
  # splitted = split(gct4plot, gct4plot$Pair) %>% lapply(., function(i) {
  #   i %>% pull(id)
  # })
  # shared = Reduce(intersect, splitted)
  # shared2 = Reduce(intersect, splitted[c(1,2)])
  # unique1 = setdiff(splitted[[1]], c(shared,shared2))
  # unique2 = setdiff(splitted[[2]], c(shared,shared2))
  # unique3 = setdiff(splitted[[3]], c(shared,shared2))
  # gct4plot$id = factor(gct4plot$id, levels=unique(rev(c(shared,shared2,unique1,unique2,unique3))))
  # gct4plot$Category = ifelse(gct4plot$id%in%shared,"Shared",
  #                          ifelse(gct4plot$id%in%unique1,"Uniq.PNLvsPL", 
  #                                 ifelse(gct4plot$id%in%unique2, "Uniq.HSvsPL", "Uniq.HSvsPNL")))
  # library(ggplot2)
  plot = ggplot2::ggplot(gct4plot, ggplot2::aes(x=NES, y=id)) +
    ggplot2::geom_point(shape=21, ggplot2::aes(size=fdr.pvalue, fill=Pair)) +
    ggplot2::facet_wrap(~Pair, ncol=3) + 
    ggplot2::ggtitle("PTM-SEA_PTMsigDB2") +
    ggplot2::xlab("NES") +
    ggplot2::ylab("") + ggplot2::labs(size="fdr.pvalue") +
    # scale_fill_manual(values=c("purple","orange","grey")) +
    ggplot2::theme_bw(base_size = 15) + 
    ggplot2::scale_size(trans="reverse") +
    ggplot2::guides(size=ggplot2::guide_legend(override.aes=list(shape=21))) +
    ggplot2::theme(strip.text.x = ggplot2::element_text(size = 14)) +
    ggplot2::labs(caption=paste0("fdr.pvalueCutoff = ", sig.thresh, ", pAdjustMethod = fdr"))

  # Prepare signficant kinase names for mapping onto the Manning kinase tree (http://kinhub.org/kinmap/index.html)
  fdrs = grep("fdr.pvalue", colnames(gct_data), value=T)
  for (fdr in fdrs) {
	  N = gsub("fdr.pvalue.", "", fdr)
	  kinaseList = gct_data %>% filter(gct_data[,fdr,drop=T]<sig.thresh) %>% pull(id) %>% 
		  grep("KINASE", ., value=T) %>% gsub("KINASE-PSP_", "", .) %>% 
		  gsub("KINASE-iKiP_", "", .) %>% stringr::str_split(., "/|-|[.]") %>% 
		  unlist() %>% paste0(., collapse=",")
    dir.create(PTMSEA_OUTDIR, showWarnings = FALSE, recursive = TRUE)
    dir.create(paste0(PTMSEA_OUTDIR, "/KinaseGroup/", N), 
               showWarnings = FALSE, recursive = TRUE)
	  writeLines(kinaseList, 
		  paste0(PTMSEA_OUTDIR, "/KinaseGroup/", N, "/SigKinaseList_", N, ".txt"))
  }

  return(list(ptmsea_rslt=ptmsea_rslt, gct4plot=gct4plot, plot=plot))
}


#' @import dplyr
#' @description Prepare the input for PhosNetVis
prepInput4PhosNetVis = function(pair, limma_output, PTM.FlankingRegion4PTMSEAanalysis, 
                                PTMSEA_output, significance_cutoff=0.05, 
                                significance_statistic="fdr.pvalue", 
                                PTMsigDB_collection_file) {
  limma_output$PTM.FlankingRegion2 = paste0(limma_output$PTM.FlankingRegion, "-p")
  limma_output = limma_output %>%
      dplyr::filter(PTM.FlankingRegion2%in%PTM.FlankingRegion4PTMSEAanalysis)
  limma_output$uniprotID = limma_output$Phosphosite %>% sub("_.*", "", .)
  limma_output$PhosLocation = limma_output$Phosphosite %>% sub(".*_", "", .)
  limma_output = limma_output%>% 
    dplyr::select(uniprotID, PhosLocation, PTM.FlankingRegion2, everything()) 
  # get significant PTMsigDB terms
  sigIDs = PTMSEA_output[PTMSEA_output[,significance_statistic,drop=T]<significance_cutoff, ] %>% 
    pull(id) %>% .[!is.na(.)] 
  # import PTMsigDB.v2 collection
  ptmsigdb = GSEABase::getGmt(PTMsigDB_collection_file)
  ptmsigdb0 = lapply(ptmsigdb, function(K) {
    data.frame(Term=K@setName, Phosphosite=K@geneIds)
  }) %>% Reduce(rbind, .)
  ptmsigdb2 = as.data.frame(stringr::str_split_fixed(ptmsigdb0$Phosphosite,";",2)) %>%
    setNames(c("Phosphosite", "Regulation"))
  ptmsigdb3 = data.frame(Term=ptmsigdb0$Term,
                               PTM.FlankingRegion2=ptmsigdb2$Phosphosite,
                               Regulation=ptmsigdb2$Regulation) 
  dir.create(paste0(OUTDIR_PTMSEA, "/CollectiveKinaseTargetNetwork/"), showWarnings = F, recursive = T)
  dplyr::inner_join(limma_output, ptmsigdb3, by="PTM.FlankingRegion2", relationship = "many-to-many") %>%
    dplyr::select(Term, uniprotID, PhosLocation, logFC, P.Value) %>%
    setNames(c("KinaseID", "TargetID", "PhosphoSiteID", "log2FC", "pValue")) %>%
    write.csv(paste0(OUTDIR_PTMSEA, "/CollectiveKinaseTargetNetwork/", 
    pair, "_input4PhosNetVis.csv"), row.names = F)
}





#' @import dplyr
#' @param PTM.FlankingRegion4PTMSEAanalysis a character vector of PTM flanking regions that have been used in PTMSEA.
#' @param limma_rslt a list of data frames containing the limma results.
#' 
processLimmaResult = function(limma_rslt, PTM.FlankingRegion4PTMSEAanalysis) {
  maps = lapply(limma_rslt, function(rslt) {
    rslt$uniprotID = rslt$Phosphosite %>% sub("_.*", "", .)
    rslt$PTM.FlankingRegion2 = paste0(rslt$PTM.FlankingRegion, "-p")
    rslt$PhosLocation = rslt$Phosphosite %>% sub(".*_", "", .)
    rslt = filter(rslt, PTM.FlankingRegion2%in%PTM.FlankingRegion4PTMSEAanalysis)
    mapUniprotPhosLocation2FlankingRegion2 = 
      setNames(paste0(rslt$uniprotID,"\n",rslt$PhosLocation), rslt$PTM.FlankingRegion2)
    mapUniprotPhosLocationFlankingRegion2FlankingRegion2 = 
      setNames(paste0(rslt$uniprotID,"\n",rslt$PhosLocation,
                    "\n",rslt$PTM.FlankingRegion), rslt$PTM.FlankingRegion2)
    mapLogFC2FlankingRegion2 = setNames(rslt$logFC, rslt$PTM.FlankingRegion2) 
    list(mapUniprotPhosLocation2FlankingRegion2 = 
            mapUniprotPhosLocation2FlankingRegion2,
         mapUniprotPhosLocationFlankingRegion2FlankingRegion2 = 
            mapUniprotPhosLocationFlankingRegion2FlankingRegion2,
         mapLogFC2FlankingRegion2 = mapLogFC2FlankingRegion2)
  })

  return(maps) 
}








#' @param significance_statistic this can be "pvalue" or "fdr.pvalue".
#' @param mapping_ID a named character vector mapping phosphosite FlankingRegion to uniprot ID/PhosLocation/FlankingRegion
#' @param PTMSEA_OUTDIR output directory for PTMSEA results
#' @param mapping_regulation a named numeric vector mapping phosphosite to regulation (usually logFC or t-statistics)
#' @param PTMsigDB_collection_file the path to the PTMsigDB collection file in GMT format
#' @param output_file_suffix a suffix for the output file names, default is empty string. If you want to add FlankingRegion to the plot, set output_file_suffix = "wPhosphosites", otherwise set it to "".
#' @import igraph
#' @import dplyr
KinaseNetwork4substrates = function(pair, PTM.FlankingRegion4PTMSEAanalysis, limma_output, PTMSEA_output, significance_cutoff=0.05, 
                                 significance_statistic="fdr.pvalue", mapping_ID, PTMSEA_OUTDIR,
                                 mapping_regulation, output_file_suffix="", PTMsigDB_collection_file) {
  limma_output$PTM.FlankingRegion2 = paste0(limma_output$PTM.FlankingRegion, "-p")
  limma_output = limma_output %>%
      filter(PTM.FlankingRegion2%in%PTM.FlankingRegion4PTMSEAanalysis)
  limma_output$uniprotID = limma_output$Phosphosite %>% sub("_.*", "", .)
  limma_output$PhosLocation = limma_output$Phosphosite %>% sub(".*_", "", .)
  limma_output = limma_output%>% dplyr::select(uniprotID, PhosLocation, PTM.FlankingRegion2, everything()) 
  # get significant PTMsigDB terms
  sigIDs = PTMSEA_output[PTMSEA_output[,significance_statistic,drop=T]<significance_cutoff, ] %>% 
    pull(id) %>% .[!is.na(.)] 

  # import PTMsigDB.v2 collection
  ptmsigdb = GSEABase::getGmt(PTMsigDB_collection_file)
  ptmsigdb0 = lapply(ptmsigdb, function(K) {
    data.frame(Term=K@setName, Phosphosite=K@geneIds)
  }) %>% Reduce(rbind, .)
  ptmsigdb2 = as.data.frame(stringr::str_split_fixed(ptmsigdb0$Phosphosite,";",2)) %>%
    setNames(c("Phosphosite", "Regulation"))
  ptmsigdb3 = data.frame(Term=ptmsigdb0$Term,
                               Phosphosite=ptmsigdb2$Phosphosite,
                               Regulation=ptmsigdb2$Regulation)  
  lapply(sigIDs, function(ID) {
    # print(ID)
    substrates = ptmsigdb3 %>% filter(Term==ID) %>% pull(Phosphosite)
    shared = intersect(substrates, limma_output$PTM.FlankingRegion2) # flanking regions
    regulation = ptmsigdb3 %>% filter(Term==ID) %>% 
      filter(Phosphosite%in%shared) %>% pull(Regulation)
    edge = data.frame(From=rep(ID, length(shared)), To=shared)    
    # library(igraph)
    g = igraph::graph_from_data_frame(d = edge, directed = FALSE)
    igraph::V(g)$label = c(igraph::V(g)[[1]]$name,   #sub(".*_", "", V(g)[[1]]$name), 
                   mapping_ID[igraph::V(g)[2:length(igraph::V(g))]$name] %>% as.character())   
    mapping_regulation = mapping_regulation[shared] 
    igraph::E(g)$color = ifelse(regulation=="u", scales::alpha("orange",0.25), scales::alpha("darkturquoise",0.25))
    igraph::V(g)$color = c(scales::alpha("black", 0.25), ifelse(mapping_regulation<0, 
                                                scales::alpha("darkturquoise",0.5), scales::alpha("orange",0.5)))
    layout = igraph::layout_with_fr(g)
    # layout = layout_in_circle(g)
    dir.create(paste0(PTMSEA_OUTDIR, "/Network/kinase_substrates/", pair,"/"), 
               showWarnings = FALSE, recursive = TRUE)
    L = length(igraph::V(g))
    ID2 = gsub("/", "_", ID)
    if (output_file_suffix!="") {
      suffix = paste0("_", output_file_suffix)
      pdf(paste0(PTMSEA_OUTDIR, "/Network/kinase_substrates/", pair,"/", ID2, suffix, ".pdf"),
               h=(L/2.5)+3, w=(L/2.5)+3)
    } else {
      suffix = ""
      pdf(paste0(PTMSEA_OUTDIR, "/Network/kinase_substrates/", pair,"/", ID2, suffix, ".pdf"),
               h=(L/4)+2.5, w=(L/4)+2.8)
    }
    print(plot(g, vertex.label.dist = 0, 
               layout = layout, 
               # vertex.label.degree = pi/2, 
               vertex.label.cex = 0.8,
               vertex.label = igraph::V(g)$label,
               vertex.color = igraph::V(g)$color,
               vertex.frame.color = igraph::V(g)$color,
               # vertex size = logFC*4
               vertex.size = abs(round(c(6,as.numeric(mapping_regulation*4)))), 
               vertex.label.color = "black",
               edge.color = igraph::E(g)$color))
    dev.off()
    print(plot(g, vertex.label.dist = 0, 
               layout = layout, 
               # vertex.label.degree = pi/2, 
               vertex.label.cex = 0.8,
               vertex.label = igraph::V(g)$label,
               vertex.color = igraph::V(g)$color,
               vertex.frame.color = igraph::V(g)$color,
               # vertex size = logFC*4
               vertex.size = abs(round(c(6,as.numeric(mapping_regulation*4)))),
               vertex.label.color = "black",
               edge.color = igraph::E(g)$color))
  })
}





#' @import igraph
#' @import dplyr
#' @param significance_statistic this can be "pvalue", "qvalue" or "p.adjust".
#' @param uniprot_provided if uniprot ID is included in mapping_ID, and you want uniprot ID to be shown in lable, set this to TRUE.
#' @param PTMsigDB_collection_file the path to the PTMsigDB collection file in GMT format.
#' @param PTMSEA_output the output from PTMSEA, which is a data frame with columns: id, Signature.set.description, Signature.set.size, and significance_statistic etc.
#' @param mapping_ID a named character vector mapping phosphosite FlankingRegion to uniprot ID/PhosLocation.
#' @param mapping_regulation a named numeric vector mapping phosphosite to regulation (usually logFC or t-statistics).
#' @param proteomics a character vector of proteomics protein uniprot IDs that are aquired from the same study system as the phosphoproteomics data.
#' @param outdir_ppi the output directory for the PPI network.
#' 
ppiNetwork4substrates = function(limma_output, PTM.FlankingRegion4PTMSEAanalysis, PTMSEA_output,       
                                 PTMSEA_OUTPUT, significance_cutoff=0.05, 
                                 significance_statistic="fdr.pvalue", mapping_ID,
                                 mapping_regulation, proteomics, outdir_ppi,
                                 PTMsigDB_collection_file,
                                 uniprot_provided=T) {
  dir.create(outdir_ppi, recursive = T, showWarnings = F)
  limma_output$PTM.FlankingRegion2 = paste0(limma_output$PTM.FlankingRegion, "-p")
  limma_output = limma_output %>%
      filter(PTM.FlankingRegion2%in%PTM.FlankingRegion4PTMSEAanalysis)
  limma_output$uniprotID = limma_output$Phosphosite %>% sub("_.*", "", .)
  limma_output$PhosLocation = limma_output$Phosphosite %>% sub(".*_", "", .)
  limma_output = limma_output%>% dplyr::select(uniprotID, PhosLocation, PTM.FlankingRegion2, everything()) 
  
  # get significant PTMsigDB terms
  sigIDs = 
    PTMSEA_output[PTMSEA_output[,significance_statistic,drop=T]<significance_cutoff, ] %>% 
    pull(id) %>% .[!is.na(.)]  
  
  # import PTMsigDB.v2 collection
  ptmsigdb = GSEABase::getGmt(PTMsigDB_collection_file)
  ptmsigdb0 = lapply(ptmsigdb, function(K) {
    data.frame(Term=K@setName, Phosphosite=K@geneIds)
  }) %>% Reduce(rbind, .)
  ptmsigdb2 = as.data.frame(stringr::str_split_fixed(ptmsigdb0$Phosphosite,";",2)) %>%
    setNames(c("Phosphosite", "Regulation"))
  ptmsigdb3 = data.frame(Term=ptmsigdb0$Term,
                               Phosphosite=ptmsigdb2$Phosphosite,
                               Regulation=ptmsigdb2$Regulation)  
  
  # Get STRING IDs for the proteomics data
  # Load STRING database
  string_db = STRINGdb::STRINGdb$new(version = "11.5", species = 9606, score_threshold = 700)
  # Map proteomics protein names to STRING IDs
  mapped_proteomics = string_db$map(data.frame(protein=proteomics), 
                            "protein", removeUnmappedRows=TRUE)  
  # Get STRING IDs for our significant kianse substrates all together
  # Get the substrates from PTMsigDB for the current ID/term (mostly kinase-related terms)
  substrates_all = ptmsigdb3 %>% filter(Term%in%sigIDs) %>% pull(Phosphosite)
  # Get the substrate phosphosite flanking regions that are with the input dataset
  shared_all = intersect(substrates_all, limma_output$PTM.FlankingRegion2)
  # Get the uniprot ID for the substrate phosphosite flanking region
  shared_uniprotID_all = mapping_ID[shared_all] %>% as.character() %>% sub("\n.*", "", .)
  # Map shared_uniprotID (shared between PTMsigDB substrates and the input dataset) to STRING IDs
  mapping_substrates = string_db$map(data.frame(protein=shared_uniprotID_all), 
                                  "protein", removeUnmappedRows=TRUE) 
  print("When one uniprot ID mapped to multiple STRING IDs, only keep the 1st one!")
  mapping_substrates =  mapping_substrates[!duplicated(mapping_substrates$protein), ] ###  
  
  # Map STRING ID to hgnc_symbol
  # Get immediate neighbors for those shared substrates between PTMsigDB and limma output
  neighbors_all = string_db$get_neighbors(mapping_substrates$STRING_id)
  # Intersect those immediate neighbors with the proteomics STRING IDs
  proteomics_neighbors_all = intersect(neighbors_all, mapped_proteomics$STRING_id)
  # Collect all substrates and their immediate neighbors that are also in the proteomics dataset
  all_proteins_all = c(mapping_substrates$STRING_id, proteomics_neighbors_all)
  # Map the STRING IDs of all substrates and their immediate neighbors back to gene names
  mapped_back_all_proteins_all = 
    string_db$add_proteins_description(data.frame(STRING_id=all_proteins_all)) %>%
    mutate(hgnc_symbol=preferred_name)
  # Prepare for ID conversion
  if (!file.exists(paste0(PTMSEA_OUTDIR, "/Network/PPI/ensembl.rds"))) {
    ensembl = biomaRt::useEnsembl(biomart="genes", 
                         dataset="hsapiens_gene_ensembl", mirror="useast")
    saveRDS(ensembl, paste0(PTMSEA_OUTDIR, "/Network/PPI/ensembl.rds"))
    # ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")    
  } else {
    ensembl = readRDS(paste0(PTMSEA_OUTDIR, "/Network/PPI/ensembl.rds"))
  }
  temp_all_proteins_all = biomaRt::getBM(
        attributes = c("hgnc_symbol", "uniprotswissprot", "uniprot_gn_symbol"),
        filters = "hgnc_symbol",
        values = mapped_back_all_proteins_all$preferred_name,
        mart = ensembl
  ) %>% filter(uniprotswissprot!="") %>% as.data.frame() %>% 
        dplyr::select(hgnc_symbol, uniprotswissprot) %>% unique() 
  # There are duplicated hgnc_symbol? if there are, it means some hgnc_symbol mapped to multiple uniprotswissprot IDs
  duplicated_hgnc_symbol = temp_all_proteins_all$hgnc_symbol[duplicated(temp_all_proteins_all$hgnc_symbol)] 
  if (length(duplicated_hgnc_symbol)>0) {
    uniprot4duplicated_hgnc_symbol = temp_all_proteins_all %>% 
      filter(temp_all_proteins_all$hgnc_symbol%in%duplicated_hgnc_symbol) %>% pull(uniprotswissprot)
    # Check whether those uniprotswissprot for duplicated_hgnc_symbol are in our data
    all_IDs_from_dat = sapply(mapping_ID, function(I) {strsplit(I, "\n")}) %>% 
      unlist() %>% as.character()
    uniprot4duplicated_hgnc_symbol_in_our_dat = uniprot4duplicated_hgnc_symbol %>% 
      .[.%in% all_IDs_from_dat]
    # If the duplicated hgnc_symbol are in our data, we will only keep those in our data, otherwise, we will keep the first one
    if (length(uniprot4duplicated_hgnc_symbol_in_our_dat)>0) {
      dup_hgnc_symbol_inOurData = 
        filter(temp_all_proteins_all, uniprotswissprot%in%uniprot4duplicated_hgnc_symbol_in_our_dat) %>%
        pull(hgnc_symbol)
      rest = filter(temp_all_proteins_all, !hgnc_symbol%in%dup_hgnc_symbol_inOurData)
      hit = filter(temp_all_proteins_all, hgnc_symbol%in%dup_hgnc_symbol_inOurData)
      temp_all_proteins_all = rbind(rest[!duplicated(rest$hgnc_symbol),], 
                   filter(hit,uniprotswissprot%in%uniprot4duplicated_hgnc_symbol_in_our_dat))
    } else {
      temp_all_proteins_all = temp_all_proteins_all[!duplicated(temp_all_proteins_all$hgnc_symbol),]
    }
  }
  # Combine STRING ID, uniprotswissprot and hgnc_symbol together
  temp_all_proteins_all = full_join(temp_all_proteins_all, mapped_back_all_proteins_all, by="hgnc_symbol") #%>% na.omit() # They share hgnc_symbol
  mapping_substrates = setNames(mapping_substrates$protein, mapping_substrates$STRING_id)
  which_substrates = temp_all_proteins_all$STRING_id%in%names(mapping_substrates)
  temp_all_proteins_all$uniprotswissprot[which_substrates] = 
    mapping_substrates[temp_all_proteins_all$STRING_id[which_substrates]]
  temp_all_proteins_all = na.omit(temp_all_proteins_all)
  # Retrieve PPI interactions from string_db for all substrates and their immediate neighbors that are also in the reference dataset
  interactions_all = string_db$get_interactions(temp_all_proteins_all$STRING_id)   
  
  # On Mac and Windows, multisession runs tasks in separate R sessions (processes)
  future::plan("multisession", workers = 4) 
  future.apply::future_lapply(sigIDs, function(ID) {
    substrates = ptmsigdb3 %>% filter(Term==ID) %>% pull(Phosphosite)
    shared = intersect(substrates, limma_output$PTM.FlankingRegion2)
    shared_uniprotID = mapping_ID[shared] %>% as.character() %>% sub("\n.*", "", .)
    mapped_interest = mapping_substrates[mapping_substrates%in%shared_uniprotID]   
    # Get immediate neighbors for those shared substrates between PTMsigDB and limma output
    neighbors = string_db$get_neighbors(names(mapped_interest))
    # Intersect those immediate neighbors with the proteomics STRING IDs
    proteomics_neighbors = intersect(neighbors, mapped_proteomics$STRING_id)
    # Collect all substrates and their immediate neighbors that are also in the proteomics dataset
    all_proteins = c(names(mapped_interest), proteomics_neighbors)
    # Map the STRING IDs of all substrates and their immediate neighbors back to gene names
    mapped_back_all_proteins = mapped_back_all_proteins_all %>% 
      filter(STRING_id %in% all_proteins)  
    if (uniprot_provided) {      
      # Convert hgnc_symbol/gene names to uniprotswissprot
      temp_all_proteins = 
        filter(temp_all_proteins_all, hgnc_symbol%in%mapped_back_all_proteins$preferred_name)
      # Map the STRING ID back to protein uniprot ID and gene name
      mapping_back_all_proteins = 
        setNames(paste0(temp_all_proteins$hgnc_symbol, "\n(", temp_all_proteins$uniprotswissprot,")"), 
                              temp_all_proteins$STRING_id)
      mapping_back_all_proteins_forKinaseSubstrates = 
        setNames(paste0(temp_all_proteins$hgnc_symbol, "\n(", temp_all_proteins$uniprotswissprot, ")"), 
                                                  temp_all_proteins$uniprotswissprot)
    } else {
      mapping_back_all_proteins = setNames(mapped_back_all_proteins$preferred_name, mapped_back_all_proteins$STRING_id)
      mapping_back_all_proteins_forKinaseSubstrates = setNames(paste0(temp_all_proteins$hgnc_symbol), 
                                                  temp_all_proteins$uniprotswissprot)}    
    # Add in the network between PTMsigDB terms (e.g. kinases) and their members (e.g. substrates)
    kinase_substrate_ineractions = data.frame(from=rep(ID, length(shared_uniprotID)), 
                                              to=shared_uniprotID, combined_score="NA")
    # Only keep interactions with To being String IDs that are included in mapping_back_all_proteins_forKinaseSubstrates
    kinase_substrate_ineractions = 
      kinase_substrate_ineractions[kinase_substrate_ineractions$to%in%names(mapping_back_all_proteins_forKinaseSubstrates),] 
    if (uniprot_provided) {
      # Change the substrate names to uniprot IDs and gene names
      kinase_substrate_ineractions$to = 
        mapping_back_all_proteins_forKinaseSubstrates[kinase_substrate_ineractions$to] 
    }
    # Remove duplicated rows
    kinase_substrate_ineractions = unique(kinase_substrate_ineractions)  
    # Only keep interactions with both From and To being String IDs that are included in mapping_back_all_proteins (Map the STRING IDs of all substrates and their immediate neighbors back to gene names/uniprot IDs)
    interactions = 
      interactions_all[interactions_all$from%in%names(mapping_back_all_proteins)&
        interactions_all$to%in%names(mapping_back_all_proteins), ] %>%
        unique()
    # Change the From and To names from String IDs to gene names/uniprot IDs
    interactions$from = mapping_back_all_proteins[as.character(interactions$from)]
    interactions$to = mapping_back_all_proteins[as.character(interactions$to)]
    # Combine the kinase-substrate interactions with the PPI interactions
    interactions = rbind(kinase_substrate_ineractions, interactions)    
    
    # Define network edge colors
    edge.colors = rep("snow3", nrow(interactions))
    # mapping_back_all_proteins (map the STRING IDs of all substrates and their immediate neighbors back to gene names/uniprot IDs)
    # mapped_interest: map shared_uniprotID (shared between PTMsigDB substrates and the input dataset) to STRING IDs
    # Set the colors for the edges that connect the kinase substrates to their proteomics immediate neighbors to orange.
    edge.colors[interactions$from%in%as.character(mapping_back_all_proteins[names(mapped_interest)])|
                  interactions$to%in%as.character(mapping_back_all_proteins[names(mapped_interest)])] = 
                  scales::alpha("orange", 0.5)    
    # Add in the regulation (up or down) of the kinase-substrate interactions from PTMsigDB
    # Subset the PTMsigDB to only keep those that match the ID
    t1 = ptmsigdb3 %>% filter(Term==ID)
    # Mapping of substrate phosphosite flankingRegions to Uniprot IDs
    t2 = mapping_ID[t1$Phosphosite] %>% na.omit()
    # Add in Regulation into the interactions
    interactions2 = filter(t1, t1$Phosphosite%in%names(t2)) %>% 
      left_join(data.frame(Phosphosite=names(t2),Names=t2,stringsAsFactors = FALSE)) %>%
      mutate(uniprotID = gsub("\n.*", "", Names)) %>% right_join(
        filter(interactions, from%in%ID) %>% 
          mutate(uniprotID = gsub(".*\n", "", to)) %>%
          mutate(uniprotID = gsub("[(]|[)]", "", uniprotID))
      )  %>% left_join(
        data.frame(Phosphosite = names(mapping_regulation),
                   effect = mapping_regulation, stringsAsFactors = FALSE),
      )
    mapping_TO_regualtion = setNames(interactions2$Regulation, interactions2$to)
    regulations_ID = mapping_TO_regualtion[filter(interactions,from%in%ID) %>% pull(to)]    
    # Set the edge colors for the edges that connect the kinase to their substrates according their regulation (up:orange or down:turqoise) as stated in PTMsigDB
    edge.colors[interactions$from%in%ID] = 
      ifelse(regulations_ID=="u", scales::alpha("purple",0.5), scales::alpha("darkturquoise",0.5))

    # Set the edge widths for the edges from or to the kinase substrates to 0.75, and all other edges to 0.1
    edge.widths=rep(0.1, nrow(interactions))
    network_vertexIDs_4substrates = as.character(mapping_back_all_proteins[names(mapped_interest)])%>%na.omit()
    edge.widths[interactions$from%in%network_vertexIDs_4substrates%>%na.omit()|
                  interactions$to%in%network_vertexIDs_4substrates%>%na.omit()] = 0.75    
    # Create igraph object
    # library(igraph)
    ppi_network = igraph::graph_from_data_frame(interactions[, 1:2], directed = FALSE)    
    #col_gradients = colorRampPalette(c("orange", "purple"))(length(mapped_interest$STRING_id))
    # Set label color of the kinase to black, set the label colors of kinase substrate to dark grey, and set the rest as grey:
    vertex.label.colors = rep(scales::alpha("black", 0.5), length(igraph::V(ppi_network)$name))
    vertex.label.colors[igraph::V(ppi_network)$name %in%ID] = scales::alpha("black", 1)
    vertex.label.colors[igraph::V(ppi_network)$name%in%network_vertexIDs_4substrates] = 
                scales::alpha("black", 0.8)    
    # Set the vertex colors of the kinase substrates based on the effect of the kinase substrates: orange: up-regulated, darkturquoise: down-regulated
    mapping_TO_effect = setNames(interactions2$effect, interactions2$to)
    # Here effect could be logFC or the t statistics from Limma analysis
    effect = mapping_TO_effect[igraph::V(ppi_network)$name %>% 
                                 .[.%in%network_vertexIDs_4substrates]]

    vertex.colors = rep("snow3", length(igraph::V(ppi_network)$name))
    vertex.colors[igraph::V(ppi_network)$name %in%ID] = scales::alpha("purple", 0.5)
    vertex.colors[igraph::V(ppi_network)$name%in%network_vertexIDs_4substrates] =
        ifelse(effect<0, scales::alpha("darkturquoise",0.5), scales::alpha("orange",0.5))    
    # Set the vertex sizes of the kinase substrates based on their effect*2, set the rest to 1, and set the kinase to 6. 
    vertex.sizes = rep(1, length(igraph::V(ppi_network)$name))
    vertex.sizes[igraph::V(ppi_network)$name%in%network_vertexIDs_4substrates] =
        abs(round(as.numeric(effect))*2)
    vertex.sizes[igraph::V(ppi_network)$name %in%ID] = 10    
    coords = igraph::layout_in_circle(ppi_network)

    L = length(igraph::V(ppi_network))/60
    ID2 = sub("/", "_", ID)
    pdf(paste0(outdir_ppi, "/", ID2, "_substrate_PPI.pdf"), h=7*L+2, w=10*L+2)
      print(plot(ppi_network,
         vertex.shape="circle",
         vertex.label.cex=0.7, 
         vertex.label.color=vertex.label.colors,
         edge.color = edge.colors,
         vertex.label.dist = 0,   # distance outward,
         edge.width=edge.widths,
         layout = coords,
         main=ID, asp = 0.7,
         vertex.color = vertex.colors,
         vertex.frame.color = "white", #vertex.colors,
         vertex.size = vertex.sizes,
         edge.curved=0))
    dev.off()
    print(plot(ppi_network,
         vertex.shape="circle",
         vertex.label.cex=0.7, 
         vertex.label.color=vertex.label.colors,
         edge.color = edge.colors,
         vertex.label.dist = 0,   # distance outward,
         edge.width=edge.widths,
         layout = coords,
         main=ID, asp = 0.7,
         vertex.color = vertex.colors,
         vertex.frame.color = "white", # vertex.colors,
         vertex.size = vertex.sizes,
         edge.curved=0))
  })
}
