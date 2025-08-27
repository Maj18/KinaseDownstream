
#' Run PTM-SEA on the limma results logFC
#' @import dplyr
#' @param species This can be "human", "mouse" or "rat"
#' @param inputtype This can be "flanking", "sitegrpid", or "uniprot"
#' @param sample.norm.type Recommend "none" for logFC, and "rank" for gene expression data
#' @param weight This will decide how much the top ranked feature been emphasized, the bigger the value is, the more weight the top features will get, recommend 0.75 for PTM data
#' @param correl.type Recommend "rank", more robust
#' @param statistic Recommend "area.under.RES" for PTM data
#' 
runPTMSEA = function(limma_rslt, PTMSEA_OUTDIR, species="human", 
  inputtype="flanking", sample.norm.type = "none", 
  weight = 0.75, correl.type = "rank", statistic = "area.under.RES", spare.cores = 4,
  output.score.type = "NES", nperm = 1000, min.overlap = 5, extended.output = TRUE,
  global.fdr = FALSE, export.signat.gct = T, param.file = T) {
  limma_rslt_4gct = lapply(seq_along(limma_rslt), function(i){
    rslt = limma_rslt[[i]] %>% arrange(-logFC)
    rownames(rslt) = NULL
    temp = rslt[, c("Feature", "logFC")]
    colnames(temp) = c("Feature", names(limma_rslt)[i])
    temp$Feature = paste0(temp$Feature, "-p")
    temp
  }) %>% Reduce(full_join, .) %>%
    .[!duplicated(.$Feature), ] 
  rownames(limma_rslt_4gct) = NULL
  limma_rslt_4gct = 
    tibble::column_to_rownames(limma_rslt_4gct, var="Feature")
  cmapR::write_gct(as.matrix(limma_rslt_4gct)%>%
    cmapR::GCT(.),
  #cmapR::new("GCT",mat=.), 
    paste0(PTMSEA_OUTDIR ,"/PTM-SEA"), precision=2)
  # run PTM-SEA
  input_gct_file = list.files(path=PTMSEA_OUTDIR, pattern = "PTM-SEA", full.names = TRUE)

  # set up the PTMsigDB file:
  if (species=="human") {
    if (inputtype=="flanking") {
      gmt_file = system.file("extdata", "ptm.sig.db.all.flanking.human.v2.0.0.gmt", 
        package="KinaseDownstream")
    } else if (inputtype=="sitegrpid") {
      gmt_file = system.file("extdata", "ptm.sig.db.all.sitegrpid.human.v2.0.0.gmt", 
        package="KinaseDownstream")
    } else if (inputtype=="uniprot") {
      gmt_file = system.file("extdata", "ptm.sig.db.all.uniprot.human.v2.0.0.gmt", 
        package="KinaseDownstream")
    }
  } else if (species=="mouse") {
    if (inputtype=="flanking") {
      gmt_file = system.file("extdata", "ptm.sig.db.all.flanking.mouse.v2.0.0.gmt", 
        package="KinaseDownstream")
    } else if (inputtype=="sitegrpid") {
      gmt_file = system.file("extdata", "ptm.sig.db.all.sitegrpid.mouse.v2.0.0.gmt", 
        package="KinaseDownstream")
    } else if (inputtype=="uniprot") {
      gmt_file = system.file("extdata", "ptm.sig.db.all.uniprot.mouse.v2.0.0.gmt", 
        package="KinaseDownstream")
    }    
  } else if (species=="rat") {
    if (inputtype=="flanking") {
      gmt_file = system.file("extdata", "ptm.sig.db.all.flanking.rat.v2.0.0.gmt", 
        package="KinaseDownstream")
    } else if (inputtype=="sitegrpid") {
      gmt_file = system.file("extdata", "ptm.sig.db.all.sitegrpid.rat.v2.0.0.gmt", 
        package="KinaseDownstream")
    } else if (inputtype=="uniprot") {
      gmt_file = system.file("extdata", "ptm.sig.db.all.uniprot.rat.v2.0.0.gmt", 
        package="KinaseDownstream")
    }    
  }
  
  set.seed(123)
  invisible(capture.output(ssGSEA2::run_ssGSEA2(input_gct_file,
        output.prefix = "PTMSEA_OUTPUT",
        gene.set.databases = gmt_file,
        output.directory = PTMSEA_OUTDIR,
        sample.norm.type = sample.norm.type,
        weight = weight, 
        correl.type = correl.type, 
        statistic = statistic,
        spare.cores = spare.cores,
        output.score.type = output.score.type, 
        nperm = nperm, 
        min.overlap = min.overlap, 
        extended.output = extended.output, 
        global.fdr = global.fdr,
        export.signat.gct = export.signat.gct ,
        param.file = param.file,
        log.file = paste0(PTMSEA_OUTDIR, "/run.log")))) 
                    
  return(limma_rslt_4gct)
}



# TODO: Add a function to get the substrate phosphosites ID that has been use in PTMSEA analysis, modify the rest of the code accordingly
getPTMsubstrates4PTMSEAanalysis = function(limma_rslt) {
  limma_rslt_4gct = lapply(seq_along(limma_rslt), function(i){
    rslt = limma_rslt[[i]] %>% arrange(-logFC)
    rslt
  }) %>% Reduce(full_join, .) %>%
    .[!duplicated(.$PTM.FlankingRegion), ] 
  if (!grepl("-p$", limma_rslt_4gct$Phosphosite[1])) {
  	limma_rslt_4gct$Phosphosite = paste0(limma_rslt_4gct$Phosphosite, "-p")
  }
  limma_rslt_4gct$Phosphosite
}

getPTMsubstrates4PTMSEAanalysis_uniprot = function(limma_rslt) {
  limma_rslt_4gct = lapply(seq_along(limma_rslt), function(i){
    rslt = limma_rslt[[i]] %>% arrange(-logFC)
    rslt
  }) %>% Reduce(full_join, .) %>%
    .[!duplicated(.$Feature), ] 
  if (!grepl("-p$", limma_rslt_4gct$Feature[1])) {
  	limma_rslt_4gct$Feature = paste0(limma_rslt_4gct$Feature, "-p")
  }
  limma_rslt_4gct$Feature
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
  ## Only significant kinases
  fdrs = grep("fdr.pvalue", colnames(gct_data), value=T)
  for (fdr in fdrs) {
	  N = gsub("fdr.pvalue.", "", fdr)
	  kinaseList = gct_data %>% filter(gct_data[,fdr,drop=T]<sig.thresh) %>% pull(id) %>% 
		  grep("KINASE", ., value=T) %>% gsub("KINASE-PSP_", "", .) %>% 
		  gsub("KINASE-iKiP_", "", .) %>% stringr::str_split(., "/|[.]") %>% 
		  unlist() %>% gsub("-.*", "", .) %>% paste0(., collapse=",")
    dir.create(PTMSEA_OUTDIR, showWarnings = FALSE, recursive = TRUE)
    dir.create(paste0(PTMSEA_OUTDIR, "/KinaseGroup/", N), 
                  showWarnings = FALSE, recursive = TRUE)
	  writeLines(kinaseList, 
		  paste0(PTMSEA_OUTDIR, "/KinaseGroup/", N, "/SigKinaseList_", N, "_onlySigKinases.txt"))
  }

  # Prepare signficant kinase names for mapping onto the Manning kinase tree (http://kinhub.org/kinmap/index.html)
  ## All kinases
  fdrs = grep("fdr.pvalue", colnames(gct_data), value=T)
  for (fdr in fdrs) {
	  N = gsub("fdr.pvalue.", "", fdr)
	  kinaseList = gct_data %>% pull(id) %>% 
		  grep("KINASE", ., value=T) %>% gsub("KINASE-PSP_", "", .) %>% 
		  gsub("KINASE-iKiP_", "", .) %>% stringr::str_split(., "/|[.]") %>% 
      unlist() %>% gsub("-.*", "", .) %>% paste0(., collapse=",")
    dir.create(PTMSEA_OUTDIR, showWarnings = FALSE, recursive = TRUE)
    dir.create(paste0(PTMSEA_OUTDIR, "/KinaseGroup/", N), 
                  showWarnings = FALSE, recursive = TRUE)
	  writeLines(kinaseList, 
		  paste0(PTMSEA_OUTDIR, "/KinaseGroup/", N, "/SigKinaseList_", N, "_allKinases.txt"))
  }

  return(list(ptmsea_rslt=ptmsea_rslt, gct4plot=gct4plot, plot=plot))
}


#' @import dplyr
#' @description Prepare the input for PhosNetVis
prepInput4PhosNetVis = function(pair, limma_output, PTMsubstrates4PTMSEAanalysis, 
                                PTMSEA_output, significance_cutoff=0.05, 
                                significance_statistic="fdr.pvalue", 
                                PTMsigDB_collection_file) {
  limma_output$PTM.FlankingRegion2 = paste0(limma_output$PTM.FlankingRegion, "-p")
  limma_output = limma_output %>%
      .[limma_output$Phosphosite%in%PTMsubstrates4PTMSEAanalysis, ]
  limma_output$uniprotID = limma_output$Phosphosite %>% sub(":.*", "", .)
  limma_output$PhosLocation = limma_output$Phosphosite %>% sub(".*:", "", .)
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
#' @param PTMsubstrates4PTMSEAanalysis a character vector of PTM flanking regions that have been used in PTMSEA.
#' @param limma_rslt a list of data frames containing the limma results.
#' 
processLimmaResult = function(limma_rslt, PTMsubstrates4PTMSEAanalysis, 
                              significance_statistic="adj.P.Val") {
  maps = lapply(limma_rslt, function(rslt) {
    rslt$uniprotID = rslt$Phosphosite %>% sub(":.*", "", .)
    rslt$PTM.FlankingRegion2 = paste0(rslt$PTM.FlankingRegion, "-p")
    rslt$PhosLocation = rslt$Phosphosite %>% sub(".*:", "", .)
    rslt = filter(rslt, Phosphosite%in%PTMsubstrates4PTMSEAanalysis)
    mapUniprotPhosLocation2FlankingRegion2 = 
      setNames(paste0(rslt$uniprotID,"\n",rslt$PhosLocation),rslt$PTM.FlankingRegion2)
    mapUniprotPhosLocationFlankingRegion2FlankingRegion2 = 
      setNames(paste0(rslt$uniprotID,"\n",rslt$PhosLocation,
                    "\n",rslt$PTM.FlankingRegion), rslt$PTM.FlankingRegion2)
    mapLogFC2FlankingRegion2 = 
    	setNames(paste0(rslt$logFC,";",rslt[,significance_statistic,drop=T]), 
    		rslt$PTM.FlankingRegion2) 
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
KinaseNetwork4substrates_uniprot = function(pair, PTMsubstrates4PTMSEAanalysis, limma_output, PTMSEA_output, 
                                significance_cutoff4PTMSEA=1, species="mouse",
                                significance_statistic4PTMSEA="fdr.pvalue", mapping_ID,
                                significance_cutoff4limma=0.05, logFCcutoff4limma=0.5, PTMSEA_OUTDIR,
                                mapping_regulation, output_file_suffix="") {
  limma_output$Feature = paste0(limma_output$Feature, "-p")
  limma_output = limma_output %>%
      filter(Feature %in% PTMsubstrates4PTMSEAanalysis)
  limma_output$uniprotID = limma_output$Feature %>% gsub(";.*", "", .)
  limma_output$PhosLocation = limma_output$Feature %>% gsub(".*;", "", .)
  limma_output = limma_output%>% dplyr::select(uniprotID, PhosLocation, Feature, everything()) 
  # get significant PTMsigDB terms
  sigIDs = PTMSEA_output[PTMSEA_output[,significance_statistic,drop=T]<significance_cutoff,] %>% 
    pull(id) %>% .[!is.na(.)] 

  # import PTMsigDB.v2 collection
  if (species=="moust") {
  	PTMsigDB_collection_file = system.file("extdata", "ptm.sig.db.all.uniprot.mouse.v2.0.0.gmt", 
        package="KinaseDownstream")
  } else if (species=="human") {
  	PTMsigDB_collection_file = system.file("extdata", "ptm.sig.db.all.uniprot.human.v2.0.0.gmt", 
        package="KinaseDownstream")
  	} else if (species=="rat") {
  	PTMsigDB_collection_file = system.file("extdata", "ptm.sig.db.all.uniprot.rat.v2.0.0.gmt", 
        package="KinaseDownstream")
  }
  ptmsigdb = GSEABase::getGmt(PTMsigDB_collection_file)
  ptmsigdb0 = lapply(ptmsigdb, function(K) {
    data.frame(Term=K@setName, Phosphosite=K@geneIds)
  }) %>% Reduce(rbind, .)
  ptmsigdb2 = as.data.frame(stringr::str_split_fixed(ptmsigdb0$Phosphosite,";",2)) %>%
    setNames(c("Uniprot", "PhosLocation", "Regulation")) #####
  ptmsigdb3 = data.frame(Term=ptmsigdb0$Term,
                               Phosphosite=paste0(ptmsigdb2$Uniprot,";",ptmsigdb2$PhosLocation),
                               Regulation=ptmsigdb2$Regulation)  
  lapply(sigIDs, function(ID) {
    # print(ID)
    substrates = ptmsigdb3 %>% filter(Term==ID) %>% pull(Phosphosite)
    shared = intersect(substrates, limma_output$Feature) # flanking regions
    mapping_regulation_db = ptmsigdb3 %>% filter(Term==ID) %>% 
      filter(Phosphosite%in%shared)
    mapping_regulation_db = setNames(mapping_regulation_db$Regulation, mapping_regulation_db$Phosphosite) 
    regulation = mapping_regulation_db[shared]
    edge = data.frame(From=rep(ID,length(shared)), To=shared) 
    # library(igraph)
    g = igraph::graph_from_data_frame(d = edge, directed = FALSE)
    mapping_all = edge%>%mutate(Feature=To) %>% left_join(limma_output) %>%
    	mutate(effect=logFC, significance=adj.P.Val) %>% 
    	tibble::column_to_rownames(var="Feature")
   	mapping_regulation = mapping_all[shared, "effect", drop=T] #############
    regulation_effect = abs(round(c(6,mapping_all[igraph::V(g)$label[-1],"effect"])))
    regulation_significance = (mapping_all[igraph::V(g)$label[-1],"significance"]<significance_cutoff4limma)&
                              (abs(mapping_all[igraph::V(g)$label[-1],"effect"])>logFCcutoff4limma)
    igraph::V(g)$label[-1][regulation_significance] = paste0(igraph::V(g)$label[-1][regulation_significance], "*")
    igraph::E(g)$color = ifelse(regulation=="u", scales::alpha("orange",0.25), scales::alpha("darkturquoise",0.25))
    igraph::V(g)$color = c(scales::alpha("black", 0.25), ifelse(mapping_regulation<0, 
                                                scales::alpha("darkturquoise",0.5), scales::alpha("orange",0.5)))

    igraph::V(g)$label = c(igraph::V(g)[[1]]$name,    
        gsub("-p","",igraph::V(g)[2:length(igraph::V(g))]$name)%>%gsub(";","\n",.))

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
               vertex.size = regulation_effect*4, 
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
               vertex.size = regulation_effect*4,
               vertex.label.color = "black",
               edge.color = igraph::E(g)$color))
  })
}





#' @param significance_statistic this can be "pvalue" or "fdr.pvalue".
#' @param mapping_ID a named character vector mapping phosphosite FlankingRegion to uniprot ID/PhosLocation/FlankingRegion
#' @param PTMSEA_OUTDIR output directory for PTMSEA results
#' @param mapping_regulation a named numeric vector mapping phosphosite to regulation (usually logFC or t-statistics)
#' @param PTMsigDB_collection_file the path to the PTMsigDB collection file in GMT format
#' @param output_file_suffix a suffix for the output file names, default is empty string. If you want to add FlankingRegion to the plot, set output_file_suffix = "wPhosphosites", otherwise set it to "".
#' @import igraph
#' @import dplyr
KinaseNetwork4substrates = function(pair, PTMsubstrates4PTMSEAanalysis, limma_output, PTMSEA_output, 
                                significance_cutoff4PTMSEA=1, 
                                significance_statistic4PTMSEA="fdr.pvalue", mapping_ID,
                                significance_cutoff4limma=0.05, logFCcutoff4limma=0.5, PTMSEA_OUTDIR,
                                mapping_regulation, output_file_suffix="", PTMsigDB_collection_file) {
  limma_output$PTM.FlankingRegion2 = paste0(limma_output$PTM.FlankingRegion, "-p")
  limma_output = limma_output %>%
      filter(Phosphosite%in%PTMsubstrates4PTMSEAanalysis)
  limma_output$uniprotID = limma_output$Phosphosite %>% sub(":.*", "", .)
  limma_output$PhosLocation = limma_output$Phosphosite %>% sub(".*:", "", .)
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

    mapping_all = edge%>%mutate(FlankingRegion=To) %>% 
      left_join(data.frame(FlankingRegion=names(mapping_ID),Name=as.character(mapping_ID))) %>% 
      left_join(data.frame(
        FlankingRegion=names(mapping_regulation), 
        effect = stringr::str_split_fixed(mapping_regulation,";",2)%>%
            as.data.frame()%>%pull(V1) %>% as.numeric(), 
        significance = stringr::str_split_fixed(mapping_regulation,";",2)%>%
            as.data.frame()%>%pull(V2) %>% as.numeric())) %>%
        tibble::column_to_rownames(var="Name")
    regulation_effect = abs(round(c(6,mapping_all[igraph::V(g)$label[-1],"effect"])))
    regulation_significance = (mapping_all[igraph::V(g)$label[-1],"significance"]<significance_cutoff4limma)&
                              (abs(mapping_all[igraph::V(g)$label[-1],"effect"])>logFCcutoff4limma)
    igraph::V(g)$label[-1][regulation_significance] = paste0(igraph::V(g)$label[-1][regulation_significance], "*")

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
               vertex.size = regulation_effect*4, 
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
               vertex.size = regulation_effect*4,
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
ppiNetwork4substrates_STRING = function(limma_output, PTMsubstrates4PTMSEAanalysis, PTMSEA_output,       
                                 PTMSEA_OUTDIR, significance_cutoff4PTMSEA=1, 
                                 significance_statistic4PTMSEA="fdr.pvalue", mapping_ID,
                                 significance_cutoff4limma=0.05, logFCcutoff4limma=0.5, 
                                 mapping_regulation, proteomics, outdir_ppi,
                                 PTMsigDB_collection_file,
                                 uniprot_provided=T) {
  dir.create(outdir_ppi, recursive = T, showWarnings = F)
  limma_output$PTM.FlankingRegion2 = paste0(limma_output$PTM.FlankingRegion, "-p")
  limma_output = limma_output %>%
      filter(Phosphosite%in%PTMsubstrates4PTMSEAanalysis)
  limma_output$uniprotID = limma_output$Phosphosite %>% sub(":.*", "", .)
  limma_output$PhosLocation = limma_output$Phosphosite %>% sub(".*:", "", .)
  limma_output = limma_output%>% dplyr::select(uniprotID, PhosLocation, PTM.FlankingRegion2, everything()) 
  
  # get significant PTMsigDB terms
  sigIDs = 
    PTMSEA_output[PTMSEA_output[,significance_statistic4PTMSEA,drop=T]<significance_cutoff, ] %>% 
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
  shared_uniprotID_all = mapping_ID[shared_all] %>% as.character() %>% 
    sub("\n.*", "", .) %>% unique() #Here, we will get duplcates for phosphosites of the same proteins
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
  
  # Convert genes names to uniprot.ID
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
  ## mapped_back_all_proteins_all: STRING ID + gene symbol
  ## temp_all_proteins_all: gene symbol + uniprot ID
  temp_all_proteins_all = full_join(temp_all_proteins_all, mapped_back_all_proteins_all, by="hgnc_symbol") #%>% na.omit() # They share hgnc_symbol
  ### The mapping between uniprot.ID and gene symbol missed some uniprot.ID, so we will get some back from mapping_substrates (protein+STRING_id)
  mapping_substrates = setNames(mapping_substrates$protein, mapping_substrates$STRING_id)
  which_substrates = temp_all_proteins_all$STRING_id%in%names(mapping_substrates)
  temp_all_proteins_all$uniprotswissprot[which_substrates] = 
    mapping_substrates[temp_all_proteins_all$STRING_id[which_substrates]]
  temp_all_proteins_all = na.omit(temp_all_proteins_all)

  # Retrieve PPI interactions from string_db for all substrates and their immediate neighbors that are also in the reference dataset
  interactions_all = string_db$get_interactions(temp_all_proteins_all$STRING_id) %>% unique() 
  # The interactions are STRING_id based, later on we need to convert them to uniprot ID
  
  # On Mac and Windows, multisession runs tasks in separate R sessions (processes)
  future::plan("multisession", workers = 4) 
  future.apply::future_lapply(sigIDs, function(ID) {
    substrates = ptmsigdb3 %>% filter(Term==ID) %>% pull(Phosphosite)
    shared = intersect(substrates, limma_output$PTM.FlankingRegion2)
    shared_uniprotID = mapping_ID[shared] %>% as.character() %>% sub("\n.*", "", .) %>% unique()
    mapped_interest = mapping_substrates[mapping_substrates%in%shared_uniprotID]  #mapping_substrates: uniprotID+STRING ID 
    # Get immediate neighbors for those shared substrates between PTMsigDB and limma output
    neighbors = string_db$get_neighbors(names(mapped_interest)) # Now need STRING IDs here
    # Intersect those immediate neighbors with the proteomics STRING IDs
    proteomics_neighbors = intersect(neighbors, mapped_proteomics$STRING_id)
    # Collect all substrates and their immediate neighbors that are also in the proteomics dataset
    all_proteins = c(names(mapped_interest), proteomics_neighbors) #Uniprot ID
    # Map the STRING IDs of all substrates and their immediate neighbors back to gene names
    mapped_back_all_proteins = mapped_back_all_proteins_all %>% ## mapped_back_all_proteins_all: STRING ID + gene symbol
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
    interactions = rbind(kinase_substrate_ineractions, interactions) %>% unique()  
    
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
    # interactions2 = filter(t1, t1$Phosphosite%in%names(t2)) %>% 
    #   left_join(data.frame(Phosphosite=names(t2),Names=t2,stringsAsFactors = FALSE)) %>%
    interactions2 = data.frame(Phosphosite=names(t2),Names=t2,stringsAsFactors=FALSE) %>%
      mutate(uniprotID = gsub("\n.*", "", Names)) %>% right_join(
        filter(interactions, from%in%ID) %>% 
          mutate(uniprotID = gsub(".*\n", "", to)) %>%
          mutate(uniprotID = gsub("[(]|[)]", "", uniprotID))
      )  %>% left_join(
        data.frame(Phosphosite = names(mapping_regulation),
                   effect = stringr::str_split_fixed(mapping_regulation,";",2)%>%
                      as.data.frame()%>%pull(V1) %>% as.numeric(),
                   sig.stat = stringr::str_split_fixed(mapping_regulation,";",2)%>%
                      as.data.frame()%>%pull(V2) %>% as.numeric(),
                   stringsAsFactors = FALSE),
      )
    # mapping_TO_regualtion = setNames(interactions2$Regulation, interactions2$to)
    # regulations_ID = mapping_TO_regualtion[filter(interactions,from%in%ID) %>% pull(to)]    
    # Set the edge colors for the edges that connect the kinase to their substrates according their regulation (up:orange or down:turqoise) as stated in PTMsigDB
    edge.colors[interactions$from%in%ID] = scales::alpha("purple",0.5)
    #  ifelse(regulations_ID=="u", scales::alpha("purple",0.5), scales::alpha("darkturquoise",0.5))

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
    vertex.label.colors[igraph::V(ppi_network)$name%in%ID] = scales::alpha("black", 1)
    vertex.label.colors[igraph::V(ppi_network)$name%in%network_vertexIDs_4substrates] = 
                scales::alpha("black", 0.8)    
    
    # Set the vertex colors of the kinase substrates based on the effect of the kinase substrates: orange: up-regulated, darkturquoise: down-regulated
    mapping_TO_effect = setNames(interactions2$effect, interactions2$to)
    # There are effects (logFC/t, for sites) share the same name (Uniprot+Gene), we will take the mean of those effects., 
    mapping_TO_effect = tapply(mapping_TO_effect, names(mapping_TO_effect), max)
    mapping_TO_effect = na.omit(mapping_TO_effect) %>% .[!is.na(names(.))]
    # Here effect could be logFC or the t statistics from Limma analysis
    effect = mapping_TO_effect[igraph::V(ppi_network)$name %>% 
                                 .[.%in%network_vertexIDs_4substrates]]
    mapping_TO_sigStat = setNames(interactions2$sig.stat, interactions2$to) %>%
      tapply(., names(.), min)
    mmapping_TO_sigStat = na.omit(mapping_TO_sigStat) %>% .[!is.na(names(.))]
    sig.stat = mapping_TO_sigStat[igraph::V(ppi_network)$name %>% 
                                 .[.%in%network_vertexIDs_4substrates]]

    vertex.colors = rep("snow3", length(igraph::V(ppi_network)$name))
    vertex.colors[igraph::V(ppi_network)$name %in%ID] = scales::alpha("purple", 0.5)
    vertex.colors[igraph::V(ppi_network)$name%in%network_vertexIDs_4substrates] =
        ifelse(effect<0, scales::alpha("darkturquoise",0.5), scales::alpha("orange",0.5))    
    # Set the vertex sizes of the kinase substrates based on their effect*4, set the rest to 1, and set the kinase to 6. 
    vertex.sizes = rep(1, length(igraph::V(ppi_network)$name))
    vertex.sizes[igraph::V(ppi_network)$name%in%network_vertexIDs_4substrates] =
        abs(round(as.numeric(effect))*4)
    vertex.sizes[igraph::V(ppi_network)$name %in%ID] = 10    
    coords = igraph::layout_in_circle(ppi_network)

    # Add signifcance info
    igraph::V(ppi_network)$name[igraph::V(ppi_network)$name%in%network_vertexIDs_4substrates][(abs(effect)>logFCcutoff4limma)&(sig.stat<significance_cutoff4limma)] =
       paste0(igraph::V(ppi_network)$name[igraph::V(ppi_network)$name%in%network_vertexIDs_4substrates][(abs(effect)>logFCcutoff4limma)&(sig.stat<significance_cutoff4limma)], "*")

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






#' @import igraph
#' @import dplyr
#' @param significance_statistic this can be "pvalue", "qvalue" or "p.adjust".
#' @param PTMsigDB_collection_file the path to the PTMsigDB collection file in GMT format.
#' @param PTMSEA_output the output from PTMSEA, which is a data frame with columns: id, Signature.set.description, Signature.set.size, and significance_statistic etc.
#' @param proteomics a character vector of proteomics protein uniprot IDs that are aquired from the same study system as the phosphoproteomics data.
#' @param outdir_ppi the output directory for the PPI network.
#' 
ppiNetwork4substrates_OmniPath = function(limma_output, PTMsubstrates4PTMSEAanalysis, PTMSEA_output,       
                                 PTMSEA_OUTDIR, significance_cutoff=1, 
                                 significance_statistic="fdr.pvalue", 
                                 proteomics, outdir_ppi,
                                 PTMsigDB_collection_file,
                                 omniPath_db_file) {
  dir.create(outdir_ppi, recursive = T, showWarnings = F)
  limma_output$PTM.FlankingRegion2 = paste0(limma_output$PTM.FlankingRegion, "-p")
  limma_output = limma_output %>%
      filter(Phosphosite%in%PTMsubstrates4PTMSEAanalysis)
  limma_output$uniprotID = limma_output$Phosphosite %>% sub(":.*", "", .)
  limma_output$PhosLocation = limma_output$Phosphosite %>% sub(".*:", "", .)
  limma_output = limma_output%>% dplyr::select(uniprotID, PhosLocation, Phosphosite, everything()) 
  
  # get significant PTMsigDB terms
  sigIDs = 
    PTMSEA_output[PTMSEA_output[,significance_statistic,drop=T]<significance_cutoff, ] %>% 
    pull(id) %>% .[!is.na(.)]  
  
  # import PTMsigDB.v2 collection
  ptmsigdb = GSEABase::getGmt(PTMsigDB_collection_file)
  ptmsigdb0 = lapply(ptmsigdb, function(K) {
    data.frame(Term=K@setName, Phosphosite=K@geneIds)
  }) %>% Reduce(rbind, .)
  ptmsigdb2 = as.data.frame(stringr::str_split_fixed(ptmsigdb0$Phosphosite,";",3)) %>%
    setNames(c("uniprotID", "PTM", "Regulation")) %>%
    mutate(PTM = sub("-p", "", PTM))
  ptmsigdb3 = data.frame(Term=ptmsigdb0$Term,
                         uniprotID=ptmsigdb2$uniprotID,
                         PTM=ptmsigdb2$PTM,
                         Phosphosite=paste0(ptmsigdb2$uniprotID, ":", ptmsigdb2$PTM),
                         Regulation=ptmsigdb2$Regulation)  
  
  # Load omniPath database
  omniPath_db = lapply(omniPath_db_file, function(file) {
    read.table(file, sep="\t", header=T) %>%
    filter(is_directed==1) # only keeop highly confident interactions
  }) %>% Reduce(rbind, .) %>% unique()
  #  Expand interactions involving complexes into interactions between all component proteins and the other partner
  omniPath_db$source = gsub("COMPLEX:","",omniPath_db$source) %>%
          stringr::str_split(., "_") 
  omniPath_db = omniPath_db %>% tidyr::unnest_longer(source) %>% unique()

  # Get the omniPath subnetwork for the proteomics dataset
  omniPath_sub_proteomics = omniPath_db %>% 
    filter(source%in%proteomics|target%in%proteomics)
  
  # On Mac and Windows, multisession runs tasks in separate R sessions (processes)
  future::plan("multisession", workers = 4) 
  future.apply::future_lapply(sigIDs, function(ID) {
    substrates = ptmsigdb3 %>% filter(Term==ID) %>% pull(Phosphosite)
    # Get the kinase substrates from PTMsigDB that are shared with the input dataset and the omniPath_sub_proteomics:
    shared_uniprotID = intersect(substrates, limma_output$Phosphosite) %>%
      gsub(":.*", "", .) %>% unique() %>% intersect(.,
      omniPath_sub_proteomics%>%dplyr::select(source,target)%>%
      unlist()%>%unique()) # uniprot IDs
    
    if (length(shared_uniprotID)>0) {
      print(ID)
      # # Get direct neighbors from the omniPath proteomics subnetwork for those shared substrates between PTMsigDB and limma output
      # neighbors = omniPath_sub_proteomics %>% 
      #   filter(source%in%shared_uniprotID|target%in%shared_uniprotID) %>%
      #   dplyr::select(source, target) %>% unlist() %>% unique()
      # # Get the omniPath proteomics subnetwork for the substrates, their direct neighbors, and the direct neighbors of the direct neighbors
      # interactions = omniPath_sub_proteomics %>% 
      #   filter(source%in%neighbors|target%in%neighbors) %>%
      interactions = omniPath_sub_proteomics %>% 
        filter(source%in%shared_uniprotID|target%in%shared_uniprotID) %>%
        unique() 
      names(interactions)[1:2] = c("from","to")
      # Add in the network between PTMsigDB terms (e.g. kinases) and their members (e.g. substrates)
      kinase_substrate_ineractions = data.frame(from=rep(ID, length(shared_uniprotID)), 
                                                to=shared_uniprotID,
                                                source_genesymbol="NA",
                                                target_genesymbol="NA",
                                                is_directed="NA",
                                                is_stimulation="NA",
                                                is_inhibition="NA",
                                                consensus_direction="NA", 
                                                consensus_stimulation="NA",
                                                consensus_inhibition="NA") %>%
                                                unique()
      interactions = rbind(kinase_substrate_ineractions, interactions) 
      #print(head(interactions))
      # Define network edge colors
      edge.colors = rep("snow3", nrow(interactions))
      # Set the colors for the edges that connect the kinase substrates to their proteomics immediate neighbors to orange.
      edge.colors[interactions$is_directed==1]=
        ifelse(interactions[interactions$is_directed==1,]$is_stimulation==1,
          scales::alpha("orange", 0.5),  
            scales::alpha("darkturquoise",0.5))
      ## Set the colors of the edges that connect the kinase and it substrates to purple
      edge.colors[interactions$from%in%ID] = scales::alpha("purple", 0.5) 

      # Set the edge widths for the edges from or to the kinase substrates to 0.75, and all other edges to 0.1
      edge.widths = rep(0.85, nrow(interactions))
      edge.widths[interactions$consensus_direction!=1] = 0.3
      edge.widths[interactions$from==ID] = 0.75 
      
      # Create igraph object
      # library(igraph)
      ppi_network = igraph::graph_from_data_frame(interactions[, 1:2], directed = TRUE) 

      #col_gradients = colorRampPalette(c("orange", "purple"))(length(mapped_interest$STRING_id))
      # Set label color of the kinase to black, set the label colors of kinase substrate to dark grey, and set the rest as grey:
      vertex.label.colors = rep(scales::alpha("black", 0.5), length(igraph::V(ppi_network)))
      vertex.label.colors[igraph::V(ppi_network)$name %in%ID] = scales::alpha("black", 1)
      vertex.label.colors[igraph::V(ppi_network)$name%in%shared_uniprotID] = 
                  scales::alpha("black", 0.8)    
      # Set the vertex colors of the kinase substrates to orange:
      vertex.colors = rep(scales::alpha("snow3", 0.25), length(igraph::V(ppi_network)))
      #****
      vertex.colors[igraph::V(ppi_network)$name%in%ID] = scales::alpha("purple", 0.5)
      vertex.colors[igraph::V(ppi_network)$name%in%shared_uniprotID] = scales::alpha("orange",0.65)   
      # Set the vertex sizes of the kinase substrates based on their effect*2, set the rest to 1, and set the kinase to 6. 

      vertex.sizes = rep(7, length(igraph::V(ppi_network)))
      # vertex.sizes[igraph::V(ppi_network)$name%in%shared_uniprotID] =
      #     abs(round(as.numeric(effect))*2)
      print(shared_uniprotID)
      vertex.sizes[igraph::V(ppi_network)$name%in%shared_uniprotID] = 7
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
          edge.curved=0,
          edge.arrow.size = 0.5))
      dev.off()
      # print(plot(ppi_network,
      #     vertex.shape="circle",
      #     vertex.label.cex=0.7, 
      #     vertex.label.color=vertex.label.colors,
      #     edge.color = edge.colors,
      #     vertex.label.dist = 0,   # distance outward,
      #     edge.width=edge.widths,
      #     layout = coords,
      #     main=ID, asp = 0.7,
      #     vertex.color = vertex.colors,
      #     vertex.frame.color = "white", # vertex.colors,
      #     vertex.size = vertex.sizes,
      #     edge.curved=0, edge.arrow.size = 0.5))
    }
  })
}