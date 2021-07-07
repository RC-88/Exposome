

output$portal_main_page <- renderUI({
  
  req(input$portal_id, input$selected_subtab, input$selected_chemical_tab)
  
  fname <- isolate({ input$portal_id }); 
  subtab <- isolate({ input$selected_subtab });
  chemical_tab <- isolate({ input$selected_chemical_tab });
  
  #get input options####
  geneset_collection <- projectlist$GS_Collection[which(projectlist$Portal == fname)]
  gs_collection(geneset_collection)
  
  geneset_collection_link <- projectlist$GS_Collection_Link[which(projectlist$Portal == fname)]
  gs_collection_link(geneset_collection_link)
  
  phenotype_test <- unlist(strsplit(as.character(projectlist$Exposure_Phenotype_Test[which(projectlist$Portal == fname)]), ",", fixed=TRUE)) %>% trimws()
  exposure_phenotype_test(phenotype_test)
  
  phenotype <- unlist(strsplit(as.character(projectlist$Exposure_Phenotype[which(projectlist$Portal == fname)]), ",", fixed=TRUE)) %>% trimws()
  exposure_phenotype(phenotype)
  
  #get landmark option
  landmark(projectlist$Landmark_Gene[which(projectlist$Portal == fname)])

  ##Getting the helpext for different gene set enrichment
  if(geneset_collection %in% "Default"){
    
    ##the default gs enrichment version####
    gs_enrichment_version <- 7
    
    ##Getting the gene set scores for diffrent gsva methods
    gs_enrichment_list <- list(
      Hallmark=paste0("gsscores_h.all.v", gs_enrichment_version, ".0"),
      C2=paste0("gsscores_c2.cp.reactome.v", gs_enrichment_version, ".0"),
      NURSA=paste0("gsscores_nursa_consensome_Cbyfdrvalue_0.01")
    )
    
    dsmap(gs_enrichment_list)
    
    ##Getting the helptext####
    helptext_geneset(
      paste0(
        "<a href=\"https://www.gsea-msigdb.org/gsea/msigdb\">MSigDB Hallmark Pathways (v", gs_enrichment_version, ")</a><br>",
        "<a href=\"https://www.gsea-msigdb.org/gsea/msigdb\">MSigDB C2 reactome Pathways (v", gs_enrichment_version, ")</a><br>",
        "<a href=\"https://signalingpathways.org\">NURSA: Nuclear Receptor Signaling Atlas, consensome data for human</a><br>"
      )
    )
    
  }else{
    
    ##Getting the gene set scores for diffrent gsva methods
    gs_enrichment_list <- paste0("gsscores_", geneset_collection)
    names(gs_enrichment_list) <- paste0("gsscores_", geneset_collection)
    
    dsmap(gs_enrichment_list)
    
    ##Getting the helptext####
    helptext_geneset(paste0("<a href=\"", geneset_collection_link, "\">", geneset_collection, "</a><br>"))
    
  }
  
  tabsetPanel(
    id="main_page", type="pills", selected = subtab,
    
    ###Annotation#####
    tabPanel(
      title="Annotation", value="annotation",
      #"Hello World!"
      source("ui_annotation.R", local=TRUE)$value
    ),
    
    ###Chemical Explorer#####
    tabPanel(
      title = "Chemical Explorer", value="chemical_explorer",
      #"Hello World!"
      source("ui_chemical.R", local=TRUE)$value
    ),
    
    ###Marker Explorer####
    tabPanel(
      title = "Marker Explorer", value="marker_explorer",
      #"Hello World!"
      source("ui_marker.R", local=TRUE)$value
    ),
    
    ###Heatmap Explorer####
    tabPanel(
      title = "Heatmap Explorer", value="heatmap_explorer",
      #"Hello World!"
      source("ui_heatmap.R", local=TRUE)$value
    ),
    
    ###Taxonomic Clustering####
    tabPanel(
      title = "K2 Taxonomer Results", value="k2_taxonomer_results",
      #"Hello World!"
      source("ui_taxonomic_clustering.R", local=TRUE)$value
    ),
    
    tabPanel(
      title = "Compare Multiple", value="compare_multiple",
      #"Hello World!"
      source("ui_compare_multiple.R", local=TRUE)$value
    )
  )

})

observeEvent(input$portal_id, {
  fname <- isolate({ input$portal_id })
  session$sendCustomMessage("SelectedPortal", fname)
})

observeEvent(input$main_page, {
  subtab <- isolate({ input$main_page })
  session$sendCustomMessage("ChangedSubtab", subtab)
})

observeEvent(input$chemical_tab, {
  chemical_tab <- isolate({ input$chemical_tab })
  session$sendCustomMessage("ChangedChemicalTab", chemical_tab)
})

observeEvent(input$chem, {
  req(input$chem)
  chemical_id <- isolate({ input$chem })
  session$sendCustomMessage("ChangedChemicalId", chemical_id)
})

## Observe the portal tab selected ####
observeEvent({
  input$portal_tab
}, {

  if(toupper(input$portal_tab) %in% toupper(c("home", "overview", "sign_in"))){

    updateQueryString(paste0("?page=", tolower(input$portal_tab)), mode="push")

  }else{
    
    if(is.null(input$main_page)){
      
      if(fname == "chemical_explorer"){
        updateQueryString(paste0("?page=", fname, "&tab=", subtab, "&chemical_id=", chemical_id, "&stat=", chemical_tab), mode="push")
      }else{
        updateQueryString(paste0("?page=", fname, "&tab=", subtab), mode="push")
      }
      
    }else{
      
      if(input$main_page == "chemical_explorer"){
        updateQueryString(paste0("?page=", input$portal_id, "&tab=", input$main_page, "&chemical_id=", input$chem, "&stat=", input$chemical_tab), mode="push")
      }else{
        updateQueryString(paste0("?page=", input$portal_id, "&tab=", input$main_page), mode="push")
      }
      
    }

  }

})

## Observe the portal id selected ####
observeEvent({
  input$portal_id
}, {
  
  req(input$main_page)
  
  if(input$main_page == "chemical_explorer"){
    updateQueryString(paste0("?page=", input$portal_id, "&tab=", input$main_page, "&chemical_id=", input$chem, "&stat=", input$chemical_tab), mode="push")
  }else{
    updateQueryString(paste0("?page=", input$portal_id, "&tab=", input$main_page), mode="push")
  }
  
})

## Observe the portal tab selected ####
observeEvent({
  input$main_page
}, {

  if(input$main_page == "chemical_explorer"){
    updateQueryString(paste0("?page=", input$portal_id, "&tab=", input$main_page, "&chemical_id=", input$chem, "&stat=", input$chemical_tab), mode="push")
  }else{
    updateQueryString(paste0("?page=", input$portal_id, "&tab=", input$main_page), mode="push")
  }

  if(input$main_page == "k2_taxonomer_results"){
    session$sendCustomMessage("ResizeK2Table", "clusteringTable")
  }
  
})

## Observe the portal tab selected ####
observeEvent({
  input$chemical_tab
}, {
  
  if(input$main_page == "chemical_explorer"){
    updateQueryString(paste0("?page=", input$portal_id, "&tab=", input$main_page, "&chemical_id=", input$chem, "&stat=", input$chemical_tab), mode="push")
  }else{
    updateQueryString(paste0("?page=", input$portal_id, "&tab=", input$main_page), mode="push")
  }

})

## Observe the portal tab selected ####
observeEvent({
  input$chem
}, {

  if(toupper(input$portal_tab) %in% toupper(c("home", "overview", "contact", "sign_in"))){
    updateQueryString(paste0("?page=", tolower(input$portal_tab)), mode="push")
  }else{
    if(input$main_page == "chemical_explorer"){
      updateQueryString(paste0("?page=", input$portal_id, "&tab=", input$main_page, "&chemical_id=", input$chem, "&stat=", input$chemical_tab), mode="push")
    }else{
      updateQueryString(paste0("?page=", input$portal_id, "&tab=", input$main_page), mode="push")
    }
  }

})
