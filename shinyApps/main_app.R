
observeEvent(input$portal, {
  
  print(input$portal)
  
})

## Go back to home page when the logo link is clicked on ####
observeEvent(input$main_navbar, {

  if(input$main_navbar %in% c("home", "about", "contact", "sign_in")){
    
    updateQueryString(paste0("?page=", input$main_navbar), mode="push")
    
  }else{
    
    if(input$main_page == "chemical_explorer"){
      if(is.null(input$chem) | input$chem == ""){
        updateQueryString(paste0("?page=", input$portal_id, "&tab=", input$main_page), mode="push")
      }else{
        updateQueryString(paste0("?page=", input$portal_id, "&tab=", input$main_page, "&chemical_id=", input$chem, "&stat=", input$chemical_tab), mode="push")
      }
    }else{
      updateQueryString(paste0("?page=", input$portal_id, "&tab=", input$main_page), mode="push")
    }
    
  }
  
})

observeEvent(input$main_page, {
  
  if(input$main_navbar %in% c("home", "about", "contact", "sign_in")){
    
    updateQueryString(paste0("?page=", input$main_navbar), mode="push")
    
  }else{
    
    if(input$main_page == "chemical_explorer"){
      if(is.null(input$chem) | input$chem == ""){
        updateQueryString(paste0("?page=", input$portal_id, "&tab=", input$main_page), mode="push")
      }else{
        updateQueryString(paste0("?page=", input$portal_id, "&tab=", input$main_page, "&chemical_id=", input$chem, "&stat=", input$chemical_tab), mode="push")
      }
    }else{
      updateQueryString(paste0("?page=", input$portal_id, "&tab=", input$main_page), mode="push")
    }
    
    if(input$main_page == "k2_taxanomer_results"){
      session$sendCustomMessage("ResizeK2Table", "clusteringTable")
    }

  }
  
})

## Go back to home page when the logo link is clicked on ####
observeEvent(input$portal_id, {
  
  if(input$main_navbar %in% c("home", "about", "contact", "sign_in")){
    
    updateQueryString(paste0("?page=", input$main_navbar), mode="push")
    
  }else{
    
    if(is.null(input$chem) | input$chem == ""){
      updateQueryString(paste0("?page=", input$portal_id, "&tab=", input$main_page), mode="push")
    }else{
      updateQueryString(paste0("?page=", input$portal_id, "&tab=", input$main_page, "&chemical_id=", input$chem, "&stat=", input$chemical_tab), mode="push")
    }
    
  }
  
  portal_id <- input$portal_id
  selected_portal(portal_id)
  
}, ignoreInit = FALSE)

output$download_portal <- downloadHandler(
  
  filename = function() {
    paste0(input$portal_id, "-Dataset.zip")
  },
  
  content = function(file) {
    
    withProgress(message = "Downloading: ", value = 0, {
      
      fs <- c()
      tmpdir <- tempdir()
      
      datasets <- listEntities("PortalDataset", portal=input$portal_id)
      
      # Sort by timestamp and extract most recent dataset to convenience object
      datasets <- datasets[order(sapply(datasets, slot, "timestamp"))]
      dataset <- datasets[[length(datasets)]]
      
      # increment the progress bar
      incProgress(1/6, detail = "profile annotation")
      
      # Read in the profile data ####
      profile_dat <- getWorkFileAsObject(
        hiveWorkFileID(dataset@ProfileAnnotationRDS)
      )
      
      saveRDS(profile_dat, file.path(tmpdir, "Profile_Annotation.RDS"))
      
      # increment the progress bar
      incProgress(2/6, detail = "chemical annotation")          
      
      # Read in the chemical data ####
      chemical_dat <- getWorkFileAsObject(
        hiveWorkFileID(dataset@ChemicalAnnotationRDS)
      )
      
      saveRDS(chemical_dat, file.path(tmpdir, "Chemical_Annotation.RDS"))
      
      # increment the progress bar
      incProgress(3/6, detail = "expression set")
      
      # Read in the expression data ####
      expression_dat <- getWorkFileAsObject(
        hiveWorkFileID(dataset@GeneExpressionRDS)
      )
      
      saveRDS(expression_dat, file.path(tmpdir, "Gene_Expression.RDS"))
      
      # increment the progress bar
      incProgress(4/6, detail = "connectivity map")
      
      # Read in the connectivity data ####
      connectivity_dat <- getWorkFileAsObject(
        hiveWorkFileID(dataset@ConnectivityRDS)
      )
      
      saveRDS(connectivity_dat, file.path(tmpdir, "Connectivity.RDS"))
      
      # increment the progress bar
      incProgress(5/6, detail = "gene set enrichment")
      
      # Read in the gs enrichment data ####
      gs_enrichment_dat <- getWorkFileAsObject(
        hiveWorkFileID(dataset@GeneSetEnrichmentRDS)
      )
      
      saveRDS(gs_enrichment_dat, file.path(tmpdir, "Gene_Set_Enrichment.RDS"))
      
      # increment the progress bar
      incProgress(6/6, detail = "K2-taxonomer") 
      
      K2summary <- getWorkFileAsObject(
        hiveWorkFileID(dataset@K2TaxonomerResultsRDS)
      )
      
      saveRDS(K2summary, file.path(tmpdir, "K2Taxonomer.RDS"))
      
      # zip the files
      file_names <- c("Profile_Annotation", "Chemical_Annotation", "Gene_Expression", "Gene_Set_Enrichment", "Connectivity", "K2Taxonomer") 
      
      for(names in file_names){
        path <- file.path(tmpdir, paste0(names, ".RDS"))
        fs <- c(fs, path)
      } 
      
      zip(zipfile=file, files=fs, flags = "-r9Xj")
      
    })
  },
  
  contentType = "application/zip"
)
