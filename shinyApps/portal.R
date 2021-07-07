
## portal page ####
output$portal_page <- renderUI({

  req(input$portal_tab == "portal")
  
  selected_portal <- isolate({ input$selected_portal });
  
  div(
    class="portal-page",
    uiOutput("portal_title"),
    br(),
    fluidRow(
      column(
        width=2,
        selectInput(inputId="portal_id", label=NULL, choices=projectlist$Portal, selected=selected_portal, width="auto")
      ),
      column(
        width=8,
        downloadButton(outputId="download_portal", label="Download Portal Data", style="padding: 5px 10px 5px 10px", class="my-btn btn btn-primary")
      )
    ),
    h4(class="highlight-title", "Description"),
    uiOutput(outputId = "portal_description"),
    br(),
    uiOutput(outputId="portal_main_page") %>% withSpinner(id="loading-page", type=4, color="#0dc5c1")
  )  
  
})

output$portal_title <- renderUI({
  
  req(input$portal_id)
  
  selected_portal <- isolate({ input$portal_id })
  
  h4(class="highlight-title", HTML(paste0("<b style='color: #3182bd;'>Project</b> - ", projectlist$Project[which(projectlist$Portal == selected_portal)], ", <b style='color: #3182bd;'>Cell-Line</b> - ", projectlist$Cell_Line[which(projectlist$Portal == selected_portal)])))

})

output$portal_description <- renderUI({
  
  req(input$portal_id)
  
  selected_portal <- isolate({ input$portal_id })
  
  HTML(paste0("<p>", projectlist$Description[which(projectlist$Portal == selected_portal)], "</p>"))
  
})

# Download the portal data
output$download_portal <- downloadHandler(

  filename = function() {
    paste0(input$portal_id, ".zip")
  },

  content = function(file) {

    withProgress(message = "Downloading: ", value = 0, {

      fs <- c()
      tmpdir <- tempdir()
      file_names <- c()
      
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
      
      # zip the files
      file_names <- c(file_names, "Profile_Annotation")
      
      saveRDS(profile_dat, file.path(tmpdir, "Profile_Annotation.RDS"))

      # increment the progress bar
      incProgress(2/6, detail = "chemical annotation")

      # Read in the chemical data ####
      chemical_dat <- getWorkFileAsObject(
        hiveWorkFileID(dataset@ChemicalAnnotationRDS)
      )

      # zip the files
      file_names <- c(file_names, "Chemical_Annotation")
      
      saveRDS(chemical_dat, file.path(tmpdir, "Chemical_Annotation.RDS"))

      # increment the progress bar
      incProgress(3/6, detail = "expression set")

      # Read in the expression data ####
      expression_dat <- getWorkFileAsObject(
        hiveWorkFileID(dataset@GeneExpressionRDS)
      )

      # zip the files
      file_names <- c(file_names, "Gene_Expression")
      
      saveRDS(expression_dat, file.path(tmpdir, "Gene_Expression.RDS"))

      # increment the progress bar
      incProgress(4/6, detail = "connectivity map")

      # Read in the connectivity data ####
      connectivity_dat <- tryCatch({

        getWorkFileAsObject(
          hiveWorkFileID(dataset@ConnectivityRDS)
        )      

      }, error=function(e){
        
        return(NULL)
        
      })
      
      if(!is.null(connectivity_dat)){
        # zip the files
        file_names <- c(file_names, "Connectivity")
        saveRDS(connectivity_dat, file.path(tmpdir, "Connectivity.RDS"))
      }
      
      # increment the progress bar
      incProgress(5/6, detail = "gene set enrichment")

      # Read in the gs enrichment data ####
      gs_enrichment_dat <- getWorkFileAsObject(
        hiveWorkFileID(dataset@GeneSetEnrichmentRDS)
      )

      # zip the files
      file_names <- c(file_names, "Gene_Set_Enrichment")
      
      saveRDS(gs_enrichment_dat, file.path(tmpdir, "Gene_Set_Enrichment.RDS"))

      # increment the progress bar
      incProgress(6/6, detail = "K2-taxonomer")

      K2summary <- getWorkFileAsObject(
        hiveWorkFileID(dataset@K2TaxonomerResultsRDS)
      )
      
      # zip the files
      file_names <- c(file_names, "K2Taxonomer")

      saveRDS(K2summary, file.path(tmpdir, "K2Taxonomer.RDS"))

      for(names in file_names){
        path <- file.path(tmpdir, paste0(names, ".RDS"))
        fs <- c(fs, path)
      }

      zip(zipfile=file, files=fs, flags = "-r9Xj")

    })
  },

  contentType = "application/zip"
  
)



