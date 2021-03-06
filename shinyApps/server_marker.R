

# Output the chemical selection ####
output$TAS_view <- renderUI({
  
  profile_dat() %...>% {
    
    prof_dat <- .
    
    if("TAS" %in% colnames(prof_dat)){
      max <- floor(max(prof_dat$TAS)/0.1)*0.1    
      sliderInput(inputId="marker_tas", label="TAS:", min=0, max=max, step=0.1, value=c(0, max))
    }  
    
  }
  
})

# update the marker selection if there is no connectivity map ####
output$marker_option <- renderUI({
  
  req(input$portal_id)  
  
  fname <- isolate({ input$portal_id }); 
  
  conn_dat <- projectlist$Connectivity_Variable[which(projectlist$Portal == fname)]
  
  if(conn_dat){
    selectizeInput(
      inputId = "marker",
      label = "Select a marker set:",
      choices = c("Please select an option below" = "", "Genes", "Gene Sets", "CMap Connectivity")
    )
  }else{
    selectizeInput(
      inputId = "marker",
      label = "Select a marker set:",
      choices = c("Please select an option below" = "", "Genes", "Gene Sets")
    )
  }
  
})

output$marker_gene_options <- renderUI({
  
  ##Getting the gene list####
  expression_dat() %...>% {
    eset <- .
    genelist <- sort(rownames(eset))
    
    div(
      selectizeInput(inputId = "marker_gene", label = "Select a gene:", choices = genelist),
      actionButton(inputId = "de_generate", label = "Generate plot", icon=icon("fas fa-arrow-circle-right"), class="mybuttons")
    )
  } 

})
 
# Output the gene set selection ###
output$marker_geneset_options <- renderUI({
  
  req(input$marker_gsname, input$marker_gsmethod)
  
  ## Get gene set list
  dsmap=dsmap(); gsname=input$marker_gsname; gsmethod=input$marker_gsmethod;
  
  ##Getting the gene list####
  gs_enrichment_dat() %...>% {
    eset <- .[[paste0(dsmap[[gsname]], "_", gsmethod)]]
    genesetlist <- sort(rownames(eset))
    selectizeInput(inputId="marker_gs", label="Select a gene set:", choices = genesetlist)
  }
  
})

# Output the gene set selection ###
output$marker_geneset_btn <- renderUI({
  
  req(input$marker_gs)
  
  actionButton(inputId = "gs_generate", label = "Generate plot", icon=icon("fas fa-arrow-circle-right"), class="mybuttons")
  
})

output$marker_conn_options <- renderUI({

  req(input$marker_conn_name)
  
  ## Get gene set list
  conn_name=input$marker_conn_name; 
  
  ##Getting the gene list####
  connectivity_dat() %...>% {
    eset <- .[[conn_name]]
    genesetlist <- sort(rownames(eset))
    selectizeInput(inputId = "marker_conn", label = "Select a gene set:", choices = genesetlist)
  }
  
})

output$marker_conn_btn <- renderUI({
  
  req(input$marker_conn)
  
  actionButton(inputId = "conn_generate", label = "Generate plot", icon=icon("fas fa-arrow-circle-right"), class="mybuttons")
  
})

##Create reactive values####
marker_header <- reactiveVal();
marker_id <- reactiveVal();
create_marker_plot <- reactiveVal();
running <- reactiveVal(FALSE);

# Clear out the plot after users selected different portals
observeEvent(input$portal_id, {

  fname <- isolate({ input$portal_id });

  phenotype <- unlist(strsplit(as.character(projectlist$Exposure_Phenotype[which(projectlist$Portal == fname)]), ",", fixed=TRUE)) %>% trimws()
  exposure_phenotype(phenotype)

  create_marker_plot(NULL)

}, ignoreInit = TRUE)

##Observe event when a button is clicked####
observeEvent(input$de_generate, {

  req(expression_dat(), profile_dat(), annot_var(), exposure_phenotype(), input$marker_gene, input$marker_tas, input$marker_view)

  #Don't do anything if in the middle of a run
  if(running()){ return(NULL) }else{ running(TRUE) }
  
  print("Starting Run")
  
  ##regenerate plots again
  create_marker_plot(NULL)

  ##Get marker header
  marker_header("Mod-Zscores")
  marker_id(input$marker_gene)

  #Get input values
  exposure_phenotype <- exposure_phenotype(); 
  header <- marker_header();
  marker_id <- marker_id();
  tas <- input$marker_tas;
  view <- input$marker_view;
  width <- input$dimension[1];

  ##Create new progress bar
  progress <- AsyncProgress$new(message=paste0("Generating ", ifelse(view %in% "Density", paste0("Density Plot"), view), "..."))

  results <<- promise_all(annot_var=annot_var(), profile_dat=profile_dat(), expression_dat=expression_dat()) %...>% 
    with({
      
      ##Create a empty list to store figures####
      marker_fig <- list(); n=length(exposure_phenotype)+2;
      
      ##Create the overall plot####
      fig1 <- get_marker_plot(
        expression_dat = expression_dat,
        profile_dat = profile_dat,
        annot_var = annot_var,
        marker_id = marker_id,
        col_id = NA,
        fill_name = "Genes",
        header = header,
        tas = tas,
        view = view
      )  %>% ggplotly(width = width)
      
      marker_fig[["Overall"]] <- fig1
      
      progress$inc(1/n)
      
      ##Color palettes for qualitative data####
      col_palette <- c("Set3", "Set1", "Paired", "Accent", "Pastel1", "Dark2", "Pastel2", "Set2")
      
      #Create the exposure phenotype plots####
      for(s in seq_along(exposure_phenotype)){
        #s=1;
        col_id=exposure_phenotype[s]
        variable=profile_dat %>% select(!!exposure_phenotype[s]) %>% distinct()
        col_colors=brewer.pal(nrow(variable), col_palette[s])
        col_names=unique(variable[,exposure_phenotype[s]])
        
        fig <- get_marker_plot(
          expression_dat = expression_dat,
          profile_dat = profile_dat,
          annot_var = annot_var,
          marker_id = marker_id,
          col_id = col_id,
          col_names = col_names,
          col_colors = col_colors,
          fill_name = col_id,
          header = header,
          tas = tas,
          view = view
        )  %>% ggplotly(width = width)
        
        marker_fig[[exposure_phenotype[s]]] <- fig
        
        progress$inc((s+1)/n)
        
      }
      
      ##Create the summary of table output####
      table <- get_de_by_gene_table(
        expression_dat = expression_dat,
        profile_dat = profile_dat,
        annot_var = annot_var,
        marker_id = marker_id,
        header = header,
        tas = tas
      ) %>% data.table.round()
      
      marker_fig[["Table"]] <- table
      
      progress$inc(n/n)
      
      progress$close()
      
      return(marker_fig)
      
    }) %...>% create_marker_plot
  
  ## Show notification on error or user interrupt
  results <- catch(
    results,
    function(e){
      create_marker_plot(NULL)
      print(e$message)
      showNotification("Task Stopped")
    })
  
  ## When done with analysis, remove progress bar
  results <- finally(results, function(){
    print("Done")
    running(FALSE) #declare done with run
  })
  
  print(paste0("Generating ", ifelse(view %in% "Density", paste0("Density Plot"), view), "..."))
  
}, ignoreInit = TRUE)

# ##Observe event when a button is clicked####
observeEvent(input$gs_generate, {

  req(gs_enrichment_dat(), profile_dat(), annot_var(), exposure_phenotype(), input$marker_gs, input$marker_gsname, input$marker_gsmethod, input$marker_tas, input$marker_view)

  #Don't do anything if in the middle of a run
  if(running()){ return(NULL) }else{ running(TRUE) }
  
  print("Starting Run")
  
  ##regenerate plots again
  create_marker_plot(NULL)
  
  ##Get marker header
  marker_header("Gene Set Scores")
  marker_id(input$marker_gs)
  
  #Get input values
  dsmap <- dsmap();
  exposure_phenotype <- exposure_phenotype(); 
  gsname <- input$marker_gsname; 
  gsmethod <- input$marker_gsmethod;  
  marker_id <- marker_id();
  header <- marker_header();
  tas <- input$marker_tas;
  view <- input$marker_view;
  width <- input$dimension[1];

  ##Create new progress bar
  progress <- AsyncProgress$new(message=paste0("Generating ", ifelse(view %in% "Density", paste0("Density Plot"), view), "..."))
  
  results <- promise_all(annot_var=annot_var(), profile_dat=profile_dat(), gs_enrichment_dat=gs_enrichment_dat()) %...>% 
    with({
      
      ##Create a empty list to store figures####
      marker_fig <- list(); n=length(exposure_phenotype)+2;
      
      ##Create the overall plot####
      fig1 <- get_marker_plot(
        expression_dat = gs_enrichment_dat[[paste0(dsmap[[gsname]], "_", gsmethod)]],
        profile_dat = profile_dat,
        annot_var = annot_var,
        marker_id = marker_id,
        col_id = NA,
        fill_name = "Gene_Sets",
        header = header,
        tas = tas,
        view = view
      )  %>% ggplotly(width = width)
      
      marker_fig[["Overall"]] <- fig1
      
      progress$inc(1/n)
      
      ##Color palettes for qualitative data####
      col_palette <- c("Set3", "Set1", "Paired", "Accent", "Pastel1", "Dark2", "Pastel2", "Set2")
      
      #Create the exposure phenotype plots####
      for(s in seq_along(exposure_phenotype)){
        #s=1;
        col_id=exposure_phenotype[s]
        variable=profile_dat %>% select(!!exposure_phenotype[s]) %>% distinct()
        col_colors=brewer.pal(nrow(variable), col_palette[s])
        col_names=unique(variable[,exposure_phenotype[s]])
        
        fig <- get_marker_plot(
          expression_dat = gs_enrichment_dat[[paste0(dsmap[[gsname]], "_", gsmethod)]],
          profile_dat = profile_dat,
          annot_var = annot_var,
          marker_id = marker_id,
          col_id = col_id,
          col_names = col_names,
          col_colors = col_colors,
          fill_name = col_id,
          header = header,
          tas = tas,
          view = view
        )  %>% ggplotly(width = width)
        
        marker_fig[[exposure_phenotype[s]]] <- fig
        
        progress$inc((s+1)/n)
        
      }
      
      ##Create the summary of table output####
      table <- get_de_by_gene_table(
        expression_dat = gs_enrichment_dat[[paste0(dsmap[[gsname]], "_", gsmethod)]],
        profile_dat = profile_dat,
        annot_var = annot_var,
        marker_id = marker_id,
        header = header,
        tas = tas
      ) %>% data.table.round()
      
      marker_fig[["Table"]] <- table
      
      progress$inc(n/n)
      
      progress$close()
      
      return(marker_fig)
      
    }) %...>% create_marker_plot()
  
  ## Show notification on error or user interrupt
  results <- catch(
    results,
    function(e){
      create_marker_plot(NULL)
      print(e$message)
      showNotification("Task Stopped")
    })
  
  ## When done with analysis, remove progress bar
  results <- finally(results, function(){
    print("Done")
    running(FALSE) #declare done with run
  })
  
  print(paste0("Generating ", ifelse(view %in% "Density", paste0("Density Plot"), view), "..."))
  
}, ignoreInit = TRUE)

##Observe event when a button is clicked####
observeEvent(input$conn_generate, {

  req(connectivity_dat(), profile_dat(), annot_var(), exposure_phenotype(), input$marker_conn, input$marker_conn_name, input$marker_tas, input$marker_view)

  ##Disable the generate button
  #shinyjs::disable(id="conn_generate")
  
  ##regenerate plots again
  create_marker_plot(NULL)
  
  ##Get marker header
  marker_header("Connectivity Score (Percentile)")
  marker_id(input$marker_conn)
  
  #Get input values
  exposure_phenotype <- exposure_phenotype(); 
  conn_name <- input$marker_conn_name; 
  marker_id <- marker_id();
  header <- marker_header();
  tas <- input$marker_tas;
  view <- input$marker_view;
  width <- input$dimension[1];

  ##Create new progress bar
  progress <- AsyncProgress$new(message=paste0("Generating ", ifelse(view %in% "Density", paste0("Density Plot"), view), "..."))
  
  results <- promise_all(annot_var=annot_var(), profile_dat=profile_dat(), connectivity_dat=connectivity_dat()) %...>% 
    with({
      
      ##Create a empty list to store figures####
      marker_fig <- list(); n=length(exposure_phenotype)+2;
      
      # throw errors that were signal (if Cancel was clicked)
      interruptor$execInterrupts()
      
      ##Create the overall plot####
      fig1 <- get_marker_plot(
        expression_dat = connectivity_dat[[conn_name]],
        profile_dat = profile_dat,
        annot_var = annot_var,
        marker_id = marker_id,
        col_id = NA,
        fill_name = "CMap_Connectivity",
        header = header,
        tas = tas,
        view = view
      )  %>% ggplotly(width = width)
      
      marker_fig[["Overall"]] <- fig1
      
      progress$inc(1/n)
      
      ##Color palettes for qualitative data####
      col_palette <- c("Set3", "Set1", "Paired", "Accent", "Pastel1", "Dark2", "Pastel2", "Set2")
      
      #Create the exposure phenotype plots####
      for(s in seq_along(exposure_phenotype)){
        #s=1;
        col_id=exposure_phenotype[s]
        variable=profile_dat %>% select(!!exposure_phenotype[s]) %>% distinct()
        col_colors=brewer.pal(nrow(variable), col_palette[s])
        col_names=unique(variable[,exposure_phenotype[s]])
        
        fig <- get_marker_plot(
          expression_dat = connectivity_dat[[conn_name]],
          profile_dat = profile_dat,
          annot_var = annot_var,
          marker_id = marker_id,
          col_id = col_id,
          col_names = col_names,
          col_colors = col_colors,
          fill_name = col_id,
          header = header,
          tas = tas,
          view = view
        )  %>% ggplotly(width = width)
        
        marker_fig[[exposure_phenotype[s]]] <- fig
        
        progress$inc((s+1)/n)
        
      }
      
      ##Create the summary of table output####
      table <- get_de_by_gene_table(
        expression_dat = connectivity_dat[[conn_name]],
        profile_dat = profile_dat,
        annot_var = annot_var,
        marker_id = marker_id,
        header = header,
        tas = tas
      ) %>% data.table.round()
      
      marker_fig[["Table"]] <- table
      
      progress$inc(n/n)
      
      progress$close()
      
      return(marker_fig)

    }) %...>% create_marker_plot()
  
  ## Show notification on error or user interrupt
  results <- catch(
    results,
    function(e){
      create_marker_plot(NULL)
      print(e$message)
      showNotification("Task Stopped")
    })
  
  ## When done with analysis, remove progress bar
  results <- finally(results, function(){
    print("Done")
    running(FALSE) #declare done with run
  })
  
  print(paste0("Generating ", ifelse(view %in% "Density", paste0("Density Plot"), view), "..."))
  
}, ignoreInit = TRUE)

##Styling the plots####
l <- function(title){
  list(
    orientation = 'v',
    title=list(
      text=paste0('<b>', title, ':</b></br>'),
      size=10,
      color="lightgray"
    ),
    font = list(
      family = "sans-serif",
      size = 10
    )
  )
}

##Output the overall plot ####
output$marker_plot <- renderPlotly({
  
  req(create_marker_plot())
  
  marker = isolate({ input$marker })
  
  fig <- create_marker_plot()[["Overall"]]
  fig <- fig %>% layout(height=400, showlegend = TRUE, legend = l(marker), margin = list(b=100), hoverlabel = list(bgcolor="white"))
  fig
  
})

##Output the exposure phenotype plots####
output$exposure_phenotype_plot <- renderPlotly({
  
  req(create_marker_plot())
  
  marker_plot <- create_marker_plot(); 
  exposure_phenotype <- exposure_phenotype(); 
  header <- marker_header(); 
  marker_id <- marker_id();
  
  if(length(exposure_phenotype) == 1){
    
    fig <- marker_plot[[exposure_phenotype]]
    fig <- fig %>% layout(showlegend = TRUE, legend = l(exposure_phenotype), margin = list(b=100), hoverlabel = list(bgcolor="white"))
    
  }else{
    
    p1 <- marker_plot[[exposure_phenotype[1]]] 
    p2 <- marker_plot[[exposure_phenotype[2]]]
    plist <- marker_plot[2:(length(exposure_phenotype)+1)]
    y_val <- seq(1/length(exposure_phenotype), 1, by=1/length(exposure_phenotype)) %>% sort(decreasing=T)

    plot_annotation <- lapply(seq_along(exposure_phenotype), function(i){ 
      
      col_id = exposure_phenotype[i]
      
      list(
        x = 0.5, 
        y = y_val[i], 
        font = list(size = 16), 
        text = paste("Distribution of ", header, " across profiles for ", marker_id, " (by ", col_id, ")\n", sep = ""), 
        xref = "paper", 
        yref = "paper", 
        xanchor = "center", 
        yanchor = "bottom", 
        showarrow = FALSE
      )
      
    })
    
    fig <- subplot(plist, nrows=length(exposure_phenotype), margin=0.05, shareX=T, shareY=T) %>% 
      layout(
        title = "",
        showlegend = TRUE, 
        legend = l("Exposure"),
        height = 400*length(exposure_phenotype),
        hoverlabel = list(bgcolor="white"),
        annotations = plot_annotation
      )
    
  }
  
  fig
  
})

##Output the marker table header####
output$marker_table_header <- renderUI({

  req(create_marker_plot())

  h3(paste0("Table of Profiles Ranked by ", marker_header()))

})

##Output the marker table####
output$marker_table <-  DT::renderDataTable({

  req(create_marker_plot())

  create_marker_plot()[["Table"]]

}, escape = FALSE, extensions = 'Buttons', server = TRUE, rownames=FALSE, selection = "none",
options = list(
  columnDefs = list(list(className = 'dt-left', targets = "_all")),
  deferRender = FALSE,
  paging = TRUE,
  searching = TRUE,
  ordering = TRUE,
  pageLength = 20,
  scrollX = TRUE,
  scrollY = 400,
  scrollCollapse = TRUE,
  dom = 'T<"clear">Blfrtip',
  buttons=c('copy','csv','print'))
)
