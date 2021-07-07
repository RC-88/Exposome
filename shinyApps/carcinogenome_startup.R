
# Defaults for gene expression subtab ####
defaults <- list(
  landmark_de = FALSE, 
  summarizefunc_de = "median", 
  filterbyinput_de = c("score", "number"),
  range_de = c(-2, 2), 
  numberthresleft_de = 10, 
  numberthresright_de = 10
)

# Create a list of gene enrichment methods ####
dsmap_method <- list("gsva"="gsva")

# The gsva helptext method####
helptext_method <- paste(
  "gsva, ssgea, zscore: from R Bioconductor package GSVA",
  sep = "<br>"
)

# Create classes for connectivity mapping #####
connmap <- list("Perturbagen Classes" = "pcl", "Perturbagens" = "pert")

# Link gene expression to genecards.org####
get_genecard_link <- function(genesymbol){
  sprintf('<a href="http://www.genecards.org/cgi-bin/carddisp.pl?gene=%s&keywords=%s" target="_blank" style="text-decoration:none;">%s</a>', genesymbol, genesymbol, genesymbol)
}

# Get the gene set enrichment hyperlink####
get_geneset_link <- function(geneset){
  label <- gsub("_", " ", geneset)
  geneset <- gsub(" ", "_", geneset)
  sprintf('<a href="http://software.broadinstitute.org/gsea/msigdb/cards/%s" target="_blank" style="text-decoration:none;">%s</a>', geneset, label)
}

# Get the 25th percentile ####
Q1 <- function(x){ 
  quantile(x, 0.25, na.rm = T)
}

# Get the 75th percentile ####
Q3 <- function(x){
  quantile(x, 0.75, na.rm = T)
}

# Get the chemical name####
get_chem_description <- function(chemical_dat, chem, chemical_id=FALSE){
  
  pos <- lapply(c("Chemical_Name", "BUID", "CAS"), function(x){
    w <- which(chemical_dat[,x] %in% chem)
    if(length(w) > 0){
      return(w)
    }
  }) %>% unlist()
  
  if(chemical_id){
    return(unique(chemical_dat$Chemical_Id[pos]))
  }else{
    return(pos)
  }
  
}

# Function to round the data table values ####
data.table.round <- function(dt, digits = 3){
  
  cols <- sapply(colnames(dt), function(i) is.numeric(dt[,i]))
  cols <- names(which(cols))
  
  for(i in cols)
    dt[,i] <- round(dt[,i], digits)
  
  dt <- data.table(dt)
  
  return(dt)
  
}

# ### Read in the project list ####
# fname = "TGGates";
# profile_dat <- readRDS(paste0("data/", fname, "/Profile_Annotation.RDS"))
# chemical_dat <- readRDS(paste0("data/", fname, "/Chemical_Annotation.RDS"))
# expression_dat <- readRDS(paste0("data/", fname, "/Expression_Set.RDS"))
# annot_var = "Chemical_Id"
# chem="acetamidofluorene"
# header = "ModZScore"
# summarize.func = "median"
# landmark = FALSE
# do.nmarkers = FALSE
# nmarkers = c(100, 100)
# do.scorecutoff = TRUE
# scorecutoff = c(-2, 2)

##Get gene expression####
get_de <- function(
  chem, annot_var,
  profile_dat, chemical_dat, expression_dat,
  header = "ModZScore",
  summarize.func = "mean",
  landmark = TRUE, 
  do.nmarkers = TRUE, nmarkers = c(100, 100),
  do.scorecutoff = TRUE, scorecutoff = c(-2, 2)){

  #getting chemical BUID####
  chemical_id <- get_chem_description(chemical_dat=chemical_dat, chem=chem, chemical_id=TRUE)

  #getting signature id ####
  profile_dat <- profile_dat[match(colnames(expression_dat), profile_dat[, annot_var]),]
  exposure <- sort(unique(profile_dat$unique_ID_by_chem[which(profile_dat$Chemical_Id %in% chemical_id)]))
  
  #getting expression data####
  eset <- expression_dat
  
  #getting feature data####
  fdat <- fData(eset)
  
  if(nrow(fdat) > 0){
    if(!"Gene" %in% colnames(fdat)){
      fdat <- data.frame(Gene=rownames(eset), fdat)
    }
  }else{
    fdat <- data.frame(Gene=rownames(eset))
  }
  
  #check it landmark is selected
  if(landmark) {
    if("Landmark_Gene" %in% colnames(fdat)){
      eset <- eset[which(toupper(fdat[,"Landmark_Gene"]) %in% "YES"),] 
    }
  }

  #summarise the expression set
  eSet <- exprs(eset); mat <- matrix(NA, nrow=nrow(eSet), ncol=length(exposure), byrow=T, dimnames=list(rownames(eSet), paste0(header, " ", exposure)))
  
  for(i in seq_along(exposure)){
    #i=1;
    sig <- unique(profile_dat[which(profile_dat$unique_ID_by_chem %in% exposure[i] & profile_dat$Chemical_Id %in% chemical_id), annot_var])
    if(length(sig) > 1){
      mat[,i] <- rowMeans(eSet[,as.character(sig)], na.rm=T)
    }else{
      mat[,i] <- eSet[,as.character(sig)]
    }
  }
  
  x <- apply(mat, 1, match.fun(summarize.func))
  x <- as.numeric(x)
  n <- length(x)
  
  if(do.nmarkers){
    
    ind_down <- which(x < 0)
    ind_up <- which(x > 0) 
    
    ord_down <- data.frame(
      pos=ind_down,
      x=x[ind_down]
    ) %>% arrange(desc(x))
    
    ord_up <- data.frame(
      pos=ind_up,
      x=x[ind_up]
    ) %>% arrange(desc(x))
    
    n1 = ifelse(nmarkers[1] > nrow(ord_down), nrow(ord_down), nmarkers[1]); 
    n2 = ifelse(nmarkers[2] > nrow(ord_up), nrow(ord_up), nmarkers[2]);
    
    if(n1 == 0 & n2 == 0){
      x.ind.nmarkers <- 1:n
    }else if(n2 == 0){
      x.ind.nmarkers <- ord_down$pos[1:n1]
    }else{
      x.ind.nmarkers <- c(ord_down$pos[1:n1], ord_up$pos[1:n2])
    }
    
  } else { 
    
    x.ind.nmarkers <- 1:n 
    
  }
  
  if(do.scorecutoff){
    #TODO: rank by score here too
    x.ind.scorecutoff <- which(x >= scorecutoff[1] & x <= scorecutoff[2])
  }else {
    x.ind.scorecutoff <- 1:n
  }
  
  inds <- intersect(x.ind.nmarkers, x.ind.scorecutoff)
  inds <- inds[order(x[inds], decreasing = TRUE)]
  
  #get expression results
  res.ind <- inds; res.scores <- x[inds]
  
  #determine whether gene is up or down regulated
  direction <- sapply(res.scores, function(i){
    if(i > 0) return("Up") else return("Down")
  })

  #creata the summary table 
  table <- cbind(fdat[res.ind,, drop = FALSE], Direction = direction, SummaryScore=res.scores, mat[res.ind,, drop = FALSE])
  colnames(table)[which(colnames(table) %in% "SummaryScore")] <- "Summary Score"
  
  #return hyperlink from genecard.org
  table$Gene <- sapply(as.character(table$Gene), get_genecard_link)

  return(table)

}

# ##Getting the gene set scores for diffrent gsva methods
# fname = "TGGates"; gs_enrichment_version=7;
# gsname="Hallmark"; gsmethod="gsva";
# dsmap <- list(
#   Hallmark=paste0("gsscores_h.all.v", gs_enrichment_version, ".0"),
#   C2=paste0("gsscores_c2.cp.reactome.v", gs_enrichment_version, ".0"),
#   NURSA=paste0("gsscores_nursa_consensome_Cbyfdrvalue_0.01")
# )
# ## Read in the project list ####
# profile_dat <- readRDS(paste0("data/", fname, "/Profile_Annotation.RDS"))
# chemical_dat <- readRDS(paste0("data/", fname, "/Chemical_Annotation.RDS"))
# assign(paste0("gs_enrichment_dat"), readRDS(paste0("data/", fname, "/GS_Enrichment.RDS")))
# expression_dat <- gs_enrichment_dat[[paste0(dsmap[[gsname]], "_", gsmethod)]]
# annot_var = "Chemical_Id"
# chem = "acetamidofluorene"
# header = "GS Score"
# summarize.func = "median"

##Summarize gene set enrichment#####
get_gsenrichment <- function(
  chem, annot_var,
  profile_dat, chemical_dat, expression_dat,
  gsname = "Hallmark",
  header = "GS Score",
  summarize.func = "mean"){
  
  #getting chemical BUID####
  chemical_id <- get_chem_description(chemical_dat=chemical_dat, chem=chem, chemical_id=TRUE)
  
  #getting signature id ####
  profile_dat <- profile_dat[match(colnames(expression_dat), profile_dat[, annot_var]),]
  exposure <- sort(unique(profile_dat$unique_ID_by_chem[which(profile_dat$Chemical_Id %in% chemical_id)]))
  
  #getting expression data####
  eset <- expression_dat
  
  #getting feature data####
  fdat <- fData(eset)
  
  if(nrow(fdat) > 0){
    if(!"Geneset" %in% colnames(fdat)){
      fdat <- data.frame(Geneset=rownames(eset), fdat)
    }
  }else{
    fdat <- data.frame(Geneset=rownames(eset))
  }
  
  #summarise the expression set
  eSet <- exprs(eset); mat <- matrix(NA, nrow=nrow(eSet), ncol=length(exposure), byrow=T, dimnames=list(rownames(eSet), paste0(header, " ", exposure)))
  
  for(i in seq_along(exposure)){
    #i=1;
    sig <- unique(profile_dat[which(profile_dat$unique_ID_by_chem %in% exposure[i] & profile_dat$Chemical_Id %in% chemical_id), annot_var])
    if(length(sig) > 1){
      mat[,i] <- rowMeans(eSet[,as.character(sig)], na.rm=T)
    }else{
      mat[,i] <- eSet[,as.character(sig)]
    }
  }
  
  res <- apply(mat, 1, match.fun(summarize.func))
  res <- as.numeric(res)
  
  res <- cbind(fdat, score = res, mat)
  res <- res[order(res$score, decreasing = TRUE),, drop = FALSE]
  colnames(res)[which(colnames(res) %in% "score")] <- "Summary Score"
  
  #return hyperlink to MSigDB genesets
  if(gsname %in% c("Hallmark", "C2")) res$Geneset <- sapply(as.character(res$Geneset), get_geneset_link)
  
  return(res)
  
}


##Summarize connectivity map#####
get_connectivity <- function(
  chem, annot_var,
  profile_dat, chemical_dat, expression_dat,
  header = "Connectivity Score",
  summarize.func = "mean"){
  
  #getting chemical BUID####
  chemical_id <- get_chem_description(chemical_dat=chemical_dat, chem=chem, chemical_id=TRUE)
  
  #getting signature id ####
  profile_dat <- profile_dat[match(colnames(expression_dat), profile_dat[, annot_var]),]
  exposure <- sort(unique(profile_dat$unique_ID_by_chem[which(profile_dat$Chemical_Id %in% chemical_id)]))
  
  #getting expression data####
  eset <- expression_dat
  
  #getting feature data####
  fdat <- fData(eset)
  
  if(nrow(fdat) > 0){
    if(!"Connectivity_Id" %in% colnames(fdat)){
      fdat <- data.frame(Connectivity_Id=rownames(eset), fdat)
    }
  }else{
    fdat <- data.frame(Connectivity_Id=rownames(eset))
  }
  
  #summarise the expression set
  eSet <- exprs(eset); mat <- matrix(NA, nrow=nrow(eSet), ncol=length(exposure), byrow=T, dimnames=list(rownames(eSet), paste0(header, " ", exposure)))
  
  for(i in seq_along(exposure)){
    #i=1;
    sig <- unique(profile_dat[which(profile_dat$unique_ID_by_chem %in% exposure[i] & profile_dat$Chemical_Id %in% chemical_id), annot_var])
    if(length(sig) > 1){
      mat[,i] <- rowMeans(eSet[,as.character(sig)], na.rm=T)
    }else{
      mat[,i] <- eSet[,as.character(sig)]
    }
  }
  
  res <- apply(mat, 1, match.fun(summarize.func))
  res <- as.numeric(res)
  
  res <- cbind(fdat, score = res, mat)
  res <- res[order(res$score, decreasing = TRUE),, drop = FALSE]
  colnames(res)[colnames(res) %in% "score"]<- "Summary Score"
  
  return(res)
  
}

#plot wrapper###
plot_wrapper <- function(df, fill_name, view = "Density", ...){
  
  if(view %in% "Density") {
    
    res <- ggplot(data=df, aes_string(x = "x", fill = fill_name)) + 
      geom_density(position = "identity", alpha = 0.5)
    
  }else if(view %in% "Boxplot"){
    
    res <- ggplot(data=df, aes_string(x = fill_name, y = "x", fill = fill_name)) + 
      geom_boxplot(
        position = "identity", 
        width = 0.2,
        alpha = 0.5, 
        outlier.fill = NULL,
        outlier.alpha = NULL
      )
    
  }
  
  return(res)
  
}

# library(tidyverse)
# library(RColorBrewer)
# 
# fname="TGGates"; gs_enrichment_version=7;
# 
# projectlist <- read_csv(paste0("data/Project_List.csv"))
# exposure <- unlist(strsplit(as.character(projectlist$Exposure[which(projectlist$Portal == fname)]), ",", fixed=TRUE)) %>% trimws()
# 
# ##Getting the gene set scores for diffrent gsva methods
# dsmap <- list(
#   Hallmark=paste0("gsscores_h.all.v", gs_enrichment_version, ".0"),
#   C2=paste0("gsscores_c2.cp.reactome.v", gs_enrichment_version, ".0"),
#   NURSA=paste0("gsscores_nursa_consensome_Cbyfdrvalue_0.01")
# )
# marker_gsname="Hallmark"; marker_gsmethod="gsva";
# ## Read in the project list ####
# profile_dat <- readRDS(paste0("data/", fname, "/Profile_Annotation.RDS"))
# #expression_dat <- readRDS(paste0("data/", fname, "/Expression_Set.RDS"))
# assign(paste0("gs_enrichment_dat"), readRDS(paste0("data/", fname, "/GS_Enrichment.RDS")))
# assign(paste0("connectivity_dat"), readRDS(paste0("data/", fname, "/Connectivity.RDS")))
# expression_dat <- gs_enrichment_dat[[paste0(dsmap[[marker_gsname]], "_", marker_gsmethod)]]
# annot_var = "Sig_Id"
# marker_id = "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
# # col_id = "Carcinogenicity"
# # col_colors = c("grey", "green", "orange")
# # col_names = c("N/A", "-", "+")
# header = "Mod Z-scores"
# tas = NULL;
# view = "Density";
# 
# ##Color palettes for qualitative data####
# col_palette <- c("Set3", "Set1", "Paired", "Accent", "Pastel1", "Dark2", "Pastel2", "Set2")
# 
# for(s in seq_along(exposure)){
#   #s=1;
#   col_id=exposure[s]
#   variable=profile_dat %>% select(!!exposure[s]) %>% distinct()
#   col_colors=brewer.pal(nrow(variable), col_palette[s])
#   col_names=unique(variable[,exposure[s]])
# }

##Function to create density and boxplot for marker explorer#####
get_marker_plot <- function(
  expression_dat,
  profile_dat,
  annot_var = "Sig_Id",
  marker_id,
  col_id = "Carcinogenicity",
  col_names = c("N/A", "-", "+"), 
  col_colors = c("grey", "green", "orange"),
  fill_name = "Genes",
  header = "Mod Z-scores",
  tas = NULL,
  view = "Density"){
  
  if(is.null(tas)){
    profile_ann <- profile_dat 
  }else{
    profile_ann <- profile_dat %>% filter(TAS >= (min(tas)-0.01) & TAS <= (max(tas)+0.01))
  }
  
  if(is.na(col_id)){
    profile_ann <- profile_ann %>% transmute(Id=(!!!syms(annot_var)), carcinogencity=marker_id) %>% distinct(Id, .keep_all=TRUE)
  }else{
    profile_ann <- profile_ann %>% transmute(Id=(!!!syms(annot_var)), carcinogencity=(!!!syms(col_id))) %>% distinct(Id, .keep_all=TRUE)
  }
  
  rowid <- which(rownames(expression_dat) %in% marker_id)
  eset <- exprs(expression_dat)[rowid,]
  query <- profile_ann %>% mutate(Id=as.character(Id)) %>% left_join(data.frame(Id=as.character(names(eset)), x=as.numeric(eset))) %>% select(x, carcinogencity)
  
  if(is.na(col_id)){
    
    background <- data.frame(x=as.numeric(exprs(expression_dat)), carcinogencity="All")
    df <- query %>% rbind(background)
    df$carcinogencity <- factor(df$carcinogencity, levels = c("All", marker_id), ordered=is.ordered(c("All", marker_id)))
    cols_match <- c("grey", "red")
    names(cols_match) <- c("All", marker_id)
    p.title <- paste("Distribution of ", header, " Across Profiles for ", marker_id, "\n", sep = "")

    df <- df %>% rename(!!paste0(fill_name) := carcinogencity)
    
    p <- plot_wrapper(df=df, fill_name=fill_name, view=view) +
      scale_fill_manual(
        name = NULL, 
        values = cols_match,
        breaks = names(cols_match),
        labels = names(cols_match)
      ) +
      xlab(header) + 
      ylab("Density") + 
      ggtitle(p.title) +
      theme_bw() +
      theme(
        plot.margin = margin(5, 5, 0, 0),
        plot.title = element_text(hjust = 0.5)
      ) +
      theme(
        legend.direction = "vertical", 
        legend.box = "vertical"
      )
    
  }else{
    
    df <- query
    cols_match <- col_colors
    names(cols_match) <- col_names
    df$carcinogencity <- factor(df$carcinogencity, levels = col_names, ordered=is.ordered(col_names))
    p.title <- paste("Distribution of ", header, " across profiles for ", marker_id, " (by ", col_id, ")\n", sep = "")
    
    df <- df %>% rename(!!paste0(fill_name) := carcinogencity)
    
    p <- plot_wrapper(df=df, fill_name=fill_name, view=view) +
      scale_fill_manual(
        name = NULL, 
        values = cols_match,
        breaks = names(cols_match),
        labels = names(cols_match)
      ) +
      xlab(header) + 
      ylab("Density") + 
      ggtitle(p.title) +
      theme_bw() +
      theme(
        plot.margin = margin(5, 5, 0, 0),
        plot.title = element_text(hjust = 0.5)
      ) +
      theme(
        legend.direction = "vertical", 
        legend.box = "vertical"
      )
    
  }
  
  return(p)
  
}


##Get the table profile ranked by mod Z-scores for marker explorer####
get_de_by_gene_table <- function(
  marker_id,
  expression_dat,
  profile_dat,
  annot_var = "Sig_Id",
  header = "Mod Z-scores",
  tas = NULL){
  
  profile_dat <- profile_dat[which(profile_dat[, annot_var] %in% colnames(expression_dat)),]
  
  if(is.null(tas)){
    sig <- unique(profile_dat[, annot_var])
  }else{
    sig <- unique(profile_dat[which(profile_dat$TAS >= (min(tas)-0.01) & profile_dat$TAS <= (max(tas)+0.01)), annot_var])
  }
  
  eset <- expression_dat[,as.character(sig)]
  rowid <- which(rownames(eset) %in% marker_id)
  x <- as.numeric(exprs(eset)[rowid,])
  
  pdat <- profile_dat[which(!duplicated(profile_dat[which(profile_dat[,annot_var] %in% as.character(sig)), annot_var])),]

  df <- cbind(value = x, pdat)
  df <- df[order(x, decreasing = TRUE),]
  
  colnames(df)[colnames(df) %in% "value"] <- header
  
  return(df)
  
}


