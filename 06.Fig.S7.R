# corralation heatmap -----
ann <- data.frame(Group = c(rep(c(rep('MC+PBS',6), rep('MO+FGF21',6), rep('MO+PBS',6)),2)),
                  Gender = c(rep('Female',18),rep('Male',18)))

rownames(ann) <- samp.name
samp.cols <- pal_npg(palette = c("nrc"), alpha = 1)(8)
group_list <- factor(c(paste0(ann$Group,'_',ann$Gender)))
group_list_large <- factor(c(ann$Group))

tpm_cor <- cor(hl.tpm.clean,method = 'spearman',use = 'pairwise.complete.ob')
ha_left.col <- list(Group = c('#0B996F','#7C6DAF','#D6570D'),
                    Gender = c('#f47068','#0e606b'))
names(ha_left.col$Group) <- factor(unique(ann$Group))
names(ha_left.col$Gender) <- factor(unique(ann$Gender))
ComplexHeatmap::pheatmap(tpm_cor,
                         name = "Spearman's coeff",
                         border_color = NA, 
                         clustering_method = 'ward.D',
                         show_colnames = F,
                         # color = scales::alpha(colorRampPalette(colors = c('#00509d','gray80','#e63946'),alpha=T,bias=1)(256),alpha = 1),
                         color = scales::alpha(colorRampPalette(colors = c('#033270','white','#cb1b16'),alpha=T,bias=1)(256),alpha = 1),
                         angle_col = '45',
                         annotation_col = ann,
                         annotation_colors = ha_left.col)

# OXI ------
library(DOSE)
library(GOSemSim)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(dplyr)

get_GO_data <- function(OrgDb, ont, keytype) {
  GO_Env <- get_GO_Env()
  use_cached <- FALSE
  
  if (exists("organism", envir=GO_Env, inherits=FALSE) &&
      exists("keytype", envir=GO_Env, inherits=FALSE)) {
    
    org <- get("organism", envir=GO_Env)
    kt <- get("keytype", envir=GO_Env)
    
    if (org == DOSE:::get_organism(OrgDb) &&
        keytype == kt &&
        exists("goAnno", envir=GO_Env, inherits=FALSE)) {
      ## https://github.com/GuangchuangYu/clusterProfiler/issues/182
      ## && exists("GO2TERM", envir=GO_Env, inherits=FALSE)){
      
      use_cached <- TRUE
    }
  }
  
  if (use_cached) {
    goAnno <- get("goAnno", envir=GO_Env)
  } else {
    OrgDb <- GOSemSim:::load_OrgDb(OrgDb)
    kt <- keytypes(OrgDb)
    if (! keytype %in% kt) {
      stop("keytype is not supported...")
    }
    
    kk <- keys(OrgDb, keytype=keytype)
    goAnno <- suppressMessages(
      AnnotationDbi::select(OrgDb, keys=kk, keytype=keytype,
                            columns=c("GOALL", "ONTOLOGYALL")))
    
    goAnno <- unique(goAnno[!is.na(goAnno$GOALL), ])
    
    assign("goAnno", goAnno, envir=GO_Env)
    assign("keytype", keytype, envir=GO_Env)
    assign("organism", DOSE:::get_organism(OrgDb), envir=GO_Env)
  }
  
  if (ont == "ALL") {
    GO2GENE <- unique(goAnno[, c(2,1)])
  } else {
    GO2GENE <- unique(goAnno[goAnno$ONTOLOGYALL == ont, c(2,1)])
  }
  
  GO_DATA <- DOSE:::build_Anno(GO2GENE, get_GO2TERM_table())
  
  goOnt.df <- goAnno[, c("GOALL", "ONTOLOGYALL")] %>% unique
  goOnt <- goOnt.df[,2]
  names(goOnt) <- goOnt.df[,1]
  assign("GO2ONT", goOnt, envir=GO_DATA)
  return(GO_DATA)
}

get_GO_Env <- function () {
  if (!exists(".GO_clusterProfiler_Env", envir = .GlobalEnv)) {
    pos <- 1
    envir <- as.environment(pos)
    assign(".GO_clusterProfiler_Env", new.env(), envir=envir)
  }
  get(".GO_clusterProfiler_Env", envir = .GlobalEnv)
}

get_GO2TERM_table <- function() {
  GOTERM.df <- get_GOTERM()
  GOTERM.df[, c("go_id", "Term")] %>% unique
}

get_GOTERM <- function() {
  pos <- 1
  envir <- as.environment(pos)
  if (!exists(".GOTERM_Env", envir=envir)) {
    assign(".GOTERM_Env", new.env(), envir)
  }
  GOTERM_Env <- get(".GOTERM_Env", envir = envir)
  if (exists("GOTERM.df", envir = GOTERM_Env)) {
    GOTERM.df <- get("GOTERM.df", envir=GOTERM_Env)
  } else {
    GOTERM.df <- toTable(GOTERM)
    assign("GOTERM.df", GOTERM.df, envir = GOTERM_Env)
  }
  return(GOTERM.df)
}

findGO <- function(pattern, method = "key"){
  
  if(!exists("GO_DATA"))
    load("GO_DATA.RData")
  if(method == "key"){
    pathways = cbind(GO_DATA$PATHID2NAME[grep(pattern, GO_DATA$PATHID2NAME)])
  } else if(method == "gene"){
    pathways = cbind(GO_DATA$PATHID2NAME[GO_DATA$EXTID2PATHID[[pattern]]])
  }
  
  colnames(pathways) = "pathway"
  
  if(length(pathways) == 0){
    cat("No results!\n")
  } else{
    return(pathways)
  }
}

getGO <- function(ID){
  
  if(!exists("GO_DATA"))
    load("GO_DATA.RData")
  allNAME = names(GO_DATA$PATHID2EXTID)
  if(ID %in% allNAME){
    geneSet = GO_DATA$PATHID2EXTID[ID]
    names(geneSet) = GO_DATA$PATHID2NAME[ID]
    return(geneSet)     
  } else{
    cat("No results!\n")
  }
}
GO_DATA <- get_GO_data("org.Mm.eg.db", "ALL", "SYMBOL")

FAT.gene <- getGO("GO:0019395")
d_p(na.omit(hl.tpm.clean[FAT.gene$`fatty acid oxidation`,]),group_list_large,'FAO.TPM')





