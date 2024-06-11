getwd()
library(Seurat)
library(ggplot2)
library(tidyverse)
library(ggsci)
library(rio)
library(scDblFinder)
library(future)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(reshape2)
plan(sequential)
nbrOfWorkers()
plan(multisession, workers=40)
options(future.globals.maxSize = 220000 * 1024^2)
options(scipen = 10)

hl.col <- c("#A6CFE5","#2078B4",'#B2DF8F',"#33A02D","#FA9A99","#E6181E","#FCC06F",'#FF8000',"#CCB2D3",'#6A3D9C',"#FFFF99","#AF5A29","#82C782","#BCAED1","#FEC182","#FFFF9A","#376BB4")
new.col <- c("#E7959B","#DB5E67","#CFD99F",'#C5EC41','#E5CD94',"#FC0C00","#7C4374","#339E55","#000376","#2A82C6","#8C6D36","#CB70C0",
             "#EBB854",'#FC8D37',"#63753A","#6D2811","#DD9AD2","#68AADE","#3B397C","#9D9AE5","#B8CF6E","#949494","#BF4F8F","#844346")
simpson.col <- pal_simpsons("springfield")(16)
cr.col <- c('#DA87B6','#A0CC58','#8F4C9A','#E47C7B','#658DC9','#D6251F','#1E2D5B','#EF7D1C','#4EA74A','#F5D23A','#A15528','#BA1D7C','#0D7C7C')
lancent.col <- pal_lancet("lanonc")(9)
futurama.col <-  pal_futurama("planetexpress")(16)
chicago.col <- pal_uchicago("default")(9)


# data import ----
for(file in c('HDF_F1','HDF_F2','HDF_M1','HDF_M2','ND_F1','ND_F2','ND_M1','ND_M2')){
  seu <- Read10X(paste0('../',file))
  seu_obj <- CreateSeuratObject(seu,
                                min.features = 200,
                                project = file)
  assign(file, seu_obj)
}
dim(ND_M2)
dim(seu_obj)
test <-  CreateSeuratObject(seu,
                            min.features = 200,names.delim="-",
                            project = 'ND_M2')
dim(test)

# double cells discard ----

seu.list <- list(HDF_F1 = HDF_F1,
                 HDF_F2 = HDF_F2,
                 HDF_M1 = HDF_M1,
                 HDF_M2 = HDF_M2,
                 ND_F1 = ND_F1,
                 ND_F2 = ND_F2,
                 ND_M1 = ND_M1,
                 ND_M2 = ND_M2)

seu.list <- lapply(X = seu.list, FUN = function(x){
  x <- as.SingleCellExperiment(x)
  x <- scDblFinder(x) 
  x$doublet_logic <- ifelse(x$scDblFinder.class == "doublet", TRUE, FALSE)
  x <- Seurat::as.Seurat(x)
  x <- subset(x, doublet_logic == 'FALSE')
})
dim(seu.list$HDF_F1)

merge_seu <- merge(x = seu.list$HDF_F1,
                   y = c(seu.list$HDF_F2,
                         seu.list$HDF_M1,
                         seu.list$HDF_M2,
                         seu.list$ND_F1,
                         seu.list$ND_F2,
                         seu.list$ND_M1,
                         seu.list$ND_M2),
                   add.cell.ids = c('HDF_F1','HDF_F2','HDF_M1','HDF_M2','ND_F1','ND_F2','ND_M1','ND_M2'))
dim(merge_seu)

#### QC 

# meta-----
merge_seu$mitoRatio <- PercentageFeatureSet(object = merge_seu, pattern = "^mt-")
merge_seu$mitoRatio <- merge_seu@meta.data$mitoRatio / 100

VlnPlot(merge_seu, features = 'mitoRatio',split.by = 'samples', pt.size = 0) + geom_hline(yintercept = 0.1)
metadata <- merge_seu@meta.data
metadata$cells <- rownames(metadata)
metadata <- metadata %>%
  dplyr::rename(samples = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
metadata$log10GenesPerUMI <- log10(metadata$nGene) / log10(metadata$nUMI)
merge_seu@meta.data <- metadata

metadata %>% 
  ggplot(aes(x=samples, fill=samples)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
metadata %>% 
  ggplot(aes(color=samples, x=nUMI, fill= samples)) + 
  geom_density(alpha = 0.5) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

metadata %>% 
  ggplot(aes(color=samples, x=nGene, fill= samples)) + 
  geom_density(alpha = 0.8) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 850) + scale_color_aaas()

metadata %>% 
  ggplot(aes(x=samples, y=log10(nGene), fill=samples)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "lightgray", high = "red") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~samples)
metadata %>% 
  ggplot(aes(color=samples, x=mitoRatio, fill=samples)) + 
  geom_density(alpha = 0.5) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = samples, fill=samples)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

filtered_seurat <- subset(x = merge_seu, 
                          subset= (nUMI >= 500) &
                            (nGene <= 5000) & 
                            (nGene >= 300) & 
                            (mitoRatio < 0.10))
dim(filtered_seurat)
table(filtered_seurat$samples)
VlnPlot(filtered_seurat, features = c('nUMI','nGene','mitoRatio'), group.by = 'samples',pt.size = 0)


counts <- GetAssayData(object = filtered_seurat, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 30
filtered_counts <- counts[keep_genes, ]

filtered_seurat <- CreateSeuratObject(counts = filtered_counts,meta.data = filtered_seurat@meta.data)
rm(counts,keep_genes,filtered_counts,metadata,nonzero)

save.image()

