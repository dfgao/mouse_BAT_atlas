plan(multisession, workers=1)
seurat_phase <- NormalizeData(filtered_seurat)
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m_genes$SYMBOL, 
                                 s.features = s_genes$SYMBOL)

# cca -----
split_seurat <- SplitObject(seurat_phase,split.by = "samples")
split_seurat <- lapply(X = split_seurat, FUN = function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = 'vst',nfeatures = 1500)
})
features <- SelectIntegrationFeatures(split_seurat,nfeatures = 1500)
split_seurat <- lapply(X = split_seurat, FUN = function(x){
  x <- ScaleData(x, vars.to.regress = c("mitoRatio",'G2M.Socre','S.Score'),features = features)
  x <- RunPCA(x,features = features)
})

anchors <- FindIntegrationAnchors(object.list = split_seurat,
                                  # normalization.method = "SCT",
                                  anchor.features = features,
                                  reduction = 'rpca',k.anchor = 20)
combined.sct <- IntegrateData(anchorset = anchors)
DefaultAssay(combined.sct)


combined.sct <- combined.sct %>% 
  ScaleData() %>%
  RunPCA( verbose = FALSE) %>% 
  RunUMAP( reduction = "pca", dims = 1:30, verbose = FALSE) %>% 
  FindNeighbors( reduction = "pca", dims = 1:30) %>% 
  FindClusters( resolution = c(0.4, 0.6, 0.8, 1.0,1.2))

p1 <- DimPlot(combined.sct,group.by = 'samples',cols = cr.col,raster = T) # 4.8*4
p1$layers[[1]]$mapping$alpha <- 0.7
p1 <- p1 + scale_alpha_continuous(range = 0.7, guide = F)
p1

combined.sct$samples <- factor(combined.sct$samples, levels = c('ND_M1','ND_M2',
                                                                'HFD_M1','HFD_M2',
                                                                'ND_F1','ND_F2',
                                                                'HFD_F1','HFD_F2'))
DimPlot(combined.sct,group.by = 'samples',split.by = 'samples',raster = F,ncol = 4) + scale_color_jco()# 15*4
# DimPlot(seurat_phase,group.by = 'samples',cols = cr.col,) # 15*4

rm(seurat_phase,anchors,p1)

Idents(combined.sct) <- 'integrated_snn_res.1'
DimPlot(combined.sct,raster = T,label = T,group.by = 'integrated_snn_res.1') + scale_color_igv()

#  regress cell cycle cell types -------

deg.res.1.re <- FindAllMarkers(combined.sct, only.pos = T)
deg.res.1.re.top20 <- deg.res.1.re %>% group_by(cluster) %>% top_n(n = 20,wt = avg_log2FC)

combined.sct <- combined.sct %>% RunUMAP( reduction = "pca", dims = 1:30,seed.use = 42 ,verbose = FALSE)
DimPlot(combined.sct,raster = T,label = T) +scale_color_igv()

# 2024.4.7
combined.sct <- RenameIdents(combined.sct, # not regress cell cycle celltypes
                             "0" = 'Adipocytes',
                             "1" = 'Adipocytes',
                             "2" = 'Adipocytes',
                             "3" = 'Adipocytes',
                             "4" = 'Adipocytes',
                             "6" = 'Adipocytes',
                             "8" = 'Adipocytes',
                             "9" = 'Adipocytes',
                             "16" = 'Adipocytes',
                             "19" = 'Adipocytes',
                             '18' = 'Adipocytes',
                             '12' = 'Adipocytes',
                             '13' = 'Adipocytes',
                             '14' = 'Fibroblast',
                             "20" = 'Adipocytes',
                             "22" = 'Adipocytes',
                             "41" = 'Adipocytes',
                             '25' = 'Pre-adipocytes',
                             '5' = 'FAPs',
                             '7' = 'FAPs',
                             '17' = 'FAPs',
                             '10' = 'FAPs',
                             '29' = 'FAPs',
                             '31' = 'FAPs',
                             '32' = 'FAPs',
                             '15' = 'ECs',
                             '16' = 'ECs',
                             '28' = 'ECs',
                             '43' = 'ECs',
                             '30' = 'ECs',
                             # '19' = 'Myofibroblast',
                             '42' = 'LECs',
                             '21' = 'Pericytes',
                             '39' = 'Pericytes',
                             '38' = 'Pericytes',
                             '23' = 'Adipocytes',
                             '24' = 'Adipocytes',
                             '35' = 'Adipocytes',
                             '40' = 'Schwann cells',
                             # '18' = 'KCs',
                             '26' = 'Erythoid cells',
                             '37' = 'Granulocytes',
                             # '38' = 'Granulocytes',
                             '21' = 'Monocytes',
                             '34' = 'DCs',
                             '33' = 'Macrophages',
                             '11' = 'Macrophages',
                             '27' = 'Macrophages',
                             '36' = 'Macrophages'
                             
)
DimPlot(combined.sct,raster = T,label = T) + scale_color_igv()

# 2024.4.11
combined.sct$large_ct[which(str_detect(combined.sct@meta.data$large_ct, "Pre-adipocytes"))] <- "Adipocytes"
combined.sct$large_ct[which(str_detect(combined.sct@meta.data$large_ct, "FAPs"))] <- "Pre-adipocytes"
combined.sct$large_ct[which(str_detect(combined.sct@meta.data$large_ct, "Prolif.macs"))] <- "Macrophages"
combined.sct$large_ct[which(str_detect(combined.sct@meta.data$large_ct, "Myocytes"))] <- "SMCs"
combined.sct$samples[which(str_detect(combined.sct@meta.data$samples, "HDF_F1"))] <- 'HFD_F1'
combined.sct$samples[which(str_detect(combined.sct@meta.data$samples, "HDF_F2"))] <- 'HFD_F2'
combined.sct$samples[which(str_detect(combined.sct@meta.data$samples, "HDF_M1"))] <- 'HFD_M1'
combined.sct$samples[which(str_detect(combined.sct@meta.data$samples, "HDF_M2"))] <- 'HFD_M2'

# 2024.4.13
combined.sct$large_ct <- as.character(combined.sct$large_ct_raw)
combined.sct$large_ct[which(str_detect(combined.sct@meta.data$large_ct_raw, "Pre-adipocytes"))] <- "Adipocytes"
combined.sct$large_ct[which(str_detect(combined.sct@meta.data$large_ct, "Fibroblast"))] <- "White adipocytes"
combined.sct$large_ct[which(str_detect(combined.sct@meta.data$large_ct, "Prolif.macs"))] <- "Macrophages"

combined.sct$large_ct[which(str_detect(combined.sct@meta.data$large_ct, "Myocytes"))] <- "SMCs"
combined.sct$large_ct[which(str_detect(combined.sct@meta.data$large_ct, "Adipocytes"))] <- "Brown adipocytes"
combined.sct$large_ct[which(str_detect(combined.sct@meta.data$large_ct_raw, "Adipocytes"))] <- "Brown adipocytes"


combined.sct$large_ct <- factor(combined.sct$large_ct, levels = c('Brown adipocytes','White adipocytes','Prolif.adipocytes','FAPs','ECs','LECs','Pericytes','SMCs','Schwann cells','Erythoid cells','Granulocytes','Monocytes','DCs','Macrophages'))

combined.sct <- combined.sct %>% RunUMAP( reduction = "pca", dims = 1:40,seed.use = 42 ,verbose = FALSE)
DimPlot(combined.sct,raster = T,label = F, cols = c(cr.col,npg.col), group.by = 'large_ct',raster.dpi = c(4096,4096),pt.size = 6) + ggtitle('Cell types') +
  NoAxes() + 
  # scale_color_igv() + 
  theme(panel.background=element_rect(fill='transparent', color='black'), 
        title = element_text(size = 10),
        legend.text = element_text(size = 10), legend.key.height=unit(.8,"line"))

# cell type markers
DefaultAssay(combined.sct) <- 'RNA'
Idents(combined.sct) <- combined.sct$large_ct
ct.markers <- FindAllMarkers(combined.sct,only.pos = T,) %>% filter(p_val_adj < 0.05)
export(ct.markers, file = '../01.analysis/01.all.samples/all.cell.type.markers.xlsx')
ct.markers.top20 <- ct.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 20, wt= avg_log2FC)
table(ct.markers$cluster)

DotPlot(combined.sct,group.by = 'large_ct', features = c('Pparg','Ppargc1b','Ppara',
                                                         'Ghr','Plcb1','Tshr',
                                                         'Mki67','Top2a','Cenpp',
                                                         'Pdgfra','Fbn1','Col3a1',
                                                         'Pecam1','Ptprb','Kdr',
                                                         'Prox1','Mmrn1','Reln',
                                                         'Abcc9','Kcnq5','Notch3',
                                                         'Myh8','Mybpc1','Myh3',
                                                         'Csmd1','Grid2','Nkain2',
                                                         'Hba-a1','Hba-a2','Hbb-bs',
                                                         'Cpa3','Hdc','Tpsb2',
                                                         'Csf3r','S100a9','S100a8',
                                                         'Arhgap15','Dock10','Dock2','Clec9a','Sirpa','Klrd1','Ms4a1','Il7r',
                                                         'Ptprc','Adgre1','Apoe'),
        cols = c('gray90','#e63946')) + 
  RotatedAxis() + NoGrid() + theme(axis.title = element_blank())
# Sankey plot -----
library(ggalluvial)
plot.data.num <- FetchData(combined.sct, 
                       vars = c("ident", "large_ct")) %>%
  dplyr::count(ident, large_ct) %>% 
  tidyr::spread(ident, n) 
plot.data.num[is.na(plot.data.num)] <- 0

# export(plot.data,file = '../cell.info.xlsx')

plot.data <- FetchData(combined.sct, 
                       vars = c("ident", "large_ct")) %>%
  dplyr::count(ident, large_ct) %>% 
  group_by(ident) %>% 
  summarise(Prop = n / sum(n))

temp <- data.frame(ident = 'HDF_F2',Prop = 0)
plot.data <- rbind( plot.data[1:19,],temp,plot.data[20:111,])

plot.data$large_ct <- factor(rep(c(levels(combined.sct$large_ct)), 8), levels = levels(combined.sct$large_ct))

plot.data$levels <- c(rep(1:14,8))
plot.data$ident <- factor(plot.data$ident, levels = c('ND_M1','ND_M2','HDF_M1','HDF_M2','ND_F1','ND_F2','HDF_F1','HDF_F2'))

ggplot(plot.data,
       aes(x = ident, stratum = large_ct, alluvium = levels,
           y = Prop,
           fill = large_ct, label = large_ct)) +
  # scale_fill_manual(values = c(lancent.col[c(2,1,3:6)])) +
  scale_fill_manual(values = c(cr.col,npg.col)) +  
  geom_flow(stat = "alluvium", 
            lode.guidance = "frontback",
            width = 0.3,
            color = "darkgray") +
  geom_stratum(alpha = .8) +
  guides(fill = guide_legend(title = 'Cell Type')) +
  theme_bw() +
  theme(legend.position = "right",
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15))

plot.data.1 <- plot.data[c(29:56,85:112),]
plot.data.1 <- plot.data[c(1:28,57:84),]

ggplot(plot.data.1,
       aes(x = ident, stratum = large_ct, alluvium = levels,
           y = Prop,
           fill = large_ct, label = large_ct)) +
  # scale_fill_manual(values = c(lancent.col[c(2,1,3:6)])) +
  scale_fill_manual(values = c(cr.col,npg.col)) +  
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "darkgray") +
  geom_stratum(alpha = .7) +
  guides(fill = guide_legend(title = 'Cell Type')) +
  theme_bw() +
  theme(legend.position = "right",
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15))

# plot.data$Prop <- round(plot.data$Prop, 4)
plot.data.wide <- plot.data %>% tidyr::spread(ident, Prop) 
cell.ratio <- list(cell.number = plot.data.num, cell.ratio = plot.data.wide)
export(cell.ratio,file = '../01.analysis/01.all.samples/cell.info.xlsx')

DimPlot(combined.sct,raster = F,label = F,split.by = 'samples', group.by = 'samples') + scale_color_jco()

combined.sct$condition <- NA
combined.sct$condition[which(str_detect(combined.sct@meta.data$cells, "^HDF_"))] <- "MO"
combined.sct$condition[which(str_detect(combined.sct@meta.data$cells, "^ND_"))] <- "MC"
combined.sct$con_gan <- paste0(combined.sct$condition, "_",combined.sct$gander)

save.image()
