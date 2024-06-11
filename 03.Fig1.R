# fig 1D ------
library(plot1cell)
combined.sct$large_ct <- Idents(combined.sct)
combined.sct$large_ct <- as.character(combined.sct$large_ct)
combined.sct$large_ct[which(str_detect(combined.sct@meta.data$large_ct, "DCs"))] <- "Lymphocytes"
combined.sct$large_ct[which(str_detect(combined.sct@meta.data$large_ct, "Prolif.adipocytes"))] <- "Brown adipocytes"
combined.sct$large_ct[which(str_detect(combined.sct@meta.data$large_ct, "Granulocytes"))] <- "Myeloid cells"
combined.sct$large_ct[which(str_detect(combined.sct@meta.data$large_ct, "Monocytes"))] <- "Myeloid cells"
combined.sct$large_ct[which(str_detect(combined.sct@meta.data$large_ct, "Macrophages"))] <- "Myeloid cells"

combined.sct$large_ct <- factor(combined.sct$large_ct, levels = c('Brown adipocytes','White adipocytes','FAPs','ECs','LECs','Pericytes','SMCs','Schwann cells','Erythoid cells','Myeloid cells','Lymphocytes'))
combined.sct$large_ct_v2 <- Idents(combined.sct)

plot.cir.test <- function (data_plot, do.label = T, contour.levels = c(0.2, 0.4, 0.6), 
                           pt.size = 0.5, kde2d.n = 1000, contour.nlevels = 100, bg.color = "#F9F2E4", 
                           col.use = NULL, label.cex = 0.5, repel = FALSE) 
{
  centers <- data_plot %>% dplyr::group_by(Cluster) %>% summarise(x = median(x = x), 
                                                                  y = median(x = y))
  z <- MASS::kde2d(data_plot$x, data_plot$y, n = kde2d.n)
  celltypes <- names(table(data_plot$Cluster))
  cell_colors <- (scales::hue_pal())(length(celltypes))
  if (!is.null(col.use)) {
    cell_colors = col.use
    col_df <- data.frame(Cluster = celltypes, color2 = col.use)
    cells_order <- rownames(data_plot)
    data_plot <- merge(data_plot, col_df, by = "Cluster")
    rownames(data_plot) <- data_plot$cells
    data_plot <- data_plot[cells_order, ]
    data_plot$Colors <- data_plot$color2
  }
  circos.clear()
  par(bg = bg.color)
  circos.par(cell.padding = c(0, 0, 0, 0), track.margin = c(0.01, 0), track.height = 0.01, gap.degree = c(rep(2, (length(celltypes) -  1)), 12), points.overflow.warning = FALSE)
  circos.initialize(sectors = data_plot$Cluster, x = data_plot$x_polar2)
  circos.track(data_plot$Cluster, data_plot$x_polar2, y = data_plot$dim2, 
               bg.border = NA, panel.fun = function(x, y) {
                 circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + 
                               mm_y(4), CELL_META$sector.index, cex = 0.5, 
                             col = "black", facing = "bending.inside", niceFacing = T)
                 breaks = seq(0, 100, by = 50)
                 circos.axis(labels.cex = 0.3, col = "black", labels.col = "black",major.at = breaks,
                             labels = paste0(breaks, "%"))
               })
  for (i in 1:length(celltypes)) {
    dd <- data_plot[data_plot$Cluster == celltypes[i], ]
    circos.segments(x0 = min(dd$x_polar2), y0 = 0, x1 = max(dd$x_polar2), 
                    y1 = 0, col = cell_colors[i], lwd = 3, sector.index = celltypes[i])
  }
  text(x = 1, y = 0.1, labels = "Cluster", cex = 0.4, col = "black", 
       srt = -90)
  points(data_plot$x, data_plot$y, pch = 19, col = alpha(data_plot$Colors, 
                                                         0.2), cex = pt.size)
  contour(z, drawlabels = F, nlevels = 100, levels = contour.levels, 
          col = "#ae9c76", add = TRUE)
  if (do.label) {
    if (repel) {
      textplot(x = centers$x, y = centers$y, words = centers$Cluster, 
               cex = label.cex, new = F, show.lines = F)
    }
    else {
      text(centers$x, centers$y, labels = centers$Cluster, 
           cex = label.cex, col = "black")
    }
  }
}

cell_order_test <- function (dat) 
{
  celltypes <- names(table(dat$Cluster))
  new_dat <- list()
  for (i in 1:length(celltypes)) {
    dat$Cluster <- as.character(dat$Cluster)
    dat1 <- dat[dat$Cluster == celltypes[i], ]
    dat1$x_polar <- (1:nrow(dat1))/nrow(dat1)
    new_dat[[i]] <- dat1
  }
  new_dat <- do.call("rbind", new_dat)
  new_dat
}

prepare_circlize_data_test <- function (seu_obj, scale = 0.8) 
{
  celltypes <- levels(seu_obj)
  cell_colors <- (scales::hue_pal())(length(celltypes))
  data_plot <- get_metadata(seu_obj, color = cell_colors, 
                            coord_scale = scale)
  data_plot <- cell_order_test(data_plot)
  data_plot$x_polar2 <- log10(data_plot$x_polar)
  data_plot
}

Idents(combined.sct) <- combined.sct$large_ct
circ_data <- prepare_circlize_data_test(combined.sct, scale = 0.8 )
circ_data$Cluster <- factor(circ_data$Cluster, levels = c('Brown adipocytes','White adipocytes','FAPs','ECs','LECs','Pericytes','SMCs','Schwann cells','Erythoid cells','Myeloid cells','Lymphocytes'))
set.seed(1234)
circ_data$x_polar2 <- circ_data$x_polar*100

cluster_colors<- c(cr.col,npg.col)[1:11]
con_colors<-c('#D6570D','#0B996F') # H N
sex_colors<-c('#f47068','#0e606b') #fe ma
con_gan_colors <-  c('#8ecae6','#ffb703', "#219ebc",'#fb8500')

plot.cir.test(circ_data,do.label = T, pt.size = 0.01, col.use = cluster_colors ,bg.color = 'white', kde2d.n = 1000, repel = T, label.cex = 1)

add_track(circ_data, group = "condition", colors = con_colors, track_num = 2) 
add_track(circ_data, group = "gander",colors = sex_colors, track_num = 3) 

# fig 1E -----
Idents(combined.sct) <- combined.sct$large_ct
com.dot.new <- function (seu_obj, features, celltypes = NULL, groups, color.palette = NULL, 
                         strip.color = NULL) 
{
  pb <- progress_bar$new(format = "  Ploting [:bar] :percent eta: :eta", 
                         clear = FALSE, total = length(features), width = 100)
  plot_list <- list()
  for (i in 1:length(features)) {
    pp <- invisible(complex_dotplot_single(seu_obj = seu_obj, 
                                           feature = features[i], groups = groups, celltypes = celltypes))
    pp <- pp$data
    pp$gene <- features[i]
    plot_list[[i]] <- pp
    pb$tick()
    Sys.sleep(1/length(features))
  }
  all_data <- do.call("rbind", plot_list)
  all_data$gene <- factor(all_data$gene, levels = rev(features))
  all_data$celltype <- factor(all_data$celltype, levels = levels(seu_obj))
  if (is.null(color.palette)) {
    color.palette <- colorRampPalette(c("grey80", "lemonchiffon1",
                                        "indianred1", "darkred"))(255)
  }
  p <- invisible(ggplot(all_data, aes(x = groups, y = gene)) + 
                   geom_tile(fill = "white", color = "white") + 
                   geom_point(aes(colour = avg.exp, size = pct.exp), alpha = 0.9) + 
                   scale_color_gradientn(colours = color.palette) + 
                   scale_size(range = c(0, 5)) +
                   theme(
                     panel.background = element_rect(fill = "white", colour = "black"), 
                     axis.text.x = element_text(angle = 45,hjust = 1),
                     plot.title = element_text(size = 10, hjust = 0.5,face = "bold"), 
                     axis.text = element_text(size = 12), 
                     axis.title = element_text(size = 8), 
                     legend.text = element_text(size = 8), 
                     legend.title = element_text(size = 12),
                     legend.position = "right", 
                     strip.text = element_text(size = 8, colour = "black",face = "bold")) + 
                   ylab("") + xlab("") + ggtitle("") + 
                   facet_wrap(~celltype, ncol = length(levels(seu_obj))))
  g <- change_strip_background(p, type = "top", strip.color = strip.color)
  print(grid.draw(g))
}
  
com.dot.new(combined.sct, feature = c('Pparg','Ppargc1b','Ppara',
                                                 'Ghr','Plcb1','Tshr',
                                                 'Pdgfra','Fbn1','Col3a1',
                                                 'Pecam1','Ptprb','Kdr',
                                                 'Prox1','Mmrn1','Reln',
                                                 'Abcc9','Kcnq5','Notch3',
                                                 'Myh8','Mybpc1','Myh3',
                                                 'Csmd1','Grid2','Nkain2',
                                                 'Hba-a1','Hba-a2','Hbb-bs',
                                                 'Adgre1','Apoe','Ptprc',
                                                 'Arhgap15','Ikzf1','Dock2')
                       ,groups = "con_gan",strip.color = cluster_colors)

# fig 1F ------
Idents(combined.sct) <- combined.sct$large_ct
adipo.only <- subset(combined.sct, cells = WhichCells(combined.sct, idents = c('Brown adipocytes')))
dim(adipo.only)
DimPlot(adipo.only,raster = T)

adipo.only <- CreateSeuratObject(counts = GetAssayData(object = adipo.only, slot = "counts"),meta.data = adipo.only@meta.data)
adipo.split <- SplitObject(adipo.only,split.by = 'samples')
plan(multisession, workers=1)
adipo.split <- lapply(X = adipo.split, FUN = function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = 'vst',nfeatures = 1500)
})
features <- SelectIntegrationFeatures(adipo.split,nfeatures = 1500)

plan(multisession, workers=50)
adipo.split <- lapply(X = adipo.split, FUN = function(x){
  x <- ScaleData(x, vars.to.regress = c("mitoRatio",'G2M.Score','S.Score'),features = features)
  x <- RunPCA(x,features = features)
})
save.image()
anchors <- FindIntegrationAnchors(object.list = adipo.split,
                                  # normalization.method = "SCT",
                                  anchor.features = features,
                                  reduction = 'rpca',k.anchor = 20)
adipo.only.inte.2 <- IntegrateData(anchorset = anchors)

DefaultAssay(adipo.only.inte.2) <- 'integrated'
adipo.only.inte.2 <- adipo.only.inte.2 %>% 
  ScaleData() %>%
  RunPCA( verbose = FALSE) %>% 
  RunUMAP( reduction = "pca", dims = 1:40, verbose = FALSE) 





