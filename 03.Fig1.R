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

# fig 1F -----
CellTypeCompositionAnalysis.R code path:https://github.com/Teichlab/sctkr/blob/2a024cafef1aae192bf9656349449c5a84d1c6ed/R/CellTypeCompositionAnalysis.R
source('CellTypeCompositionAnalysis.R')
plot.data.num <- FetchData(combined.sct, 
                       vars = c("ident", "large_ct")) %>%
  dplyr::count(ident, large_ct) %>% 
  tidyr::spread(ident, n) 
plot.data.num[is.na(plot.data.num)] <- 0
export(plot.data,file = '../cell.number.xlsx')

cell_types <- cell.number$large_ct
n_cells_per_sample <- colSums(cell.number)
n_var_cats <- 2
sample_cats <- tibble(
  Sample_ID = sample_ids,
  Treatment = c(rep('MO',4),rep('MC',4)),
  Gender = c(rep('Female',2),rep('Male',2),rep('Female',2),rep('Male',2)),
  cell.num = n_cells_per_sample
)
sample_num1_values <- rnorm(length(sample_ids), mean = 10, sd = 1) 
obs_tbl <- data.frame(
  Sample_ID = rep(sample_ids, c(n_cells_per_sample)),
  Treatment = rep(sample_cats$Treatment, c(n_cells_per_sample)),
  Gender = rep(sample_cats$Gender, c(n_cells_per_sample)),
  #Var_Num1 = rnorm(n_cells_per_sample * length(sample_ids), mean = 10, sd = 1)
  Var_Num1 = rep(sample_num1_values, c(n_cells_per_sample))
)
obs_tbl$Cell_type <- c(rep(cell.number$large_ct,c(cell.number$HDF_F1)),
                       rep(cell.number$large_ct,c(cell.number$HDF_F2)),
                       rep(cell.number$large_ct,c(cell.number$HDF_M1)),
                       rep(cell.number$large_ct,c(cell.number$HDF_M2)),
                       rep(cell.number$large_ct,c(cell.number$ND_F1)),
                       rep(cell.number$large_ct,c(cell.number$ND_F2)),
                       rep(cell.number$large_ct,c(cell.number$ND_M1)),
                       rep(cell.number$large_ct,c(cell.number$ND_M2))
                       )
results <- CellTypeCompositionAnalysis(obs_tbl, "Sample_ID", "Cell_type", c("Treatment", "Gender"), "Var_Num1")
vars3 <- list(Treatment = c('MC','MO'), Gender = c('Female','Male'))
ranef_plot <- plot_ranef(ranef_tbl, vars = vars3, celltypes = cell_types, celltype_order = rev(cell_types),
                         maxFC = 2, LTSR2p = FALSE)
sdse_plot <- plot_sdse(sdse_tbl, "Sample_ID", ci = 0.95, xlim = c(0, .5))

# fig 1G & 1H------
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

adipo.only.inte.2$con_gan <- factor(adipo.only.inte.2$con_gan, levels = c('MC_male',
                                                                          'MO_male',
                                                                          'MC_female',
                                                                          'MO_female'))

DimPlot(adipo.only.inte.2,raster = F,label = F, cols = simpson.col, group.by = 'con_gan',raster.dpi = c(4096,4096),split.by = 'con_gan') + ggtitle('') +
  NoAxes() + 
  # scale_color_igv() + 
  theme(panel.background=element_rect(fill='transparent', color='black'), 
        title = element_text(size = 10),
        legend.text = element_text(size = 10), legend.key.height=unit(.8,"line"))

DefaultAssay(adipo.only.inte.2) <-'RNA'
adipo.only.inte.2 <- NormalizeData(adipo.only.inte.2)
FeaturePlot(adipo.only.inte.2, features = c('Peg3','Ucp1'),ncol = 2,
            split.by = 'con_gan',
            # blend = T,
            cols =c("grey80", "lemonchiffon1","indianred1", "darkred"),
            order = T,
            min.cutoff = 0.5,raster = F) 

p1 <- FeaturePlot(adipo.only.inte.2, features = c('Klb'),
                  pt.size = .05,
                  cols = c("grey80", "lemonchiffon1","indianred1", "darkred"),
                  combine = F,raster = F,order = T,
                  split.by = 'con_gan',min.cutoff = 0.5) 
for(i in 1:length(p1)) {
  p1[[i]] <- p1[[i]] + NoLegend() + NoAxes() + theme(panel.background=element_rect(fill='transparent', color='black'), title = element_text(size = 8))
}
patchwork::wrap_plots(c(p1))

col_fun = colorRamp2(c(0, 0.25, 0.75, 1), c("grey80", "lemonchiffon1","indianred1", "darkred"))
lgd = Legend(col_fun = col_fun, title = "Relative \n expression", at = c(0, 1), 
             labels = c("min", "max"))
draw(lgd)

for (ct in unique(adipo.only.inte.2$gander)) {
  print(ct)
  Idents(adipo.only.inte.2) <- 'gander'
  sub.seu.1 <- subset(adipo.only.inte.2, idents = ct)
  for (gp in nrow(comlist)) {
    age1 <- comlist[gp,1]
    age2 <- comlist[gp,2]
    DEGs <- FindMarkers(sub.seu.1,
                        ident.1 = age1,
                        ident.2 = age2, 
                        logfc.threshold = 0.01,
                        group.by = 'condition') %>% dplyr::filter(p_val_adj < 0.05)
    DEGs.name <- paste('b.adi.all',ct,age2,age1,sep = '_')
    assign(DEGs.name,DEGs)
  }
}

b.adi.all.deg <- list(b.adi.all_male_MO_MC = b.adi.all_male_MO_MC,
                      b.adi.all_female_MO_MC = b.adi.all_female_MO_MC)
intersect(rownames(b.adi.all_male_MO_MC), rownames(b.adi.all_female_MO_MC))
rio::export(b.adi.all.deg, file = '../01.analysis/063fig.final//Table4.brwon adipocytes wilcox.test.lfc0.01.xlsx',row.names = T)

b.adi.male.lfc.0.2 <- b.adi.all_male_MO_MC[abs(b.adi.all_male_MO_MC$avg_log2FC) > 0.2,]
b.adi.male.lfc.0.2$gene<- rownames(b.adi.male.lfc.0.2)
b.adi.male.lfc.0.2$cluster <- ifelse(b.adi.male.lfc.0.2$avg_log2FC> 0.2,
                                      'MC_male','MO_male')

b.adi.female.lfc.0.2 <- b.adi.all_female_MO_MC[abs(b.adi.all_female_MO_MC$avg_log2FC) > 0.2,]
b.adi.female.lfc.0.2$gene<- rownames(b.adi.female.lfc.0.2)
b.adi.female.lfc.0.2$cluster <- ifelse(b.adi.female.lfc.0.2$avg_log2FC> 0.2,
                                      'MC_female','MO_female')
tmp <- rbind(b.adi.male.lfc.0.2,
             b.adi.female.lfc.0.2)

Idents(adipo.only.inte.2) <- adipo.only.inte.2$con_gan
b.adi.data <- prepareDataFromScRNA(adipo.only.inte.2,
                                diffData = tmp,
                                showAverage = T,
                                group.by = 'con_gan',keep.uniqGene = F,
                                scale.data = T)
saveRDS(b.adi.data, file = '../01.analysis/063fig.final/b.adi.data.RDS',compress = F)

library(EnhancedVolcano)

keyvals.colour <- ifelse(
  b.adi.all_male_MO_MC$avg_log2FC < -0.25, 'royalblue',
  ifelse(b.adi.all_male_MO_MC$avg_log2FC > 0.25, 'red3',
         'gray80'))
keyvals.colour[is.na(keyvals.colour)] <- 'gray80'
names(keyvals.colour)[keyvals.colour == 'red3'] <- 'mND high'
names(keyvals.colour)[keyvals.colour == 'gray80'] <- 'No sig'
names(keyvals.colour)[keyvals.colour == 'royalblue'] <- 'mHFD high'

EnhancedVolcano(b.adi.all_male_MO_MC,
                lab = rownames(b.adi.all_male_MO_MC),
                selectLab = 'Ucp1',
                x = 'avg_log2FC',
                y = 'p_val_adj',
                FCcutoff = 0.25,
                pCutoff = 0.05,
                ylab = bquote(~-Log[10]~ 'padj'),
                pointSize = 1.0,
                labSize = 4,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = T,
                # col = c('gray80', 'gray80', 'gray80', 'red3'),
                colCustom = keyvals.colour,
                colAlpha = 2/5,
                legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = .5,
                max.overlaps = 300,
                colConnectors = '#14213d',
                title = 'Male',
                subtitle = '',
                xlim = c(-1.5,1.5)
) 


keyvals.colour <- ifelse(
  b.adi.all_female_MO_MC$avg_log2FC < -0.25, 'royalblue',
  ifelse(b.adi.all_female_MO_MC$avg_log2FC > 0.25, 'red3',
         'gray80'))
keyvals.colour[is.na(keyvals.colour)] <- 'gray80'
names(keyvals.colour)[keyvals.colour == 'red3'] <- 'mND high'
names(keyvals.colour)[keyvals.colour == 'gray80'] <- 'No sig'
names(keyvals.colour)[keyvals.colour == 'royalblue'] <- 'mHFD high'

EnhancedVolcano(b.adi.all_female_MO_MC,
                lab = rownames(b.adi.all_female_MO_MC),
                selectLab = 'Ucp1',
                x = 'avg_log2FC',
                y = 'p_val_adj',
                FCcutoff = 0.25,
                pCutoff = 0.05,
                ylab = bquote(~-Log[10]~ 'padj'),
                pointSize = 1.0,
                labSize = 4,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = T,
                # col = c('gray80', 'gray80', 'gray80', 'red3'),
                colCustom = keyvals.colour,
                colAlpha = 2/5,
                legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = .5,
                max.overlaps = 300,
                colConnectors = '#14213d',
                title = 'Female',
                subtitle = '',
                xlim = c(-1.5,1.5)
) 

# fig 1H & I ------
library(ggplot2)
# install.packages('cols4all')
library(cols4all)

go.plot <- read.table('clipboard',header = T,sep = '\t')
male.go <- go.plot[1:16,]
male.go$'-log10pvalue' <- ifelse(male.go$group == 'MC_male', -log10(male.go$pvalue), -(-log10(male.go$pvalue)))
male.go$Description <- factor(male.go$Description,
                          levels = rev(male.go$Description))

female.go <- go.plot[17:32,]
female.go$'-log10pvalue' <- ifelse(female.go$group == 'MC_female', -log10(female.go$pvalue), -(-log10(female.go$pvalue)))
female.go$Description <- factor(female.go$Description,
                              levels = rev(female.go$Description))

mytheme3 <- theme(legend.position = 'none',
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  axis.line.x = element_line(color = 'grey60', size = 1.1),
                  axis.text = element_text(size = 12),
                  plot.margin = margin(t = 5.5, r = 10, l = 5.5, b = 5.5))

up <- female.go[which(female.go$`-log10pvalue` > 0),]
down <- female.go[which(female.go$`-log10pvalue` < 0),]

p3 <- ggplot(female.go,aes(x = `-log10pvalue`, y = Description)) +
  geom_col(aes(fill = `-log10pvalue`), width = 0.1) + #添加条形图，收窄柱子为一条线
  geom_point(aes(size = ratio, color = `-log10pvalue`)) +
  scale_size_continuous(range = c(2, 7)) +
  scale_color_continuous_c4a_div('div_gn_wh_rd', mid = 0, reverse = F) +
  scale_fill_continuous_c4a_div('div_gn_wh_rd', mid = 0, reverse = F) +
  scale_x_continuous(breaks = seq(-20, 20, by = 5),
                     labels = abs(seq(-20, 20, by = 5))) +
  ylab('')

p5 <- p3 + theme_bw() + mytheme3 +
  geom_text(data = up,aes(x = -0.5, y = Description, label = Description),
            size = 3.5, hjust = 1)+ 
  geom_text(data = down, aes(x = 0.5, y = Description, label = Description),
            size = 3.5, hjust = 0) + 
  geom_text(x = 15, y = 15, label = "Up", size = 6, color = '#EE4E4A') +
  geom_text(x = -16, y = 5, label = "Down", size = 6, color = '#419A18')
p5 

# fig 1J ------
p1 <- DimPlot(adipo.only.inte.2,group.by = 'con_gan',cols = cr.col,raster = T) # 4.8*4
p1$layers[[1]]$mapping$alpha <- 0.7
p1 <- p1 + scale_alpha_continuous(range = 0.7, guide = F)
p1

# fig 1K -----
FeaturePlot(adipo.only.inte.2, features = c('Ucp1'), cols = c('#a8dadc','#f1faee',"#e63946"),split.by = 'con_gan',order = T)
