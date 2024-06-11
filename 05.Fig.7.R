# data prepare -------

library(tidyverse)
# gene info 
mus.gene.info <- read.table('../mouse.110.gtf',header = F,sep = '\t')
colnames(mus.gene.info) <- c('geneid','symbol','biotype')

# counts 
hl.counts <- read.table('../../02.quant/all_counts (2).txt',header = T,sep = '\t',row.names = NULL)
colnames(hl.counts)

samp.name <- c(paste0('MC+PBS_F',seq(1,6)), 
               paste0('MO+FGF21_F',seq(1,6)),
               paste0('MO+PBS_F',seq(1,6)),
               paste0('MC+PBS_M',seq(1,6)),
               paste0('MO+FGF21_M',seq(1,6)),
               paste0('MO+PBS_M',seq(1,6)))

colnames(hl.counts) <- c('geneid',samp.name)
hl.counts <- hl.counts %>% 
  left_join(y = mus.gene.info, by = 'geneid') %>% 
  na.omit() %>% 
  filter(biotype == 'protein_coding') %>% 
  group_by(`symbol`) %>% slice_sample(n = 1) %>% 
  column_to_rownames(var = 'symbol') %>% 
  dplyr::select(!c('geneid','biotype')) %>% 
  round(0)

# TPM 
hl.tpm <- read.table('../../02.quant/all_quant_TPM (2).txt',header = T,sep = '\t',row.names = NULL)
colnames(hl.tpm) <- c('geneid',samp.name)
hl.tpm <- hl.tpm %>% 
  left_join(y = mus.gene.info, by = 'geneid') %>% 
  na.omit() %>% 
  filter(biotype == 'protein_coding') %>% 
  group_by(`symbol`) %>% slice_sample(n = 1) %>% 
  column_to_rownames(var = 'symbol') %>% 
  dplyr::select(!c('geneid','biotype')) %>% round(2)

# pick genes 
pick.genes <- unique(c(rownames(hl.tpm[matrixStats::rowMedians(as.matrix(hl.tpm),cols = 1:6) > 0.1, ]),
                       rownames(hl.tpm[matrixStats::rowMedians(as.matrix(hl.tpm),cols = 7:12) > 0.1, ]),
                       rownames(hl.tpm[matrixStats::rowMedians(as.matrix(hl.tpm),cols = 13:18) > 0.1, ]),
                       rownames(hl.tpm[matrixStats::rowMedians(as.matrix(hl.tpm),cols = 19:24) > 0.1, ]),
                       rownames(hl.tpm[matrixStats::rowMedians(as.matrix(hl.tpm),cols = 25:30) > 0.1, ]),
                       rownames(hl.tpm[matrixStats::rowMedians(as.matrix(hl.tpm),cols = 30:36) > 0.1, ])
))

smallestGroupSize <- 3
keep <- rowSums(hl.tpm >= 0.1) >= smallestGroupSize
test <- hl.tpm[keep,]
Vennerable::Venn(list(a = pick.genes, b = rownames(test)))
rm(smallestGroupSize,keep,test)

# clean dataset 
hl.counts.clean <- hl.counts[pick.genes,]
hl.counts.clean$Gene <- rownames(hl.counts.clean)
hl.counts.clean <- hl.counts.clean[,c(37,1:36)]
hl.tpm.clean <- hl.tpm[pick.genes,]
hl.tpm.clean$Gene <- rownames(hl.tpm.clean)
hl.tpm.clean <- hl.tpm.clean[,c(37,1:36)]

# fig 7H -----

library(scales)
library(FactoMineR)
library(factoextra)
library(ggsci)
library(reshape2)
library(GGally)
library(grid)
library(DESeq2)
library(EDASeq)

ann <- data.frame(Group = c(rep(c(rep('MC+PBS',6), rep('MO+FGF21',6), rep('MO+PBS',6)),2)),
                  Gender = c(rep('Female',18),rep('Male',18)))

rownames(ann) <- samp.name
samp.cols <- pal_npg(palette = c("nrc"), alpha = 1)(8)
group_list <- factor(c(paste0(ann$Group,'_',ann$Gender)))
group_list_large <- factor(c(ann$Group))

d_p <- function(tpm,group_list,name){
  tpm.t = as.data.frame(t(tpm))
  tpm.t = cbind(tpm.t,group_list)
  tpm.pca <- PCA(tpm.t[,-ncol(tpm.t)],graph = FALSE)
  fviz_screeplot(tpm.pca, addlabels = TRUE, ylim = c(0, 40),title='Dim choose')
  fviz_pca_ind(tpm.pca,
               mean.point=F,
               addEllipses = T,
               ellipse.type = "t", #convex norm confidence t euclid
               # repel = T,
               habillage = group_list,
               label = "none",
               geom.ind = c("point"),
               fill.ind = tpm.t$group_list,
               palette = c("#0B996F",'#7C6DAF','#D6570D','#D6570D','#e29578'),
               # palette = c("#7C6DAF",'#C7AED4','#0B996F','#83c5be','#D6570D','#e29578'),
               legend.title = "Groups",
               pointsize = 3,
               pointshape = 21,
               col.ind = 'dark',
               title = name
  ) + theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15))
}
d_p(hl.tpm.clean,group_list_large,'tpm')

# fig 7I ------
d_p_pca <- function(tpm,group_list){
  tpm.t = as.data.frame(t(tpm))
  tpm.t = cbind(tpm.t,group_list)
  tpm.pca <- PCA(tpm.t[,-ncol(tpm.t)],graph = FALSE)
}
tpm.pca <- d_p_pca(na.omit(hl.tpm.clean[FAT.gene$`fatty acid oxidation`,]),group_list_large)

dist.pca <- tpm.pca$ind$coord[,1:2] %>% as.data.frame()
dist.pca <- dist.pca[c(1:6,19:24,7:12,25:30,13:18,31:36),]
dist.pca.euc <- dist(dist.pca[1:36,],method = 'euclidean') %>% as.matrix()
dist.ND_FGF <- dist.pca.euc[13:24,1:12] %>% as.data.frame()
dist.FGF_HFD <- dist.pca.euc[13:24,25:36] %>% as.data.frame()

dist.ND_FGF.long <- melt(dist.ND_FGF)
dist.FGF_HFD.long <- melt(dist.FGF_HFD)

dist.group <- data.frame(dist = c(dist.ND_FGF.long$value,dist.FGF_HFD.long$value),
                         group = c(rep('ND_FGF',144),rep('FGF_HFD',144)))

p1 <- ggpubr::ggboxplot(dist.group, x="group", y="dist", width = 0.5, 
                color = "black"
                fill="group"
                palette = "npg",
                xlab = F, 
                bxp.errorbar=T,
                bxp.errorbar.width=0.5,
                size=1, 
                outlier.shape=NA,
                legend = "right",
                alpha = 0.8) + 
  geom_jitter(color="black", size=0.4, alpha=1,width = 0.1, height = 0.5)+
  ylab('Euclidean distance')  + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),axis.title.x = element_blank(),legend.position = 'none')

my_comparisons <- list(c("ND_FGF", "FGF_HFD"))
p1+ggpubr::stat_compare_means(comparisons = my_comparisons,
                      method = "wilcox.test")

# fig 7J -----
n.mito.gene <- import('E:/GSEA/线粒体数据库/mmu.mito.genes.xlsx')
m.mito.gene <- mus.gene.info[mus.gene.info$biotype == 'protein_coding',]$symbol[which(str_detect(mus.gene.info[mus.gene.info$biotype == 'protein_coding',]$symbol, "^mt-"))]

mito.tpm <- hl.tpm.clean[unique(c(n.mito.gene$Symbol,m.mito.gene)),] %>% na.omit()
mito.rlog <- rld_mat[c(n.mito.gene$Symbol,m.mito.gene),] %>% na.omit() # 02.131 
mito.vst <-  rld_mat[c(n.mito.gene$Symbol,m.mito.gene),] %>% na.omit() # 02.131 
d_p(mito.tpm,group_list_large,'mito.gene.TPM')

# fig 7K&L ------
options(scipen = 10)

library(IHW)
library(ashr)
library(tidyverse)

comp_fat <- hl.counts.clean[,c(13:18,31:36,1:6,19:24)]
colData_fat <- data.frame(group_list = colData[c(13:18,31:36,1:6,19:24),])
colData_fat$group_list_large <- c(rep('Obesity',12),rep('Normal',12))
rownames(colData_fat) <- colnames(comp_fat)
dds <- DESeqDataSetFromMatrix(countData = comp_fat,colData = colData_fat,design = ~group_list_large)
dds <- DESeq(dds,minReplicatesForReplace = 7)
fat_res_ashr <- lfcShrink(dds, contrast = c('group_list_large','Normal','Obesity'), type = 'ashr') # adjust LFC
summary(fat_res_ashr)
fat_res_ashr <- fat_res_ashr %>% 
  data.frame() %>%
  rownames_to_column(var="geneid") %>%
  as_tibble() %>%
  arrange('padj')


comp_efct <- hl.counts.clean[,c(13:18,31:36,7:12,25:30)]
colData_efct <- data.frame(group_list = colData[c(13:18,31:36,7:12,25:30),])
colData_efct$group_list_large <- c(rep('Obesity',12),rep('FGF',12))
rownames(colData_efct) <- colnames(comp_efct)
dds <- DESeqDataSetFromMatrix(countData = comp_efct,colData = colData_efct,design = ~group_list_large)
dds <- DESeq(dds,minReplicatesForReplace = 7)
efct_res_ashr <- lfcShrink(dds, contrast = c('group_list_large','FGF','Obesity'), type = 'ashr') # adjust LFC
summary(efct_res_ashr)
efct_res_ashr <- efct_res_ashr %>% 
  data.frame() %>%
  rownames_to_column(var="geneid") %>%
  as_tibble() %>%
  arrange('padj')

library(clusterProfiler)
library(enrichplot)
library(GSEABase)
library(org.Mm.eg.db)

gene_for_gsea <- fat_res_ashr$log2FoldChange 
names(gene_for_gsea) <- fat_res_ashr$geneid
geneList=sort(gene_for_gsea,decreasing = T)
mus.gmt <- read.gmt('E:/GSEA/v2023.mouse/msigdb_v2023.2.Mm_files_to_download_locally/msigdb_v2023.2.Mm_GMTs/m5.go.bp.v2023.2.Mm.symbols.gmt')

length(unique(mus.gmt$term))
egmt <- GSEA(geneList, TERM2GENE=mus.gmt, 
             minGSSize = 1,
             pvalueCutoff = 0.1,
             verbose=FALSE)

egmt.order <- egmt[order(egmt$enrichmentScore, decreasing = T),]

gseaplot2(egmt,geneSetID = 'GOBP_REGULATION_OF_FATTY_ACID_BETA_OXIDATION', color = "orange", rel_heights=c(1, .2, .6),pvalue_table = T)
gseaplot2(egmt,geneSetID = 'GOBP_ADAPTIVE_THERMOGENESIS', color = "orange", rel_heights=c(1, .2, .6),pvalue_table = T)
