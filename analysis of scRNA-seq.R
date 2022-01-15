options(stringsAsFactors = FALSE)
library(Seurat)
library(MAESTRO)
library(dplyr)
library(plyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(GGally)
library(reticulate)
CT1 <- CreateSeuratObject(Read10X('/data/hxl/con1/filtered_feature_bc_matrix'), project = "CT1")
CT2 <- CreateSeuratObject(Read10X('/data/hxl/con2/filtered_feature_bc_matrix'), project = "CT2")
AI1 <- CreateSeuratObject(Read10X('/data/hxl/AI1/filtered_feature_bc_matrix'), project = "AI1")
AI2 <- CreateSeuratObject(Read10X('/data/hxl/AI2/filtered_feature_bc_matrix'), project = "AI2")
CI1 <- CreateSeuratObject(Read10X('/data/hxl/CI1/filtered_feature_bc_matrix'), project = "CI1")
CI2 <- CreateSeuratObject(Read10X('/data/hxl/CI2/filtered_feature_bc_matrix'), project = "CI2")
all.set <- merge(x = CT1, y = c(CT2, AI1, AI2, CI1, CI2),add.cell.ids = c('CT1', 'CT2', 'AI1', 'AI2', 'CI1', 'CI2'), project = 'ALLset')
all.set[["percent.mt"]] <- PercentageFeatureSet(all.set, pattern = "^mt-")
all.set <- subset(all.set, percent.mt < 10 & nFeature_RNA > 200)
VlnPlot(all.set, features = c("nFeature_RNA", "percent.mt"))
#batch effect
all <- NormalizeData(all.set, verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = all.set@var.genes, npcs = 50, verbose = FALSE)
p1 <- DimPlot(object = all, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = all, features = "PC_1", group.by = "orig.ident", pt.size = .1)
p1 + p2
all <- all %>% RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(all, 'harmony')
harmony_embeddings[1:5, 1:5]
DimPlot(object = all, reduction = "harmony", pt.size = .1, group.by = "orig.ident") +
  VlnPlot(object = all, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
all <- all %>% 
  RunUMAP(reduction = "harmony", dims = 1:35) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:35) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
DimPlot(all, group.by = 'orig.ident') + scale_color_brewer(palette = 'Paired') + theme_bw()
DimPlot(all, label = TRUE)
save(all, file = "/data/hxl/all.Rdata")
p <- DimPlot(all, label = T) + theme_bw() + theme(legend.position = 'none')
ggsave('UMAP_all_clusters.png', p, width=7, height=5)
p <- DimPlot(all, group.by = 'orig.ident') + scale_color_brewer(palette = 'Paired') + theme_bw()
ggsave('UMAP_all_samples.png', p, width=7, height=5)
#find markers
cluster_de <- FindAllMarkersMAESTRO(all, only.pos = T)
write.csv(cluster_de, file='cluster_SdInt_de.csv', quote = F, row.names = F)

all <- RNAAnnotateCelltype(all, genes = cluster_de)
DimPlot(all, group.by = 'assign.ident', label = T) + NoLegend()
#remove red cells and re-run
no_red <- all[,!all@meta.data[["seurat_clusters"]] %in% c ("13") ]
no_red <- NormalizeData(no_red, verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = no_red@var.genes, npcs = 50, verbose = FALSE)
p1 <- DimPlot(object = no_red, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p1
all_no_red <- no_red %>% RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(all_no_red, 'harmony')
harmony_embeddings[1:5, 1:5]
DimPlot(object = no_red, reduction = "harmony", pt.size = .1, group.by = "orig.ident") +
  VlnPlot(object = no_red, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
all_no_red <- all_no_red %>% 
  RunUMAP(reduction = "harmony", dims = 1:35) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:35) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
DimPlot(all_no_red, group.by = 'orig.ident') + scale_color_brewer(palette = 'Paired') + theme_bw()
DimPlot(all_no_red, label = TRUE)
save(all_no_red, file = "/data/hxl/all_no_red.Rdata")

#find markers
cluster_nored <- FindAllMarkersMAESTRO(all_no_red, only.pos = TRUE)
write.csv(cluster_nored, file='cluster_no_red.csv', quote = FALSE, row.names = FALSE)
all_no_red <- RNAAnnotateCelltype(all_no_red, genes = cluster_nored)
DimPlot(all_no_red, group.by = 'assign.ident', label = TRUE) + NoLegend()
#only immune cells
immune <- all_no_red[,!all_no_red@meta.data[["assign.ident"]] %in% c ("Fibroblasts") ]
immune <- NormalizeData(immune, verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = immune@var.genes, npcs = 50, verbose = FALSE)
p1 <- DimPlot(object = immune, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p1
immune <- immune %>% RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(immune, 'harmony')
harmony_embeddings[1:5, 1:5]
DimPlot(object = immune, reduction = "harmony", pt.size = .1, group.by = "orig.ident") +
  VlnPlot(object = immune, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
immune <- immune %>% 
  RunUMAP(reduction = "harmony", dims = 1:35) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:35) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
DimPlot(immune, group.by = 'orig.ident') + scale_color_brewer(palette = 'Paired') + theme_bw()
DimPlot(immune, label = TRUE)
save(immune, file = "/data/hxl/Rdata/immune.Rdata")
#find markers
cluster_immune <- FindAllMarkersMAESTRO(immune, only.pos = TRUE)
cluster_immune_label <- FindAllMarkers(immune, only.pos = TRUE, group.by = 'Label')
write.csv(cluster_immune_label, file='cluster_immune_label.csv', quote = FALSE, row.names = FALSE)
write.csv(cluster_immune, file='cluster_immune.csv', quote = FALSE, row.names = FALSE)
cluster_immune <- read.csv('/data/hxl/csv/cluster_immune.csv')
cluster_immune_label <- read.csv('/data/hxl/csv/cluster_immune_label.csv')
immune <- RNAAnnotateCelltype(immune, genes = cluster_immune)
DimPlot(immune, group.by = 'assign.ident', label = TRUE) + NoLegend()
DimPlot(immune, group.by = 'seurat_clusters', label = TRUE)
immune <- immune[,!immune@meta.data[["assign.ident"]] %in% c ("Fibroblasts") ]
immune <- immune[,!immune@meta.data[["seurat_clusters"]] %in% c ("22") ]
immune_dotplot <- DotPlot(immune, features = c("Cd3g", 'Cd3e', "Cd3d", "Cd4", "Cd8a", "Klrb1c", "Klrk1", "Klrg1", "Ncr1", "Klrc2",
                                               "Cd79a", "Cd79b", "Cd19", "Fcmr", "Ms4a1", "Jchain", "Sdc1", "Cd81",
                                               "Ly6g", "S100a8", "Mmp8", "Mmp9", 
                                               "Gngt2", "Cd68", "Lyz2", 
                                               "Cpa3", "Mcpt8", "Il18r1",
                                               "Mki67", "Stmn1", "Cdc20", "Birc5", "Ube2c"), group.by = 'Label', cols = c('gray', 'Tomato')) +
  RotatedAxis()+
  theme(axis.text.x = element_text(family = 'Helvetica', size=18, color = 'black', angle=45, hjust=1),
        axis.text.y = element_text(family = 'Helvetica', size = 18, color = 'black'),
        axis.title.x = element_blank(), axis.title.y.left = element_blank())+
  theme(panel.grid.major=element_line(colour='gray'))+
  guides(size = guide_legend(title = "% of Exp")) + guides(color = guide_colorbar(title = "Avg. Exp"))+
  theme(legend.title = element_text(family = 'Helvetica', size=18, color = 'black'))+
  annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 1)
ggsave("immune_dotplot_new_1.pdf", immune_dotplot, width = 15, height = 5)
FeaturePlot(immune, c("Cd3e", "Cd79a", "Jchain", "Klrb1c", "Mki67", "Ly6g", "Hdc", "Cd68", "Cpa3"))
DimPlot(immune, group.by = "seurat_clusters", label = T)  
FeaturePlot(immune, c("Cd3d", 'Cd4', 'Cd8a', "Klrb1c", "Ly6g", "Cd79a", "Mki67", "Cd68", "Jchain", "Cpa3"),
            cols = c("Lightgray", "blue")) + DimPlot(immune, group.by = "seurat_clusters", label = T) 
ggsave("immune_featureplot.pdf", immune_featureplot, width = 10, height = 5)

immune <- RenameIdents(immune, '0' = 'CD8+ T', '1' = 'CD8+ T', '2' = 'NK', '3' = 'CD4+ T', '4' = 'Neutrophils', '5' = 'B', 
                       '6' = 'Neutrophils', '8' = 'Neutrophils', '7' = 'CD4+ T', '9' = 'Proliferating cells', '10' = 'Monocytes',
                       '11' = 'Monocytes', '12' = 'Monocytes', '13' = 'Monocytes', '14'= 'Neutrophils', '15' = 'CD4+ T', '16' = 'CD8+ T', 
                       '17' = 'B', '18' = 'NKT', '19' = 'Plasma', '20' = 'Basophils', '21' = 'Neutrophils')
immune <- RenameIdents(immune, 'Proliferating cells' = 'Proliferating')
immune$sample <- immune$orig.ident
table(immune$sample)
str(immune@meta.data)
immune$sample <- as.factor(immune$orig.ident)
str(immune@meta.data)
levels(immune$sample)
levels(immune$sample) <- c('Acute1', 'Acute2', 'Chronic1', 'Chronic2', 'Ctrl1','Ctrl2')
str(immune@meta.data)
immune$sample <- factor(immune$sample, levels = c('Ctrl1', 'Ctrl2', 'Acute1', 'Acute2', 'Chronic1', 'Chronic2'))
colourCount = length(unique(immune$Label))
getPalette = colorRampPalette(brewer.pal(8, "Set2"))
immune_umap <- DimPlot(immune, group.by = 'Label', label = TRUE, label.size = 6.5, repel =T, cols = brewer.pal(12, "Set3")[c(1,3,4,5,6,7,8,9,10,11,12)])+
  annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(legend.position = "none")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+#移除网格线
  xlab("UMAP 1") + ylab("UMAP 2")+
  theme(axis.text.x = element_text(size = 24, family = "Helvetica", color = 'black'), axis.text.y = element_text(size = 24, family = "Helvetica", color = 'black'))+
  theme(axis.title.x = element_text(size = 24, family = "Helvetica", color = 'black'), axis.title.y = element_text(size = 24, family = "Helvetica", color = 'black'))+
  theme(plot.title = element_blank())

ggsave("immune_umap.pdf", immune_umap, width = 5.6, height = 5)
ggsave("immune_umap.jpg", immune_umap, width = 5.6, height = 5, dpi = 300)

immune_umap_nolabel <- DimPlot(immune, group.by = 'Label', label = F, label.size = 6.5, repel =T, cols = brewer.pal(12, "Set3")[c(1,3,4,5,6,7,8,9,10,11,12)]) +
  annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 1)+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(legend.position = "none")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+#移除网格线
  xlab("UMAP 1") + ylab("UMAP 2")+
  theme(axis.text.x = element_text(size = 26, family = "Helvetica", color = 'black'), axis.text.y = element_text(size = 26, family = "Helvetica", color = 'black'))+
  theme(axis.title.x = element_text(size = 26, family = "Helvetica", color = 'black'), axis.title.y = element_text(size = 26, family = "Helvetica", color = 'black'))+
  theme(plot.title = element_blank())
ggsave("immune_umap_nolabel.pdf", immune_umap_nolabel, width = 5.5, height = 5) 
ggsave("immune_umap_nolabel.jpg", immune_umap_nolabel, width = 5.5, height = 5, dpi = 300)
immune_sample <- DimPlot(immune, group.by = "sample", label = F, cols = brewer.pal(12, "Paired")[c(3,4,9,10,11,12)])+
  annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 1)+
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.text.x = element_text(size = 24, family = "Helvetica", color = "black"), axis.text.y = element_text(size = 24, family = "Helvetica", color = "black"), title = element_blank())+
  theme(axis.title.x = element_text(size = 24, family = "Helvetica"), axis.title.y = element_text(size = 24, family = "Helvetica"))+
  theme(legend.text = element_text(size = 20, family = "Helvetica"))+
  theme(plot.title = element_blank())+
  guides(color=guide_legend(override.aes = list(size=6)))
ggsave("immune_sample.pdf", immune_sample, width = 7.25, height = 5)
ggsave("immune_sample.jpg", immune_sample, width = 7.25, height = 5, dpi = 300)
immune_umap + immune_sample
ggsave("immune_sample_umap.jpg", immune_umap + immune_sample, width =17, height = 7)
immune_group <- DimPlot(immune, group.by = "group", split.by = "group", ncol = 3) +
  scale_color_manual(values = c('Olive Drab', 'Pink 1', 'Steel Blue 4'))+
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.text.x = element_text(size = 20, family = "Helvetica", face = "bold"), axis.text.y = element_text(size = 20, family = "Helvetica", face = "bold"), title = element_blank())+
  theme(axis.title.x = element_text(size = 20, family = "Helvetica"), axis.title.y = element_text(size = 20, family = "Helvetica"))+
  theme(legend.text = element_text(size = 18, family = "Helvetica"))+
  theme(plot.title = element_blank())
ggsave("immune_group.pdf", immune_group, width = 14, height = 5)
ggsave("immune_group.jpg", immune_group, width = 14, height = 5, dpi = 300)

DimPlot(immune,label = T, label.size = 5.5, repel = T, split.by = 'group') + scale_color_brewer(palette = "Set3") +
  annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 1)+
  theme_bw(base_size = 20, base_family = "Helvetica")+
  theme(legend.position = "none")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+#移除网格线
  xlab("UMAP 1") + ylab("UMAP 2")+
  theme(axis.text.x = element_text(size = 20, family = "Helvetica", face = "plain", color = 'black'), axis.text.y = element_text(size = 20, family = "Helvetica", face = "plain", color = 'black'))+
  theme(axis.title.x = element_text(size = 20, family = "Helvetica", color = 'black'), axis.title.y = element_text(size = 20, family = "Helvetica", color = 'black'))

ggsave("immune_group_umap_nolabel.pdf", width = 12, height = 4.5)
ggsave("immune_group_umap_nolabel.jpg", width = 12, height = 4.5, dpi = 300)
save(immune, file = "/data/hxl/Rdata/immune.Rdata")

#cell proportion
immune$Label <- Idents(immune)
immune$Label <- factor(immune$Label, levels = c('CD4+ T', 'CD8+ T', 'NK', 'NKT', 'B', 'Plasma',  'Neutrophils', "Monocytes", 'Basophils', 'Proliferating'))
head(immune@meta.data)
table(immune$orig.ident, immune$Label)
t.immune=as.data.frame.matrix(table(immune$orig.ident, immune$Label))
n_sample_immune <- ddply(immune@meta.data, .(orig.ident, Label), nrow)
n_sample_immune <- ddply(n_sample_immune, .(orig.ident), transform, percent=V1/sum(V1)*100)
write.csv(n_sample_immune, file = "/data/hxl/csv/n_sample_immune.csv")
lv <- c("CT1", "CT2", "AI1", "AI2", "CI1", "CI2")
n_sample_immune$orig.ident <- factor(as.character(n_sample_immune$orig.ident), levels = lv)

pn_sample_immune <- ggplot(n_sample_immune, aes(x=orig.ident, y=percent, fill=Label)) +
  labs(x = '', y='Percentage', fill='Cell Type')+
  geom_bar(stat="identity", colour="black") +
  scale_fill_manual(values = brewer.pal(12, "Set3")[c(1,3,4,5,6,7,8,9,10,11,12)])+
  theme_bw(base_size = 22, base_family = "Helvetica")+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5,color="black", family = "Helvetica")) +
  scale_x_discrete(labels = c("Ctrl1", "Ctrl2", "Acute1", "Acute2", "Chronic1", "Chronic2"))+
  theme(axis.text.y = element_text(family = "Helvetica", size = 22, color = "black"))+
  theme(legend.title = element_text(family = "Helvetica", size = 22, face = "plain", color = "black"))
ggsave('immune_cell_percentage.pdf', pn_sample_immune, width = 7.2, height = 5)
save(immune, file = "/data/hxl/Rdata/immune.Rdata")

#reanalyse T and NK
T_NK <- immune[,immune$seurat_clusters %in% c("0", "1","2", "3","8", "9", "14", "16", "18")]
T_NK<- NormalizeData(T_NK, verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = T_NK@var.features, npcs = 50, verbose = FALSE)
T_NK <- T_NK %>% RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(T_NK, 'harmony')
harmony_embeddings[1:5, 1:5]
T_NK <- T_NK %>% 
  RunUMAP(reduction = "harmony", dims = 1:35) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:35) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
DimPlot(T_NK, group.by = 'orig.ident') + scale_color_brewer(palette = 'Paired') + theme_bw()
DimPlot(T_NK, label = TRUE)
save(T_NK, file = "/data/hxl/T_NK.Rdata")
cluster_T_NK <- FindAllMarkersMAESTRO(T_NK, only.pos = TRUE)
write.csv(cluster_NK, file='cluster_T_NK.csv', quote = FALSE, row.names = FALSE)
cluster_T_NK_label <- FindAllMarkers(T_NK, only.pos = TRUE, group.by = 'Label')
write.csv(cluster_T_NK_label, file='/data/hxl/csv/cluster_T_NK_label.csv', quote = FALSE, row.names = FALSE)
FeaturePlot(T_NK, c("Cd3e", "Cd4", "Cd8a", "Cxcr3", "Cxcr6", "Rora", "Gzmb", "Sell", "Ccr7", "Foxp3")) + DimPlot(T_NK, group.by = "seurat_clusters", label = TRUE)
T_NK <- RenameIdents(T_NK, "0" = "Te", "1" = "NK", "2" = "CD8+ Naive", "3" = "CD4+ Naive", "4" = "Trm", "5" = "Mixed T", "6" = "Treg", "7" = "Proliferating CD8+ T", "8" = "Proliferating CD8+ T", "9" = "ILC1",
                     '10' = 'ILC3', '11' = 'DP T', '12' = 'Te', '13' = 'Esm+ T', '14' = 'NKT', "15" = 'NK',
                     '16' = 'Il1b+ T', '17' = 'ILC2', '18' = 'MHCII+ lymphocytes')
T_NK <- RenameIdents(T_NK, 'CD8+ Naive' = 'CD8+ Naive T', 'CD4+ Naive' = 'CD4+ Naive T', 'MHCII+ lymphocytes' = 'MHCII+ Lymphocytes', 'Te' = 'Teff')
T_NK <- RenameIdents(T_NK, 'Proliferating CD8+ T' = 'Proliferating T')
DimPlot(T_NK, label = TRUE, label.size = 3)
colourCount = length(unique(T_NK$Label))
getPalette = colorRampPalette(brewer.pal(12, "Set2"))
T_NK_umap <- DimPlot(T_NK,label = TRUE, label.size = 5.5, repel = T) + scale_color_manual(values = getPalette(colourCount)) +
  annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 1)+
  theme_bw(base_size = 20, base_family = "Helvetica")+
  theme(legend.position = "none")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+#移除网格线
  xlab("UMAP 1") + ylab("UMAP 2")+
  theme(axis.text.x = element_text(size = 20, family = "Helvetica", face = "plain", color = 'black'), axis.text.y = element_text(size = 20, family = "Helvetica", face = "plain", color = 'black'))+
  theme(axis.title.x = element_text(size = 20, family = "Helvetica", color = 'black'), axis.title.y = element_text(size = 20, family = "Helvetica", color = 'black'))
ggsave("T_NK_umap.pdf", T_NK_umap, width = 5.5, height = 5)
T_NK_umap_s <- DimPlot(T_NK,label = TRUE, label.size = 5.5, repel = T) + scale_color_manual(values = getPalette(colourCount)) +
  theme(legend.position = "none")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+#移除网格线
  xlab("UMAP 1") + ylab("UMAP 2")+
  theme(axis.text.x = element_text(size = 20, family = "Helvetica", face = "plain", color = 'black'), axis.text.y = element_text(size = 20, family = "Helvetica", face = "plain", color = 'black'))+
  theme(axis.title.x = element_text(size = 20, family = "Helvetica", color = 'black'), axis.title.y = element_text(size = 20, family = "Helvetica", color = 'black'))
ggsave("T_NK_umap_s.pdf", T_NK_umap_s, width = 5.5, height = 5)
T_NK_group <- DimPlot(T_NK, group.by = "group", label = F)+
  scale_color_manual(values = c("Orange", "Yellow Green", "Steel Blue")) + 
  annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 1)+
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.text.x = element_text(size = 20, family = "Helvetica"), axis.text.y = element_text(size = 20, family = "Helvetica"), title = element_blank())+
  theme(axis.title.x = element_text(size = 20, family = "Helvetica"), axis.title.y = element_text(size = 20, family = "Helvetica"))+
  theme(legend.text = element_text(size = 18, family = "Helvetica"))+
  theme(plot.title = element_blank())
ggsave("T_NK_group.pdf", T_NK_group, width = 6.8, height = 5)
T_NK_group_split <- DimPlot(T_NK, group.by = "group", split.by = 'group', label = F)+
  scale_color_manual(values = c("Orange", "Yellow Green", "Steel Blue")) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.text.x = element_text(size = 20, family = "Helvetica"), axis.text.y = element_text(size = 20, family = "Helvetica"), title = element_blank())+
  theme(axis.title.x = element_text(size = 20, family = "Helvetica"), axis.title.y = element_text(size = 20, family = "Helvetica"))+
  theme(legend.text = element_text(size = 18, family = "Helvetica"))+
  theme(plot.title = element_blank())
ggsave("T_NK_group_split.pdf", T_NK_group_split, width = 14, height = 5)
DimPlot(T_NK,label = TRUE) + scale_color_brewer(palette = "Set3")
FeaturePlot(T_NK, c("Cd4", "Cd8a", "Ccr7", "Sell")) 
FeaturePlot(T_NK, c("Cd4", "Cd8a", "Cd7", "Klrb1c", "Cxcr3", "Il17a", "Il17f", "Ahr", "Rora", "Gata3", "Arg1"))
FeaturePlot(T_NK, c("Cd4", "Cd8a", "Mki67", "Stmn1", "Foxp3", "Ctla4", "Nrp1")) 
FeaturePlot(T_NK, c("Cd4", "Cd8a", "Gzmb", "Klrg1", "Cd44", "Cxcr3", "Cxcr6", "Cd69", "Ccl5"))
FeaturePlot(T_NK, c("Cd4", "Cd8a", "Il1b", "Esm1", "Cd79a", "H2-Aa"))
FeaturePlot(T_NK, c("Cd4", "Cd8a", "Id3", "Il10"))
T_NK_Feature <- FeaturePlot(T_NK, c( "Cd3e","Cd4", "Cd8a", "Ccr7", "Sell", "Cd7", "Klrb1c", "Cxcr3",
                                     "Il17a", "Il17f", "Ahr", "Rora", "Gata3", "Arg1", "Mki67", 
                                     "Stmn1", "Foxp3", "Ctla4", "Nrp1", "Gzmb", "Klrg1", "Cd44",
                                     "Cxcr6", "Cd69", "Il1b", "Esm1", "Cd79a", "H2-Aa", "Id3", "Il10"),
                            cols = c("lightgray", "Orange red"), ncol = 5)
ggsave('t_nk_featureplot.pdf',T_NK_Feature, width = 12, height = 15)
ggsave('t_nk_featureplot.png',T_NK_Feature, width = 15.5, height = 15)
DotPlot(T_NK, features = c("Cd4", "Cd8a", "Ccr7", "Sell", "Cd7", "Klrb1c", "Cxcr3", "Il17a", "Il17f", "Ahr", "Rora", "Gata3", "Arg1",
                           "Mki67", "Stmn1", "Foxp3", "Ctla4", "Nrp1", "Gzmb", "Klrg1", "Cd44", "Cxcr6", "Cd69", 
                           "Il1b", "Esm1", "Cd79a", "H2-Aa", "Id3", "Il10")) + RotatedAxis()
dir.create("./T_NK_Features")
setwd("./T_NK_Features")
pdf("T_NK_Dotplot.pdf", width = 12.5, height = 5.5)
T_NK_dotplot <- DotPlot(T_NK, features = c("Cd3e", "Cd8a", "Mki67", "Stmn1", "Ccr7", "Sell", "Cd4", "Cd79a", "H2-Aa",
                                           "Gzmb", "Klrg1", "Klrb1c", "Ncr1", "Cxcr3", "Cxcr6", "Cd69", "Foxp3", "Ctla4", "Nrp1",
                                           "Cd7", "Il17a", "Il17f", "Ahr", "Esm1", "Il1b", "Rora", "Gata3", "Arg1"),
                        cols = c("lightgray", "salmon")) + RotatedAxis() +
  theme(
    axis.text.x = element_text(family = 'Helvetica', size = 18, angle=45, hjust=1, color = 'black'),
    axis.text.y = element_text(family = 'Helvetica', size = 18, face = 'plain'))+
  annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 1)+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
  theme(panel.grid.major=element_line(colour='gray'))+
  guides(size = guide_legend(title = "% of Exp")) + guides(color = guide_colorbar(title = "Avg. Exp"))+
  theme(legend.title = element_text(family = 'Helvetica', size=18, color = 'black'))
garbage <- dev.off()
ggsave("T_NK_Dotplot.pdf", T_NK_dotplot, width = 12.5, height = 5.5)
DimPlot(T_NK, group.by = "seurat_clusters", label = T) + NoLegend()
DimPlot(T_NK, label = T) + scale_color_manual(values = c('Light Coral', 'Orange', 'Olive Drab', 'Coral','Dark Cyan',
                                                         'Light Green', 'Sky Blue', 'Pink', 'Steel Blue',
                                                         'Plum', 'Dark Slate Blue', 'Tomato 1', 'Slate Blue',
                                                         'Gray', 'Indian Red', 'Slate Gray 1')) 
#cell proportion
T_NK$Label <- Idents(T_NK)
head(T_NK@meta.data)
table(T_NK$orig.ident, T_NK$Label)
t.T_NK=as.data.frame.matrix(table(T_NK$orig.ident, T_NK$Label))
n_sample_T_NK <- ddply(T_NK@meta.data, .(orig.ident, Label), nrow)
n_sample_T_NK <- ddply(n_sample_T_NK, .(orig.ident), transform, percent=V1/sum(V1)*100)
write.csv(n_sample_T_NK, file = "/data/hxl/csv/n_sample_T_NK.csv")
lv <- c("CT1", "CT2", "AI1", "AI2", "CI1", "CI2")
n_sample_T_NK$orig.ident <- factor(as.character(n_sample_T_NK$orig.ident), levels = lv)
colourCount = length(unique(T_NK$Label))
getPalette = colorRampPalette(brewer.pal(12, "Set3"))
ggplot(n_sample_T_NK, aes(x=orig.ident, y=percent, fill=Label)) +
  labs(x='', y='Percentage', fill='Cell Type')+
  geom_bar(stat="identity", colour="black") +
  scale_fill_manual(values = getPalette(colourCount))+
  theme_bw(base_size = 20, base_family = "Helvetica")+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5,color="black")) + 
  scale_x_discrete(labels = c("Ctrl1", "Ctrl2", "Acute1", "Acute2", "Chronic1", "Chronic2"))+
  theme(axis.text.y = element_text(family = "Helvetica", size = 20, color = "black"))+
  theme(legend.title = element_text(family = "Helvetica", size = 20, color = "black")) 
ggsave('T_NK_cell_percentage.pdf', width = 8, height = 6)
save(T_NK, file = "/data/hxl/T_NK.Rdata")

#reanalyse B
B <- immune[,immune$seurat_clusters %in% c("5","16", "19", "18")]
B <- NormalizeData(B, verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = B@var.features, npcs = 50, verbose = FALSE)
B <- B %>% RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(B, 'harmony')
harmony_embeddings[1:5, 1:5]
B <- B %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.2) %>% 
  identity()
DimPlot(B, group.by = 'orig.ident') + scale_color_brewer(palette = 'Paired') + theme_bw()
DimPlot(B, label = TRUE)
save(B, file = "/data/hxl/B.Rdata")
cluster_B <- FindAllMarkersMAESTRO(B, only.pos = TRUE)
write.csv(cluster_B, file='cluster_B.csv', quote = FALSE, row.names = FALSE)
B <- RNAAnnotateCelltype(B, genes = cluster_B)
DimPlot(B, group.by = 'assign.ident', label = TRUE) + NoLegend()


cluster_B2 <- FindAllMarkers(B, only.pos = TRUE)
write.csv(cluster_B2, file='cluster_B2.csv', quote = FALSE, row.names = FALSE)
B <- B[,!B@meta.data$seurat_clusters %in% c("6")]
B <- RenameIdents(B, "0" = "Cxcl2+ B", "1" = "Marginal zone B", "2" = "Memory B", "3" = "Gzmk+ B", "4" = "Plasma B", "5" = "Il7r+ B", "7" = "Follicular B")
B <- RenameIdents(B, "Marginal zone B" = "MZ B")
B <- RenameIdents(B, "Plasma B" = "Plasma")
DimPlot(B, label = TRUE)
DimPlot(B, group.by = "orig.ident", label = F)
B_umap <- DimPlot(B,label = TRUE, repel = T, label.size = 7.5) + scale_color_brewer(palette = "Set2") +
  annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 1)+
  theme(legend.position = "none")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+#移除网格线
  xlab("UMAP 1") + ylab("UMAP 2")+
  theme(axis.text.x = element_text(size = 24, family = "Helvetica", face = "plain", color = 'black'), axis.text.y = element_text(size = 24, family = "Helvetica", face = "plain", color = 'black'))+
  theme(axis.title.x = element_text(size = 24, family = "Helvetica", color = 'black'), axis.title.y = element_text(size = 24, family = "Helvetica", color = 'black'))
ggsave("B_umap.pdf", B_umap, width = 5.5, height = 5)
B_group <- DimPlot(B,group.by = "group", label = T, repel = T, label.size = 7) + scale_color_manual(values = c('Yellow Green', 'Pink 3', 'Steel Blue 4')) +
  annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)+
  theme_bw(base_size = 20, base_family = "Helvetica")+
  theme(legend.position = "none")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+#移除网格线
  xlab("UMAP 1") + ylab("UMAP 2")+
  theme(axis.text.x = element_text(size = 20, family = "Helvetica", face = "bold", color = 'black'), axis.text.y = element_text(size = 20, family = "Helvetica", face = "bold", color = 'black'))+
  theme(axis.title.x = element_text(size = 20, family = "Helvetica", color = 'black'), axis.title.y = element_text(size = 20, family = "Helvetica", color = 'black'))+
  theme(plot.title = element_blank())
ggsave("B_group.pdf", B_group, width = 5.5, height = 5)
ggsave('b_featureplot.pdf', width = 11, height = 4.5)
ggsave('b_featureplot.png', width = 11, height = 4.5, dpi = 600)
DimPlot(B, group.by = 'group', label = F, cols = c("Olive Drab", "Pink 3", "Steel Blue 4"))
VlnPlot(B, features = c("Cxcl2", "Cr2", "Spib", "Gzmk", "Jchain", "Il7r", "Fcer1g"), group.by = "Label", pt.size = 0) 
B$sample <- as.factor(B$orig.ident)
str(B@meta.data)
levels(B$sample)
levels(B$sample) <- c('Acute1', 'Acute2', 'Chronic1', 'Chronic2', 'Ctrl1','Ctrl2')
str(B@meta.data)
B$sample <- factor(B$sample, levels = c('Ctrl1', 'Ctrl2', 'Acute1', 'Acute2', 'Chronic1', 'Chronic2'))
#cell proportion
B$Label <- Idents(B)
head(B@meta.data)
table(B$orig.ident, B$Label)
t.B=as.data.frame.matrix(table(B$orig.ident, B$Label))
n_sample_B <- ddply(B@meta.data, .(orig.ident, Label), nrow)
n_sample_B <- ddply(n_sample_B, .(orig.ident), transform, percent=V1/sum(V1)*100)
write.csv(n_sample_B, file = "/data/hxl/csv/n_sample_B.csv")
lv <- c("CT1", "CT2", "AI1", "AI2", "CI1", "CI2")
n_sample_B$orig.ident <- factor(as.character(n_sample_B$orig.ident), levels = lv)
ggplot(n_sample_B, aes(x=orig.ident, y=percent, fill=Label)) +
  labs(x = '', y='Percentage', fill='Cell Type')+
  geom_bar(stat="identity", colour="black") +
  scale_fill_brewer(palette="Set2")+
  theme_bw(base_size = 22, base_family = "Helvetica")+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5,color="black", family = "Helvetica")) +
  scale_x_discrete(labels = c("Ctrl1", "Ctrl2", "Acute1", "Acute2", "Chronic1", "Chronic2"))+
  theme(axis.text.y = element_text(family = "Helvetica", size = 22, color = "black"))+
  theme(legend.title = element_text(family = "Helvetica", size = 22, face = "plain", color = "black"))
ggsave('B_cell_percentage.pdf', width = 7, height = 5.5)
save(B, file = "/data/hxl/Rdata/B.Rdata")
#reanalyse myeloid
myeloid <- immune[,immune$seurat_clusters %in% c("10","11", "12", "7", "4", "6", "13", "20", "21", "23")]
myeloid <- NormalizeData(myeloid, verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = myeloid@var.features, npcs = 50, verbose = FALSE)
myeloid <- myeloid %>% RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(myeloid, 'harmony')
harmony_embeddings[1:5, 1:5]
myeloid <- myeloid %>% 
  RunUMAP(reduction = "harmony", dims = 1:35) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:35) %>% 
  FindClusters(resolution = 1.0) %>% 
  identity()
DimPlot(myeloid, group.by = 'orig.ident') + scale_color_brewer(palette = 'Paired') + theme_bw()
DimPlot(myeloid, label = TRUE)
save(myeloid, file = "/data/hxl/myeloid.Rdata")
cluster_myeloid <- FindAllMarkersMAESTRO(myeloid, only.pos = TRUE)
write.csv(cluster_myeloid, file='cluster_myeloid.csv', quote = FALSE, row.names = FALSE)
myeloid <- RNAAnnotateCelltype(myeloid, genes = cluster_myeloid)
DimPlot(myeloid, group.by = 'assign.ident', label = TRUE) + NoLegend()
#remove high "Rp-" and "Hbb-" cluster
myeloid[["percent.Hbb"]] <- PercentageFeatureSet(myeloid, pattern = "^Hbb-")
VlnPlot(myeloid, features = c("percent.Hbb"), group.by = "seurat_clusters")
myeloid[["percent.Rp"]] <- PercentageFeatureSet(myeloid, pattern = "Rp-")
VlnPlot(myeloid, features = c("percent.Rp"), group.by = "seurat_clusters")
myeloid <- myeloid[, !myeloid$seurat_clusters %in% c("9")]
myeloid <- NormalizeData(myeloid, verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = myeloid@var.features, npcs = 50, verbose = FALSE)
myeloid <- myeloid %>% 
  RunUMAP(reduction = "harmony", dims = 1:35) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:35) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
cluster_myeloid <- FindAllMarkersMAESTRO(myeloid, only.pos = TRUE)
write.csv(cluster_myeloid, file='cluster_myeloid.csv', quote = FALSE, row.names = FALSE)
FeaturePlot(myeloid, c("Mrc1", "Ear1", "S100a8", "Stmn1", "Cd14", "C1qa", "H2-Aa",
                       "Itgax", "Itgam", "Mertk", "F13a1", "Adgre4")) + 
  DimPlot(myeloid, group.by = "seurat_clusters", label = TRUE)
FeaturePlot(myeloid, c("Nkg7", "Gzmb", "Ms4a4b", "Ebf1", "Ly6d", "Cpa3", "Mcpt8")) 
FeaturePlot(myeloid, c("Itgax", "H2-Aa", "C1qa", "C1qc", "Itgam", "Mertk"))
FeaturePlot(myeloid, c("Itgax", "Siglecf", "Stmn1", "Ccr2", "Adgre4", "Itgam", "Bst2", "Siglech"))
FeaturePlot(myeloid, c("Sell", "Ly6c2", "Ccr2", "Spn", "Itgax", "Siglecf", "Stmn1", "Adgre4"))

myeloid$Label <- Idents(myeloid)

myeloid$Label <- factor(myeloid$Label, levels = c("Neutrophils", "Basophils", "Stmn+ AMs", "Stmn- AMs",
                                                  "IM1", "IM2", "cDCs", "Inflammatory monocytes",
                                                  "Nkg+ Monocytes", "Ly6d+ Monocytes"))
myeloid_dotplot <- DotPlot(myeloid, features = c("Mrc1", "Ear1", "Ear2", "Sell", "Ly6c2", "Ccr2",
                                                 "Ly6g", "Mmp8", 'Mmp9', "Adgre4", "Nkg7", 
                                                 "Stmn1", "Birc5", "C1qa", "C1qb", "H2-Aa","Mertk", "Ly6d",
                                                 "Cpa3", "Mcpt8", "Gata2"),dot.scale = 8) +
  RotatedAxis() + theme(text = element_text(family = 'Helvetica', size=18),
                        axis.text.x = element_text(family = 'Helvetica', size=18, angle=45, hjust=1),
                        axis.text.y = element_text(family = 'Helvetica', size = 18))+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
  annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 1)+
  theme(panel.grid.major=element_line(colour='gray'))+
  guides(size = guide_legend(title = "% of Exp")) + guides(color = guide_colorbar(title = "Avg. Exp"))+
  theme(legend.title = element_text(family = 'Helvetica', size=18, color = 'black'))+
  scale_colour_viridis_c(
    alpha = 1,
    begin = 0.3,
    end = 0.9,
    direction = 1,
    option ="D",
    aesthetics = "colour",
    na.value = "grey50",
    guide = "colourbar"
  )
ggsave("myeloid_dotplot_new_2.pdf", myeloid_dotplot, width = 12.5, height = 5)
myeloid_featureplot <- FeaturePlot(myeloid, c("Mrc1", "Ear1", "Sell", "Ly6c2", "Ccr2",
                                              "Ly6g", "Hdc", 'Mmp9', "Adgre4", "Nkg7", 
                                              "Stmn1", "Birc5", "C1qa", "C1qb", "H2-Aa","Mertk", "Ly6d",
                                              "Cpa3", "Mcpt8", "Gata2"), cols = c("lightgray", "Tomato"))
ggsave(myeloid_featureplot, filename = 'myeloid_featureplot.pdf', width = 10, height = 9.3)
myeloid <- RenameIdents(myeloid, '0' = 'Neutrophils', '1' = 'Neutrophils', '2' = 'Neutrophils', '3' = 'AMs',
                        '4' = 'pDCs', '5' = 'Neutrophils', '6' = 'cDCs', '7' = 'Nkg7+ Monocytes', '8' = 'Neutrophils', '9' = 'Stmn+ AMs', '10' = 'IM2',
                        '11' = 'IM1', '12' = 'Ly6d+ Monocytes', '13' = 'Nkg7+ Monocytes', '14' = 'Basophils')
myeloid <- RenameIdents(myeloid, 'AMs(alveolar macrophages)' = 'AMs', 'IM1(Type 1 interstitial macrophages)' = 'IM1', 'IM2(Type 2 interstitial macrophages)' = 'IM2',
                        'PDCs' = 'pDCs', 'Nkg7+ monocytes' = 'Nkg7+ Monocytes', 'Stmn+ AMs(alveolar macrophages)' = 'Stmn+ AMs',
                        'Ly6d+ monocytes' = 'Ly6d+ Monocytes')
myeloid <- RenameIdents(myeloid, 'pDCs' = 'Inflammatory Monocytes')
myeloid <- RenameIdents(myeloid, 'AMs' = 'Stmn- AMs')
myeloid_umap <- DimPlot(myeloid, label = T, label.size = 6, repel = T) + scale_color_brewer(palette = "Set3") +
  annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)+
  theme_bw(base_size = 24, base_family = "Helvetica")+
  theme(legend.position = "none")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+#移除网格线
  xlab("UMAP 1") + ylab("UMAP 2")+
  theme(axis.text.x = element_text(size = 24, family = "Helvetica", face = "plain", color = 'black'), axis.text.y = element_text(size = 24, family = "Helvetica", face = "plain", color = 'black'))+
  theme(axis.title.x = element_text(size = 24, family = "Helvetica", color = 'black'), axis.title.y = element_text(size = 24, family = "Helvetica", color = 'black'))+
  theme(plot.title = element_blank())
ggsave("myeloid_umap.pdf", myeloid_umap, width = 5.5, height = 5)
myeloid_group <- DimPlot(myeloid,group.by = "group", label = F, repel = T, label.size = 7) + scale_color_manual(values = c('Yellow Green', 'Pink 3', 'Steel Blue 4')) +
  annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)+
  theme_bw(base_size = 22, base_family = "Helvetica")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+#移除网格线
  xlab("UMAP 1") + ylab("UMAP 2")+
  theme(axis.text.x = element_text(size = 24, family = "Helvetica", face = "plain", color = 'black'), axis.text.y = element_text(size = 24, family = "Helvetica", face = "plain", color = 'black'))+
  theme(axis.title.x = element_text(size = 24, family = "Helvetica", color = 'black'), axis.title.y = element_text(size = 24, family = "Helvetica", color = 'black'))+
  theme(plot.title = element_blank())+
  guides(color=guide_legend(override.aes = list(size=6)))

ggsave("myeloid_group.pdf", myeloid_group, width = 7.3, height = 5)
save(myeloid, file = "/data/hxl/myeloid.Rdata")
VlnPlot(neutrophils, features = c("Nfkbiz", "Il1r2", "Cxcl2"), cols = c("Olive Drab", "Pink", "Lightblue"), pt.size = 0)
#cell proportion
myeloid$Label <- Idents(myeloid)
myeloid$Label <- factor(myeloid$Label, levels = c("Neutrophils", "Basophils", "Stmn+ AMs", "Stmn- AMs",
                                                  "IM1", "IM2", "cDCs", "Inflammatory Monocytes",
                                                  "Nkg7+ Monocytes", "Ly6d+ Monocytes"))
head(myeloid@meta.data)
table(myeloid$Label)
table(myeloid$sample, myeloid$Label)
t.myeloid=as.data.frame.matrix(table(myeloid$sample, myeloid$Label))
n_sample_myeloid <- ddply(myeloid@meta.data, .(sample, Label), nrow)
n_sample_myeloid <- ddply(n_sample_myeloid, .(sample), transform, percent=V1/sum(V1)*100)
write.csv(n_sample_myeloid, file = "/data/hxl/csv/n_sample_myeloid.csv")
n_sample_myeloid$sample <- factor(as.character(n_sample_myeloid$sample), levels = c("Ctrl1", "Ctrl2", "Acute1", "Acute2", "Chronic1", "Chronic2"))
ggplot(n_sample_myeloid, aes(x=sample, y=percent, fill=Label)) +
  labs(x='', y='Percentage', fill='Cell Type')+
  geom_bar(stat="identity", colour="black") +
  scale_fill_brewer(palette="Set3")+
  theme_bw(base_size = 22)+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5,color="black"))+ 
  scale_x_discrete(labels = c("Ctrl1", "Ctrl2", "Acute1", "Acute2", "Chronic1", "Chronic2"))+
  theme(axis.text.y = element_text(family = "Helvetica", size = 22, color = "black"))
ggsave('myeloid_cell_percentage.pdf', width = 7, height = 4.5)
myeloid$sample <- myeloid$orig.ident
table(myeloid$sample)
str(myeloid@meta.data)
myeloid$sample <- as.factor(myeloid$orig.ident)
str(myeloid@meta.data)
levels(myeloid$sample)
levels(myeloid$sample) <- c('Acute1', 'Acute2', 'Chronic1', 'Chronic2', 'Ctrl1','Ctrl2')
str(myeloid@meta.data)
myeloid$sample <- factor(myeloid$sample, levels = c('Ctrl1', 'Ctrl2', 'Acute1', 'Acute2', 'Chronic1', 'Chronic2'))
DimPlot(myeloid, group.by = "sample", split.by = "sample", ncol = 2) +
  scale_color_manual(values = c('Yellow Green', 'Olive Drab', 'Pink 1', 'Pink 3', 'Steel Blue 2', 'Steel Blue 4'))
save(myeloid, file = "/data/hxl/Rdata/myeloid.Rdata")
#cell proportion(no neutrophils)

myeloid1 <- myeloid[, !myeloid$Label %in% c("Neutrophils")]
myeloid1$Label <- Idents(myeloid)
head(myeloid1@meta.data)
table(myeloid1$Label)
table(myeloid1$orig.ident, myeloid1$Label)
t.myeloid1=as.data.frame.matrix(table(myeloid1$orig.ident, myeloid1$Label))
n_sample_myeloid1 <- ddply(myeloid1@meta.data, .(orig.ident, Label), nrow)
n_sample_myeloid1 <- ddply(n_sample_myeloid1, .(orig.ident), transform, percent=V1/sum(V1)*100)
write.csv(n_sample_myeloid1, file = "/data/hxl/csv/n_sample_myeloid1.csv")
n_sample_myeloid1$orig.ident <- factor(as.character(n_sample_myeloid1$orig.ident), levels = lv)
ggplot(n_sample_myeloid1, aes(x=orig.ident, y=percent, fill=Label)) +
  labs(x='', y='Percentage', fill='Cell Type')+
  geom_bar(stat="identity", colour="black") +
  scale_fill_brewer(palette = ("Set3")) +
  theme_bw(base_size = 19, base_family = "Helvetica")+
  theme(axis.text.x = element_text(family = "Helvetica", size = 19, angle=90, hjust=1, vjust=.5,color="black"))+ 
  scale_x_discrete(labels = c("Ctrl1", "Ctrl2", "Acute1", "Acute2", "Chronic1", "Chronic2"))+
  theme(axis.text.y = element_text(family = "Helvetica", size = 20, color = "black"))
ggsave('myeloid1_cell_percentage.pdf', width = 7.5, height = 5)
save(myeloid1, file = "/data/hxl/Rdata/myeloid1.Rdata")