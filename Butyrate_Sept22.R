#Removed header info from RSEC raw file using TextEdit
library(tidyverse)
library(ggplot2)
library(Seurat)
library(monocle3)
library(ggplot2)
library(ggrepel)



# Read in clean Treg cells RDS into Seurat object called SubsetNoFresh

SubsetNoFresh <- readRDS("ExpandedTregsClean.rds")
SubsetNoFresh <- subset(SubsetNoFresh, idents = c("Control", "Butyrate"))

SubsetNoFresh$Treatment <- factor(SubsetNoFresh$Treatment, levels = c("Control", "Butyrate"))

SubsetNoFresh <- FindVariableFeatures(SubsetNoFresh)
SubsetNoFresh <- ScaleData(SubsetNoFresh)
SubsetNoFresh <- RunPCA(SubsetNoFresh, verbose = FALSE)
SubsetNoFresh <- FindNeighbors(SubsetNoFresh, dims = 1:30)
SubsetNoFresh <- FindClusters(SubsetNoFresh, resolution = 1.2, verbose = FALSE)
SubsetNoFresh <- RunUMAP(SubsetNoFresh, dims = 1:30)
DimPlot(SubsetNoFresh, reduction = "umap", label = TRUE, repel = TRUE) 

DefaultAssay(SubsetNoFresh) <- "ADT"
SubsetNoFresh <- NormalizeData(SubsetNoFresh, normalization.method = "CLR", margin = 2)



### RNA DE analysis
DefaultAssay(SubsetNoFresh) <- "RNA"
Subset3_markers <- FindAllMarkers(SubsetNoFresh, only.pos = TRUE, 
                                  min.pct = 0.25 , logfc.threshold = 0.25)
top10 <- Subset3_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top10_2 <- top10[which(top10$p_val_adj < 0.05),]
### RNA DE analysis
top10_2
write.csv(top10_2, "transcript.csv", sep = ",")

### ADT DE analysis
DefaultAssay(SubsetNoFresh) <- "ADT"
Subset3_markers <- FindAllMarkers(SubsetNoFresh, only.pos = TRUE, 
                                  min.pct = 0.25, logfc.threshold = 0.9)
top10 <- Subset3_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top10_2 <- top10[which(top10$p_val_adj < 0.05),]
write.csv(top10_2, "top10.csv")
### ADT DE analysis



######Saving RDS object and exploratory analysis######
saveRDS(SubsetNoFresh, "Expanded_Tregs.rds")

Idents(SubsetNoFresh) <- "Treatment"
VlnPlot(SubsetNoFresh, c("TBX21", "GATA3", "FOXP3", "KLRB1", "RORA", "RORC"), ncol=3)

DefaultAssay(SubsetNoFresh) <- "ADT"
VlnPlot(SubsetNoFresh, c("CD25.M.A251.IL2RA.AHS0166.pAbO", 
                         "HLA.DR.CD74.AHS0035.pAbO",
                         "CD279.EH12.1.PDCD1.AHS0014.pAbO",
                         "CD161.DX12.KLRB1.AHS0002.pAbO","CD194.CCR4.AHS0038.pAbO", "CD162.SELPLG.AHS0139.pAbO"), ncol=3)


#######Correlations between RNA and ADT#######


FeatureScatter(object = SubsetNoFresh, feature1 = 'IL2RA', feature2 = 'CD25.M.A251.IL2RA.AHS0166.pAbO',
      pt.size = 0.1, plot.cor = T, group.by = "seurat_clusters")

FeatureScatter(object = SubsetNoFresh, feature1 = 'PDCD1', feature2 = 'CD279.EH12.1.PDCD1.AHS0014.pAbO',
               pt.size = 0.1, plot.cor = T, group.by = "seurat_clusters")


FeatureScatter(object = SubsetNoFresh, feature1 = 'KLRB1', feature2 = 'CD161.DX12.KLRB1.AHS0002.pAbO',
               pt.size = 0.1, plot.cor = T, group.by = "seurat_clusters")

FeatureScatter(object = SubsetNoFresh, feature1 = 'CD74', feature2 = 'HLA.DR.CD74.AHS0035.pAbO',
               pt.size = 0.1, plot.cor = T, group.by = "seurat_clusters")


########PROGENY based on Tutorials by Saez Lab: https://github.com/saezlab/transcriptutorial##########

BiocManager::install("progeny")
library(progeny)

## We create a data frame with the specification of the cells that belong to 
## each cluster to match with the Progeny scores.

Idents(SubsetNoFresh) <- "seurat_clusters"

CellsClusters <- data.frame(Cell = names(Idents(SubsetNoFresh)), 
                            CellType = as.character(Idents(SubsetNoFresh)),
                            stringsAsFactors = FALSE)



## We compute the Progeny activity scores and add them to our Seurat object
## as a new assay called Progeny. 
SubsetNoFresh <- progeny(SubsetNoFresh, scale=FALSE, organism="Human", top=500, perm=1, 
                         return_assay = TRUE)

## We can now directly apply Seurat functions in our Progeny scores. 
## For instance, we scale the pathway activity scores. 
SubsetNoFresh <- Seurat::ScaleData(SubsetNoFresh, assay = "progeny") 

## We transform Progeny scores into a data frame to better handling the results
progeny_scores_df <- 
  as.data.frame(t(GetAssayData(SubsetNoFresh, slot = "scale.data", 
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell) 

## We match Progeny scores with the cell clusters.
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

## We summarize the Progeny scores by cellpopulation
summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))



## We prepare the data for the plot
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 

paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)

progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0, 
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength, 
                      max(summarized_progeny_scores_df), 
                      length.out=floor(paletteLength/2)))

library(pheatmap)
progeny_hmap = pheatmap(t(summarized_progeny_scores_df[,-1]),fontsize=14, 
                        fontsize_row = 10, 
                        color=myColor, breaks = progenyBreaks, 
                        main = "Heterogenous Signalling Pathways in Treg cells", angle_col = 45,
                        treeheight_col = 0,  border_color = NA)


######################DOROTHEA based on Tutorials from Saez Lab:https://github.com/saezlab/transcriptutorial#########


BiocManager::install("dorothea")
library(dorothea)
library(tibble)
library(pheatmap)
library(tidyr)
library(viper)
library(Seurat)

## We read Dorothea Regulons for Human:
dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))

## We obtain the regulons based on interactions with confidence level A, B and C
regulon <- dorothea_regulon_human %>%
  dplyr::filter(confidence %in% c("A"))

DimPlot(Subset2_stat, reduction = 'umap', label = TRUE, repel = TRUE)
levels(SubsetNoFresh)
## We compute Viper Scores 
Subset2_stat <- run_viper(SubsetNoFresh, regulon,
                          options = list(method = "scale", minsize = 4, 
                                         eset.filter = FALSE, cores = 4, 
                                         verbose = FALSE))
save.image()

## We compute the Nearest Neighbours to perform cluster
levels(Subset2_stat)

DefaultAssay(object = Subset2_stat) <- "dorothea"
Subset2_stat <- ScaleData(Subset2_stat)
Subset2_stat <- RunPCA(Subset2_stat, features = rownames(Subset2_stat), verbose = FALSE)
Subset2_stat <- FindNeighbors(Subset2_stat, dims = 1:30, verbose = FALSE)
Subset2_stat <- FindClusters(Subset2_stat, resolution = 0.8, verbose = FALSE)

Idents(Subset2_stat) <- "Treatment"
## We transform Viper scores, scaled by seurat, into a data frame to better 
## handling the results
viper_scores_df <- GetAssayData(Subset2_stat, slot = "scale.data", 
                                assay = "dorothea") %>%
  data.frame(check.names = F) %>%
  t()

## We create a data frame containing the cells and their clusters
CellsClusters <- data.frame(cell = names(Idents(Subset2_stat)), 
                            cell_type = as.character(Idents(Subset2_stat)),
                            check.names = F)

## We create a data frame with the Viper score per cell and its clusters
viper_scores_clusters <- viper_scores_df  %>%
  data.frame() %>% 
  rownames_to_column("cell") %>%
  gather(tf, activity, -cell) %>%
  inner_join(CellsClusters)

## We summarize the Viper scores by cellpopulation
summarized_viper_scores <- viper_scores_clusters %>% 
  group_by(tf, cell_type) %>%
  summarise(avg = mean(activity),
            std = sd(activity))

## We select the 20 most variable TFs. (20*9 populations = 180)
highly_variable_tfs <- summarized_viper_scores %>%
  group_by(tf) %>%
  mutate(var = var(avg))  %>%
  ungroup() %>%
  top_n(180, var) %>%
  distinct(tf)

## We prepare the data for the plot
summarized_viper_scores_df <- summarized_viper_scores %>%
  semi_join(highly_variable_tfs, by = "tf") %>%
  dplyr::select(-std) %>%   
  spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE) 

palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(min(summarized_viper_scores_df), 0, 
                   length.out=ceiling(palette_length/2) + 1),
               seq(max(summarized_viper_scores_df)/palette_length, 
                   max(summarized_viper_scores_df), 
                   length.out=floor(palette_length/2)))

viper_hmap <- pheatmap(t(summarized_viper_scores_df),fontsize=14, 
                       fontsize_row = 10, 
                       color=my_color, breaks = my_breaks, 
                       main = "Heterogenous Regulons in Treg cells", angle_col = 0,
                       treeheight_col = 0,  border_color = NA) 

DimPlot(SubsetNoFresh)
SubsetNoFresh

VlnPlot(SubsetNoFresh, features = c("FFAR1", "FFAR2", "FFAR3", "FFAR4"), log= F, same.y.lims = F)

#### Nebulosa for kernel-density-based visualisation ###
#### Nebulosa for kernel-density-based visualisation ###
#### Nebulosa for kernel-density-based visualisation ###

BiocManager::install("Nebulosa")
library(Nebulosa)
library(Seurat)
library("BiocFileCache")

DefaultAssay(SubsetNoFresh) <- "RNA"

SubsetNoFresh

DimPlot(SubsetNoFresh)
plot_density(SubsetNoFresh, reduction = "umap", c("FOXP3", "IL2RA", "CTLA4","PDCD1", "PTMA", "HMGB2", "TUBA1B", 
                                                  "LAG3", "SELPLG", "TIGIT", "LGALS3", "MAL"), 
             joint = FALSE, pal = "inferno" )

DefaultAssay(SubsetNoFresh) <- "ADT"
plot_density(SubsetNoFresh, reduction = "umap", c("CD25.M.A251.IL2RA.AHS0166.pAbO", 
                                                  "CD69.CD69.AHS0010.pAbO",
                                                  "HLA.DR.CD74.AHS0035.pAbO",
                                                  "CD279.EH12.1.PDCD1.AHS0014.pAbO",
                                                  "CXCR5.CXCR5.AHS0039.pAbO",
                                                  "CD194.CCR4.AHS0038.pAbO",
                                                  "CD161.DX12.KLRB1.AHS0002.pAbO", 
                                                  "HLA.A-B-C.HLA.A-B-C.AHS0066.pAbO",
                                                  "CD31.WM59.PECAM1.AHS0170.pAbO",
                                                  "CD45RA.HI100.PTPRC.AHS0009.pAbO",
                                                  "CD134.ACT35.TNFRSF4.AHS0013.pAbO",
                                                  "CD137.TNFRSF9.AHS0003.pAbO",
                                                  "CD154.CD40LG.AHS0077.pAbO",
                                                  "CD49b.ITGA2.AHS0093.pAbO", 
                                                  "CD162.SELPLG.AHS0139.pAbO",
                                                  "CD44.515.CD44.AHS0140.pAbO",
                                                  "CD2.CD2.AHS0029.pAbO"), pal = "inferno") + NoLegend()
#New Figure 24/3/22
FeaturePlot(SubsetNoFresh, reduction = "umap", c("CD25.M.A251.IL2RA.AHS0166.pAbO", 
                                                  "CD69.CD69.AHS0010.pAbO",
                                                  "HLA.DR.CD74.AHS0035.pAbO",
                                                  "CD279.EH12.1.PDCD1.AHS0014.pAbO",
                                                 "HLA.A-B-C.HLA.A-B-C.AHS0066.pAbO",
                                                 "CXCR5.CXCR5.AHS0039.pAbO",
                                                  "CD194.CCR4.AHS0038.pAbO",
                                                  "CD161.DX12.KLRB1.AHS0002.pAbO", 
                                                  "CD134.ACT35.TNFRSF4.AHS0013.pAbO",
                                                  "CD137.TNFRSF9.AHS0003.pAbO",
                                                  "CD49b.ITGA2.AHS0093.pAbO", 
                                                  "CD162.SELPLG.AHS0139.pAbO",
                                                  "CD44.515.CD44.AHS0140.pAbO",
                                                 "CD2.CD2.AHS0029.pAbO"), ncol = 7)
#New Figure 24/3/22
DefaultAssay(SubsetNoFresh) <- "RNA"
FeaturePlot(SubsetNoFresh, c("FOXP3", "IL2RA", "CTLA4","PDCD1", "CD74", "KLRB1", "PTMA", "HMGB2", "TUBA1B", 
"LAG3", "SELPLG", "TIGIT", "LGALS3", "MAL"), ncol=7)


#### Nebulosa ###
#### Nebulosa ###
#### Nebulosa ###


### Cell frequency graph ###
library(tidyverse)
library(RColorBrewer)


pt <- table(Idents(SubsetNoFresh), SubsetNoFresh$orig.ident)
pt <- as.data.frame(pt)
View(pt)
pt$Var1 <- as.character(pt$Var1)
colnames(pt) <- c("Cluster", "Data", "Frequency")

ggplot(pt, aes(x = Data, y = Frequency, fill = Cluster)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) + ggtitle("Proportions")
  ylab("Proportions") +
  scale_fill_manual(values = my_color_palette) +
  theme(legend.title = element_blank())

save.image()



#Exploratory analyses for Ridgeplots to visualise Abseqs# 

First6Abs <- c("CD45.PTPRC.AHS0040.pAbO", "CD3.UCHT1.CD3E.AHS0231.pAbO", 
               "CD4.SK3.CD4.AHS0032.pAbO", "CD8.SK1.CD8A.AHS0228.pAbO", 
               "CD31.WM59.PECAM1.AHS0170.pAbO","CD45RA.HI100.PTPRC.AHS0009.pAbO")

Second6Abs <- c("CD69.CD69.AHS0010.pAbO", "CD25.M.A251.IL2RA.AHS0166.pAbO",
                "CD134.ACT35.TNFRSF4.AHS0013.pAbO","CD137.TNFRSF9.AHS0003.pAbO",
                "CD154.CD40LG.AHS0077.pAbO", "CD161.DX12.KLRB1.AHS0002.pAbO")

Third6Abs <- c("HLA.A-B-C.HLA.A-B-C.AHS0066.pAbO","HLA.DR.CD74.AHS0035.pAbO",
               "CD279.EH12.1.PDCD1.AHS0014.pAbO","CXCR5.CXCR5.AHS0039.pAbO",
               "CD162.SELPLG.AHS0139.pAbO","CD18.ITGB2.AHS0091.pAbO") 


Fourth6Abs <- c("CD11a.ITGAL.AHS0081.pAbO", "CD194.CCR4.AHS0038.pAbO",
                "CD44.515.CD44.AHS0140.pAbO","CD46.CD46.AHS0071.pAbO", 
                "CD49b.ITGA2.AHS0093.pAbO", "CD59.CD59.AHS0224.pAbO")

Fifth6Abs <- c("CD278.ICOS.AHS0012.pAbO", "CD2.CD2.AHS0029.pAbO", "CD5.UCHT2.CD5.AHS0047.pAbO")

RidgePlot(SubsetNoFresh, group.by = "Treatment", features = First6Abs)
RidgePlot(SubsetNoFresh, group.by = "Treatment", features = Second6Abs)
RidgePlot(SubsetNoFresh, group.by = "Treatment", features = Third6Abs)
RidgePlot(SubsetNoFresh, group.by = "Treatment", features = Fourth6Abs)
RidgePlot(SubsetNoFresh, group.by = "Treatment", features = Fifth6Abs)


#FeaturePlots to visualise transcriptome and AbSeq as well as correlations between them


p1 <- FeaturePlot(SubsetNoFresh, "CD45.PTPRC.AHS0040.pAbO", 
                  cols = c("lightgrey", "darkgreen")) +
  ggtitle("CD45 Protein")
p2 <- FeaturePlot(SubsetNoFresh, "PTPRC") + ggtitle("PTPRC RNA")
p1 | p2
FeatureScatter(SubsetNoFresh, 
               feature1 = "CD45.PTPRC.AHS0040.pAbO", 
               feature2 = "PTPRC")

p1 <- FeaturePlot(SubsetNoFresh, "CD3.UCHT1.CD3E.AHS0231.pAbO", 
                  cols = c("lightgrey", "darkgreen")) +
  ggtitle("CD3 Protein")
p2 <- FeaturePlot(SubsetNoFresh, "CD3E") + ggtitle("CD3E RNA")
p1 | p2
FeatureScatter(SubsetNoFresh, 
               feature1 = "CD3.UCHT1.CD3E.AHS0231.pAbO", 
               feature2 = "CD3E")

p1 <- FeaturePlot(SubsetNoFresh, "CD4.SK3.CD4.AHS0032.pAbO", 
                  cols = c("lightgrey", "darkgreen")) +
  ggtitle("CD4 Protein")
p2 <- FeaturePlot(SubsetNoFresh, "CD4") + ggtitle("CD4 RNA")
p1 | p2
FeatureScatter(SubsetNoFresh, 
               feature1 = "CD4.SK3.CD4.AHS0032.pAbO", 
               feature2 = "CD4")

p1 <- FeaturePlot(SubsetNoFresh, "CD8.SK1.CD8A.AHS0228.pAbO", 
                  cols = c("lightgrey", "darkgreen")) +
  ggtitle("CD8 Protein")
p2 <- FeaturePlot(SubsetNoFresh, "CD8A") + ggtitle("CD8 RNA")
p1 | p2
FeatureScatter(SubsetNoFresh, 
               feature1 = "CD8.SK1.CD8A.AHS0228.pAbO", 
               feature2 = "CD8A")




p1 <- FeaturePlot(SubsetNoFresh, "CD31.WM59.PECAM1.AHS0170.pAbO", 
                  cols = c("lightgrey", "darkgreen")) +
  ggtitle("CD31 Protein")
p2 <- FeaturePlot(SubsetNoFresh, "PECAM1") + ggtitle("PECAM1 RNA")
p1 | p2
FeatureScatter(SubsetNoFresh, 
               feature1 = "CD31.WM59.PECAM1.AHS0170.pAbO", 
               feature2 = "PECAM1")

p1 <- FeaturePlot(SubsetNoFresh, "CD45RA.HI100.PTPRC.AHS0009.pAbO", 
                  cols = c("lightgrey", "darkgreen")) +
  ggtitle("CD45RA Protein")
p2 <- FeaturePlot(SubsetNoFresh, "PTPRC") + ggtitle("PTPRC RNA")
p1 | p2
FeatureScatter(SubsetNoFresh, 
               feature1 = "CD45RA.HI100.PTPRC.AHS0009.pAbO", 
               feature2 = "PTPRC")

p1 <- FeaturePlot(SubsetNoFresh, "CD69.CD69.AHS0010.pAbO", 
                  cols = c("lightgrey", "darkgreen")) +
  ggtitle("CD69 Protein")
p2 <- FeaturePlot(SubsetNoFresh, "CD69") + ggtitle("CD69 RNA")
p1 | p2
FeatureScatter(SubsetNoFresh, 
               feature1 = "CD69.CD69.AHS0010.pAbO", 
               feature2 = "CD69")

p1 <- FeaturePlot(SubsetNoFresh, "CD25.M.A251.IL2RA.AHS0166.pAbO", 
                  cols = c("lightgrey", "darkgreen")) +
  ggtitle("CD25 Protein")
p2 <- FeaturePlot(SubsetNoFresh, "IL2RA") + ggtitle("IL2RA RNA")
p1 | p2
FeatureScatter(SubsetNoFresh, 
               feature1 = "CD25.M.A251.IL2RA.AHS0166.pAbO", 
               feature2 = "IL2RA")

p1 <- FeaturePlot(SubsetNoFresh, "CD134.ACT35.TNFRSF4.AHS0013.pAbO", 
                  cols = c("lightgrey", "darkgreen")) +
  ggtitle("OX40 Protein")
p2 <- FeaturePlot(SubsetNoFresh, "TNFRSF4") + ggtitle("OX40 RNA")
p1 | p2
FeatureScatter(SubsetNoFresh, 
               feature1 = "CD134.ACT35.TNFRSF4.AHS0013.pAbO", 
               feature2 = "TNFRSF4")

p1 <- FeaturePlot(SubsetNoFresh, "CD137.TNFRSF9.AHS0003.pAbO", 
                  cols = c("lightgrey", "darkgreen")) +
  ggtitle("4-1BB Protein")
p2 <- FeaturePlot(SubsetNoFresh, "TNFRSF9") + ggtitle("4-1BB RNA")
p1 | p2
FeatureScatter(SubsetNoFresh, 
               feature1 = "CD137.TNFRSF9.AHS0003.pAbO", 
               feature2 = "TNFRSF9")

p1 <- FeaturePlot(SubsetNoFresh, "CD154.CD40LG.AHS0077.pAbO", 
                  cols = c("lightgrey", "darkgreen")) +
  ggtitle("CD40L Protein")
p2 <- FeaturePlot(SubsetNoFresh, "CD40LG") + ggtitle("CD40LG RNA")
p1 | p2
FeatureScatter(SubsetNoFresh, 
               feature1 = "CD154.CD40LG.AHS0077.pAbO", 
               feature2 = "CD40LG")

p1 <- FeaturePlot(SubsetNoFresh, "CD161.DX12.KLRB1.AHS0002.pAbO", 
                  cols = c("lightgrey", "darkgreen")) +
  ggtitle("CD161 Protein")
p2 <- FeaturePlot(SubsetNoFresh, "KLRB1") + ggtitle("KLRB1 RNA")
p1 | p2
FeatureScatter(SubsetNoFresh, 
               feature1 = "CD161.DX12.KLRB1.AHS0002.pAbO", 
               feature2 = "KLRB1")

p1 <- FeaturePlot(SubsetNoFresh, "HLA.A-B-C.HLA.A-B-C.AHS0066.pAbO", 
                  cols = c("lightgrey", "darkgreen")) +
  ggtitle("HLA-ABC Protein")
p1
# p2 <- FeaturePlot(SubsetNoFresh, "HLA.A-B-C") + ggtitle("HLA.A-B-C RNA")
# p1 | p2
# FeatureScatter(SubsetNoFresh, 
#               feature1 = "HLA.A-B-C.HLA.A-B-C.AHS0066.pAbO", 
#              feature2 = "HLA.A-B-C")

p1 <- FeaturePlot(SubsetNoFresh, "HLA.DR.CD74.AHS0035.pAbO", 
                  cols = c("lightgrey", "darkgreen")) +
  ggtitle("HLA-DR Protein")
p2 <- FeaturePlot(SubsetNoFresh, "CD74") + ggtitle("CD74 RNA")
p1 | p2
FeatureScatter(SubsetNoFresh, 
               feature1 = "HLA.DR.CD74.AHS0035.pAbO", 
               feature2 = "CD74")

p1 <- FeaturePlot(SubsetNoFresh, "CD279.EH12.1.PDCD1.AHS0014.pAbO", 
                  cols = c("lightgrey", "darkgreen")) +
  ggtitle("PD1 Protein")
p2 <- FeaturePlot(SubsetNoFresh, "PDCD1") + ggtitle("PD1 RNA")
p1 | p2
FeatureScatter(SubsetNoFresh, 
               feature1 = "CD279.EH12.1.PDCD1.AHS0014.pAbO", 
               feature2 = "PDCD1")

p1 <- FeaturePlot(SubsetNoFresh, "CXCR5.CXCR5.AHS0039.pAbO", 
                  cols = c("lightgrey", "darkgreen")) +
  ggtitle("CXCR5 Protein")
p2 <- FeaturePlot(SubsetNoFresh, "CXCR5") + ggtitle("CXCR5 RNA")
p1 | p2
FeatureScatter(SubsetNoFresh, 
               feature1 = "CXCR5.CXCR5.AHS0039.pAbO", 
               feature2 = "CXCR5")

p1 <- FeaturePlot(SubsetNoFresh, "CD162.SELPLG.AHS0139.pAbO", 
                  cols = c("lightgrey", "darkgreen")) +
  ggtitle("P-selectin Glycoprotein Ligand-1")
p2 <- FeaturePlot(SubsetNoFresh, "SELPLG") + ggtitle("SELPLG RNA")
p1 | p2
FeatureScatter(SubsetNoFresh, 
               feature1 = "CD162.SELPLG.AHS0139.pAbO", 
               feature2 = "SELPLG")

p1 <- FeaturePlot(SubsetNoFresh, "CD18.ITGB2.AHS0091.pAbO", 
                  cols = c("lightgrey", "darkgreen")) +
  ggtitle("Integrin Beta-2 Protein")
p2 <- FeaturePlot(SubsetNoFresh, "ITGB2") + ggtitle("ITGB2 RNA")
p1 | p2
FeatureScatter(SubsetNoFresh, 
               feature1 = "CD18.ITGB2.AHS0091.pAbO", 
               feature2 = "ITGB2")

p1 <- FeaturePlot(SubsetNoFresh, "CD11a.ITGAL.AHS0081.pAbO", 
                  cols = c("lightgrey", "darkgreen")) +
  ggtitle("CD11a Protein")
p2 <- FeaturePlot(SubsetNoFresh, "ITGAL") + ggtitle("ITGAL RNA")
p1 | p2
FeatureScatter(SubsetNoFresh, 
               feature1 = "CD11a.ITGAL.AHS0081.pAbO", 
               feature2 = "ITGAL")

p1 <- FeaturePlot(SubsetNoFresh, "CD194.CCR4.AHS0038.pAbO", 
                  cols = c("lightgrey", "darkgreen")) +
  ggtitle("CCR4 Protein")
p2 <- FeaturePlot(SubsetNoFresh, "CCR4") + ggtitle("CCR4 RNA")
p1 | p2
FeatureScatter(SubsetNoFresh, 
               feature1 = "CD194.CCR4.AHS0038.pAbO", 
               feature2 = "CCR4")

p1 <- FeaturePlot(SubsetNoFresh, "CD44.515.CD44.AHS0140.pAbO", 
                  cols = c("lightgrey", "darkgreen")) +
  ggtitle("CD44 Protein")
p2 <- FeaturePlot(SubsetNoFresh, "CD44") + ggtitle("CD44 RNA")
p1 | p2
FeatureScatter(SubsetNoFresh, 
               feature1 = "CD44.515.CD44.AHS0140.pAbO", 
               feature2 = "CD44")

p1 <- FeaturePlot(SubsetNoFresh, "CD46.CD46.AHS0071.pAbO", 
                  cols = c("lightgrey", "darkgreen")) +
  ggtitle("CD46 Protein")
p2 <- FeaturePlot(SubsetNoFresh, "CD46") + ggtitle("CD46 RNA")
p1 | p2
FeatureScatter(SubsetNoFresh, 
               feature1 = "CD46.CD46.AHS0071.pAbO", 
               feature2 = "CD46")

p1 <- FeaturePlot(SubsetNoFresh, "CD49b.ITGA2.AHS0093.pAbO", 
                  cols = c("lightgrey", "darkgreen")) +
  ggtitle("CD49b Protein")
p2 <- FeaturePlot(SubsetNoFresh, "ITGA2") + ggtitle("ITGA2 RNA")
p1 | p2
FeatureScatter(SubsetNoFresh, 
               feature1 = "CD49b.ITGA2.AHS0093.pAbO", 
               feature2 = "ITGA2")

p1 <- FeaturePlot(SubsetNoFresh, "CD59.CD59.AHS0224.pAbO", 
                  cols = c("lightgrey", "darkgreen")) +
  ggtitle("CD59 Protein")
p2 <- FeaturePlot(SubsetNoFresh, "CD59") + ggtitle("CD59 RNA")
p1 | p2
FeatureScatter(SubsetNoFresh, 
               feature1 = "CD59.CD59.AHS0224.pAbO", 
               feature2 = "CD59")

p1 <- FeaturePlot(SubsetNoFresh, "CD2.CD2.AHS0029.pAbO", 
                  cols = c("lightgrey", "darkgreen")) +
  ggtitle("CD2 Protein")
p2 <- FeaturePlot(SubsetNoFresh, "CD2") + ggtitle("CD2 RNA")
p1 | p2
FeatureScatter(SubsetNoFresh, 
               feature1 = "CD2.CD2.AHS0029.pAbO", 
               feature2 = "CD2")

p1 <- FeaturePlot(SubsetNoFresh, "CD5.UCHT2.CD5.AHS0047.pAbO", 
                  cols = c("lightgrey", "darkgreen")) +
  ggtitle("CD5 Protein")
p2 <- FeaturePlot(SubsetNoFresh, "CD5") + ggtitle("CD5 RNA")
p1 | p2
FeatureScatter(SubsetNoFresh, 
               feature1 = "CD5.UCHT2.CD5.AHS0047.pAbO", 
               feature2 = "CD5")
save.image()




