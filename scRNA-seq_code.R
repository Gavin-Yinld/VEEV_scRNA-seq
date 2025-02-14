options(stringAsFactors=F)
library(dplyr)
library(Seurat)
library(patchwork)
library(stringr)
library(harmony)
#### read 10X genomics scRNAseq data
MOCK_106_DPI_Spike <- Read10X("MOCK-90-DPI_results_FRP/sample_filtered_feature_bc_matrix/sample_filtered_feature_bc_matrix/")
VEEV_106_DPI_Spike <- Read10X("VEEV-90-DPI_results_FRP/sample_filtered_feature_bc_matrix/sample_filtered_feature_bc_matrix/")
VEEV_neuro_106_DPI_Spike <- Read10X("VEEV-neuro-90-DPI_results_FRP/sample_filtered_feature_bc_matrix/sample_filtered_feature_bc_matrix/")
MOCK_7_DPI_Spike <- Read10X("MOCK-7-DPI_results_FRP/sample_filtered_feature_bc_matrix/sample_filtered_feature_bc_matrix/")
VEEV_7_DPI_Spike <- Read10X("VEEV-7-DPI_results_FRP/sample_filtered_feature_bc_matrix/sample_filtered_feature_bc_matrix/")
VEEV_neuro_106_DPI_v2_Spike <- Read10X("VEEV-neuro-90-DPI_results_FRP_additional/sample_filtered_feature_bc_matrix/")

MOCK_106_DPI_Spike.obj <- CreateSeuratObject(counts =  MOCK_106_DPI_Spike, project = "MOCK_106_DPI_Spike", min.cells = 0, min.features = 0)
VEEV_106_DPI_Spike.obj <- CreateSeuratObject(counts =  VEEV_106_DPI_Spike, project = "VEEV_106_DPI_Spike", min.cells = 0, min.features = 0)
VEEV_neuro_106_DPI_Spike.obj <- CreateSeuratObject(counts =  VEEV_neuro_106_DPI_Spike, project = "VEEV_neuro_106_DPI_Spike", min.cells = 0, min.features = 0)
MOCK_7_DPI_Spike.obj <- CreateSeuratObject(counts =  MOCK_7_DPI_Spike, project = "MOCK_7_DPI_Spike", min.cells = 0, min.features = 0)
VEEV_7_DPI_Spike.obj <- CreateSeuratObject(counts =  VEEV_7_DPI_Spike, project = "VEEV_7_DPI_Spike", min.cells = 0, min.features = 0)
VEEV_neuro_106_DPI_v2_Spike.obj <- CreateSeuratObject(counts =  VEEV_neuro_106_DPI_v2_Spike, project = "VEEV_neuro_106_DPI_v2_Spike", min.cells = 0, min.features = 0)


merge_obj <- merge(x=MOCK_106_DPI_Spike.obj, 
                   y = c(VEEV_106_DPI_Spike.obj,
                         VEEV_neuro_106_DPI_Spike.obj,
                         MOCK_7_DPI_Spike.obj,
                         VEEV_7_DPI_Spike.obj,
                         VEEV_neuro_106_DPI_v2_Spike.obj
                         ),
                   add.cell.ids = c("Study-2-MOCK-90-DPI",
                                    "Study-2-VEEV-90-DPI",
                                    "Study-2-VEEV-neuro-90-DPI",
                                    "Study-3-MOCK-7-DPI",
                                    "Study-3-VEEV-7-DPI",
                                    "VEEV_neuro_106_DPI_add"
                                    ) )
merge_obj <- JoinLayers(merge_obj)

#################################################
######################## merge two batches of VEEV_neuro_106_DPI

length(unique(c(colnames(VEEV_neuro_106_DPI_Spike) , colnames(VEEV_neuro_106_DPI_v2_Spike))))
all(rownames(VEEV_neuro_106_DPI_Spike) == rownames(VEEV_neuro_106_DPI_v2_Spike))
VEEV_neuro_106_DPI_merge_Spike <- VEEV_neuro_106_DPI_v2_Spike

cells <- intersect(colnames(VEEV_neuro_106_DPI_Spike) , colnames(VEEV_neuro_106_DPI_v2_Spike))
VEEV_neuro_106_DPI_merge_Spike <- VEEV_neuro_106_DPI_v2_Spike
VEEV_neuro_106_DPI_merge_Spike[,cells] <- VEEV_neuro_106_DPI_merge_Spike[,cells] + VEEV_neuro_106_DPI_Spike[,cells]
VEEV_neuro_106_DPI_merge_Spike <- cbind(VEEV_neuro_106_DPI_merge_Spike,
                                                    VEEV_neuro_106_DPI_Spike[,-which(colnames(VEEV_neuro_106_DPI_Spike) %in% cells)]
                                                    )

VEEV_neuro_106_DPI_merge_Spike.obj <- CreateSeuratObject(counts =  VEEV_neuro_106_DPI_merge_Spike, project = "VEEV_neuro_106_DPI_merge_Spike", min.cells = 0, min.features = 0)

merge_obj <- merge(x=MOCK_106_DPI_Spike.obj, 
                   y = c(VEEV_106_DPI_Spike.obj,
                         VEEV_neuro_106_DPI_merge_Spike.obj,
                         MOCK_7_DPI_Spike.obj,
                         VEEV_7_DPI_Spike.obj
                   ),
                   add.cell.ids = c("Study-2-MOCK-90-DPI",
                                    "Study-2-VEEV-90-DPI",
                                    "Study-2-VEEV-neuro-90-DPI",
                                    "Study-3-MOCK-7-DPI",
                                    "Study-3-VEEV-7-DPI"
                   ) )
merge_obj <- JoinLayers(merge_obj)
#######################################################################
##### summary parameters after merge two batches
table(merge_obj$orig.ident)
tapply(merge_obj$nCount_RNA,merge_obj$orig.ident,median)
tapply(merge_obj$nFeature_RNA,merge_obj$orig.ident,median)

length(which(rowSums(MOCK_106_DPI_Spike)>0))
length(which(rowSums(VEEV_106_DPI_Spike)>0))
length(which(rowSums(VEEV_neuro_106_DPI_merge_Spike)>0))
length(which(rowSums(MOCK_7_DPI_Spike)>0))
length(which(rowSums(VEEV_7_DPI_Spike)>0))

###########################################
##  summary infected cells

merge_obj$probe1=merge_obj[["RNA"]]$counts["Probe1",]
merge_obj$probe2=merge_obj[["RNA"]]$counts["Probe2",]
merge_obj$probe3=merge_obj[["RNA"]]$counts["Probe3",]
merge_obj$probe4=merge_obj[["RNA"]]$counts["Probe4",]
merge_obj$probe5=merge_obj[["RNA"]]$counts["Probe5",]
merge_obj$probe6=merge_obj[["RNA"]]$counts["Probe6",]
merge_obj$probe7=merge_obj[["RNA"]]$counts["Probe7",]

merge_obj$probe9=merge_obj[["RNA"]]$counts["Probe9",]
merge_obj$probe10=merge_obj[["RNA"]]$counts["Probe10",]
merge_obj$probe11=merge_obj[["RNA"]]$counts["Probe10",]
merge_obj$probe12=merge_obj[["RNA"]]$counts["Probe12",]


probe_mat <- data.frame(probe1=merge_obj[["RNA"]]$counts["Probe1",],
                        probe2=merge_obj[["RNA"]]$counts["Probe2",],
                        probe3=merge_obj[["RNA"]]$counts["Probe3",],
                        probe4=merge_obj[["RNA"]]$counts["Probe4",],
                        probe5=merge_obj[["RNA"]]$counts["Probe5",],
                        probe6=merge_obj[["RNA"]]$counts["Probe6",],
                        probe7=merge_obj[["RNA"]]$counts["Probe7",],
                        probe9=merge_obj[["RNA"]]$counts["Probe9",],
                        probe10=merge_obj[["RNA"]]$counts["Probe10",],
                        probe11=merge_obj[["RNA"]]$counts["Probe11",],
                        probe12=merge_obj[["RNA"]]$counts["Probe12",]
)


index <- apply(probe_mat,1,function(x){return(all(x>0))}   ) 

merge_obj$infected2 <- as.numeric(0)
merge_obj$infected2[which(index == TRUE )] <- 1

##############################################
######### quality control

merge_obj[["percent.mt"]] <- PercentageFeatureSet(merge_obj, pattern = "mt-")
# Visualize QC metrics as a violin plot
VlnPlot(merge_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(merge_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(merge_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

####################################
######### filter low quality cells
mt.cutoff <- 10
merge_obj <- merge_obj[,which(merge_obj$percent.mt <= mt.cutoff)]
merge_obj <- merge_obj[,which(merge_obj$nFeature_RNA >= 200)]


#######################################################################
##### summary parameters after quality control
table(merge_obj$orig.ident)
tapply(merge_obj$nCount_RNA,merge_obj$orig.ident,median)
tapply(merge_obj$nFeature_RNA,merge_obj$orig.ident,median)

#######################################################################
## clustering

merge_obj <- NormalizeData(merge_obj)
merge_obj <- FindVariableFeatures(merge_obj, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(merge_obj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(merge_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(merge_obj)

merge_obj <- ScaleData(merge_obj, features = all.genes)
merge_obj <- RunPCA(merge_obj, features = VariableFeatures(object = merge_obj))
merge_obj <- RunHarmony(merge_obj, group.by.vars = "orig.ident")

DimPlot(merge_obj, reduction = "pca") + NoLegend()
DimHeatmap(merge_obj, dims = 1:3, cells = 500, balanced = TRUE)
ElbowPlot(merge_obj,ndims =30)
merge_obj <- FindNeighbors(merge_obj, dims = 1:20,reduction='harmony')
merge_obj <- FindClusters(merge_obj, resolution = 0.2)
merge_obj <- RunUMAP(merge_obj, dims = 1:20,reduction = "harmony")

DimPlot(merge_obj, reduction = "umap",label = T,pt.size=0.5)
save(merge_obj,file=merge_obj,cluster.Rdata)

#################################################################
############################## cell type annotation

merge_obj$celltype <- ""
merge_obj$celltype[which(merge_obj$seurat_clusters %in% c(0,11,14))] = "Inhibitory_neurons"
merge_obj$celltype[which(merge_obj$seurat_clusters %in% c(4:6,13))] = "Excitatory_neurons"
merge_obj$celltype[which(merge_obj$seurat_clusters %in% c(1))] = "Oligodendrocytes"
merge_obj$celltype[which(merge_obj$seurat_clusters %in% c(2))] = "Microglia"
merge_obj$celltype[which(merge_obj$seurat_clusters %in% c(3))] = "Astrocytes"
merge_obj$celltype[which(merge_obj$seurat_clusters %in% c(7))] = "Pericytes "
merge_obj$celltype[which(merge_obj$seurat_clusters %in% c(8,16))] = "Immune_cells"
merge_obj$celltype[which(merge_obj$seurat_clusters %in% c(9))] = "Endothelial_cells"
merge_obj$celltype[which(merge_obj$seurat_clusters %in% c(10))] = "OPCs"
merge_obj$celltype[which(merge_obj$seurat_clusters %in% c(12))] = "Neural_progenitor_cell"
merge_obj$celltype[which(merge_obj$seurat_clusters %in% c(15))] = "Ependymal_cells"

###### plot cell type markers and UMAP
library(ggplot2)
colors <- c(  "#FB8072",  "#45BFC0",  "#96B85D", "#C790E0","#DEB068", "#4EAEE3","#FF8DC6",  "#52B77C",  "#7DACD8", "#FFB090",   "#B89FC8"   
)

DimPlot(merge_obj,group.by = "celltype", cols = colors,label=F )

selected_features = c(
  "Aqp4", "Slc1a3",           ## Astrocytes
  "Slc17a7","Slc17a6",        ## EXC neuron
  "Gad1", "Gad2",             ## INH neuron 
  "Ttr","Ecrg4",              ## NPC
  "Itgam","Aif1",             ## Microglia
  "Mog","Mag","Ermn",         ## Oligodendrocytes
  "Cspg4", "Vcan",            ## OPC
  "Cd3e" , "Nkg7" ,           ## Immune cells
  "Pecam1", "Cldn5" ,         ## Endothelial
  "Col1a2","Pdgfrb",          ## Pericytes
  "Foxj1", "Ccdc153","Cfap65","Pax6","Hes5"       ## Ependymal cells
)


merge_obj$celltype <- factor(merge_obj$celltype,
                             levels = c("Astrocytes",
                                        "Excitatory_neurons",
                                        "Inhibitory_neurons",
                                        "NPC",
                                        "Microglia",
                                        "Oligodendrocytes",
                                        "OPCs",
                                        "Immune_cells",
                                        "Endothelial_cells",
                                        "Pericytes",
                                        "Ependymal_cells"
                             ))


p1 <- DoHeatmap(merge_obj, features = selected_features,group.by="celltype",
                group.colors = colors
)+ 
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))

VlnPlot(merge_obj,
        features = selected_features,group.by="celltype",raster=T,
        flip = F,stack = T,fill.by ="ident", cols =  rev(colors))

p2 <- DotPlot(object = merge_obj, 
              features = selected_features, cols =c("white","firebrick3"),
              group.by = 'celltype',scale.min=30) + RotatedAxis() 

cowplot::plot_grid(p1,p2,ncol=1)





