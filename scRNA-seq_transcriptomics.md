## Project: Comparative transcriptomics profile of FS HGG and MS
#### Running title: "Single cell RNA study identifies an analogous genetic profile of fetal skin and gingiva fibroblast involved in tissue regeneration"
#### Author: "Rajneesh Srivastava"
#### Date: "23/01/2022"

### 1. Setup the directory for data analysis in R (version 4.0.4)
```setwd ("/path/public_dataset/Result/")```

### 2. Software installation
```
install.packages("Seurat")   # version 4.0.2
install.packages("hdf5r")
install.packages("ggplot2")  # version 3.3.4
install.packages("future")
```
### 3. Load libraries in R
```
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(hdf5r)
library(patchwork)
library(future)
library(sctransform)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 15000 * 1024^2)
```

### 4. Upload Data
Data were downloaded from Gene Expression Omnibus (GEO) repository with following GSE IDs
```
GSE-ID		Groups	Human Tissue
GSE156972	FS	Fetal Skin
GSE164241	HGG	Gingiva
GSE153760	MS	Biopsy -skin
GSE158924	MS	arm -skin
GSE147944	MS	Biopsy -skin
In-house	MS	Biopsy -skin
```
#### Fetal_Skin
```
FET.dir = "C:/path/public_dataset/Fetal_GSE156972/"

samples <- read.table("/path/public_dataset/Fetal_GSE156972/LibraryID.txt", stringsAsFactors=F, sep="\t",header=T)

FET_data <- Read10X_h5(file="/path/public_dataset/Fetal_GSE156972/GSE156972_raw_gene_bc_matrices_h5.h5")

cells <- new("seurat", raw.data = FET_data)
cellcodes <- as.data.frame(cells@raw.data@Dimnames[[2]])
colnames(cellcodes) <- "barcodes"
rownames(cellcodes) <- cellcodes$barcodes
cellcodes$libcodes <- as.factor(gsub(pattern=".+-", replacement="", cellcodes$barcodes))
cellcodes$samples <- as.vector(samples$library_id[cellcodes$libcodes])
sampleidentity <- cellcodes["samples"]

GSE156972 <- CreateSeuratObject(FET_data,
				    meta.data=sampleidentity,
                                    min.cells=3,
                                    min.features = 200,
                                    project = "GSE156972")

```
#### Gingiva
```
HGG.dir = "/path/public_dataset/GSE164241_HGG/healthy/"

GSE164241.list=list("GSM5005048","GSM5005049","GSM5005050","GSM5005051","GSM5005052","GSM5005053","GSM5005054","GSM5005055","GSM5005056","GSM5005057","GSM5177039","GSM5177040","GSM5177041")


for (file in GSE164241.list){
               OC_data <- Read10X(data.dir =    
                                    paste0(OC.dir, file))
               
               OC_obj <- CreateSeuratObject(counts = 
                                    OC_data,
                                    min.cells=3,
                                    min.features = 200,
                                    project = file)
               assign(file, OC_obj)
                            }
```
#### Mature Skin 1 - GSE158924
```
MS1.dir = "/path/public_dataset/HC_GSE158924/"

GSE158924.list=list("GSM4815804","GSM4815805")

for (file in GSE158924.list){
               MS1_data <- read.table(file = 
                                    paste0(MS1.dir, file,".txt.gz"))
               
               MS1_obj <- CreateSeuratObject(counts = 
                                    MS1_data,
                                    min.cells=3,
                                    min.features = 200,
                                    project = file)
               assign(file, MS1_obj)
                            }
```
#### Mature Skin 2 - GSE153760
```
MS2.dir = "/path/public_dataset/MS_GSE153760/"

GSE153760.list=list("GSM4653868","GSM4653869")

for (file in GSE153760.list){
               MS2_data <- Read10X(data.dir =    
                                    paste0(MS2.dir, file))
               
               MS2_obj <- CreateSeuratObject(counts = 
                                    MS2_data,
                                    min.cells=3,
                                    min.features = 200,
                                    project = file)
               assign(file, MS2_obj)
                            }
```
#### Mature Skin 3 - GSE147944
```
MS3.dir = "/path/public_dataset/MS_GSE147944/"

GSE147944.list=list("GSM4450726","GSM4450727","GSM4450728","GSM4450729")

for (file in GSE147944.list){
               MS3_data <- read.csv(file = 
                                    paste0(MS3.dir, file,"_raw.csv.gz"), header = TRUE, row.names = 1)
               
               MS3_obj <- CreateSeuratObject(counts = 
                                    MS3_data,
                                    min.cells=3,
                                    min.features = 200,
                                    project = file)
               assign(file, MS3_obj)
                            }
```
#### Mature Skin 4 - In-house
```
MS4.dir = "/path/public_dataset/MS_Sen-Lab-data/dataset/"

Sen.list=list("SCNS1","SCNS2","SCNS8","SCNS9","SCNS10")

for (file in Sen.list){
               MS4_data <- Read10X(data.dir =  
                                    paste0(MS4.dir, file))
               
               MS4_obj <- CreateSeuratObject(counts = 
                                    MS4_data,
                                    min.cells=3,
                                    min.features = 200,
                                    project = file)
               assign(file, MS4_obj)
                            }

```
### 5. Compute the mitochondrial content and filter with additional cutoffs
```
sample.list=list(GSE156972,GSM5005048,GSM5005049,GSM5005050,GSM5005051,GSM5005052,GSM5005053,GSM5005054,GSM5005055,GSM5005056,GSM5005057,GSM5177039,GSM5177040,GSM5177041,GSM4815804,GSM4815805,GSM4653868,GSM4653869,GSM4450726,GSM4450727,GSM4450728,GSM4450729,SCNS1,SCNS2,SCNS8,SCNS9,SCNS10)

for (i in 1:length(sample.list)) {
    sample.list[[i]][["percent.mt"]] <-
            PercentageFeatureSet(sample.list[[i]], pattern = "^MT-")
    sample.list[[i]] <- subset(sample.list[[i]], 
                            subset = nFeature_RNA > 200 & 
                            nFeature_RNA < 5000 & percent.mt < 15 & 
                            nCount_RNA < 25000 & nCount_RNA > 2000)
                                }
```
### 6. Data transformation and scaling using 'SCT' module in Seurat

#### SCT transformation
```
for (i in 1:length(sample.list)) {
     sample.list[[i]] <- SCTransform(sample.list[[i]],
                                    vars.to.regress = "percent.mt", 
                                    return.only.var.genes = FALSE,
                                    verbose = FALSE)
                                 }
```
### 7. Check if the sample meets the analysis criteria (Please see methods)
```
VlnPlot(sample.list[[1]],                #all files can be checked from 1-29 or can be looped
	features = c("nCount_Spatial",
			"nFeature_Spatial",     
			"percent.mt"), 
			 ncol = 3, 
			 pt.size=0)
```
### 8. Data integration
```
sample.features <- SelectIntegrationFeatures(object.list = sample.list,    
                   nfeatures = 3000)

sample.list <- lapply(X = sample.list, FUN = function(x) {
                   x <- RunPCA(x, features = sample.features)
                                                         } )

sample.list <- PrepSCTIntegration(object.list = sample.list, 
                   anchor.features = sample.features,
                   verbose = FALSE)
                   
sample.anchors <- FindIntegrationAnchors(object.list = sample.list, 
                   normalization.method = "SCT", 
                   anchor.features = sample.features, 
                   reference = c(1, 2), 
                   reduction = "rpca", 
                   verbose = TRUE)

sample.integrated <- IntegrateData(anchorset = sample.anchors, 
                   normalization.method = "SCT",
                   verbose = TRUE)
```
### 9. Clustering analysis

```
sample.integrated <- RunPCA(object = sample.integrated, verbose = FALSE) 
```
###### Test the %variance of principle components (ranked)
```
ElbowPlot(sample.integrated, ndims = 50)
ElbowPlot(sample.integrated, ndims = 20) #after manual judgement
```
###### Explore the top PCs
```
DimHeatmap(sample.integrated, dims = 1, cells = 500, balanced=TRUE) #1st PC
DimHeatmap(sample.integrated, dims = 2, cells = 500, balanced=TRUE) #2nd PC
```
###### Run tSNE/UMAP
```
sample.integrated = RunUMAP(sample.integrated, dims = 1:30) # optional
sample.integrated = RunTSNE(sample.integrated, dims = 1:30)
```
### 10. Find cluster with defined resolution
```
sample.integrated <- FindNeighbors(sample.integrated) # dim=1:10
sample25 <- FindClusters(sample.integrated, resolution = 0.25)
#sample50 <- FindClusters(sample.integrated, resolution = 0.50)
```
### 11. Find average expression of integrated and clustered object [sample25]
```
AvgExpS25 = AverageExpression(sample25, return.seurat = FALSE, verbose = TRUE)
#head (AvgExpS25$integrated)
#head (AvgExpS25$RNA)
#head(AvgExpS25$SCT)
```
### 12. Save RDS file
```
#setwd ("/path/")
#write.table(AvgExpS25$SCT,"AvgExpS25_SCT.txt)
#saveRDS (sample.integrated, file = "sample.integrated.rds")
#saveRDS (sample25, file = "FS_HGG_MS_25.rds")

```
### 13. Transcriptomic profile of FS, HGG, MS groups and respective signature genes for cell clusters.
###### DOT Plot for Markers
```
DefaultAssay(sample25)="SCT"
DimPlot(sample25,reduction = "tsne")
DimPlot(sample25,reduction = "tsne", split.by="Type")
#Idents(sample25)="seurat_clusters"
sample25_markers=FindAllMarkers(sample25,only.pos = T,logfc.threshold = 0.30, min.pct = 0.10)
#write.table(sample25_markers,"All_marker_sample25.txt") # Explore the table
top10 <- sample25_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#head(top10)
top10genes=unique(top10$gene)
DotPlot(FB_MYE_re25, features = top10genes) + 
    theme(axis.text.x=element_text(size=7.5, angle=45, hjust=1)) + 
    theme(axis.text.y=element_text(size=7.5, face="italic"))
```
###### Signature genes for cell clusters
```
FigS1=c("COL1A1","SELE","IL7R","KRT1","ACTA2","LYZ","KRT14","KRT6A","IGFBP3","TYRP1","TFF3","IGHG1","CTSG")
FeaturePlot(sample25,features=FigS1,reduction="tsne",cols=c("grey","red","black"))

fb=c("DCN","COL1A1","CFD","FBLN1")
FeaturePlot(sample25, reduction="tsne", features=fb, cols=c("grey","red","black"))

mye=c("LYZ","CXCL8","HLA-DQB1","IL1B")
FeaturePlot(sample25, reduction="tsne", features=mye, cols=c("grey","red","black"))
```
## Cell-Chat-Analysis
```
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
future::plan("multiprocess", workers = 4)
```
###### INPUT Seurat object
```
Idents(sample25)="Type"
FS_s25=subset(sample25,subset=Type=="FS")
HGG_s25=subset(sample25,subset=Type=="HGG")
MS_s25=subset(sample25,subset=Type=="MS")
```
###### Create CellChat objects
```
FS_cell_chat=createCellChat(object=FS_s25,meta=FS_s25@meta.data,group.by="seurat_clusters")
HGG_cell_chat=createCellChat(object=HGG_s25,meta=HGG_s25@meta.data,group.by="seurat_clusters")
MS_cell_chat=createCellChat(object=MS_s25,meta=MS_s25@meta.data,group.by="seurat_clusters")
```
###### Assign Idents
```
FS_cell_chat <-setIdent(FS_cell_chat, ident.use="seurat_clusters")
HGG_cell_chat <-setIdent(HGG_cell_chat, ident.use="seurat_clusters")
MS_cell_chat <-setIdent(MS_cell_chat, ident.use="seurat_clusters")
#levels(FS_cell_chat@idents)
#FS_groupSize<-as.numeric(table(FS_cell_chat@idents)) #number of cells in each grp
```
###### Database usage
```
CellChatDB<- CellChatDB.human
#showDatabaseCategory(CellChatDB)
CellChatDB.use<- CellChatDB
FS_cell_chat@DB <-CellChatDB.use
HGG_cell_chat@DB <-CellChatDB.use
MS_cell_chat@DB <-CellChatDB.use
```
###### Subset the dataset for enhanced performance
```
FS_cell_chat <-subsetData(FS_cell_chat)
HGG_cell_chat <-subsetData(HGG_cell_chat)
MS_cell_chat <-subsetData(MS_cell_chat)
```
##### Individual Analysis
```
FS_cell_chat <- identifyOverExpressedGenes(FS_cell_chat)
FS_cell_chat <- identifyOverExpressedInteractions(FS_cell_chat)
FS_cell_chat <- projectData(FS_cell_chat, PPI.human) #(optional)
FS_cell_chat <- computeCommunProb(FS_cell_chat, raw.use = TRUE)
FS_cell_chat <- filterCommunication(FS_cell_chat,min.cells = 10)
FS_cell_chat <- computeCommunProbPathway(FS_cell_chat)
FS_cell_chat <- netAnalysis_computeCentrality(FS_cell_chat, slot.name = "netP")
FS_cell_chat <- aggregateNet(FS_cell_chat)

HGG_cell_chat <- identifyOverExpressedGenes(HGG_cell_chat)
HGG_cell_chat <- identifyOverExpressedInteractions(HGG_cell_chat)
HGG_cell_chat <- projectData(HGG_cell_chat, PPI.human) #(optional)
HGG_cell_chat <- computeCommunProb(HGG_cell_chat, raw.use = TRUE)
HGG_cell_chat <- filterCommunication(HGG_cell_chat,min.cells = 10)
HGG_cell_chat <- computeCommunProbPathway(HGG_cell_chat)
HGG_cell_chat <- netAnalysis_computeCentrality(HGG_cell_chat, slot.name = "netP")
HGG_cell_chat <- aggregateNet(HGG_cell_chat)

MS_cell_chat <- identifyOverExpressedGenes(MS_cell_chat)
MS_cell_chat <- identifyOverExpressedInteractions(MS_cell_chat)
MS_cell_chat <- projectData(MS_cell_chat, PPI.human) #(optional)
MS_cell_chat <- computeCommunProb(MS_cell_chat, raw.use = TRUE)
MS_cell_chat <- filterCommunication(MS_cell_chat,min.cells = 10)
MS_cell_chat <- computeCommunProbPathway(MS_cell_chat)
MS_cell_chat <- netAnalysis_computeCentrality(MS_cell_chat, slot.name = "netP")
MS_cell_chat <- aggregateNet(MS_cell_chat)
```
##### cell_chat_objects
###### save files
```
FS_cell_chat=readRDS(file = "FS_cell_chat.rds")
HGG_cell_chat=readRDS(file = "HGG_cell_chat.rds")
MS_cell_chat=readRDS(file = "MS_cell_chat.rds")
merged_cellchat=readRDS(file = "Merged_cell_chat.rds")
```
##### Comparison of cell chat objects
###### S2A
```
g1=compareInteractions(merged_cellchat, show.legend = F, group = c(1,2,3))
g2=compareInteractions(merged_cellchat, show.legend = F, group = c(1,2,3), measure = "weight")
g1+g2
```
###### S2B
```
object.list=c(FS_cell_chat,HGG_cell_chat,MS_cell_chat)
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,3), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[[i]]))
}
```
###### S2C
```
object.list=c(FS_cell_chat,MS_cell_chat)
group.cellType <- c(rep("0", 4), rep("5", 4))
group.cellType <- factor(group.cellType, levels = c("0", "5"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)

object.list=c(HGG_cell_chat,MS_cell_chat)
group.cellType <- c(rep("0", 4), rep("5", 4))
group.cellType <- factor(group.cellType, levels = c("0", "5"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
```
##### S2D - Information_flow
```
if1=rankNet(merged_cellchat, comparison=c(1,3),stacked = T, do.stat = TRUE)
if2=rankNet(merged_cellchat, comparison=c(1,3),stacked = T, do.stat = TRUE)
if3=rankNet(merged_cellchat, comparison=c(1,2,3),stacked = T, do.stat = TRUE)
#if1+if2
if3
```
##### Compare_interactions
###### heatmap
```
netVisual_heatmap(merged_cellchat,comparison=c(1,3))
netVisual_heatmap(merged_cellchat,comparison=c(2,3))
netVisual_heatmap(merged_cellchat,comparison=c(1,3),measure = "weight")
netVisual_heatmap(merged_cellchat,comparison=c(2,3),measure = "weight")
```
##### network
```
netVisual_diffInteraction(merged_cellchat, weight.scale = T,comparison = c(1,3))
netVisual_diffInteraction(merged_cellchat, weight.scale = T,comparison = c(2,3))
```
###### S2E-F
```
library(ComplexHeatmap)
object.list=c(FS_cell_chat,HGG_cell_chat,MS_cell_chat)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 16)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 16)
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+2], width = 5, height = 16)
draw(ht1 + ht2 + ht3, ht_gap = unit(0.5, "cm"))
```
###### SX
```
FSMS.chat.list=c(FS_cell_chat,MS_cell_chat)
chat.names=list("FS","MS")
FSvsMSchat <- mergeCellChat(FSMS.chat.list, add.names = names(chat.names))
#par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(FSvsMSchat, weight.scale = T)
netVisual_diffInteraction(FSvsMSchat, weight.scale = T, measure = "weight")

HGGMS.chat.list=c(HGG_cell_chat,MS_cell_chat)
chat2.names=list("HGG","MS")
HGGvsMSchat <- mergeCellChat(HGGMS.chat.list, add.names = names(chat2.names))
#par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(HGGvsMSchat, weight.scale = T)
netVisual_diffInteraction(HGGvsMSchat, weight.scale = T, measure = "weight")
```
##### Cluster specific analysis
```
mat <- FS_cell_chat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
```
##### Pathway specific analysis
```
pathways.show <- c("ANNEXIN") 
weight.max <- getMaxWeight(merged_cellchat, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
object.list=c(FS_cell_chat,HGG_cell_chat)#,MS_cell_chat)
for (i in 1:length(object.list)) {
              netVisual_aggregate(object.list[[i]], 
                                  signaling = pathways.show, 
                                  layout = "circle")
                                 }
#OR

netVisual_aggregate(FS_cell_chat, signaling = "ANNEXIN",layout="circle")
netVisual_aggregate(HGG_cell_chat, signaling = "ANNEXIN",layout="circle")
```
## Cell type specific Differential Expression Analysis
##### Isolate Fb and Myeloid cell clusters
```
Idents(sample25)="seurat_clusters"
FB=subset(sample25,subset=seurat_clusters=="0")
MYE=subset(sample25,subset=seurat_clusters=="5")
```
###### Cell type specific Differential Expression Analysis
```
Idents(FB)="Type"
DefaultAssay(FB)="SCT"
FB_FSvsMS <- FindMarkers(FB, ident.1 = "FS", ident.2 = "MS",logfc.threshold = 0.6, min.pct=0.30)
Idents(MYE)="Type"
DefaultAssay(MYE)="SCT"
MYE_FSvsMS <- FindMarkers(MYE, ident.1 = "FS", ident.2 = "MS",logfc.threshold = 0.6, min.pct=0.30)
```
###### Intersection of DEGs across cell types (refer to Fig. 4)
```
FB_MYE=subset(sample25,subset=seurat_clusters=="0"|seurat_clusters=="5")

g31=c("ASPN","ATP5G1","CCL2","CLEC11A","COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","FAM195B","IGF2","LAPTM4A","LHFP","LOXL2","LUM","MDK","NGFRAP1","NREP","OSTC","PHLDA1","POSTN","PTK7","PTRF","SEPP1","SEPW1","SERPINE2","SNAI2","SPARC","TMED2","VIMP","WBP5")
g35=c("AC090498.1","ALDOA","ATP5A1","ATP5B","ATP5C1","ATP5D","ATP5E","ATP5F1","ATP5G2","ATP5G3","ATP5H","ATP5I","ATP5J","ATP5J2","ATP5L","ATP5O","ATPIF1","C11orf31","C14orf166","C14orf2","C19orf43","C19orf60","C7orf73","GNB2L1","GPX1","LINC00493","MYEOV2","RPL13A","RPS17","SELK","SEP15","SHFM1","TCEB1","TCEB2","USMG5")
g25=c("BRK1","CD69","CORO1A","COX7A2","CTD-3252C9.4","CXCR4","FAM26F","FYB","GLTSCR2","PFN1","POMP","RPL21","RPL27A","RPL31","RPL36A","RPL7","RPS19","S100A8","S100A9","SEC61G","SELT","SERF2","SLC16A3","TMSB4X","UQCR11")

VlnPlot(FB_MYE,features=g25,split.by="Type",stack=T,flip=T,cols=c("#F8766D","#00BA38","#619CFF"))
VlnPlot(FB_MYE,features=g31,split.by="Type",stack=T,flip=T,cols=c("#F8766D","#00BA38","#619CFF"))
VlnPlot(FB_MYE,features=g35,split.by="Type",stack=T,flip=T,cols=c("#F8766D","#00BA38","#619CFF"))
```
##### FB_MYE clusters
###### tweak-in for assigning the celltype
```
meta=read.table("metadata.txt",sep="\t", header=T)

GSM=FB_MYE@meta.data
GSM$celltype=1
head(GSM)

for(i in 1:nrow(GSM)){
  for (j in 1:nrow(meta)){
   if (GSM[i,9] == meta[j,1])
        { GSM[i,11] = meta[j,2] }
                 }
             }

FB_MYE@meta.data=GSM
```
#### ReCluster ANALYSIS
```
#FB_MYE=readRDS(file="FB_MYE.rds")
DefaultAssay(FB_MYE)="integrated"
FB_MYE_re <- RunPCA(object = FB_MYE, verbose = FALSE)
#ElbowPlot(FB_MYE_re, ndims = 50)
FB_MYE_re = RunUMAP(FB_MYE_re, dims = 1:30)
FB_MYE_re = RunTSNE(FB_MYE_re, dims = 1:30)
FB_MYE_re <- FindNeighbors(FB_MYE_re) # dim=1:10
FB_MYE_re25 <- FindClusters(FB_MYE_re, resolution = 0.25)
#saveRDS(FB_MYE_re25,file="FB_MYE_re25.rds")
```
###### Visualization of fibroblast and myeloid cell sub clusters 
```
#FB_MYE_re25=readRDS(file="FB_MYE_re25.rds")
DefaultAssay(FB_MYE_re25)="SCT"
DimPlot(FB_MYE_re25,reduction = "tsne")
DimPlot(FB_MYE_re25,split.by = "Type",reduction = "tsne")
```
###### Get top 10 Markers
```
DefaultAssay(FB_MYE_re25)="SCT"
FB_MYE_re25_markers=FindAllMarkers(FB_MYE_re25,only.pos = T,logfc.threshold = 0.3, min.pct = 0.10)
#write.table(FB_MYE_re25_markers,"All_marker_FB_MYE_re25.txt") # Explore the table
Idents(FB_MYE_re25)="seurat_clusters"
top10 <- FB_MYE_re25_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#head(top10)
top10genes=unique(top10$gene)
DotPlot(MY_5_re50, features = top10genes) + 
    theme(axis.text.x=element_text(size=7.5, angle=45, hjust=1)) + 
    theme(axis.text.y=element_text(size=7.5, face="italic"))
```
### Pseudotime Analysis
##### Create monocle object from seurat object
```
#setwd("C:/Users/rsrivast/Desktop/scRNA/scRNA-Fetal-vs-Adult/Results_July/Saved_RDS_final/")
#sample25 <-readRDS(file = "FS_HGG_MS_25.rds")
#FB_MYE_re25=readRDS(file="FB_MYE_re25.rds")
```
## Thank you
