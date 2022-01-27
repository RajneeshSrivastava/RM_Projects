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
In-house	MS	NA
GSE158924	MS	arm -skin
GSE147944	MS	Biopsy -skin
```
#### Data upload

# Fetal_Skin
```
FET.dir = "C:/Users/rsrivast/Desktop/scRNA/Fetal_GSE156972/"

samples <- read.table("C:/Users/rsrivast/Desktop/scRNA/Fetal_GSE156972/LibraryID.txt", stringsAsFactors=F, sep="\t",header=T)

FET_data <- Read10X_h5(file="C:/Users/rsrivast/Desktop/scRNA/Fetal_GSE156972/GSE156972_raw_gene_bc_matrices_h5.h5")

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
#Gingiva

```
HGG.dir = "C:/Users/rsrivast/Desktop/scRNA/GSE164241_HGG/healthy/"

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
# Mature Skin 1 - GSE158924

```
MS1.dir = "C:/Users/rsrivast/Desktop/scRNA/HC_GSE158924/"

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
# Mature Skin 2 - GSE153760

```
MS2.dir = "C:/Users/rsrivast/Desktop/scRNA/HC_GSE153760/"

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
##### Mature Skin - GSE147944]
```
MS3.dir = "C:/Users/rsrivast/Desktop/scRNA/MS_GSE147944/"

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
##### Mature Skin 4 - In-house]
```
MS4.dir = "C:/Users/rsrivast/Desktop/scRNA/MS_Sen-Lab-data/dataset/"

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

RS.dir = "/path/"                                             
RS_data <- Load10X_Spatial(data.dir = paste0(RS.dir,
						filename = "filtered_feature_bc_matrix.h5", 
						assay = "Spatial", 
						slice = NS9, 
						filter.matrix = TRUE)
            assign (NS009, KS_data)
```
### 5. Compute mitochondrial content and filter
```
sample.list=list(NS1,NS2,IMQ1,IMQ2,SHP1,SHP2)

for (i in 1:length(sample.list)) {
    sample.list[[i]][["percent.mt"]] <-
            PercentageFeatureSet(sample.list[[i]], pattern = "^MT-")
    sample.list[[i]] <- subset(sample.list[[i]], 
                            subset = nFeature_RNA > 200 & 
                            nFeature_RNA < 5000 & percent.mt < 15 & 
                            nCount_RNA < 25000 & nCount_RNA > 2000)
#                                }
### 6. Data transformation and scaling using 'SCT' module in Seurat
#SCT transformation

#for (i in 1:length(sample.list)) {
     sample.list[[i]] <- SCTransform(sample.list[[i]],
                                    vars.to.regress = "percent.mt", 
                                    return.only.var.genes = FALSE,
                                    verbose = FALSE)
                                 }
```
### 7. Check if the sample meets the analysis criteria (Please see methods)
```
VlnPlot(NS9, 
	features = c("nCount_Spatial",
			"nFeature_Spatial",     
			"percent.mt"), 
			 ncol = 3, 
			 pt.size=0)
```

#Attempt 2 (refer to May18.rmd)
# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features

# Data integration

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

# CLUSTERING ANALYSIS - RUN PCA
sample.integrated <- RunPCA(object = sample.integrated, verbose = FALSE) 

#test the %variance of principle components (ranked)
ElbowPlot(sample.integrated, ndims = 50)
ElbowPlot(sample.integrated, ndims = 20) #after manual judgement

## Explore the top PCs
DimHeatmap(sample.integrated, dims = 1, cells = 500, balanced=TRUE) #1st PC
DimHeatmap(sample.integrated, dims = 2, cells = 500, balanced=TRUE) #2nd PC

#RUN tSNE/UMAP
sample.integrated = RunUMAP(sample.integrated, dims = 1:30)
#sample.integrated = RunTSNE(sample.integrated, dims = 1:30)

## FIND CLUSTERS WITH DEFINED RESOLUTION
sample.integrated <- FindNeighbors(sample.integrated) # dim=1:10
sample25 <- FindClusters(sample.integrated, resolution = 0.25)
#sample50 <- FindClusters(sample.integrated, resolution = 0.50)

## FIND AVERAGE EXPRESSION OF INTEGRATED OBJECT [sample25]
AvgExpS25 = AverageExpression(sample25, return.seurat = FALSE, verbose = TRUE)
head (AvgExpS25$integrated)
head (AvgExpS25$RNA)
head(AvgExpS25$SCT)

### 8. Save RDS file
```
#setwd ("/path/")
#write.table(AvgExpS25$SCT,"AvgExpS25_SCT.txt)
#saveRDS (sample.integrated, file = "sample.integrated.rds")
#saveRDS (sample25, file = "FS_HGG_MS_25.rds")

```
### 9. Spatial transcriptomic profile of signature genes for Kera1 and Kera 2 clusters.
```
Kera1=c("KRT14","KRT1")
Kera2=c("KRT19","KRT7")
SpatialFeaturePlot(NS9, features = Kera1)
SpatialFeaturePlot(NS9, features = Kera2)
```
### 10. Spatial transcriptome profile (and localization) of top transcription factors,enzymes and metabolically active genes with increased expression levels in Kera 2
```
Genes=c("GAPDH", "CITED4", "COX7C", "COX8A", "COX5B", "NDUFA4", "COX7B", "PRDX5", "COX6A1")
SpatialFeaturePlot(NS9, features = Genes)
```
## Thank you
