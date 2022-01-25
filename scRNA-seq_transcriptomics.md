## Project: Fetal skin fibroblast
#### Running title: " Comparative transcriptomics profile of FS HGG and MS"
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
```
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
                                }
```
### 6. Data transformation and scaling using 'SCT' module in Seurat
```
#SCT transformation

for (i in 1:length(sample.list)) {
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
sample.features <- SelectIntegrationFeatures(object.list = sample.list, nfeatures = 3000)
sample.list <- lapply(X = sample.list, FUN = function(x) {
                                x <- RunPCA(x, features = sample.features)
                                                         } )
sample.list <- PrepSCTIntegration(object.list = sample.list, anchor.features = sample.features, verbose = FALSE)
sample.anchors <- FindIntegrationAnchors(object.list = sample.list, normalization.method = "SCT", anchor.features = sample.features, reference = c(1, 2), reduction = "rpca", verbose = TRUE)
#reference = c(1, 2) additional option
sample.integrated <- IntegrateData(anchorset = sample.anchors, normalization.method = "SCT", verbose = TRUE)

### 8. Save RDS file
```
setwd ("/path/")
saveRDS (NS9, file = "NS9.rds")
NS9=readRDS(file="./path/NS9.rds")
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
