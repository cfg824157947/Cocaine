library(Seurat)
library(Signac)
library(rtracklayer)
library(dplyr)
library(ggplot2)


print(packageVersion("Seurat"))

filename <- "/data2/duren_lab/cham/cocain/filtered_feature_bc_matrix.h5"

inputdata.10x <- Read10X_h5(filename, use.names = TRUE, unique.features = TRUE)

                                        # extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
SeuratO<- CreateSeuratObject(counts = rna_counts)
SeuratO[["percent.mt"]] <- PercentageFeatureSet(SeuratO, pattern = "^mt:")
                                        # the 10x hdf5 file contains both data types. 
meta <- read.table("/data2/duren_lab/cham/cocain/Data/meta_data.tsv", sep= '\t', header=TRUE, row.names=1)

SeuratO$filter <- colnames(SeuratO) %in% rownames(meta)
SeuratO$sample_name <- meta$sample_name
SeuratO$cluster <- meta['cluster']
SeuratO$sample <- meta$sample

SeuratO <- subset(SeuratO, subset = filter == TRUE)

DefaultAssay(SeuratO) <- "RNA"
SeuratO <- SCTransform(SeuratO, verbose = FALSE)

saveRDS(SeuratO, "/data2/duren_lab/cham/cocain/Data/Coc_seuratObj.rds")
