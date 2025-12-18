#Set Environmental Variables
Sys.setenv(RETICULATE_PYTHON = "/usr/bin/python3")
library(cowplot)
library(ggplot2)
library(dplyr)
library(Seurat)
library(rlang)
#library(TopKLists)
#library(readxl)
library(gtools)
library(openxlsx)
library(readr)
options(future.globals.maxSize= 600000000000)
setwd("/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Analysis/Gene_expression/")
load("Input_data_wrkspc.RData")

#Read data in from filtered feature counts BC matrix. Because these are multiome data, there will two separate matrices per library and Seurat will warn you as such.
RP0457_Coc_F_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_1_8/RP0457_Coc_F_r2/outs/filtered_feature_bc_matrix/')
RP0231_Suc_M_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_1_8/RP0231_Suc_M_r2/outs/filtered_feature_bc_matrix/')
RP0609_Coc_M_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_1_8/RP0609_Coc_M_r2/outs/filtered_feature_bc_matrix/')
RP0231_Suc_F_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_1_8/RP0231_Suc_F_r2/outs/filtered_feature_bc_matrix/')
RP0446_Coc_M_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_1_8/RP0446_Coc_M_r2/outs/filtered_feature_bc_matrix/')
RP0457_Coc_M_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_1_8/RP0457_Coc_M_r2/outs/filtered_feature_bc_matrix/')
RP0897_Suc_M_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_1_8/RP0897_Suc_M_r2/outs/filtered_feature_bc_matrix/')
RP0897_Coc_M_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_1_8/RP0897_Coc_M_r1/outs/filtered_feature_bc_matrix/')
RP0231_Coc_F_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_1_8/RP0231_Coc_F_r1/outs/filtered_feature_bc_matrix/')
RP0897_Coc_F_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_09_16/RP0897_Coc_F_r2/outs/filtered_feature_bc_matrix/')
RP0446_Suc_M_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_09_16/RP0446_Suc_M_r2/outs/filtered_feature_bc_matrix/')
RP0457_Coc_F_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_09_16/RP0457_Coc_F_r1/outs/filtered_feature_bc_matrix/')
RP0457_Suc_M_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_09_16/RP0457_Suc_M_r2/outs/filtered_feature_bc_matrix/')
RP0883_Suc_M_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_09_16/RP0883_Suc_M_r1/outs/filtered_feature_bc_matrix/')
RP0457_Coc_M_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_09_16/RP0457_Coc_M_r1/outs/filtered_feature_bc_matrix/')
RP0457_Suc_F_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_09_16/RP0457_Suc_F_r2/outs/filtered_feature_bc_matrix/')
RP0609_Coc_F_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_09_16/RP0609_Coc_F_r1/outs/filtered_feature_bc_matrix/')
RP0897_Suc_F_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_17_24/RP0897_Suc_F_r1/outs/filtered_feature_bc_matrix/')
RP0231_Suc_M_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_17_24/RP0231_Suc_M_r1/outs/filtered_feature_bc_matrix/')
RP0446_Suc_F_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_17_24/RP0446_Suc_F_r2/outs/filtered_feature_bc_matrix/')
RP0446_Coc_M_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_17_24/RP0446_Coc_M_r1/outs/filtered_feature_bc_matrix/')
RP0883_Suc_F_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_17_24/RP0883_Suc_F_r1/outs/filtered_feature_bc_matrix/')
RP0231_Suc_F_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_17_24/RP0231_Suc_F_r1/outs/filtered_feature_bc_matrix/')
RP0897_Suc_M_r3_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_17_24/RP0897_Suc_M_r3/outs/filtered_feature_bc_matrix/')
RP0446_Coc_F_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_17_24/RP0446_Coc_F_r1/outs/filtered_feature_bc_matrix/')
RP0897_Suc_F_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_25_32/RP0897_Suc_F_r2/outs/filtered_feature_bc_matrix/')
RP0883_Coc_F_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_25_32/RP0883_Coc_F_r2/outs/filtered_feature_bc_matrix/')
RP0457_Suc_F_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_25_32/RP0457_Suc_F_r1/outs/filtered_feature_bc_matrix/')
RP0883_Suc_F_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_25_32/RP0883_Suc_F_r2/outs/filtered_feature_bc_matrix/')
RP0231_Coc_M_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_25_32/RP0231_Coc_M_r2/outs/filtered_feature_bc_matrix/')
RP0609_Suc_M_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_25_32/RP0609_Suc_M_r1/outs/filtered_feature_bc_matrix/')
RP0897_Coc_F_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_25_32/RP0897_Coc_F_r1/outs/filtered_feature_bc_matrix/')
RP0883_Coc_M_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_25_32/RP0883_Coc_M_r2/outs/filtered_feature_bc_matrix/')
RP0457_Suc_M_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_33_40/RP0457_Suc_M_r1/outs/filtered_feature_bc_matrix/')
RP0609_Coc_F_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_33_40/RP0609_Coc_F_r2/outs/filtered_feature_bc_matrix/')
RP0231_Coc_M_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_33_40/RP0231_Coc_M_r1/outs/filtered_feature_bc_matrix/')
RP0609_Suc_F_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_33_40/RP0609_Suc_F_r1/outs/filtered_feature_bc_matrix/')
RP0446_Suc_M_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_33_40/RP0446_Suc_M_r1/outs/filtered_feature_bc_matrix/')
RP0897_Coc_M_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_33_40/RP0897_Coc_M_r2/outs/filtered_feature_bc_matrix/')
RP0609_Suc_M_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_33_40/RP0609_Suc_M_r2/outs/filtered_feature_bc_matrix/')
RP0883_Coc_F_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_33_40/RP0883_Coc_F_r1/outs/filtered_feature_bc_matrix/')
RP0609_Suc_F_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_41_48/RP0609_Suc_F_r2/outs/filtered_feature_bc_matrix/')
RP0883_Coc_M_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_41_48/RP0883_Coc_M_r1/outs/filtered_feature_bc_matrix/')
RP0446_Suc_F_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_41_48/RP0446_Suc_F_r1/outs/filtered_feature_bc_matrix/')
RP0231_Coc_F_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_41_48/RP0231_Coc_F_r2/outs/filtered_feature_bc_matrix/')
RP0446_Coc_F_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_41_48/RP0446_Coc_F_r2/outs/filtered_feature_bc_matrix/')
RP0897_Suc_M_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_41_48/RP0897_Suc_M_r1/outs/filtered_feature_bc_matrix/')
RP0883_Suc_M_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_41_48/RP0883_Suc_M_r2/outs/filtered_feature_bc_matrix/')
RP0609_Coc_M_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_41_48/RP0609_Coc_M_r1/outs/filtered_feature_bc_matrix/')

#Creat Seurat Objects using input BC matrices only. This is done by adding "$`Gene Expression` to the object name from the previous step.
RP0457_Coc_F_r2 <- CreateSeuratObject(counts=RP0457_Coc_F_r2_Data$`Gene Expression`,project="RP0457_Coc_F",min.cells=5)
RP0231_Suc_M_r2 <- CreateSeuratObject(counts=RP0231_Suc_M_r2_Data$`Gene Expression`,project="RP0231_Suc_M",min.cells=5)
RP0609_Coc_M_r2 <- CreateSeuratObject(counts=RP0609_Coc_M_r2_Data$`Gene Expression`,project="RP0609_Coc_M",min.cells=5)
RP0231_Suc_F_r2 <- CreateSeuratObject(counts=RP0231_Suc_F_r2_Data$`Gene Expression`,project="RP0231_Suc_F",min.cells=5)
RP0446_Coc_M_r2 <- CreateSeuratObject(counts=RP0446_Coc_M_r2_Data$`Gene Expression`,project="RP0446_Coc_M",min.cells=5)
RP0457_Coc_M_r2 <- CreateSeuratObject(counts=RP0457_Coc_M_r2_Data$`Gene Expression`,project="RP0457_Coc_M",min.cells=5)
RP0897_Suc_M_r2 <- CreateSeuratObject(counts=RP0897_Suc_M_r2_Data$`Gene Expression`,project="RP0897_Suc_M",min.cells=5)
RP0897_Coc_M_r1 <- CreateSeuratObject(counts=RP0897_Coc_M_r1_Data$`Gene Expression`,project="RP0897_Coc_M",min.cells=5)
RP0231_Coc_F_r1 <- CreateSeuratObject(counts=RP0231_Coc_F_r1_Data$`Gene Expression`,project="RP0231_Coc_F",min.cells=5)
RP0897_Coc_F_r2 <- CreateSeuratObject(counts=RP0897_Coc_F_r2_Data$`Gene Expression`,project="RP0897_Coc_F",min.cells=5)
RP0446_Suc_M_r2 <- CreateSeuratObject(counts=RP0446_Suc_M_r2_Data$`Gene Expression`,project="RP0446_Suc_M",min.cells=5)
RP0457_Coc_F_r1 <- CreateSeuratObject(counts=RP0457_Coc_F_r1_Data$`Gene Expression`,project="RP0457_Coc_F",min.cells=5)
RP0457_Suc_M_r2 <- CreateSeuratObject(counts=RP0457_Suc_M_r2_Data$`Gene Expression`,project="RP0457_Suc_M",min.cells=5)
RP0883_Suc_M_r1 <- CreateSeuratObject(counts=RP0883_Suc_M_r1_Data$`Gene Expression`,project="RP0883_Suc_M",min.cells=5)
RP0457_Coc_M_r1 <- CreateSeuratObject(counts=RP0457_Coc_M_r1_Data$`Gene Expression`,project="RP0457_Coc_M",min.cells=5)
RP0457_Suc_F_r2 <- CreateSeuratObject(counts=RP0457_Suc_F_r2_Data$`Gene Expression`,project="RP0457_Suc_F",min.cells=5)
RP0609_Coc_F_r1 <- CreateSeuratObject(counts=RP0609_Coc_F_r1_Data$`Gene Expression`,project="RP0609_Coc_F",min.cells=5)
RP0897_Suc_F_r1 <- CreateSeuratObject(counts=RP0897_Suc_F_r1_Data$`Gene Expression`,project="RP0897_Suc_F",min.cells=5)
RP0231_Suc_M_r1 <- CreateSeuratObject(counts=RP0231_Suc_M_r1_Data$`Gene Expression`,project="RP0231_Suc_M",min.cells=5)
RP0446_Suc_F_r2 <- CreateSeuratObject(counts=RP0446_Suc_F_r2_Data$`Gene Expression`,project="RP0446_Suc_F",min.cells=5)
RP0446_Coc_M_r1 <- CreateSeuratObject(counts=RP0446_Coc_M_r1_Data$`Gene Expression`,project="RP0446_Coc_M",min.cells=5)
RP0883_Suc_F_r1 <- CreateSeuratObject(counts=RP0883_Suc_F_r1_Data$`Gene Expression`,project="RP0883_Suc_F",min.cells=5)
RP0231_Suc_F_r1 <- CreateSeuratObject(counts=RP0231_Suc_F_r1_Data$`Gene Expression`,project="RP0231_Suc_F",min.cells=5)
RP0897_Suc_M_r3 <- CreateSeuratObject(counts=RP0897_Suc_M_r3_Data$`Gene Expression`,project="RP0897_Suc_M",min.cells=5)
RP0446_Coc_F_r1 <- CreateSeuratObject(counts=RP0446_Coc_F_r1_Data$`Gene Expression`,project="RP0446_Coc_F",min.cells=5)
RP0897_Suc_F_r2 <- CreateSeuratObject(counts=RP0897_Suc_F_r2_Data$`Gene Expression`,project="RP0897_Suc_F",min.cells=5)
RP0883_Coc_F_r2 <- CreateSeuratObject(counts=RP0883_Coc_F_r2_Data$`Gene Expression`,project="RP0883_Coc_F",min.cells=5)
RP0457_Suc_F_r1 <- CreateSeuratObject(counts=RP0457_Suc_F_r1_Data$`Gene Expression`,project="RP0457_Suc_F",min.cells=5)
RP0883_Suc_F_r2 <- CreateSeuratObject(counts=RP0883_Suc_F_r2_Data$`Gene Expression`,project="RP0883_Suc_F",min.cells=5)
RP0231_Coc_M_r2 <- CreateSeuratObject(counts=RP0231_Coc_M_r2_Data$`Gene Expression`,project="RP0231_Coc_M",min.cells=5)
RP0609_Suc_M_r1 <- CreateSeuratObject(counts=RP0609_Suc_M_r1_Data$`Gene Expression`,project="RP0609_Suc_M",min.cells=5)
RP0897_Coc_F_r1 <- CreateSeuratObject(counts=RP0897_Coc_F_r1_Data$`Gene Expression`,project="RP0897_Coc_F",min.cells=5)
RP0883_Coc_M_r2 <- CreateSeuratObject(counts=RP0883_Coc_M_r2_Data$`Gene Expression`,project="RP0883_Coc_M",min.cells=5)
RP0457_Suc_M_r1 <- CreateSeuratObject(counts=RP0457_Suc_M_r1_Data$`Gene Expression`,project="RP0457_Suc_M",min.cells=5)
RP0609_Coc_F_r2 <- CreateSeuratObject(counts=RP0609_Coc_F_r2_Data$`Gene Expression`,project="RP0609_Coc_F",min.cells=5)
RP0231_Coc_M_r1 <- CreateSeuratObject(counts=RP0231_Coc_M_r1_Data$`Gene Expression`,project="RP0231_Coc_M",min.cells=5)
RP0609_Suc_F_r1 <- CreateSeuratObject(counts=RP0609_Suc_F_r1_Data$`Gene Expression`,project="RP0609_Suc_F",min.cells=5)
RP0446_Suc_M_r1 <- CreateSeuratObject(counts=RP0446_Suc_M_r1_Data$`Gene Expression`,project="RP0446_Suc_M",min.cells=5)
RP0897_Coc_M_r2 <- CreateSeuratObject(counts=RP0897_Coc_M_r2_Data$`Gene Expression`,project="RP0897_Coc_M",min.cells=5)
RP0609_Suc_M_r2 <- CreateSeuratObject(counts=RP0609_Suc_M_r2_Data$`Gene Expression`,project="RP0609_Suc_M",min.cells=5)
RP0883_Coc_F_r1 <- CreateSeuratObject(counts=RP0883_Coc_F_r1_Data$`Gene Expression`,project="RP0883_Coc_F",min.cells=5)
RP0609_Suc_F_r2 <- CreateSeuratObject(counts=RP0609_Suc_F_r2_Data$`Gene Expression`,project="RP0609_Suc_F",min.cells=5)
RP0883_Coc_M_r1 <- CreateSeuratObject(counts=RP0883_Coc_M_r1_Data$`Gene Expression`,project="RP0883_Coc_M",min.cells=5)
RP0446_Suc_F_r1 <- CreateSeuratObject(counts=RP0446_Suc_F_r1_Data$`Gene Expression`,project="RP0446_Suc_F",min.cells=5)
RP0231_Coc_F_r2 <- CreateSeuratObject(counts=RP0231_Coc_F_r2_Data$`Gene Expression`,project="RP0231_Coc_F",min.cells=5)
RP0446_Coc_F_r2 <- CreateSeuratObject(counts=RP0446_Coc_F_r2_Data$`Gene Expression`,project="RP0446_Coc_F",min.cells=5)
RP0897_Suc_M_r1 <- CreateSeuratObject(counts=RP0897_Suc_M_r1_Data$`Gene Expression`,project="RP0897_Suc_M",min.cells=5)
RP0883_Suc_M_r2 <- CreateSeuratObject(counts=RP0883_Suc_M_r2_Data$`Gene Expression`,project="RP0883_Suc_M",min.cells=5)
RP0609_Coc_M_r1 <- CreateSeuratObject(counts=RP0609_Coc_M_r1_Data$`Gene Expression`,project="RP0609_Coc_M",min.cells=5)

#Set treatment factor as stim in each object:
RP0457_Coc_F_r2$stim <- "Coc"
RP0231_Suc_M_r2$stim <- "Suc"
RP0609_Coc_M_r2$stim <- "Coc"
RP0231_Suc_F_r2$stim <- "Suc"
RP0446_Coc_M_r2$stim <- "Coc"
RP0457_Coc_M_r2$stim <- "Coc"
RP0897_Suc_M_r2$stim <- "Suc"
RP0897_Coc_M_r1$stim <- "Coc"
RP0231_Coc_F_r1$stim <- "Coc"
RP0897_Coc_F_r2$stim <- "Coc"
RP0446_Suc_M_r2$stim <- "Suc"
RP0457_Coc_F_r1$stim <- "Coc"
RP0457_Suc_M_r2$stim <- "Suc"
RP0883_Suc_M_r1$stim <- "Suc"
RP0457_Coc_M_r1$stim <- "Coc"
RP0457_Suc_F_r2$stim <- "Suc"
RP0609_Coc_F_r1$stim <- "Coc"
RP0897_Suc_F_r1$stim <- "Suc"
RP0231_Suc_M_r1$stim <- "Suc"
RP0446_Suc_F_r2$stim <- "Suc"
RP0446_Coc_M_r1$stim <- "Coc"
RP0883_Suc_F_r1$stim <- "Suc"
RP0231_Suc_F_r1$stim <- "Suc"
RP0897_Suc_M_r3$stim <- "Suc"
RP0446_Coc_F_r1$stim <- "Coc"
RP0897_Suc_F_r2$stim <- "Suc"
RP0883_Coc_F_r2$stim <- "Coc"
RP0457_Suc_F_r1$stim <- "Suc"
RP0883_Suc_F_r2$stim <- "Suc"
RP0231_Coc_M_r2$stim <- "Coc"
RP0609_Suc_M_r1$stim <- "Suc"
RP0897_Coc_F_r1$stim <- "Coc"
RP0883_Coc_M_r2$stim <- "Coc"
RP0457_Suc_M_r1$stim <- "Suc"
RP0609_Coc_F_r2$stim <- "Coc"
RP0231_Coc_M_r1$stim <- "Coc"
RP0609_Suc_F_r1$stim <- "Suc"
RP0446_Suc_M_r1$stim <- "Suc"
RP0897_Coc_M_r2$stim <- "Coc"
RP0609_Suc_M_r2$stim <- "Suc"
RP0883_Coc_F_r1$stim <- "Coc"
RP0609_Suc_F_r2$stim <- "Suc"
RP0883_Coc_M_r1$stim <- "Coc"
RP0446_Suc_F_r1$stim <- "Suc"
RP0231_Coc_F_r2$stim <- "Coc"
RP0446_Coc_F_r2$stim <- "Coc"
RP0897_Suc_M_r1$stim <- "Suc"
RP0883_Suc_M_r2$stim <- "Suc"
RP0609_Coc_M_r1$stim <- "Coc"

#Set gender by treatment factor as gender_stim variable.
RP0457_Coc_F_r2$gender_stim <- "Coc_F"
RP0231_Suc_M_r2$gender_stim <- "Suc_M"
RP0609_Coc_M_r2$gender_stim <- "Coc_M"
RP0231_Suc_F_r2$gender_stim <- "Suc_F"
RP0446_Coc_M_r2$gender_stim <- "Coc_M"
RP0457_Coc_M_r2$gender_stim <- "Coc_M"
RP0897_Suc_M_r2$gender_stim <- "Suc_M"
RP0897_Coc_M_r1$gender_stim <- "Coc_M"
RP0231_Coc_F_r1$gender_stim <- "Coc_F"
RP0897_Coc_F_r2$gender_stim <- "Coc_F"
RP0446_Suc_M_r2$gender_stim <- "Suc_M"
RP0457_Coc_F_r1$gender_stim <- "Coc_F"
RP0457_Suc_M_r2$gender_stim <- "Suc_M"
RP0883_Suc_M_r1$gender_stim <- "Suc_M"
RP0457_Coc_M_r1$gender_stim <- "Coc_M"
RP0457_Suc_F_r2$gender_stim <- "Suc_F"
RP0609_Coc_F_r1$gender_stim <- "Coc_F"
RP0897_Suc_F_r1$gender_stim <- "Suc_F"
RP0231_Suc_M_r1$gender_stim <- "Suc_M"
RP0446_Suc_F_r2$gender_stim <- "Suc_F"
RP0446_Coc_M_r1$gender_stim <- "Coc_M"
RP0883_Suc_F_r1$gender_stim <- "Suc_F"
RP0231_Suc_F_r1$gender_stim <- "Suc_F"
RP0897_Suc_M_r3$gender_stim <- "Suc_M"
RP0446_Coc_F_r1$gender_stim <- "Coc_F"
RP0897_Suc_F_r2$gender_stim <- "Suc_F"
RP0883_Coc_F_r2$gender_stim <- "Coc_F"
RP0457_Suc_F_r1$gender_stim <- "Suc_F"
RP0883_Suc_F_r2$gender_stim <- "Suc_F"
RP0231_Coc_M_r2$gender_stim <- "Coc_M"
RP0609_Suc_M_r1$gender_stim <- "Suc_M"
RP0897_Coc_F_r1$gender_stim <- "Coc_F"
RP0883_Coc_M_r2$gender_stim <- "Coc_M"
RP0457_Suc_M_r1$gender_stim <- "Suc_M"
RP0609_Coc_F_r2$gender_stim <- "Coc_F"
RP0231_Coc_M_r1$gender_stim <- "Coc_M"
RP0609_Suc_F_r1$gender_stim <- "Suc_F"
RP0446_Suc_M_r1$gender_stim <- "Suc_M"
RP0897_Coc_M_r2$gender_stim <- "Coc_M"
RP0609_Suc_M_r2$gender_stim <- "Suc_M"
RP0883_Coc_F_r1$gender_stim <- "Coc_F"
RP0609_Suc_F_r2$gender_stim <- "Suc_F"
RP0883_Coc_M_r1$gender_stim <- "Coc_M"
RP0446_Suc_F_r1$gender_stim <- "Suc_F"
RP0231_Coc_F_r2$gender_stim <- "Coc_F"
RP0446_Coc_F_r2$gender_stim <- "Coc_F"
RP0897_Suc_M_r1$gender_stim <- "Suc_M"
RP0883_Suc_M_r2$gender_stim <- "Suc_M"
RP0609_Coc_M_r1$gender_stim <- "Coc_M"

#Setting sample ID as a factor.
RP0457_Coc_F_r2$sample_id <- "RP0457_Coc_F_r2"
RP0231_Suc_M_r2$sample_id <- "RP0231_Suc_M_r2"
RP0609_Coc_M_r2$sample_id <- "RP0609_Coc_M_r2"
RP0231_Suc_F_r2$sample_id <- "RP0231_Suc_F_r2"
RP0446_Coc_M_r2$sample_id <- "RP0446_Coc_M_r2"
RP0457_Coc_M_r2$sample_id <- "RP0457_Coc_M_r2"
RP0897_Suc_M_r2$sample_id <- "RP0897_Suc_M_r2"
RP0897_Coc_M_r1$sample_id <- "RP0897_Coc_M_r1"
RP0231_Coc_F_r1$sample_id <- "RP0231_Coc_F_r1"
RP0897_Coc_F_r2$sample_id <- "RP0897_Coc_F_r2"
RP0446_Suc_M_r2$sample_id <- "RP0446_Suc_M_r2"
RP0457_Coc_F_r1$sample_id <- "RP0457_Coc_F_r1"
RP0457_Suc_M_r2$sample_id <- "RP0457_Suc_M_r2"
RP0883_Suc_M_r1$sample_id <- "RP0883_Suc_M_r1"
RP0457_Coc_M_r1$sample_id <- "RP0457_Coc_M_r1"
RP0457_Suc_F_r2$sample_id <- "RP0457_Suc_F_r2"
RP0609_Coc_F_r1$sample_id <- "RP0609_Coc_F_r1"
RP0897_Suc_F_r1$sample_id <- "RP0897_Suc_F_r1"
RP0231_Suc_M_r1$sample_id <- "RP0231_Suc_M_r1"
RP0446_Suc_F_r2$sample_id <- "RP0446_Suc_F_r2"
RP0446_Coc_M_r1$sample_id <- "RP0446_Coc_M_r1"
RP0883_Suc_F_r1$sample_id <- "RP0883_Suc_F_r1"
RP0231_Suc_F_r1$sample_id <- "RP0231_Suc_F_r1"
RP0897_Suc_M_r3$sample_id <- "RP0897_Suc_M_r3"
RP0446_Coc_F_r1$sample_id <- "RP0446_Coc_F_r1"
RP0897_Suc_F_r2$sample_id <- "RP0897_Suc_F_r2"
RP0883_Coc_F_r2$sample_id <- "RP0883_Coc_F_r2"
RP0457_Suc_F_r1$sample_id <- "RP0457_Suc_F_r1"
RP0883_Suc_F_r2$sample_id <- "RP0883_Suc_F_r2"
RP0231_Coc_M_r2$sample_id <- "RP0231_Coc_M_r2"
RP0609_Suc_M_r1$sample_id <- "RP0609_Suc_M_r1"
RP0897_Coc_F_r1$sample_id <- "RP0897_Coc_F_r1"
RP0883_Coc_M_r2$sample_id <- "RP0883_Coc_M_r2"
RP0457_Suc_M_r1$sample_id <- "RP0457_Suc_M_r1"
RP0609_Coc_F_r2$sample_id <- "RP0609_Coc_F_r2"
RP0231_Coc_M_r1$sample_id <- "RP0231_Coc_M_r1"
RP0609_Suc_F_r1$sample_id <- "RP0609_Suc_F_r1"
RP0446_Suc_M_r1$sample_id <- "RP0446_Suc_M_r1"
RP0897_Coc_M_r2$sample_id <- "RP0897_Coc_M_r2"
RP0609_Suc_M_r2$sample_id <- "RP0609_Suc_M_r2"
RP0883_Coc_F_r1$sample_id <- "RP0883_Coc_F_r1"
RP0609_Suc_F_r2$sample_id <- "RP0609_Suc_F_r2"
RP0883_Coc_M_r1$sample_id <- "RP0883_Coc_M_r1"
RP0446_Suc_F_r1$sample_id <- "RP0446_Suc_F_r1"
RP0231_Coc_F_r2$sample_id <- "RP0231_Coc_F_r2"
RP0446_Coc_F_r2$sample_id <- "RP0446_Coc_F_r2"
RP0897_Suc_M_r1$sample_id <- "RP0897_Suc_M_r1"
RP0883_Suc_M_r2$sample_id <- "RP0883_Suc_M_r2"
RP0609_Coc_M_r1$sample_id <- "RP0609_Coc_M_r1"

#Set genotype as a factor.
RP0457_Coc_F_r2$genotype <- "RP0457"
RP0231_Suc_M_r2$genotype <- "RP0231"
RP0609_Coc_M_r2$genotype <- "RP0609"
RP0231_Suc_F_r2$genotype <- "RP0231"
RP0446_Coc_M_r2$genotype <- "RP0446"
RP0457_Coc_M_r2$genotype <- "RP0457"
RP0897_Suc_M_r2$genotype <- "RP0897"
RP0897_Coc_M_r1$genotype <- "RP0897"
RP0231_Coc_F_r1$genotype <- "RP0231"
RP0897_Coc_F_r2$genotype <- "RP0897"
RP0446_Suc_M_r2$genotype <- "RP0446"
RP0457_Coc_F_r1$genotype <- "RP0457"
RP0457_Suc_M_r2$genotype <- "RP0457"
RP0883_Suc_M_r1$genotype <- "RP0883"
RP0457_Coc_M_r1$genotype <- "RP0457"
RP0457_Suc_F_r2$genotype <- "RP0457"
RP0609_Coc_F_r1$genotype <- "RP0609"
RP0897_Suc_F_r1$genotype <- "RP0897"
RP0231_Suc_M_r1$genotype <- "RP0231"
RP0446_Suc_F_r2$genotype <- "RP0446"
RP0446_Coc_M_r1$genotype <- "RP0446"
RP0883_Suc_F_r1$genotype <- "RP0883"
RP0231_Suc_F_r1$genotype <- "RP0231"
RP0897_Suc_M_r3$genotype <- "RP0897"
RP0446_Coc_F_r1$genotype <- "RP0446"
RP0897_Suc_F_r2$genotype <- "RP0897"
RP0883_Coc_F_r2$genotype <- "RP0883"
RP0457_Suc_F_r1$genotype <- "RP0457"
RP0883_Suc_F_r2$genotype <- "RP0883"
RP0231_Coc_M_r2$genotype <- "RP0231"
RP0609_Suc_M_r1$genotype <- "RP0609"
RP0897_Coc_F_r1$genotype <- "RP0897"
RP0883_Coc_M_r2$genotype <- "RP0883"
RP0457_Suc_M_r1$genotype <- "RP0457"
RP0609_Coc_F_r2$genotype <- "RP0609"
RP0231_Coc_M_r1$genotype <- "RP0231"
RP0609_Suc_F_r1$genotype <- "RP0609"
RP0446_Suc_M_r1$genotype <- "RP0446"
RP0897_Coc_M_r2$genotype <- "RP0897"
RP0609_Suc_M_r2$genotype <- "RP0609"
RP0883_Coc_F_r1$genotype <- "RP0883"
RP0609_Suc_F_r2$genotype <- "RP0609"
RP0883_Coc_M_r1$genotype <- "RP0883"
RP0446_Suc_F_r1$genotype <- "RP0446"
RP0231_Coc_F_r2$genotype <- "RP0231"
RP0446_Coc_F_r2$genotype <- "RP0446"
RP0897_Suc_M_r1$genotype <- "RP0897"
RP0883_Suc_M_r2$genotype <- "RP0883"
RP0609_Coc_M_r1$genotype <- "RP0609"

#Subset features based on gold standards and representative sample values (RP0457_Coc_F_r2)
RP0457_Coc_F_r2 <- subset(RP0457_Coc_F_r2,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0231_Suc_M_r2 <- subset(RP0231_Suc_M_r2,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0609_Coc_M_r2 <- subset(RP0609_Coc_M_r2,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0231_Suc_F_r2 <- subset(RP0231_Suc_F_r2,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0446_Coc_M_r2 <- subset(RP0446_Coc_M_r2,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0457_Coc_M_r2 <- subset(RP0457_Coc_M_r2,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0897_Suc_M_r2 <- subset(RP0897_Suc_M_r2,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0897_Coc_M_r1 <- subset(RP0897_Coc_M_r1,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0231_Coc_F_r1 <- subset(RP0231_Coc_F_r1,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0897_Coc_F_r2 <- subset(RP0897_Coc_F_r2,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0446_Suc_M_r2 <- subset(RP0446_Suc_M_r2,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0457_Coc_F_r1 <- subset(RP0457_Coc_F_r1,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0457_Suc_M_r2 <- subset(RP0457_Suc_M_r2,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0883_Suc_M_r1 <- subset(RP0883_Suc_M_r1,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0457_Coc_M_r1 <- subset(RP0457_Coc_M_r1,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0457_Suc_F_r2 <- subset(RP0457_Suc_F_r2,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0609_Coc_F_r1 <- subset(RP0609_Coc_F_r1,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0897_Suc_F_r1 <- subset(RP0897_Suc_F_r1,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0231_Suc_M_r1 <- subset(RP0231_Suc_M_r1,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0446_Suc_F_r2 <- subset(RP0446_Suc_F_r2,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0446_Coc_M_r1 <- subset(RP0446_Coc_M_r1,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0883_Suc_F_r1 <- subset(RP0883_Suc_F_r1,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0231_Suc_F_r1 <- subset(RP0231_Suc_F_r1,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0897_Suc_M_r3 <- subset(RP0897_Suc_M_r3,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0446_Coc_F_r1 <- subset(RP0446_Coc_F_r1,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0897_Suc_F_r2 <- subset(RP0897_Suc_F_r2,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0883_Coc_F_r2 <- subset(RP0883_Coc_F_r2,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0457_Suc_F_r1 <- subset(RP0457_Suc_F_r1,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0883_Suc_F_r2 <- subset(RP0883_Suc_F_r2,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0231_Coc_M_r2 <- subset(RP0231_Coc_M_r2,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0609_Suc_M_r1 <- subset(RP0609_Suc_M_r1,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0897_Coc_F_r1 <- subset(RP0897_Coc_F_r1,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0883_Coc_M_r2 <- subset(RP0883_Coc_M_r2,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0457_Suc_M_r1 <- subset(RP0457_Suc_M_r1,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0609_Coc_F_r2 <- subset(RP0609_Coc_F_r2,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0231_Coc_M_r1 <- subset(RP0231_Coc_M_r1,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0609_Suc_F_r1 <- subset(RP0609_Suc_F_r1,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0446_Suc_M_r1 <- subset(RP0446_Suc_M_r1,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0897_Coc_M_r2 <- subset(RP0897_Coc_M_r2,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0609_Suc_M_r2 <- subset(RP0609_Suc_M_r2,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0883_Coc_F_r1 <- subset(RP0883_Coc_F_r1,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0609_Suc_F_r2 <- subset(RP0609_Suc_F_r2,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0883_Coc_M_r1 <- subset(RP0883_Coc_M_r1,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0446_Suc_F_r1 <- subset(RP0446_Suc_F_r1,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0231_Coc_F_r2 <- subset(RP0231_Coc_F_r2,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0446_Coc_F_r2 <- subset(RP0446_Coc_F_r2,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0897_Suc_M_r1 <- subset(RP0897_Suc_M_r1,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0883_Suc_M_r2 <- subset(RP0883_Suc_M_r2,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)
RP0609_Coc_M_r1 <- subset(RP0609_Coc_M_r1,subset=nFeature_RNA > 300 & nFeature_RNA < 2000)

#SCTransform
RP0457_Coc_F_r2 <- SCTransform(RP0457_Coc_F_r2,verbose = FALSE,return.only.var.genes = FALSE)
RP0231_Suc_M_r2 <- SCTransform(RP0231_Suc_M_r2,verbose = FALSE,return.only.var.genes = FALSE)
RP0609_Coc_M_r2 <- SCTransform(RP0609_Coc_M_r2,verbose = FALSE,return.only.var.genes = FALSE)
RP0231_Suc_F_r2 <- SCTransform(RP0231_Suc_F_r2,verbose = FALSE,return.only.var.genes = FALSE)
RP0446_Coc_M_r2 <- SCTransform(RP0446_Coc_M_r2,verbose = FALSE,return.only.var.genes = FALSE)
RP0457_Coc_M_r2 <- SCTransform(RP0457_Coc_M_r2,verbose = FALSE,return.only.var.genes = FALSE)
RP0897_Suc_M_r2 <- SCTransform(RP0897_Suc_M_r2,verbose = FALSE,return.only.var.genes = FALSE)
RP0897_Coc_M_r1 <- SCTransform(RP0897_Coc_M_r1,verbose = FALSE,return.only.var.genes = FALSE)
RP0231_Coc_F_r1 <- SCTransform(RP0231_Coc_F_r1,verbose = FALSE,return.only.var.genes = FALSE)
RP0897_Coc_F_r2 <- SCTransform(RP0897_Coc_F_r2,verbose = FALSE,return.only.var.genes = FALSE)
RP0446_Suc_M_r2 <- SCTransform(RP0446_Suc_M_r2,verbose = FALSE,return.only.var.genes = FALSE)
RP0457_Coc_F_r1 <- SCTransform(RP0457_Coc_F_r1,verbose = FALSE,return.only.var.genes = FALSE)
RP0457_Suc_M_r2 <- SCTransform(RP0457_Suc_M_r2,verbose = FALSE,return.only.var.genes = FALSE)
RP0883_Suc_M_r1 <- SCTransform(RP0883_Suc_M_r1,verbose = FALSE,return.only.var.genes = FALSE)
RP0457_Coc_M_r1 <- SCTransform(RP0457_Coc_M_r1,verbose = FALSE,return.only.var.genes = FALSE)
RP0457_Suc_F_r2 <- SCTransform(RP0457_Suc_F_r2,verbose = FALSE,return.only.var.genes = FALSE)
RP0609_Coc_F_r1 <- SCTransform(RP0609_Coc_F_r1,verbose = FALSE,return.only.var.genes = FALSE)
RP0897_Suc_F_r1 <- SCTransform(RP0897_Suc_F_r1,verbose = FALSE,return.only.var.genes = FALSE)
RP0231_Suc_M_r1 <- SCTransform(RP0231_Suc_M_r1,verbose = FALSE,return.only.var.genes = FALSE)
RP0446_Suc_F_r2 <- SCTransform(RP0446_Suc_F_r2,verbose = FALSE,return.only.var.genes = FALSE)
RP0446_Coc_M_r1 <- SCTransform(RP0446_Coc_M_r1,verbose = FALSE,return.only.var.genes = FALSE)
RP0883_Suc_F_r1 <- SCTransform(RP0883_Suc_F_r1,verbose = FALSE,return.only.var.genes = FALSE)
RP0231_Suc_F_r1 <- SCTransform(RP0231_Suc_F_r1,verbose = FALSE,return.only.var.genes = FALSE)
RP0897_Suc_M_r3 <- SCTransform(RP0897_Suc_M_r3,verbose = FALSE,return.only.var.genes = FALSE)
RP0446_Coc_F_r1 <- SCTransform(RP0446_Coc_F_r1,verbose = FALSE,return.only.var.genes = FALSE)
RP0897_Suc_F_r2 <- SCTransform(RP0897_Suc_F_r2,verbose = FALSE,return.only.var.genes = FALSE)
RP0883_Coc_F_r2 <- SCTransform(RP0883_Coc_F_r2,verbose = FALSE,return.only.var.genes = FALSE)
RP0457_Suc_F_r1 <- SCTransform(RP0457_Suc_F_r1,verbose = FALSE,return.only.var.genes = FALSE)
RP0883_Suc_F_r2 <- SCTransform(RP0883_Suc_F_r2,verbose = FALSE,return.only.var.genes = FALSE)
RP0231_Coc_M_r2 <- SCTransform(RP0231_Coc_M_r2,verbose = FALSE,return.only.var.genes = FALSE)
RP0609_Suc_M_r1 <- SCTransform(RP0609_Suc_M_r1,verbose = FALSE,return.only.var.genes = FALSE)
RP0897_Coc_F_r1 <- SCTransform(RP0897_Coc_F_r1,verbose = FALSE,return.only.var.genes = FALSE)
RP0883_Coc_M_r2 <- SCTransform(RP0883_Coc_M_r2,verbose = FALSE,return.only.var.genes = FALSE)
RP0457_Suc_M_r1 <- SCTransform(RP0457_Suc_M_r1,verbose = FALSE,return.only.var.genes = FALSE)
RP0609_Coc_F_r2 <- SCTransform(RP0609_Coc_F_r2,verbose = FALSE,return.only.var.genes = FALSE)
RP0231_Coc_M_r1 <- SCTransform(RP0231_Coc_M_r1,verbose = FALSE,return.only.var.genes = FALSE)
RP0609_Suc_F_r1 <- SCTransform(RP0609_Suc_F_r1,verbose = FALSE,return.only.var.genes = FALSE)
RP0446_Suc_M_r1 <- SCTransform(RP0446_Suc_M_r1,verbose = FALSE,return.only.var.genes = FALSE)
RP0897_Coc_M_r2 <- SCTransform(RP0897_Coc_M_r2,verbose = FALSE,return.only.var.genes = FALSE)
RP0609_Suc_M_r2 <- SCTransform(RP0609_Suc_M_r2,verbose = FALSE,return.only.var.genes = FALSE)
RP0883_Coc_F_r1 <- SCTransform(RP0883_Coc_F_r1,verbose = FALSE,return.only.var.genes = FALSE)
RP0609_Suc_F_r2 <- SCTransform(RP0609_Suc_F_r2,verbose = FALSE,return.only.var.genes = FALSE)
RP0883_Coc_M_r1 <- SCTransform(RP0883_Coc_M_r1,verbose = FALSE,return.only.var.genes = FALSE)
RP0446_Suc_F_r1 <- SCTransform(RP0446_Suc_F_r1,verbose = FALSE,return.only.var.genes = FALSE)
RP0231_Coc_F_r2 <- SCTransform(RP0231_Coc_F_r2,verbose = FALSE,return.only.var.genes = FALSE)
RP0446_Coc_F_r2 <- SCTransform(RP0446_Coc_F_r2,verbose = FALSE,return.only.var.genes = FALSE)
RP0897_Suc_M_r1 <- SCTransform(RP0897_Suc_M_r1,verbose = FALSE,return.only.var.genes = FALSE)
RP0883_Suc_M_r2 <- SCTransform(RP0883_Suc_M_r2,verbose = FALSE,return.only.var.genes = FALSE)
RP0609_Coc_M_r1 <- SCTransform(RP0609_Coc_M_r1,verbose = FALSE,return.only.var.genes = FALSE)

#Start integration procedure
Integration_feature_set <-SelectIntegrationFeatures(object.list = c(RP0457_Coc_F_r2,
RP0231_Suc_M_r2,
RP0609_Coc_M_r2,
RP0231_Suc_F_r2,
RP0446_Coc_M_r2,
RP0457_Coc_M_r2,
RP0897_Suc_M_r2,
RP0897_Coc_M_r1,
RP0231_Coc_F_r1,
RP0897_Coc_F_r2,
RP0446_Suc_M_r2,
RP0457_Coc_F_r1,
RP0457_Suc_M_r2,
RP0883_Suc_M_r1,
RP0457_Coc_M_r1,
RP0457_Suc_F_r2,
RP0609_Coc_F_r1,
RP0897_Suc_F_r1,
RP0231_Suc_M_r1,
RP0446_Suc_F_r2,
RP0446_Coc_M_r1,
RP0883_Suc_F_r1,
RP0231_Suc_F_r1,
RP0897_Suc_M_r3,
RP0446_Coc_F_r1,
RP0897_Suc_F_r2,
RP0883_Coc_F_r2,
RP0457_Suc_F_r1,
RP0883_Suc_F_r2,
RP0231_Coc_M_r2,
RP0609_Suc_M_r1,
RP0897_Coc_F_r1,
RP0883_Coc_M_r2,
RP0457_Suc_M_r1,
RP0609_Coc_F_r2,
RP0231_Coc_M_r1,
RP0609_Suc_F_r1,
RP0446_Suc_M_r1,
RP0897_Coc_M_r2,
RP0609_Suc_M_r2,
RP0883_Coc_F_r1,
RP0609_Suc_F_r2,
RP0883_Coc_M_r1,
RP0446_Suc_F_r1,
RP0231_Coc_F_r2,
RP0446_Coc_F_r2,
RP0897_Suc_M_r1,
RP0883_Suc_M_r2,
RP0609_Coc_M_r1), nfeatures=1500)

Integration_list <- PrepSCTIntegration(object.list = c(RP0457_Coc_F_r2,
RP0231_Suc_M_r2,
RP0609_Coc_M_r2,
RP0231_Suc_F_r2,
RP0446_Coc_M_r2,
RP0457_Coc_M_r2,
RP0897_Suc_M_r2,
RP0897_Coc_M_r1,
RP0231_Coc_F_r1,
RP0897_Coc_F_r2,
RP0446_Suc_M_r2,
RP0457_Coc_F_r1,
RP0457_Suc_M_r2,
RP0883_Suc_M_r1,
RP0457_Coc_M_r1,
RP0457_Suc_F_r2,
RP0609_Coc_F_r1,
RP0897_Suc_F_r1,
RP0231_Suc_M_r1,
RP0446_Suc_F_r2,
RP0446_Coc_M_r1,
RP0883_Suc_F_r1,
RP0231_Suc_F_r1,
RP0897_Suc_M_r3,
RP0446_Coc_F_r1,
RP0897_Suc_F_r2,
RP0883_Coc_F_r2,
RP0457_Suc_F_r1,
RP0883_Suc_F_r2,
RP0231_Coc_M_r2,
RP0609_Suc_M_r1,
RP0897_Coc_F_r1,
RP0883_Coc_M_r2,
RP0457_Suc_M_r1,
RP0609_Coc_F_r2,
RP0231_Coc_M_r1,
RP0609_Suc_F_r1,
RP0446_Suc_M_r1,
RP0897_Coc_M_r2,
RP0609_Suc_M_r2,
RP0883_Coc_F_r1,
RP0609_Suc_F_r2,
RP0883_Coc_M_r1,
RP0446_Suc_F_r1,
RP0231_Coc_F_r2,
RP0446_Coc_F_r2,
RP0897_Suc_M_r1,
RP0883_Suc_M_r2,
RP0609_Coc_M_r1), anchor.features = Integration_feature_set, verbose=FALSE)

#Run PCA for Reciprocal PCA approach in FIA step:
RP0457_Coc_F_r2 <- RunPCA(RP0457_Coc_F_r2,features = Integration_feature_set)
RP0231_Suc_M_r2 <- RunPCA(RP0231_Suc_M_r2,features = Integration_feature_set)
RP0609_Coc_M_r2 <- RunPCA(RP0609_Coc_M_r2,features = Integration_feature_set)
RP0231_Suc_F_r2 <- RunPCA(RP0231_Suc_F_r2,features = Integration_feature_set)
RP0446_Coc_M_r2 <- RunPCA(RP0446_Coc_M_r2,features = Integration_feature_set)
RP0457_Coc_M_r2 <- RunPCA(RP0457_Coc_M_r2,features = Integration_feature_set)
RP0897_Suc_M_r2 <- RunPCA(RP0897_Suc_M_r2,features = Integration_feature_set)
RP0897_Coc_M_r1 <- RunPCA(RP0897_Coc_M_r1,features = Integration_feature_set)
RP0231_Coc_F_r1 <- RunPCA(RP0231_Coc_F_r1,features = Integration_feature_set)
RP0897_Coc_F_r2 <- RunPCA(RP0897_Coc_F_r2,features = Integration_feature_set)
RP0446_Suc_M_r2 <- RunPCA(RP0446_Suc_M_r2,features = Integration_feature_set)
RP0457_Coc_F_r1 <- RunPCA(RP0457_Coc_F_r1,features = Integration_feature_set)
RP0457_Suc_M_r2 <- RunPCA(RP0457_Suc_M_r2,features = Integration_feature_set)
RP0883_Suc_M_r1 <- RunPCA(RP0883_Suc_M_r1,features = Integration_feature_set)
RP0457_Coc_M_r1 <- RunPCA(RP0457_Coc_M_r1,features = Integration_feature_set)
RP0457_Suc_F_r2 <- RunPCA(RP0457_Suc_F_r2,features = Integration_feature_set)
RP0609_Coc_F_r1 <- RunPCA(RP0609_Coc_F_r1,features = Integration_feature_set)
RP0897_Suc_F_r1 <- RunPCA(RP0897_Suc_F_r1,features = Integration_feature_set)
RP0231_Suc_M_r1 <- RunPCA(RP0231_Suc_M_r1,features = Integration_feature_set)
RP0446_Suc_F_r2 <- RunPCA(RP0446_Suc_F_r2,features = Integration_feature_set)
RP0446_Coc_M_r1 <- RunPCA(RP0446_Coc_M_r1,features = Integration_feature_set)
RP0883_Suc_F_r1 <- RunPCA(RP0883_Suc_F_r1,features = Integration_feature_set)
RP0231_Suc_F_r1 <- RunPCA(RP0231_Suc_F_r1,features = Integration_feature_set)
RP0897_Suc_M_r3 <- RunPCA(RP0897_Suc_M_r3,features = Integration_feature_set)
RP0446_Coc_F_r1 <- RunPCA(RP0446_Coc_F_r1,features = Integration_feature_set)
RP0897_Suc_F_r2 <- RunPCA(RP0897_Suc_F_r2,features = Integration_feature_set)
RP0883_Coc_F_r2 <- RunPCA(RP0883_Coc_F_r2,features = Integration_feature_set)
RP0457_Suc_F_r1 <- RunPCA(RP0457_Suc_F_r1,features = Integration_feature_set)
RP0883_Suc_F_r2 <- RunPCA(RP0883_Suc_F_r2,features = Integration_feature_set)
RP0231_Coc_M_r2 <- RunPCA(RP0231_Coc_M_r2,features = Integration_feature_set)
RP0609_Suc_M_r1 <- RunPCA(RP0609_Suc_M_r1,features = Integration_feature_set)
RP0897_Coc_F_r1 <- RunPCA(RP0897_Coc_F_r1,features = Integration_feature_set)
RP0883_Coc_M_r2 <- RunPCA(RP0883_Coc_M_r2,features = Integration_feature_set)
RP0457_Suc_M_r1 <- RunPCA(RP0457_Suc_M_r1,features = Integration_feature_set)
RP0609_Coc_F_r2 <- RunPCA(RP0609_Coc_F_r2,features = Integration_feature_set)
RP0231_Coc_M_r1 <- RunPCA(RP0231_Coc_M_r1,features = Integration_feature_set)
RP0609_Suc_F_r1 <- RunPCA(RP0609_Suc_F_r1,features = Integration_feature_set)
RP0446_Suc_M_r1 <- RunPCA(RP0446_Suc_M_r1,features = Integration_feature_set)
RP0897_Coc_M_r2 <- RunPCA(RP0897_Coc_M_r2,features = Integration_feature_set)
RP0609_Suc_M_r2 <- RunPCA(RP0609_Suc_M_r2,features = Integration_feature_set)
RP0883_Coc_F_r1 <- RunPCA(RP0883_Coc_F_r1,features = Integration_feature_set)
RP0609_Suc_F_r2 <- RunPCA(RP0609_Suc_F_r2,features = Integration_feature_set)
RP0883_Coc_M_r1 <- RunPCA(RP0883_Coc_M_r1,features = Integration_feature_set)
RP0446_Suc_F_r1 <- RunPCA(RP0446_Suc_F_r1,features = Integration_feature_set)
RP0231_Coc_F_r2 <- RunPCA(RP0231_Coc_F_r2,features = Integration_feature_set)
RP0446_Coc_F_r2 <- RunPCA(RP0446_Coc_F_r2,features = Integration_feature_set)
RP0897_Suc_M_r1 <- RunPCA(RP0897_Suc_M_r1,features = Integration_feature_set)
RP0883_Suc_M_r2 <- RunPCA(RP0883_Suc_M_r2,features = Integration_feature_set)
RP0609_Coc_M_r1 <- RunPCA(RP0609_Coc_M_r1,features = Integration_feature_set)

#Parallelism
plan("multicore",workers=6)

##Set the mode of integration to SCT (check 2017 paper for this procedure)
Integration_anchors <- FindIntegrationAnchors(object.list = Integration_list, normalization.method = "SCT", anchor.features = Integration_feature_set, verbose = FALSE)

##Normalize using SCT
Integrated <- IntegrateData(anchorset = Integration_anchors, normalization.method = "SCT", verbose = FALSE)

##Run dimensionality reduction procedures
Integrated <- RunPCA(object = Integrated, verbose = FALSE)

Integrated <- FindNeighbors(Integrated, reduction = "pca", dims = 1:10)

Integrated <- RunUMAP(Integrated, reduction = "pca", dims= 1:10)

##Iterate with different resolutions to parametrize granularity and identify cluster number plateau. Range of 0.4 to 2.0 is from Seurat documentation
for (x in seq(from=0.4, to=2.0, by=0.1)) {
cluster <- FindClusters(Integrated, resolution = x)
}


###Version2 Using Standard procedure without SCTransform. Pickup after subset stage.

ifnb.list <- c(RP0457_Coc_F_r2,
RP0231_Suc_M_r2,
RP0609_Coc_M_r2,
RP0231_Suc_F_r2,
RP0446_Coc_M_r2,
RP0457_Coc_M_r2,
RP0897_Suc_M_r2,
RP0897_Coc_M_r1,
RP0231_Coc_F_r1,
RP0897_Coc_F_r2,
RP0446_Suc_M_r2,
RP0457_Coc_F_r1,
RP0457_Suc_M_r2,
RP0883_Suc_M_r1,
RP0457_Coc_M_r1,
RP0457_Suc_F_r2,
RP0609_Coc_F_r1,
RP0897_Suc_F_r1,
RP0231_Suc_M_r1,
RP0446_Suc_F_r2,
RP0446_Coc_M_r1,
RP0883_Suc_F_r1,
RP0231_Suc_F_r1,
RP0897_Suc_M_r3,
RP0446_Coc_F_r1,
RP0897_Suc_F_r2,
RP0883_Coc_F_r2,
RP0457_Suc_F_r1,
RP0883_Suc_F_r2,
RP0231_Coc_M_r2,
RP0609_Suc_M_r1,
RP0897_Coc_F_r1,
RP0883_Coc_M_r2,
RP0457_Suc_M_r1,
RP0609_Coc_F_r2,
RP0231_Coc_M_r1,
RP0609_Suc_F_r1,
RP0446_Suc_M_r1,
RP0897_Coc_M_r2,
RP0609_Suc_M_r2,
RP0883_Coc_F_r1,
RP0609_Suc_F_r2,
RP0883_Coc_M_r1,
RP0446_Suc_F_r1,
RP0231_Coc_F_r2,
RP0446_Coc_F_r2,
RP0897_Suc_M_r1,
RP0883_Suc_M_r2,
RP0609_Coc_M_r1)

# normalize and identify variable features for each dataset independently

ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1500)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = ifnb.list)
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features, reduction = "rpca")

combined <- IntegrateData(anchorset = anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
combined <- ScaleData(immune.combined, verbose = FALSE)
combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
combined <- FindClusters(immune.combined, resolution = 0.5)

##This failed with vector size errors:
##Trying with BPCells and Seurat v5:

library(Seurat)
library(BPCells)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
options(future.globals.maxSize= 1000000000000)
setwd("/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Analysis/Gene_expression/")


#Read data in from filtered feature counts BC matrix. Because these are multiome data, there will two separate matrices per library and Seurat will warn you as such.
RP0457_Coc_F_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_1_8/RP0457_Coc_F_r2/outs/filtered_feature_bc_matrix/')
RP0231_Suc_M_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_1_8/RP0231_Suc_M_r2/outs/filtered_feature_bc_matrix/')
RP0609_Coc_M_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_1_8/RP0609_Coc_M_r2/outs/filtered_feature_bc_matrix/')
RP0231_Suc_F_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_1_8/RP0231_Suc_F_r2/outs/filtered_feature_bc_matrix/')
RP0446_Coc_M_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_1_8/RP0446_Coc_M_r2/outs/filtered_feature_bc_matrix/')
RP0457_Coc_M_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_1_8/RP0457_Coc_M_r2/outs/filtered_feature_bc_matrix/')
RP0897_Suc_M_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_1_8/RP0897_Suc_M_r2/outs/filtered_feature_bc_matrix/')
RP0897_Coc_M_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_1_8/RP0897_Coc_M_r1/outs/filtered_feature_bc_matrix/')
RP0231_Coc_F_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_1_8/RP0231_Coc_F_r1/outs/filtered_feature_bc_matrix/')
RP0897_Coc_F_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_09_16/RP0897_Coc_F_r2/outs/filtered_feature_bc_matrix/')
RP0446_Suc_M_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_09_16/RP0446_Suc_M_r2/outs/filtered_feature_bc_matrix/')
RP0457_Coc_F_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_09_16/RP0457_Coc_F_r1/outs/filtered_feature_bc_matrix/')
RP0457_Suc_M_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_09_16/RP0457_Suc_M_r2/outs/filtered_feature_bc_matrix/')
RP0883_Suc_M_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_09_16/RP0883_Suc_M_r1/outs/filtered_feature_bc_matrix/')
RP0457_Coc_M_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_09_16/RP0457_Coc_M_r1/outs/filtered_feature_bc_matrix/')
RP0457_Suc_F_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_09_16/RP0457_Suc_F_r2/outs/filtered_feature_bc_matrix/')
RP0609_Coc_F_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_09_16/RP0609_Coc_F_r1/outs/filtered_feature_bc_matrix/')
RP0897_Suc_F_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_17_24/RP0897_Suc_F_r1/outs/filtered_feature_bc_matrix/')
RP0231_Suc_M_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_17_24/RP0231_Suc_M_r1/outs/filtered_feature_bc_matrix/')
RP0446_Suc_F_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_17_24/RP0446_Suc_F_r2/outs/filtered_feature_bc_matrix/')
RP0446_Coc_M_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_17_24/RP0446_Coc_M_r1/outs/filtered_feature_bc_matrix/')
RP0883_Suc_F_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_17_24/RP0883_Suc_F_r1/outs/filtered_feature_bc_matrix/')
RP0231_Suc_F_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_17_24/RP0231_Suc_F_r1/outs/filtered_feature_bc_matrix/')
RP0897_Suc_M_r3_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_17_24/RP0897_Suc_M_r3/outs/filtered_feature_bc_matrix/')
RP0446_Coc_F_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_17_24/RP0446_Coc_F_r1/outs/filtered_feature_bc_matrix/')
RP0897_Suc_F_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_25_32/RP0897_Suc_F_r2/outs/filtered_feature_bc_matrix/')
RP0883_Coc_F_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_25_32/RP0883_Coc_F_r2/outs/filtered_feature_bc_matrix/')
RP0457_Suc_F_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_25_32/RP0457_Suc_F_r1/outs/filtered_feature_bc_matrix/')
RP0883_Suc_F_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_25_32/RP0883_Suc_F_r2/outs/filtered_feature_bc_matrix/')
RP0231_Coc_M_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_25_32/RP0231_Coc_M_r2/outs/filtered_feature_bc_matrix/')
RP0609_Suc_M_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_25_32/RP0609_Suc_M_r1/outs/filtered_feature_bc_matrix/')
RP0897_Coc_F_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_25_32/RP0897_Coc_F_r1/outs/filtered_feature_bc_matrix/')
RP0883_Coc_M_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_25_32/RP0883_Coc_M_r2/outs/filtered_feature_bc_matrix/')
RP0457_Suc_M_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_33_40/RP0457_Suc_M_r1/outs/filtered_feature_bc_matrix/')
RP0609_Coc_F_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_33_40/RP0609_Coc_F_r2/outs/filtered_feature_bc_matrix/')
RP0231_Coc_M_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_33_40/RP0231_Coc_M_r1/outs/filtered_feature_bc_matrix/')
RP0609_Suc_F_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_33_40/RP0609_Suc_F_r1/outs/filtered_feature_bc_matrix/')
RP0446_Suc_M_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_33_40/RP0446_Suc_M_r1/outs/filtered_feature_bc_matrix/')
RP0897_Coc_M_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_33_40/RP0897_Coc_M_r2/outs/filtered_feature_bc_matrix/')
RP0609_Suc_M_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_33_40/RP0609_Suc_M_r2/outs/filtered_feature_bc_matrix/')
RP0883_Coc_F_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_33_40/RP0883_Coc_F_r1/outs/filtered_feature_bc_matrix/')
RP0609_Suc_F_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_41_48/RP0609_Suc_F_r2/outs/filtered_feature_bc_matrix/')
RP0883_Coc_M_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_41_48/RP0883_Coc_M_r1/outs/filtered_feature_bc_matrix/')
RP0446_Suc_F_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_41_48/RP0446_Suc_F_r1/outs/filtered_feature_bc_matrix/')
RP0231_Coc_F_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_41_48/RP0231_Coc_F_r2/outs/filtered_feature_bc_matrix/')
RP0446_Coc_F_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_41_48/RP0446_Coc_F_r2/outs/filtered_feature_bc_matrix/')
RP0897_Suc_M_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_41_48/RP0897_Suc_M_r1/outs/filtered_feature_bc_matrix/')
RP0883_Suc_M_r2_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_41_48/RP0883_Suc_M_r2/outs/filtered_feature_bc_matrix/')
RP0609_Coc_M_r1_Data <- Read10X(data.dir='/data/Palmetto_sync/Projects/vshanka_cocaine_scrnadna/Counts/SH10X_multiome_41_48/RP0609_Coc_M_r1/outs/filtered_feature_bc_matrix/')

#Rename cell IDs

colnames(RP0457_Coc_F_r2_Data$`Gene Expression`) <- paste0(colnames(RP0457_Coc_F_r2_Data$`Gene Expression`),"_","RP0457_Coc_F_r2")
colnames(RP0231_Suc_M_r2_Data$`Gene Expression`) <- paste0(colnames(RP0231_Suc_M_r2_Data$`Gene Expression`),"_","RP0231_Suc_M_r2")
colnames(RP0609_Coc_M_r2_Data$`Gene Expression`) <- paste0(colnames(RP0609_Coc_M_r2_Data$`Gene Expression`),"_","RP0609_Coc_M_r2")
colnames(RP0231_Suc_F_r2_Data$`Gene Expression`) <- paste0(colnames(RP0231_Suc_F_r2_Data$`Gene Expression`),"_","RP0231_Suc_F_r2")
colnames(RP0446_Coc_M_r2_Data$`Gene Expression`) <- paste0(colnames(RP0446_Coc_M_r2_Data$`Gene Expression`),"_","RP0446_Coc_M_r2")
colnames(RP0457_Coc_M_r2_Data$`Gene Expression`) <- paste0(colnames(RP0457_Coc_M_r2_Data$`Gene Expression`),"_","RP0457_Coc_M_r2")
colnames(RP0897_Suc_M_r2_Data$`Gene Expression`) <- paste0(colnames(RP0897_Suc_M_r2_Data$`Gene Expression`),"_","RP0897_Suc_M_r2")
colnames(RP0897_Coc_M_r1_Data$`Gene Expression`) <- paste0(colnames(RP0897_Coc_M_r1_Data$`Gene Expression`),"_","RP0897_Coc_M_r1")
colnames(RP0231_Coc_F_r1_Data$`Gene Expression`) <- paste0(colnames(RP0231_Coc_F_r1_Data$`Gene Expression`),"_","RP0231_Coc_F_r1")
colnames(RP0897_Coc_F_r2_Data$`Gene Expression`) <- paste0(colnames(RP0897_Coc_F_r2_Data$`Gene Expression`),"_","RP0897_Coc_F_r2")
colnames(RP0446_Suc_M_r2_Data$`Gene Expression`) <- paste0(colnames(RP0446_Suc_M_r2_Data$`Gene Expression`),"_","RP0446_Suc_M_r2")
colnames(RP0457_Coc_F_r1_Data$`Gene Expression`) <- paste0(colnames(RP0457_Coc_F_r1_Data$`Gene Expression`),"_","RP0457_Coc_F_r1")
colnames(RP0457_Suc_M_r2_Data$`Gene Expression`) <- paste0(colnames(RP0457_Suc_M_r2_Data$`Gene Expression`),"_","RP0457_Suc_M_r2")
colnames(RP0883_Suc_M_r1_Data$`Gene Expression`) <- paste0(colnames(RP0883_Suc_M_r1_Data$`Gene Expression`),"_","RP0883_Suc_M_r1")
colnames(RP0457_Coc_M_r1_Data$`Gene Expression`) <- paste0(colnames(RP0457_Coc_M_r1_Data$`Gene Expression`),"_","RP0457_Coc_M_r1")
colnames(RP0457_Suc_F_r2_Data$`Gene Expression`) <- paste0(colnames(RP0457_Suc_F_r2_Data$`Gene Expression`),"_","RP0457_Suc_F_r2")
colnames(RP0609_Coc_F_r1_Data$`Gene Expression`) <- paste0(colnames(RP0609_Coc_F_r1_Data$`Gene Expression`),"_","RP0609_Coc_F_r1")
colnames(RP0897_Suc_F_r1_Data$`Gene Expression`) <- paste0(colnames(RP0897_Suc_F_r1_Data$`Gene Expression`),"_","RP0897_Suc_F_r1")
colnames(RP0231_Suc_M_r1_Data$`Gene Expression`) <- paste0(colnames(RP0231_Suc_M_r1_Data$`Gene Expression`),"_","RP0231_Suc_M_r1")
colnames(RP0446_Suc_F_r2_Data$`Gene Expression`) <- paste0(colnames(RP0446_Suc_F_r2_Data$`Gene Expression`),"_","RP0446_Suc_F_r2")
colnames(RP0446_Coc_M_r1_Data$`Gene Expression`) <- paste0(colnames(RP0446_Coc_M_r1_Data$`Gene Expression`),"_","RP0446_Coc_M_r1")
colnames(RP0883_Suc_F_r1_Data$`Gene Expression`) <- paste0(colnames(RP0883_Suc_F_r1_Data$`Gene Expression`),"_","RP0883_Suc_F_r1")
colnames(RP0231_Suc_F_r1_Data$`Gene Expression`) <- paste0(colnames(RP0231_Suc_F_r1_Data$`Gene Expression`),"_","RP0231_Suc_F_r1")
colnames(RP0897_Suc_M_r3_Data$`Gene Expression`) <- paste0(colnames(RP0897_Suc_M_r3_Data$`Gene Expression`),"_","RP0897_Suc_M_r3")
colnames(RP0446_Coc_F_r1_Data$`Gene Expression`) <- paste0(colnames(RP0446_Coc_F_r1_Data$`Gene Expression`),"_","RP0446_Coc_F_r1")
colnames(RP0897_Suc_F_r2_Data$`Gene Expression`) <- paste0(colnames(RP0897_Suc_F_r2_Data$`Gene Expression`),"_","RP0897_Suc_F_r2")
colnames(RP0883_Coc_F_r2_Data$`Gene Expression`) <- paste0(colnames(RP0883_Coc_F_r2_Data$`Gene Expression`),"_","RP0883_Coc_F_r2")
colnames(RP0457_Suc_F_r1_Data$`Gene Expression`) <- paste0(colnames(RP0457_Suc_F_r1_Data$`Gene Expression`),"_","RP0457_Suc_F_r1")
colnames(RP0883_Suc_F_r2_Data$`Gene Expression`) <- paste0(colnames(RP0883_Suc_F_r2_Data$`Gene Expression`),"_","RP0883_Suc_F_r2")
colnames(RP0231_Coc_M_r2_Data$`Gene Expression`) <- paste0(colnames(RP0231_Coc_M_r2_Data$`Gene Expression`),"_","RP0231_Coc_M_r2")
colnames(RP0609_Suc_M_r1_Data$`Gene Expression`) <- paste0(colnames(RP0609_Suc_M_r1_Data$`Gene Expression`),"_","RP0609_Suc_M_r1")
colnames(RP0897_Coc_F_r1_Data$`Gene Expression`) <- paste0(colnames(RP0897_Coc_F_r1_Data$`Gene Expression`),"_","RP0897_Coc_F_r1")
colnames(RP0883_Coc_M_r2_Data$`Gene Expression`) <- paste0(colnames(RP0883_Coc_M_r2_Data$`Gene Expression`),"_","RP0883_Coc_M_r2")
colnames(RP0457_Suc_M_r1_Data$`Gene Expression`) <- paste0(colnames(RP0457_Suc_M_r1_Data$`Gene Expression`),"_","RP0457_Suc_M_r1")
colnames(RP0609_Coc_F_r2_Data$`Gene Expression`) <- paste0(colnames(RP0609_Coc_F_r2_Data$`Gene Expression`),"_","RP0609_Coc_F_r2")
colnames(RP0231_Coc_M_r1_Data$`Gene Expression`) <- paste0(colnames(RP0231_Coc_M_r1_Data$`Gene Expression`),"_","RP0231_Coc_M_r1")
colnames(RP0609_Suc_F_r1_Data$`Gene Expression`) <- paste0(colnames(RP0609_Suc_F_r1_Data$`Gene Expression`),"_","RP0609_Suc_F_r1")
colnames(RP0446_Suc_M_r1_Data$`Gene Expression`) <- paste0(colnames(RP0446_Suc_M_r1_Data$`Gene Expression`),"_","RP0446_Suc_M_r1")
colnames(RP0897_Coc_M_r2_Data$`Gene Expression`) <- paste0(colnames(RP0897_Coc_M_r2_Data$`Gene Expression`),"_","RP0897_Coc_M_r2")
colnames(RP0609_Suc_M_r2_Data$`Gene Expression`) <- paste0(colnames(RP0609_Suc_M_r2_Data$`Gene Expression`),"_","RP0609_Suc_M_r2")
colnames(RP0883_Coc_F_r1_Data$`Gene Expression`) <- paste0(colnames(RP0883_Coc_F_r1_Data$`Gene Expression`),"_","RP0883_Coc_F_r1")
colnames(RP0609_Suc_F_r2_Data$`Gene Expression`) <- paste0(colnames(RP0609_Suc_F_r2_Data$`Gene Expression`),"_","RP0609_Suc_F_r2")
colnames(RP0883_Coc_M_r1_Data$`Gene Expression`) <- paste0(colnames(RP0883_Coc_M_r1_Data$`Gene Expression`),"_","RP0883_Coc_M_r1")
colnames(RP0446_Suc_F_r1_Data$`Gene Expression`) <- paste0(colnames(RP0446_Suc_F_r1_Data$`Gene Expression`),"_","RP0446_Suc_F_r1")
colnames(RP0231_Coc_F_r2_Data$`Gene Expression`) <- paste0(colnames(RP0231_Coc_F_r2_Data$`Gene Expression`),"_","RP0231_Coc_F_r2")
colnames(RP0446_Coc_F_r2_Data$`Gene Expression`) <- paste0(colnames(RP0446_Coc_F_r2_Data$`Gene Expression`),"_","RP0446_Coc_F_r2")
colnames(RP0897_Suc_M_r1_Data$`Gene Expression`) <- paste0(colnames(RP0897_Suc_M_r1_Data$`Gene Expression`),"_","RP0897_Suc_M_r1")
colnames(RP0883_Suc_M_r2_Data$`Gene Expression`) <- paste0(colnames(RP0883_Suc_M_r2_Data$`Gene Expression`),"_","RP0883_Suc_M_r2")
colnames(RP0609_Coc_M_r1_Data$`Gene Expression`) <- paste0(colnames(RP0609_Coc_M_r1_Data$`Gene Expression`),"_","RP0609_Coc_M_r1")

count_list <- c(RP0457_Coc_F_r2_Data$`Gene Expression`,
RP0231_Suc_M_r2_Data$`Gene Expression`,
RP0609_Coc_M_r2_Data$`Gene Expression`,
RP0231_Suc_F_r2_Data$`Gene Expression`,
RP0446_Coc_M_r2_Data$`Gene Expression`,
RP0457_Coc_M_r2_Data$`Gene Expression`,
RP0897_Suc_M_r2_Data$`Gene Expression`,
RP0897_Coc_M_r1_Data$`Gene Expression`,
RP0231_Coc_F_r1_Data$`Gene Expression`,
RP0897_Coc_F_r2_Data$`Gene Expression`,
RP0446_Suc_M_r2_Data$`Gene Expression`,
RP0457_Coc_F_r1_Data$`Gene Expression`,
RP0457_Suc_M_r2_Data$`Gene Expression`,
RP0883_Suc_M_r1_Data$`Gene Expression`,
RP0457_Coc_M_r1_Data$`Gene Expression`,
RP0457_Suc_F_r2_Data$`Gene Expression`,
RP0609_Coc_F_r1_Data$`Gene Expression`,
RP0897_Suc_F_r1_Data$`Gene Expression`,
RP0231_Suc_M_r1_Data$`Gene Expression`,
RP0446_Suc_F_r2_Data$`Gene Expression`,
RP0446_Coc_M_r1_Data$`Gene Expression`,
RP0883_Suc_F_r1_Data$`Gene Expression`,
RP0231_Suc_F_r1_Data$`Gene Expression`,
RP0897_Suc_M_r3_Data$`Gene Expression`,
RP0446_Coc_F_r1_Data$`Gene Expression`,
RP0897_Suc_F_r2_Data$`Gene Expression`,
RP0883_Coc_F_r2_Data$`Gene Expression`,
RP0457_Suc_F_r1_Data$`Gene Expression`,
RP0883_Suc_F_r2_Data$`Gene Expression`,
RP0231_Coc_M_r2_Data$`Gene Expression`,
RP0609_Suc_M_r1_Data$`Gene Expression`,
RP0897_Coc_F_r1_Data$`Gene Expression`,
RP0883_Coc_M_r2_Data$`Gene Expression`,
RP0457_Suc_M_r1_Data$`Gene Expression`,
RP0609_Coc_F_r2_Data$`Gene Expression`,
RP0231_Coc_M_r1_Data$`Gene Expression`,
RP0609_Suc_F_r1_Data$`Gene Expression`,
RP0446_Suc_M_r1_Data$`Gene Expression`,
RP0897_Coc_M_r2_Data$`Gene Expression`,
RP0609_Suc_M_r2_Data$`Gene Expression`,
RP0883_Coc_F_r1_Data$`Gene Expression`,
RP0609_Suc_F_r2_Data$`Gene Expression`,
RP0883_Coc_M_r1_Data$`Gene Expression`,
RP0446_Suc_F_r1_Data$`Gene Expression`,
RP0231_Coc_F_r2_Data$`Gene Expression`,
RP0446_Coc_F_r2_Data$`Gene Expression`,
RP0897_Suc_M_r1_Data$`Gene Expression`,
RP0883_Suc_M_r2_Data$`Gene Expression`,
RP0609_Coc_M_r1_Data$`Gene Expression`)

names(count_list) <- c("RP0457_Coc_F_r2",
"RP0231_Suc_M_r2",
"RP0609_Coc_M_r2",
"RP0231_Suc_F_r2",
"RP0446_Coc_M_r2",
"RP0457_Coc_M_r2",
"RP0897_Suc_M_r2",
"RP0897_Coc_M_r1",
"RP0231_Coc_F_r1",
"RP0897_Coc_F_r2",
"RP0446_Suc_M_r2",
"RP0457_Coc_F_r1",
"RP0457_Suc_M_r2",
"RP0883_Suc_M_r1",
"RP0457_Coc_M_r1",
"RP0457_Suc_F_r2",
"RP0609_Coc_F_r1",
"RP0897_Suc_F_r1",
"RP0231_Suc_M_r1",
"RP0446_Suc_F_r2",
"RP0446_Coc_M_r1",
"RP0883_Suc_F_r1",
"RP0231_Suc_F_r1",
"RP0897_Suc_M_r3",
"RP0446_Coc_F_r1",
"RP0897_Suc_F_r2",
"RP0883_Coc_F_r2",
"RP0457_Suc_F_r1",
"RP0883_Suc_F_r2",
"RP0231_Coc_M_r2",
"RP0609_Suc_M_r1",
"RP0897_Coc_F_r1",
"RP0883_Coc_M_r2",
"RP0457_Suc_M_r1",
"RP0609_Coc_F_r2",
"RP0231_Coc_M_r1",
"RP0609_Suc_F_r1",
"RP0446_Suc_M_r1",
"RP0897_Coc_M_r2",
"RP0609_Suc_M_r2",
"RP0883_Coc_F_r1",
"RP0609_Suc_F_r2",
"RP0883_Coc_M_r1",
"RP0446_Suc_F_r1",
"RP0231_Coc_F_r2",
"RP0446_Coc_F_r2",
"RP0897_Suc_M_r1",
"RP0883_Suc_M_r2",
"RP0609_Coc_M_r1")

#Standard processing and object creation
object <- CreateSeuratObject(counts = count_list)
object <- subset(object, subset = nFeature_RNA > 200 & nFeature_RNA < 2000)
object <- NormalizeData(object)
object <- FindVariableFeatures(object, verbose = FALSE)
#save.image("sketch_pipeline_1.RData")
#load("sketch_pipeline.RData")
object <- SketchData(object = object, ncells = 5000, method = "LeverageScore", sketched.assay = "sketch")

#Sketch Processing (https://satijalab.org/seurat/articles/parsebio_sketch_integration)
DefaultAssay(object) <- "sketch"
object <- FindVariableFeatures(object, verbose = F)
object <- ScaleData(object, verbose = F)
object <- RunPCA(object, verbose = F)
# integrate the datasets using the first sample as the reference. I think the reference is just a starting point for convergence.
object <- IntegrateLayers(object, method = RPCAIntegration, orig = "pca", new.reduction = "integrated.rpca",
    dims = 1:30, k.anchor = 20, reference = which(Layers(object, search = "data") %in% c("data.RP0457_Coc_F_r2")),
    verbose = F)
save.image("sketch_integrated_1.RData")
# cluster the integrated data

load("sketch_integrated_1.RData")
object <- FindNeighbors(object, reduction = "integrated.rpca", dims = 1:30)
object <- FindClusters(object, resolution = 2)
object <- RunUMAP(object, reduction = "integrated.rpca", dims = 1:30, return.model = T, verbose = F)

# you can now rejoin the layers in the sketched assay this is required to perform differential
# expression
object[["sketch"]] <- JoinLayers(object[["sketch"]])
#object <- FindAllMarkers(object = object, max.cells.per.ident = 500, only.pos = TRUE)
object$celltype.manual <- Idents(object)
#save.image("sketch_integrated_2.RData")

extract_sampleID <- function(input_string) {
  parts <- unlist(strsplit(input_string, "_"))
  concatenated <- paste(parts[2:5], collapse = "_")
  return(concatenated)
}

test <- sapply(names(object$orig.ident),extract_sampleID)
object$sample <- test

extract_genotype <- function(input_string) {
  parts <- unlist(strsplit(input_string, "_"))
  concatenated <- paste(parts[2], collapse = "_")
  return(concatenated)
}

test2 <- sapply(names(object$orig.ident),extract_genotype)
object$genotype <- test2

extract_treatment <- function(input_string) {
  parts <- unlist(strsplit(input_string, "_"))
  concatenated <- paste(parts[3], collapse = "_")
  return(concatenated)
}
 test3 <- sapply(names(object$orig.ident),extract_treatment)
object$treatment <- test3

extract_sex <- function(input_string) {
  parts <- unlist(strsplit(input_string, "_"))
  concatenated <- paste(parts[4], collapse = "_")
  return(concatenated)
}

test4 <- sapply(names(object$orig.ident),extract_sex)
object$sex <- test4

# resplit the sketched cell assay into layers this is required to project the integration onto
# all cells
object[["sketch"]] <- split(object[["sketch"]], f = object$sample)

object <- ProjectIntegration(object = object, sketched.assay = "sketch", assay = "RNA", reduction = "integrated.rpca")


object <- ProjectData(object = object, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "integrated.rpca.full",
    full.reduction = "integrated.rpca.full", dims = 1:30, refdata = list(celltype.full = "celltype.manual"))

object <- RunUMAP(object, reduction = "integrated.rpca.full", dims = 1:30, reduction.name = "umap.full",
    reduction.key = "UMAP_full_")

p1 <- DimPlot(object, reduction = "umap.full", group.by = "sample", alpha = 0.1,raster=FALSE) + ggtitle("UMAP by Sample ID")
p2 <- DimPlot(object, reduction = "umap.full", group.by = "celltype.full", alpha = 0.1,raster=FALSE,order = order(levels(factor(object$celltype.full,ordered = FALSE)),decreasing=TRUE)) + ggtitle("UMAP by Celltype Cluster")
p1 + p2 + plot_layout(ncol = 1)

save.image("sketch_integrate_3.RData")
