#Presumably you have called all relevant libraries in earlier script, if not, do so.
#library(reticulate)

#make sure you're using a version of python which is less than 3.9 or fa2 won't work. Check system preferences and which environment it's pointing at. 

#Make sure you set up your C++ references to clang using: usethis::edit_r_makevars() with directory directions as detailed by Supatt

#run this to make sure you have scanorama in the console/terminal
#pip3 install scanorama 

#conda_create("r-reticulate") #only use if you haven't used a reticulate environment yet otherwise... 
#use_condaenv("r-reticulate")


#py_install("fa2")
#py_install("networkx")
#py_install("scrublet")

#install.packages("magick")
#library(magick)
library(devtools)
library(SingCellaR)
library(magick)

#Additional packages for the singleCellaR integration
library(Seurat)
library(harmony)
library(liger)
library(sva)

org_comp <- new("SingCellaR_Int")
org_comp@dir_path_SingCellR_object_files<-"/Users/khanas/Desktop/scRNAseq_2021/BM_scRNA_count/Count_RProject/SingleCellaR/SingleCellaR_BM/Integration_Files_TKI"
org_comp@SingCellR_object_files=c("VEGFA_integration.rds",
                                     "VEGFAC_integration.rds"
                                     )

preprocess_integration(org_comp)
org_comp

#Subsample cells
#Vsubsample_cells(org_comp,n.subsample = 8000,n.seed = 42)

#Call correct object - may have to come back to this as accidentally called on to VEGFAC, shouldn't make a difference but still
filter_cells_and_genes(org_comp,
                       min_UMIs=500,
                       max_UMIs=50000,
                       min_detected_genes=300,
                       max_detected_genes=6000,
                       max_percent_mito=12,
                       isRemovedDoublets = FALSE,
                       genes_with_expressing_cells = 10)


#Annotate samples: 
cell_anno.info<-get_cells_annotation(org_comp)
head(cell_anno.info)

table(cell_anno.info$sampleID)

org_comp@meta.data$status[org_comp@meta.data$sampleID=="1_VEGFA"]<-"VEGFA"
org_comp@meta.data$status[org_comp@meta.data$sampleID=="1_VEGFAC"]<-"VEGFAC"

cell_anno.info<-get_cells_annotation(org_comp)
head(cell_anno.info)

normalize_UMIs(org_comp,use.scaled.factor = T)

#There is another option for get variable genes but it isn't clear what it is
get_variable_genes_by_fitting_GLM_model(org_comp,mean_expr_cutoff = 0.1,disp_zscore_cutoff = 0.1)

tiff("Integrated_Variable_genes.tiff", units="in", width=10, height=10, res=300)
plot_variable_genes(org_comp)
dev.off()

runPCA(org_comp,use.components=50,use.regressout.data = F)

tiff("Integrated_ElbowPlot.tiff", units="in", width=10, height=10, res=300)
plot_PCA_Elbowplot(org_comp)
dev.off()

SingCellaR::runUMAP(org_comp,dim_reduction_method = "pca",n.dims.use = 50,n.neighbors = 30,
                    uwot.metric = "euclidean")


tiff("Integrated_UMAP_byStatus.tiff", units="in", width=10, height=10, res=300)
plot_umap_label_by_a_feature_of_interest(org_comp,feature = "status",point.size = 0.5)
dev.off()


#Harmony integration
runHarmony(org_comp,covariates = c("sampleID"),n.dims.use = 50,harmony.max.iter = 20,n.seed = 1)


SingCellaR::runUMAP(org_comp,useIntegrativeEmbeddings = T, integrative_method = "harmony",n.dims.use = 50,
                    n.neighbors = 30,uwot.metric = "euclidean")

tiff("Integrated_UMAP_HARMONY.tiff", units="in", width=10, height=10, res=300)
plot_umap_label_by_a_feature_of_interest(org_comp,feature = "status", mark.feature = FALSE, point.size = 0.5)
dev.off()

runFA2_ForceDirectedGraph(org_comp,n.dims.use = 20, useIntegrativeEmbeddings = T,integrative_method = "harmony",n.neighbors = 5,n.seed = 1,fa2_n_iter = 1000)

tiff("Integrated_Force_directed_unlabelled.tiff", units="in", width=10, height=10, res=300)
plot_forceDirectedGraph_label_by_clusters(org_comp,show_method = "louvain",vertex.size = 0.85, mark.clusters = FALSE,
                                          background.color = "white")
dev.off()


#IF sticking to harmony then 
identifyClusters(org_comp,useIntegrativeEmbeddings = T,integrative_method = "harmony", n.dims.use = 50,n.neighbors = 50, knn.metric = "euclidean")

tiff("Integrated_Harmony_UMAP.tiff", units="in", width=10, height=10, res=300)
plot_umap_label_by_clusters(org_comp,show_method = "louvain")
dev.off()

tiff("Integrated_Harmony_UMAP_unlabelled.tiff", units="in", width=10, height=10, res=300)
plot_umap_label_by_clusters(org_comp,show_method = "louvain", mark.clusters = FALSE)
dev.off()

tiff("Integrated_Harmony_UMAP_byStatus.tiff", units="in", width=10, height=10, res=300)
plot_umap_label_by_a_feature_of_interest(org_comp,feature = "status",point.size = 0.5, mark.feature = FALSE)
dev.off()

findMarkerGenes(org_comp,cluster.type = "louvain")

runFA2_ForceDirectedGraph(org_comp,n.dims.use = 50, useIntegrativeEmbeddings = T,integrative_method = "harmony",n.neighbors = 5,n.seed = 1,fa2_n_iter = 1000)

tiff("Integrated_Harmony_ForceDirected.tiff", units="in", width=10, height=10, res=300)
plot_forceDirectedGraph_label_by_clusters(org_comp,show_method = "louvain")
dev.off()

tiff("Integrated_Harmony_ForceDirected_unlabelled.tiff", units="in", width=10, height=10, res=300)
plot_forceDirectedGraph_label_by_clusters(org_comp,show_method = "louvain",mark.clusters = FALSE)
dev.off()

findMarkerGenes(org_comp,cluster.type = "louvain")

tiff("Integrated_Harmony_Heatmap_for_Marker_Genes.tiff", units="in", width=10, height=10, res=300)
plot_heatmap_for_marker_genes(org_comp,cluster.type = "louvain",n.TopGenes = 5,rowFont.size = 5)
dev.off()


tiff("Integrated_feature_plot_stroma.tiff", units="in", width=10, height=10, res=300)
plot_umap_label_by_genes(org_comp,gene_list = c("CXCL12","COL3A1","CDH5","ENG"))
dev.off()


tiff("Integrated_feature_plot_haemato.tiff", units="in", width=10, height=10, res=300)
plot_umap_label_by_genes(org_comp,gene_list = c("PPBP","PF4","CD14","GNLY"))
dev.off()

### ANALYSIS FOR CELLPHONE DB 11TH JAN ###

### GSEA For Combined Data Sets - org_comp ####

pre_rankedGenes_for_GSEA_org_comp<-identifyGSEAPrerankedGenes_for_all_clusters(org_comp,
                                                                                  cluster.type = "louvain")
#My own curated list
fgsea_Results_total_1<-Run_fGSEA_for_multiple_comparisons(GSEAPrerankedGenes_list = pre_rankedGenes_for_GSEA_org_comp,minSize = 10,
                                                         gmt.file = "/Users/khanas/Desktop/scRNAseq_2021/BM_scRNA_count/SingleCellaR_supp_data/GSEA_Set_1.gmt", eps = 0)
tiff("VEGFA_VEGFAC_GSEA_list_AB.tiff", units="in", width=10, height=10, res=300)
plot_heatmap_for_fGSEA_all_clusters(fgsea_Results_total_1,isApplyCutoff = TRUE,
                                    use_pvalues_for_clustering=T,
                                    show_NES_score = T,fontsize_row =7,
                                    adjusted_pval = 0.1,
                                    show_only_NES_positive_score = T,format.digits = 3,
                                    clustering_method = "ward.D",
                                    clustering_distance_rows = "euclidean",
                                    clustering_distance_cols = "euclidean",show_text_for_ns = F)
dev.off()

##### rename clusters as per GSEA #####

org_comp@sc.clusters$annotations[org_comp@sc.clusters$louvain_cluster %in% "cl1"] <- "HSPC 1"
org_comp@sc.clusters$annotations[org_comp@sc.clusters$louvain_cluster %in% "cl2"] <- "Late Erythroid"
org_comp@sc.clusters$annotations[org_comp@sc.clusters$louvain_cluster %in% "cl3"] <- "Myeloid Progenitor 1"
org_comp@sc.clusters$annotations[org_comp@sc.clusters$louvain_cluster %in% "cl4"] <- "Early Erythroid"
org_comp@sc.clusters$annotations[org_comp@sc.clusters$louvain_cluster %in% "cl5"] <- "Myeloid Progenitor 2"
org_comp@sc.clusters$annotations[org_comp@sc.clusters$louvain_cluster %in% "cl6"] <- "Megakaryocyte 1"
org_comp@sc.clusters$annotations[org_comp@sc.clusters$louvain_cluster %in% "cl7"] <- "MSC"
org_comp@sc.clusters$annotations[org_comp@sc.clusters$louvain_cluster %in% "cl8"] <- "Monocyte"
org_comp@sc.clusters$annotations[org_comp@sc.clusters$louvain_cluster %in% "cl9"] <- "Myeloid Progenitor 3"
org_comp@sc.clusters$annotations[org_comp@sc.clusters$louvain_cluster %in% "cl10"] <- "Fibroblast 1"
org_comp@sc.clusters$annotations[org_comp@sc.clusters$louvain_cluster %in% "cl11"] <- "Endothelium"
org_comp@sc.clusters$annotations[org_comp@sc.clusters$louvain_cluster %in% "cl12"] <- "Mid Erythroid"
org_comp@sc.clusters$annotations[org_comp@sc.clusters$louvain_cluster %in% "cl13"] <- "Megakaryocyte 2"
org_comp@sc.clusters$annotations[org_comp@sc.clusters$louvain_cluster %in% "cl14"] <- "Fibroblast 2"
org_comp@sc.clusters$annotations[org_comp@sc.clusters$louvain_cluster %in% "cl15"] <- "HSPC 2"

## Write integrated object
save(org_comp,file="/Users/khanas/Desktop/scRNAseq_2021/BM_scRNA_count/Integrated_VEGFA_C_R_Object.rds")


#### Make a CellPhone DB Thingy ####
### Work out What You Need from Previous Code ####

normalize_UMIs(org_comp,use.scaled.factor = T)

#sc.clusters.stromal <- VEGFAC_stroma@sc.clusters
#sc.clusters.hemat <- VEGFAC_haem@sc.clusters

#head(sc.clusters.stromal)
#head(sc.clusters.hemat)
###########################

head(org_comp@meta.data)
head(org_comp@sc.clusters)
dim(org_comp@sc.clusters)
#dim(sc.clusters.hemat)
#dim(sc.clusters.stromal)
############################

#sc.clusters <- rbind(sc.clusters.stromal,sc.clusters.hemat)
#dim(sc.clusters)


meta.data.comp <- org_comp@meta.data
org.sc.clusters  <- org_comp@sc.clusters
meta.data.comp <- merge(meta.data.comp,org.sc.clusters, by = "Cell")
dim(meta.data.comp)


head(org_comp@meta.data)
head(org_comp@sc.clusters)
head(meta.data.comp)
#Prepare for CellPhoneDB - meta data! With status for VEGF samples 
meta_data_comp_db <- meta.data.comp[,c('Cell','annotations','status')]
write.table(meta_data_comp_db, '/Users/khanas/Desktop/scRNAseq_2021/cellphonedb_integrated_VEGF_data/cellphonedb_meta.txt', sep = '\t',quote=F, row.names=F)

#Normalise and sort out counts
norm_counts_comp <- get_normalized_umi(org_comp) #this will give normalised counts from CellPhoneDB
norm_counts_comp <- norm_counts_comp[,meta_data_comp_db$Cell] #merge ,meta.data with 
# ########h ead(norm_counts) this may be where you have to come back to if you run into trouble 
norm_counts_df_comp <- data.frame(norm_counts_comp, check.names = F)
norm_counts_df_comp <- cbind('Gene' = rownames(norm_counts_df_comp),norm_counts_df_comp) #This will add a column with Gene as title from rownames which is gene values

#load gene name if you don't have the relevant data frames: 
#genes <- read.table("/Users/khanas/Desktop/scRNAseq_2021/BM_scRNA_count/VEGFAC/filtered_feature_bc_matrix/features.tsv.gz", stringsAsFactors = FALSE) # read in genes... 
#gene_names <- genes[['V2']] #create a list of gene names as a string
#gene_ensembl <- genes[['V1']] #create a list of ENSEMBL IDs as a string 

norm_counts_df_comp$Gene[norm_counts_df_comp$Gene == genes$V2] <-genes$V1[norm_counts_df_comp$Gene == genes$V2]
head(norm_counts_df_comp[35000:35010,1:10])
#Swap out contents of 'Gene'
#norm_counts_df[['Gene']] <- plyr::mapvalues(x =umap.results[['Gene']], from = gene_names, to = gene_ensembl) #plyr::mapvalues is positional
#norm_counts_df$Gene[norm_counts_df$Gene == gene_names] <- gene_ensembl
#don't work
head(norm_counts_df_comp[1:10,1:10])

#You need to make the Gene column row names: 
norm_counts_df_comp_II <- norm_counts_df_comp
row.names(norm_counts_df_comp_II) <- norm_counts_df_comp_II$Gene #convert row names to Ensemble IDs
norm_counts_df_comp_II$Gene <- NULL ## Drops the column called Gene
dim(norm_counts_df_comp_II) # Check this is the same dimension as meta 
head(norm_counts_df_comp_II[1:10,1:10])
dim(meta_data_comp_db) # Check this is the same dimension as meta 
head(meta_data_comp_db)

write.table(norm_counts_df_comp, '/Users/khanas/Desktop/scRNAseq_2021/cellphonedb_integrated_VEGF_data/cellphonedb_count.txt', sep = '\t',quote=F, row.names = F)


#use isna to find if there are any incorrect values ...

