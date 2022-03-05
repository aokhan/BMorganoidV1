#Try and separate force directed graph : 

VEGFAC_haem <- sub_clusters(
  VEGFAC,
  cluster.type = c("louvain"),
  cluster_ids = c('cl1','cl2','cl3','cl4','cl5','cl7','cl8','cl11','cl12','cl15','cl16','cl17') 
)
#Not sure about the rest 

process_cells_annotation(VEGFAC_haem,mitochondiral_genes_start_with="MT-")

plot_cells_annotation(VEGFAC_haem,type="histogram")

plot_UMIs_vs_Detected_genes(VEGFAC_haem)

filter_cells_and_genes(VEGFAC_haem,
                       min_UMIs=500,
                       max_UMIs=50000,
                       min_detected_genes=300,
                       max_detected_genes=6000,
                       max_percent_mito=12,
                       isRemovedDoublets = FALSE,
                       genes_with_expressing_cells = 10)

normalize_UMIs(VEGFAC_haem,use.scaled.factor = T)

remove_unwanted_confounders(VEGFAC_haem,residualModelFormulaStr="~UMI_count+percent_mito")


get_variable_genes_by_fitting_GLM_model(VEGFAC_haem,mean_expr_cutoff = 0.05,disp_zscore_cutoff = 0.05)

tiff("VEGFAC_haem_Variable_Gene.tiff", units="in", width=10, height=10, res=300)
plot_variable_genes(VEGFAC_haem)
dev.off()



runPCA(VEGFAC_haem,use.components=50,use.regressout.data = T)

tiff("VEGFAC_haem_Elbow_Plot", units="in", width=10, height=10, res=300)
plot_PCA_Elbowplot(VEGFAC_haem)
dev.off()

#UMAP

runUMAP(VEGFAC_haem,dim_reduction_method = "pca",n.dims.use = 40,n.neighbors = 120,
        uwot.metric = "euclidean")
#lot_umap_label_by_a_feature_of_interest(VEGFAC_haem,feature = "UMI_count",point.size = 0.1)
identifyClusters(VEGFAC_haem,n.dims.use = 40,n.neighbors = 120,knn.metric = "euclidean")

plot_umap_label_by_clusters(VEGFAC_haem,show_method = "louvain",point.size = 3)

tiff("VEGFAC_haem_UMAP.tiff", units="in", width=6, height=5, res=300)
plot_umap_label_by_clusters(VEGFAC_haem,show_method = "louvain",point.size = 3)
dev.off()

tiff("VEGFAC_haem_UMAP_no_label.tiff", units="in", width=6, height=5, res=300)
plot_umap_label_by_clusters(VEGFAC_haem,show_method = "louvain",point.size = 3, mark.clusters = FALSE)
dev.off()

runFA2_ForceDirectedGraph(VEGFAC_haem,n.dims.use = 30,
                          n.neighbors = 5,n.seed = 25,fa2_n_iter = 1000)

tiff("VEGFAC_haem_Force_directed.tiff", units="in", width=6, height=5, res=300)
plot_forceDirectedGraph_label_by_clusters(VEGFAC_haem,show_method = "louvain",vertex.size = 2,
                                          background.color = "white")
dev.off()

#VEGFAC GSEA Pre-Ranked Genes
pre_rankedGenes_for_GSEA_VEGFAC_haem<-identifyGSEAPrerankedGenes_for_all_clusters(VEGFAC_haem,
                                                                      cluster.type = "louvain")
#My own curated list
fgsea_Results_haem_1<-Run_fGSEA_for_multiple_comparisons(GSEAPrerankedGenes_list = pre_rankedGenes_for_GSEA_VEGFAC_haem,minSize = 10,
                                                    gmt.file = "/Users/khanas/Desktop/scRNAseq_2021/BM_scRNA_count/SingleCellaR_supp_data/GSEA_Set_1.gmt", eps = 0)
tiff("VEGFAC_haem_GSEA_list_AB.tiff", units="in", width=10, height=10, res=300)
plot_heatmap_for_fGSEA_all_clusters(fgsea_Results_haem_1,isApplyCutoff = TRUE,
                                    use_pvalues_for_clustering=T,
                                    show_NES_score = T,fontsize_row =7,
                                    adjusted_pval = 0.1,
                                    show_only_NES_positive_score = T,format.digits = 3,
                                    clustering_method = "ward.D",
                                    clustering_distance_rows = "euclidean",
                                    clustering_distance_cols = "euclidean",show_text_for_ns = F)
dev.off()

write.table(fgsea_Results_haem_1, '/Users/khanas/Desktop/scRNAseq_2021/fgsea_Results_haem_1.txt', sep = '\t',quote=F, row.names = F)



findMarkerGenes(
  VEGFAC_haem,
  cluster.type = c("louvain"),
  min.log2FC = 0.3,
  min.expFraction = 0.25,
)

tiff("VEGFAC_haem_HeatMap.tiff", units="in", width=10, height=10, res=300)
plot_heatmap_for_marker_genes(VEGFAC_haem,cluster.type = "louvain",n.TopGenes = 10,rowFont.size = 5)
dev.off()


tiff("VEGFAC_HAEM_COLOURED_FD_1.tiff", units="in", width=6, height=5, res=300)
plot_forceDirectedGraph_label_by_multiple_gene_sets(VEGFAC_haem,gmt.file = "/Users/khanas/Desktop/scRNAseq_2021/BM_scRNA_count/SingleCellaR_supp_data/gene_set_colour_code.gmt",
                                                    show_gene_sets = c("Erythroid","HSPC","Monocyte/Neutrophil","Megakaryocyte"),
                                                    custom_color = c("red","green","orange","blue"),
                                                    isNormalizedByHouseKeeping = T,vertex.size = 2,edge.size = 0.05,
                                                    background.color = "white")
dev.off()


#rename clusters 

VEGFAC_haem@sc.clusters$annotations[VEGFAC_haem@sc.clusters$louvain_cluster %in% "cl1"] <- "Early Erythroid 1"
VEGFAC_haem@sc.clusters$annotations[VEGFAC_haem@sc.clusters$louvain_cluster %in% "cl2"] <- "Early Erythroid 2"
VEGFAC_haem@sc.clusters$annotations[VEGFAC_haem@sc.clusters$louvain_cluster %in% "cl3"] <- "Late Erythroid 1"
VEGFAC_haem@sc.clusters$annotations[VEGFAC_haem@sc.clusters$louvain_cluster %in% "cl4"] <- "HSPC"
VEGFAC_haem@sc.clusters$annotations[VEGFAC_haem@sc.clusters$louvain_cluster %in% "cl5"] <- "Megakaryocyte 1"
VEGFAC_haem@sc.clusters$annotations[VEGFAC_haem@sc.clusters$louvain_cluster %in% "cl6"] <- "Late Erythroid"
VEGFAC_haem@sc.clusters$annotations[VEGFAC_haem@sc.clusters$louvain_cluster %in% "cl7"] <- "Mid Erythroid"
VEGFAC_haem@sc.clusters$annotations[VEGFAC_haem@sc.clusters$louvain_cluster %in% "cl8"] <- "Monocyte/Neutrophil Progenitors"
VEGFAC_haem@sc.clusters$annotations[VEGFAC_haem@sc.clusters$louvain_cluster %in% "cl9"] <- "Eosinophil/Basophil/Mast Cell Progenitors"
VEGFAC_haem@sc.clusters$annotations[VEGFAC_haem@sc.clusters$louvain_cluster %in% "cl10"] <- "Monocyte"
VEGFAC_haem@sc.clusters$annotations[VEGFAC_haem@sc.clusters$louvain_cluster %in% "cl11"] <- "Megakaryocyte 2"

VEGFAC_haem@sc.clusters$type <- "Haematopoeitic"


tiff("17_11_2021_VEGFAC_HAEM_bubble_3.tiff", units="in", width=5.5, height=5, res=300)
plot_bubble_for_genes_per_cluster(VEGFAC_haem,cluster.type = "louvain",
                                  gene_list=c("CD34","WNT5B",
                                              "PRSS57",
                                              "GYPB","HEMGN","SLC2A1","KLF1",
                                              "GP9","PPBP","PF4","KIT","TPSB2","RUNX1","GATA2","GATA1","ITGAM","CD14"),show.percent = F,buble.scale = 8,point.color1 = "WHITE",
                                  point.color2 = "#e74c3c",
                                  gene.font.size = 7,
                                  axist.x.font.size = 10
                                  )
dev.off()
																												

tiff("Erythroid _test.tiff", units="in", width=5.5, height=5, res=300)
plot_bubble_for_genes_per_cluster(VEGFAC_haem,cluster.type = "louvain",
                                  gene_list=c("HBB","HBG2",
                                              "HBZ","HBM",
                                              "ALAS2","AHSP","SLC25A37",
                                              "GYPA","GYPB","HEMGN","SNCA","PRDX2","SLC4A1","HMBS","BPGM","HBE1","MYL4","BLVRB","SLC25A39","SELENBP1","FECH","EPB42","SLC2A1","UROD","ARL4A","GYPC","MGST3"),show.percent = F,buble.scale = 8,point.color1 = "WHITE",
                                  point.color2 = "#e74c3c",
                                  gene.font.size = 7,
                                  axist.x.font.size = 10
)
dev.off()

## CHECK for MAST CELL genes ##
tiff("17_11_2021_VEGFAC_HAEM_Mast_bubble.tiff", units="in", width=5.5, height=5, res=300)
plot_bubble_for_genes_per_cluster(VEGFAC_haem,cluster.type = "louvain",
                                  gene_list=c("KIT","TPSAB1",
                                              "TPSB2","HDC","FCER1A","PRG2","ERG"),show.percent = F,buble.scale = 8,point.color1 = "WHITE",
                                  point.color2 = "#e74c3c",
                                  gene.font.size = 7,
                                  axist.x.font.size = 10
)
dev.off()


##### ########## Custom Heatmap for GSEA data set ########## #####

library(tidyverse)
library(ComplexHeatmap)
library(tidyr)
library(RColorBrewer)

#Make a DF to work with 
GSEA_Haem <- fgsea_Results_haem_1

GSEA_Haem$cluster[GSEA_Haem$cluster  %in% "cl1"] <- "Early Erythroid 1"
GSEA_Haem$cluster[GSEA_Haem$cluster  %in% "cl2"] <- "Early Erythroid 2"
GSEA_Haem$cluster[GSEA_Haem$cluster  %in% "cl3"] <- "Late Erythroid 1"
GSEA_Haem$cluster[GSEA_Haem$cluster  %in% "cl4"] <- "HSPC"
GSEA_Haem$cluster[GSEA_Haem$cluster  %in% "cl5"] <- "Megakaryocyte 1"
GSEA_Haem$cluster[GSEA_Haem$cluster  %in% "cl6"] <- "Late Erythroid"
GSEA_Haem$cluster[GSEA_Haem$cluster  %in% "cl7"] <- "Mid Erythroid"
GSEA_Haem$cluster[GSEA_Haem$cluster  %in% "cl8"] <- "Monocyte/Neutrophil Progenitors"
GSEA_Haem$cluster[GSEA_Haem$cluster  %in% "cl9"] <- "Eosinophil/Basophil/Mast Cell Progenitors"
GSEA_Haem$cluster[GSEA_Haem$cluster  %in% "cl10"] <- "Monocyte"
GSEA_Haem$cluster[GSEA_Haem$cluster  %in% "cl11"] <- "Megakaryocyte 2"



GSEA_Haem$pathway[GSEA_Haem$pathway %in% "CMP_Chen_et_al"] <- "Chen et al. CMP"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "Curated_HSPC"] <- "Curated HSPC"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "DC1_Popescu_et_al"] <- "Popescu et al. DC1"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "DC_precursor_Popescu_et_al"] <- "Popescu et al. DC Precursor"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "EPPERT_HSC_R"] <- "Eppert et al. HSC"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "EPPERT_PROGENITOR"] <- "Eppert et al. Progenitor"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "Early_Erythroid_Popescu_et_al"] <- "Popescu et al. Early Erythroid"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "Ery_Progenitor1_donor1_Velten_et_al"] <- "Velten et al. Erythroid Progenitor 1 Donor 1"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "Ery_Progenitor2_donor1_Velten_et_al"] <- "Velten et al. Erythroid Progenitor 2 Donor 1"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "GMP_Drissen_et_al"] <- "Drissen et al. GMP"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "HALLMARK_TNFA_SIGNALING_VIA_NFKB"] <- "Hallmark TNF Signaling vis NFKb"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "HAY_BONE_MARROW_CD34_POS_HSC"] <- "Hay et al. Bone Marrow CD34+ HSC"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "HAY_BONE_MARROW_CD34_POS_LMPP"] <- "Hay et al. Bone Marrow CD34+ LMPP"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "HSC/MPP_Popescu_et_al"] <- "Popescu et al. HSC/MPP"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "HSC_Chen_et_al"] <- "Chen et al. HSC"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "HSC_MPP_Drissen_et_al"] <- "Drissen et al. HSC MPP"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "HSPC_MPP"] <- "HSPC MPP Curated"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "KEGG_HEMATOPOIETIC_CELL_LINEAGE"] <- "Kegg Hematopoeitic Cell Lineage"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "Kupffer_Popescu_et_al"] <- "Popescu et al. Kupffer"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "LMPP_Drissen_et_al"] <- "Drissen et al. LMPP"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "Late_Erythroid_Popescu_et_al"] <- "Popescu et al. Late Erythroid"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "MEMP_Popescu_et_al"] <- "Popescu et al. MEMP"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "MEP_Drissen_et_al"] <- "Drissen et al. MEP"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "MK_Chen_et_al"] <- "Chen et al. MK"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "Macrophage_Nano_princeton"] <- "Princeton Macrophage Nano"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "Mast_Popescu_et_al"] <- "Popescu et al. Mast"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "Mid_Erythroid_Popescu_et_al"] <- "Popescu et al. Mid Erythroid"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "MkE_Progenitor1_donor1_Velten_et_al"] <- "Velten et al. MkE Progenitor 1 Donor"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "MkE_Progenitor2_donor1_Velten_et_al"] <- "Velten et al. MkE Progenitor 2 Donor"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "Mk_Popescu_et_al"] <- "Popescu et al. MK"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "Mk_Progenitor_donor1_Velten_et_al"] <- "Velten et al. MK Progenitor Donor 1"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "Mono/Mac_Popescu_et_al"] <- "Popescu et al. Monocyte/Macrophage"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "Monocyte/Dendritic_Progenitor_donor2_Velten_et_al"] <- "Velten et al. Monocyte/Dendritic Progenitor Donor 2"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "Monocyte/Dendritic_donor1_Velten_et_al"] <- "Velten et al. Monocyte/Dendritic Donor 1"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "Monocyte_Nano_princeton"] <- "Princeton Monocyte Nano"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "Monocyte_Popescu_et_al"] <- "Popescu et al. Monocyte"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "Monocyte_precursor_Popescu_et_al"] <- "Popescu et al. Monocyte precursor"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "Myeloid_Progenitor1_donor1_Velten_et_al"] <- "Velten et al. Myeloid Progenitor 1 Donor 1"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "Myeloid_Progenitor2_donor1_Velten_et_al"] <- "Velten et al. Myeloid Progenitor 2 Donor 1"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "Myeloid_Progenitor_donor2_Velten_et_al"] <- "Velten et al. Myeloid Progenitor Donor 2"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "NK_Popescu_et_al"] <- "Popescu et al. NK"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "Neutrophil_Nano_princeton"] <- "Princeton Neutrophil Nano"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "Neutrophil_Porgenitor3_donor1_Velten_et_al"] <- "Velten et al. Neutrophil Progenitor 3 Donor 1"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "Neutrophil_Progenitor1_donor1_Velten_et_al"] <- "Velten et al. Neutrophil Progenitor 2 Donor 1"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "Neutrophil_Progenitor2_donor1_Velten_et_al"] <- "Velten et al. Neutrophil Progenitor 1 Donor 1"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "Neutrophil_myeloid_progenitor_Popescu_et_al"] <- "Popescu et al. Neutrophil Myeloid Progenitor"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "Quiescence"] <- "Hallmark Quiescence"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "VCAM1_EI_macrophage_Popescu_et_al"] <- "Popescu et al. VCAM EI Macrophage"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "WP_HEMATOPOIETIC_STEM_CELL_DIFFERENTIATION"] <- "WP Hematopoeitic Stem Cell Differentiation"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "ZHENG_CORD_BLOOD_C6_HSC_MULTIPOTENT_PROGENITOR"] <- "Zheng et al. Cord Blood C6 HSC Multipotent Progenitor"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "pDC_precursor_Popescu_et_al"] <- "Popescu et al. pDC Precursor"
GSEA_Haem$pathway[GSEA_Haem$pathway %in% "DC2_Popescu_et_al"] <- "Popescu et al. DC2"

######## Working with NES #####
GSEA_Haem_2 <- GSEA_Haem[, c("pathway", "NES","cluster")]
GSEA_Haem_wide <- pivot_wider(GSEA_Haem_2, names_from = "cluster", values_from = "NES") #Reshape data using tidyr package - pivot wider converts from long to wide. 
GSEA_Haem_wide <- GSEA_Haem_wide %>% remove_rownames %>% column_to_rownames(var="pathway")# Make a column into rownames in this dataframe

GSEA_Haem_matrix <- data.matrix(GSEA_Haem_wide, rownames.force = T) #convert to a matrix for complex heatmap

pheatmap(GSEA_Haem_matrix)  # this is ComplexHeatmap::pheatmap with matrix of choice. 

pheatmap(GSEA_Haem_matrix, 
         color = colorRampPalette(c("#3498db", "white", "#e74c3c"))(50),
         cellwidth = 8,
         cellheight = 8,
         border_color = "white",
         fontsize = 6 
)

######## Working with padj #####

#GSEA_Haem is the original data frame: 
GSEA_Haem_Test <- GSEA_Haem
GSEA_Haem_Test <- subset(GSEA_Haem_Test,padj < adjusted_pval)


GSEA_Haem_P <- GSEA_Haem[, c("pathway", "padj","cluster")]
GSEA_Haem_wide_P <- pivot_wider(GSEA_Haem_P, names_from = "cluster", values_from = "padj") #Reshape data using tidyr package - pivot wider converts from long to wide. 
GSEA_Haem_wide_P <-GSEA_Haem_wide_P %>% remove_rownames %>% column_to_rownames(var="pathway")# Make a column into rownames in this dataframe

GSEA_Haem_matrix_p <- data.matrix(GSEA_Haem_wide_P, rownames.force = T) #convert to a matrix for complex heatmap
#GSEA_Haem_matrix_p[GSEA_Haem_matrix_p[, "three"] == 11,] #WHAT IS THIS ?
#[is.na(GSEA_Haem_matrix_p)] <- 0 #

GSEA_Haem_matrix_p
# this is ComplexHeatmap::pheatmap with matrix of choice. 
tiff("29_01_2022_VEGFAC_fGSEA.tiff", units="in", width=5.5, height=5, res=300)
pheatmap(GSEA_Haem_matrix_p, 
         color = colorRampPalette(c("white","#e74c3c"))(50),
         cellwidth = 8,
         cellheight = 8,
         border_color = "white",
         fontsize = 6 
)
dev.off()

######## PLOTTING HEATMAP based on SingleCellaR Source Code ##### THIS IS WHAT YOU USED IN THE END ######
#GSEA_Haem is the original data frame: 
library(reshape2)
GSEA_Haem_Test <- GSEA_Haem

#takes output dataframe, performs essentially the pivot function to convert to wider with value padj.
FM<- dcast(GSEA_Haem_Test, pathway~cluster,value.var="padj",na.rm=T)
#Push column pathway to rownames
rownames(FM)<-FM$pathway
#Remove column 'pathways' as is now rownames
FM<-FM[-c(1)]
FM[is.na(FM)==T]<-1
add_small_value = 0.0001 #Do we need this? 
FM <- -log10(FM+add_small_value)
#######
# Now filter out only positive NES values 
FM2<- dcast(GSEA_Haem_Test, pathway~cluster,value.var="NES",na.rm=T)
rownames(FM2)<-FM2$pathway
FM2<-FM2[-c(1)]
FM2[is.na(FM2)==T]<-0
FM2[FM2 < 0]<-0
FM2<-FM2[rowSums(FM2) > 0,]
FM<-FM[rownames(FM2),]
FM[FM2 == 0]<-0


###
GSEA_Haem_Test_Matrix <- data.matrix(FM, rownames.force = T) #convert to a matrix for complex heatmap

tiff("29_01_2022_VEGFAC_fGSEA.tiff", units="in", width=5.5, height=5, res=300)
pheatmap(GSEA_Haem_Test_Matrix, 
         color = colorRampPalette(c("white","#e74c3c"))(50),
         cellwidth = 8,
         cellheight = 8,
         border_color = "white",
         fontsize = 6 
)
dev.off()


VEGFAC_stroma <- sub_clusters(
  VEGFAC,
  cluster.type = c("louvain"),
  cluster_ids = c('cl6','cl9','cl10','cl13','cl14')
)
#Not sure about the rest 

process_cells_annotation(VEGFAC_stroma,mitochondiral_genes_start_with="MT-")

plot_cells_annotation(VEGFAC_stroma,type="histogram")

plot_UMIs_vs_Detected_genes(VEGFAC_stroma)

filter_cells_and_genes(VEGFAC_stroma,
                       min_UMIs=500,
                       max_UMIs=50000,
                       min_detected_genes=300,
                       max_detected_genes=6000,
                       max_percent_mito=12,
                       isRemovedDoublets = FALSE,
                       genes_with_expressing_cells = 10)

normalize_UMIs(VEGFAC_stroma,use.scaled.factor = T)

remove_unwanted_confounders(VEGFAC_stroma,residualModelFormulaStr="~UMI_count+percent_mito")


get_variable_genes_by_fitting_GLM_model(VEGFAC_stroma,mean_expr_cutoff = 0.05,disp_zscore_cutoff = 0.05)

tiff("VEGFAC_stroma_Variable_Gene.tiff", units="in", width=10, height=10, res=300)
plot_variable_genes(VEGFAC_stroma)
dev.off()



runPCA(VEGFAC_stroma,use.components=50,use.regressout.data = T)

tiff("VEGFAC_stroma_Elbow_Plot", units="in", width=10, height=10, res=300)
plot_PCA_Elbowplot(VEGFAC_stroma)
dev.off()

#UMAP

runUMAP(VEGFAC_stroma,dim_reduction_method = "pca",n.dims.use = 15,n.neighbors = 150,
        uwot.metric = "euclidean")
identifyClusters(VEGFAC_stroma,n.dims.use = 15,n.neighbors = 150,knn.metric = "euclidean")

#remove_clusters(VEGFAC_stroma,cluster.type = "louvain",cluster_ids = "cl5")

plot_umap_label_by_clusters(VEGFAC_stroma,show_method = "louvain",point.size = 3)

tiff("17_11_2021_VEGFAC_stroma_UMAP.tiff", units="in", width=6, height=5, res=300)
plot_umap_label_by_clusters(VEGFAC_stroma,show_method = "louvain",point.size = 3)
dev.off()

tiff("17_11_2021_VEGFAC_stroma_UMAP_no_label.tiff", units="in", width=6, height=5, res=300)
plot_umap_label_by_clusters(VEGFAC_stroma,show_method = "louvain",point.size = 3, mark.clusters = FALSE)
dev.off()



runFA2_ForceDirectedGraph(VEGFAC_stroma,n.dims.use = 20,
                          n.neighbors = 5,n.seed = 1,fa2_n_iter = 1000)

tiff("17_11_2021_VEGFAC_stroma_Force_directed.tiff", units="in", width=6, height=5, res=300)
plot_forceDirectedGraph_label_by_clusters(VEGFAC_stroma,show_method = "louvain",vertex.size = 3,
                                          background.color = "white")
dev.off()

tiff("17_11_2021_VEGFAC_stroma_Violin_Stroma.tiff", units="in", width=6, height=5, res=300)
plot_violin_for_marker_genes(VEGFAC_stroma,gene_list = c("CXCL12","COL3A1"))
dev.off()

tiff("17_11_2021_VEGFAC_stroma_Violin_Endo.tiff", units="in", width=6, height=5, res=300)
plot_violin_for_marker_genes(VEGFAC_stroma,gene_list = c("CDH5","PECAM1"))
dev.off()

tiff("17_11_2021_VEGFAC_stroma_Violin_PRDGFR.tiff", units="in", width=6, height=5, res=300)
plot_violin_for_marker_genes(VEGFAC_stroma,gene_list = c("PDGFRA","PDGFRB"))
dev.off()

tiff("17_11_2021_VEGFAC_stroma_Feature_Populations_CXCL12.tiff", units="in", width=6, height=5, res=300)
plot_umap_label_by_genes(VEGFAC_stroma,gene_list = c("CXCL12"), point.size = 2)
dev.off()

tiff("17_11_2021_VEGFAC_stroma_Feature_Populations_CDH5.tiff", units="in", width=6, height=5, res=300)
plot_umap_label_by_genes(VEGFAC_stroma,gene_list = c("CDH5"), point.size = 2)
dev.off()

tiff("17_11_2021_VEGFAC_stroma_Feature_Populations_PDGFRB.tiff", units="in", width=6, height=5, res=300)
plot_umap_label_by_genes(VEGFAC_stroma,gene_list = c("PDGFRB"), point.size = 2)
dev.off()

tiff("17_11_2021_VEGFAC_stroma_bubble_2.tiff", units="in", width=4, height=5, res=300)
plot_bubble_for_genes_per_cluster(VEGFAC_stroma,cluster.type = "louvain",
                                  gene_list=c("CXCL12","COL3A1",
                                              "PDGFRA","PDGFRB",
                                              "CSPG4","NES",
                                              "CDH5","PECAM1","FLT4","ENG","EMCN","MCAM","ITGA4"),show.percent = F,buble.scale = 6,point.color1 = "WHITE",
                                  point.color2 = "#3498db",
                                  gene.font.size = 8,
                                  axist.x.font.size = 10,
                                  #IsCustomOrder = TRUE,
                                  #IsCustomOrderByRow = TRUE,
                                  #IsCustomOrderByColumn = TRUE,
                                  #row.custom.order = "",
                                  #column.custom.order = c("cl2","cl1","cl5","cl6","cl3","cl4")
)
dev.off()


findMarkerGenes(
  VEGFAC_stroma,
  cluster.type = c("louvain"),
  min.log2FC = 0.3,
  min.expFraction = 0.25,
)

tiff("17_11_2021_VEGFAC_stroma_HeatMap.tiff", units="in", width=10, height=10, res=300)
plot_heatmap_for_marker_genes(VEGFAC_stroma,cluster.type = "louvain",n.TopGenes = 10,rowFont.size = 5)
dev.off()

#Rename cluster IDs
VEGFAC_stroma@sc.clusters$annotations[VEGFAC_stroma@sc.clusters$louvain_cluster %in% "cl1"] <- "Fibroblast 1"
VEGFAC_stroma@sc.clusters$annotations[VEGFAC_stroma@sc.clusters$louvain_cluster %in% "cl2"] <- "Endothelium"
VEGFAC_stroma@sc.clusters$annotations[VEGFAC_stroma@sc.clusters$louvain_cluster %in% "cl3"] <- "MSC 1"
VEGFAC_stroma@sc.clusters$annotations[VEGFAC_stroma@sc.clusters$louvain_cluster %in% "cl4"] <- "MSC 2"
VEGFAC_stroma@sc.clusters$annotations[VEGFAC_stroma@sc.clusters$louvain_cluster %in% "cl5"] <- "Fibroblast 2"
VEGFAC_stroma@sc.clusters$annotations[VEGFAC_stroma@sc.clusters$louvain_cluster %in% "cl6"] <- "Fibroblast 3"
VEGFAC_stroma@sc.clusters$type <- "Stromal"

tiff("VEGFAC_STROMA_COLOURED_FD_white.tiff", units="in", width=10, height=10, res=300)
plot_forceDirectedGraph_label_by_multiple_gene_sets(VEGFAC_stroma,gmt.file = "/Users/khanas/Desktop/scRNAseq_2021/BM_scRNA_count/SingleCellaR_supp_data/GSEA_Set_Stroma.gmt",
                                                    show_gene_sets = c("Fibroblast","MSC","Endothelium"),
                                                    custom_color = c("red","grey","blue"),
                                                    isNormalizedByHouseKeeping = T,vertex.size = 1.5,edge.size = 0.05,
                                                    background.color = "white")
dev.off()


#VEGFAC GSEA Pre-Ranked Genes
pre_rankedGenes_for_GSEA_VEGFAC_stroma<-identifyGSEAPrerankedGenes_for_all_clusters(VEGFAC_stroma,
                                                                                  cluster.type = "louvain")
#My own curated list
fgsea_Results_stroma_1<-Run_fGSEA_for_multiple_comparisons(GSEAPrerankedGenes_list = pre_rankedGenes_for_GSEA_VEGFAC_stroma,minSize = 10,
                                                         gmt.file = "/Users/khanas/Desktop/scRNAseq_2021/BM_scRNA_count/SingleCellaR_supp_data/GSEA_Set_Stroma.gmt", eps = 0)
tiff("VEGFAC_haem_GSEA_list_stroma.tiff", units="in", width=10, height=10, res=300)
plot_heatmap_for_fGSEA_all_clusters(fgsea_Results_stroma_1,isApplyCutoff = TRUE,
                                    use_pvalues_for_clustering=T,
                                    show_NES_score = T,fontsize_row =7,
                                    adjusted_pval = 0.1,
                                    show_only_NES_positive_score = T,format.digits = 3,
                                    clustering_method = "ward.D",
                                    clustering_distance_rows = "euclidean",
                                    clustering_distance_cols = "euclidean",show_text_for_ns = F)
dev.off()

#trying to show single plot for enodthelium 
plot_violin_for_genes_per_custom_group_of_cells(
  VEGFAC_stroma,
  custom_group_of_cells = list("cl2"),
  gene_list = c("CD34"),
  take_log2 = T,
  xlab.text.size = 5,
  point.size = 0.2,
  point.alpha = 0.1,
  #grid.ncol = 1,
  #grid.nrow = 7
)

#Plot to show stromal support 
#plot_violin_for_marker_genes(VEGFAC_stroma, gene_list = c('SCF','TPO','GCSF','IL3','IL6','IL10','IL11','EPO','GCSF','GMCSF'))

## Colour coded FDG
tiff("VEGFAC_Stroma_COLOURED_FD_white_V4.tiff", units="in", width=10, height=10, res=300)
plot_forceDirectedGraph_label_by_multiple_gene_sets(VEGFAC_stroma,gmt.file = "/Users/khanas/Desktop/scRNAseq_2021/BM_scRNA_count/SingleCellaR_supp_data/GSEA_Set_Stroma_2.gmt",
                                                    show_gene_sets = c("MSC","Fibroblast","Endothelial"),
                                                    custom_color = c("#ff3838","#0096FF","#ffaf40"),
                                                    isNormalizedByHouseKeeping = T,vertex.size = 1.5,edge.size = 0.05,
                                                    background.color = "white") #86bc86", "cl10" = "#79706e", "cl4" = '#ff7f0e',"cl13" = '#bcbd22'
dev.off()


save(VEGFAC_stroma,file="/Users/khanas/Desktop/scRNAseq_2021/BM_scRNA_count/17_11_2021_VEGFAC_Reanalysis/17_11_2021_VEGFAC_Stroma_SingleCellaR.rds")
save(VEGFAC_haem,file="/Users/khanas/Desktop/scRNAseq_2021/BM_scRNA_count/17_11_2021_VEGFAC_Reanalysis/17_11_2021_VEGFAC_Haem_SingleCellaR.rds")


#17th November 2021 - subset analysis of Stromal and Haematopoeitic compartments complete. Now integrate for CellPhoneDB etc. 
library(SingCellaR)
library(Matrix)

##################
VEGFAC_comp <- new("SingCellaR_Int")
VEGFAC_comp@dir_path_SingCellR_object_files<-"/Users/khanas/Desktop/scRNAseq_2021/BM_scRNA_count/17_11_2021_VEGFAC_Reanalysis/"
VEGFAC_comp@SingCellR_object_files=c( "17_11_2021_VEGFAC_Haem_SingleCellaR.rds",
                                "17_11_2021_VEGFAC_Stroma_SingleCellaR.rds"
                                
)

preprocess_integration(VEGFAC_comp)
VEGFAC_comp
table(VEGFAC_comp@meta.data$sampleID) #don't think I need this 

filter_cells_and_genes(VEGFAC_comp,
                       min_UMIs=500,
                       max_UMIs=50000,
                       min_detected_genes=300,
                       max_detected_genes=6000,
                       max_percent_mito=12,
                       isRemovedDoublets = FALSE,
                       genes_with_expressing_cells = 10)
# [1] ""0/9301 cells will be filtered out from the downstream analyses!"

table(VEGFAC_comp@meta.data$sampleID)
table(VEGFAC_comp_stromal_updated@meta.data$sampleID,VEGFAC_comp_stromal_updated@meta.data$source)

#Deleted batch, source etc as we don;t have this. 
#table(VEGFAC_comp@meta.data$sampleID,VEGFAC_comp@meta.data$data_set)

normalize_UMIs(VEGFAC_comp,use.scaled.factor = T)

sc.clusters.stromal <- VEGFAC_stroma@sc.clusters
sc.clusters.hemat <- VEGFAC_haem@sc.clusters

head(sc.clusters.stromal)
head(sc.clusters.hemat)
###########################

head(VEGFAC_comp@meta.data)
dim(sc.clusters.hemat)
dim(sc.clusters.stromal)
###########################

sc.clusters <- rbind(sc.clusters.stromal,sc.clusters.hemat)
dim(sc.clusters)


meta.data.full <- VEGFAC_comp@meta.data
meta.data.full <- merge(meta.data.full,sc.clusters, by = "Cell")
dim(meta.data.full)
VEGFAC_comp@sc.clusters  <- sc.clusters

head(VEGFAC_comp@meta.data)
head(VEGFAC_comp@sc.clusters)

#Prepare for CellPhoneDB
meta_data <- meta.data.full[,c('Cell','annotations','type')]
write.table(meta_data, '/Users/khanas/Desktop/scRNAseq_2021/CellphoneDB_ktplots_II/cellphonedb_meta.txt', sep = '\t',quote=F,row.names=F)
head(meta_data)
#Normalise and sort out counts
norm_counts_II <- get_normalized_umi(VEGFAC_comp) #this will give normalised counts from CellPhoneDB
norm_counts_II <- norm_counts_II[,meta_data$Cell] #merge ,meta.data with 
#head(norm_counts)
norm_counts_df_II <- data.frame(norm_counts_II, check.names = F)
norm_counts_df_II <- cbind('Gene' = rownames(norm_counts_df_II),norm_counts_df_II) #This will add a column with Gene as title from rownames which is gene values

#load gene names: 
genes <- read.table("/Users/khanas/Desktop/scRNAseq_2021/BM_scRNA_count/VEGFAC/filtered_feature_bc_matrix/features.tsv.gz", stringsAsFactors = FALSE) # read in genes... 
gene_names <- genes[['V2']] #create a list of gene names as a string
gene_ensembl <- genes[['V1']] #create a list of ENSEMBL IDs as a string 

norm_counts_df_II$Gene[norm_counts_df_II$Gene == genes$V2] <-genes$V1[norm_counts_df_II$Gene == genes$V2]
head(norm_counts_df_II[35000:35010,1:10])
#Swap out contents of 'Gene'
#norm_counts_df[['Gene']] <- plyr::mapvalues(x =umap.results[['Gene']], from = gene_names, to = gene_ensembl) #plyr::mapvalues is positional
#norm_counts_df$Gene[norm_counts_df$Gene == gene_names] <- gene_ensembl
#don't work
head(norm_counts_df_II[1:10,1:10])

#You need to make the Gene column the row name... 


write.table(norm_counts_df_II, '/Users/khanas/Desktop/scRNAseq_2021/CellphoneDB_ktplots_II/cellphonedb_count.txt', sep = '\t',quote=F, row.names = F)


save(VEGFAC_comp,file="/Users/khanas/Desktop/scRNAseq_2021/BM_scRNA_count/17_11_2021_VEGFAC_Reanalysis/07_02_2021_VEGFAC_Comp_SingleCellaR.rds")


dim(meta_data)
dim(VEGFAC_comp)
### GSEA for combined data set: 



findMarkerGenes(VEGFAC_comp,cluster.type = "louvain")

export_marker_genes_to_table(
  VEGFAC_comp,
  cluster.type = c("louvain"),
  n.TopGenes = 20,
  min.log2FC = 0.5,
  min.expFraction = 0.3,
  write.t

