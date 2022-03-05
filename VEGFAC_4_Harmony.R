 
#This is a script for analysis of a single sample BEFORE integration
library(SingCellaR)

data_matrices_dir<-"/Users/khanas/Desktop/scRNAseq_2021/BM_scRNA_count/VEGFAC/filtered_feature_bc_matrix/"
VEGFAC<-new("SingCellaR")
VEGFAC@dir_path_10x_matrix<-data_matrices_dir
VEGFAC@sample_uniq_id<-"VEGFAC"

load_matrices_from_cellranger(VEGFAC,cellranger.version = 6)
#Sparse matrix created run: VEGFAC to check should return 'An object of class SingleCellaR...
#Annotate MT genes
process_cells_annotation(VEGFAC,mitochondiral_genes_start_with="MT-")

#Produce a QC Tiff
tiff("VEGFAC_QC.tiff", units="in", width=10, height=10, res=300)
plot_cells_annotation(VEGFAC,type="histogram")
dev.off()


tiff("VEGFAC_QC_Boxplot.tiff", units="in", width=10, height=10, res=300)
plot_cells_annotation(VEGFAC,type="boxplot")
dev.off()

tiff("VEGFAC_QC_UMIvsDet.tiff", units="in", width=10, height=10, res=300)
plot_UMIs_vs_Detected_genes(VEGFAC)
dev.off()

#First version of this: 
#filter_cells_and_genes(VEGFAC,
                      # min_UMIs=500,
                      # max_UMIs=50000,
                      # min_detected_genes=300,
                      # max_detected_genes=6000,
                      # max_percent_mito=12,
                      # isRemovedDoublets = FALSE,
                      # genes_with_expressing_cells = 10)

filter_cells_and_genes(VEGFAC,
                       min_UMIs=500,
                       max_UMIs=50000,
                       min_detected_genes=300,
                       max_detected_genes=6000,
                       max_percent_mito=12,
                       isRemovedDoublets = FALSE,
                       genes_with_expressing_cells = 10)

normalize_UMIs(VEGFAC,use.scaled.factor = T)

remove_unwanted_confounders(VEGFAC,residualModelFormulaStr="~UMI_count+percent_mito")


get_variable_genes_by_fitting_GLM_model(VEGFAC,mean_expr_cutoff = 0.05,disp_zscore_cutoff = 0.05)

tiff("VEGFAC_Variable_Gene.tiff", units="in", width=10, height=10, res=300)
plot_variable_genes(VEGFAC)
dev.off()



runPCA(VEGFAC,use.components=50,use.regressout.data = T)

tiff("VEGFAC_Elbow_Plot", units="in", width=10, height=10, res=300)
plot_PCA_Elbowplot(VEGFAC)
dev.off()

#UMAP

runUMAP(VEGFAC,dim_reduction_method = "pca",n.dims.use = 50,n.neighbors = 30,
        uwot.metric = "euclidean")
#plot_umap_label_by_a_feature_of_interest(VEGFAC,feature = "UMI_count",point.size = 0.1)
identifyClusters(VEGFAC,n.dims.use = 50,n.neighbors = 30,knn.metric = "euclidean")

#tiff("VEGFAC_UMAP.tiff", units="in", width=10, height=10, res=300)
plot_umap_label_by_clusters(VEGFAC,show_method = "louvain",point.size = 3, mark.clusters =  TRUE)
#dev.off()



runFA2_ForceDirectedGraph(VEGFAC,n.dims.use = 50,
                          n.neighbors = 5,n.seed = 1,fa2_n_iter = 1000)

tiff("VEGFAC_Force_directed.tiff", units="in", width=10, height=10, res=300)
plot_forceDirectedGraph_label_by_clusters(VEGFAC,show_method = "louvain",vertex.size = 0.85,
                                          background.color = "white")
dev.off()

tiff("VEGFAC_Feature_Test.tiff", units="in", width=10, height=10, res=300)
plot_umap_label_by_genes(VEGFAC,gene_list = c("GATA2","CD14","TFRC","PPBP","CDH5","CXCL12","GNLY"))
dev.off()

tiff("VEGFAC_UMAP_unlabelled.tiff", units="in", width=10, height=10, res=300)
plot_umap_label_by_clusters(VEGFAC,show_method = "louvain",point.size = 3, mark.clusters = FALSE)
dev.off()

plot_umap_label_by_clusters(VEGFAC,show_method = "louvain",point.size = 3, mark.clusters = TRUE)

tiff("VEGFAC_Force_directed_unlabelled.tiff", units="in", width=10, height=10, res=300)
plot_forceDirectedGraph_label_by_clusters(VEGFAC,show_method = "louvain",vertex.size = 0.85, mark.clusters = FALSE,
                                          background.color = "white")
dev.off()


#Plot according to GSEA cilourise populations

#tiff("VEGFAC_UMAP_Gene_Coloured.tiff", units="in", width=10, height=10, res=300)
#plot_umap_label_by_multiple_gene_sets(VEGFAC,gmt.file = "/Users/khanas/Desktop/scRNAseq_2021/BM_scRNA_count/SingleCellaR_supp_data/human.signature.genes.v2.gmt",
#                                      show_gene_sets = c("Erythroid","Myeloid","Megakaryocyte","Monocyte","HSPC_MPP","Endothelium", "Fibroblast","MSC"),
#                                      custom_color = c("red","orange","cyan","purple","green","grey","yellow","white"),
#                                      isNormalizedByHouseKeeping = T,point.size = 1,background.color = "black")
#dev.off()


findMarkerGenes(
  VEGFAC,
  cluster.type = c("louvain"),
  min.log2FC = 0.3,
  min.expFraction = 0.25,
)

export_marker_genes_to_table(
  VEGFAC,
  cluster.type = c("louvain"),
  n.TopGenes = 50,
  min.log2FC = 0.3,
  min.expFraction = 0.25,
  write.to.file = "/Users/khanas/Desktop/scRNAseq_2021/BM_scRNA_count/Count_RProject/SingleCellaR/VEGFAC_MarkerGenes.csv"
)

tiff("VEGFAC_HeatMap.tiff", units="in", width=10, height=10, res=300)
plot_heatmap_for_marker_genes(VEGFAC,cluster.type = "louvain",n.TopGenes = 10,rowFont.size = 5)
dev.off()


#### SAVE OBJECT ####
save(VEGFAC,file="/Users/khanas/Desktop/scRNAseq_2021/BM_scRNA_count/Count_RProject/SingleCellaR/SingleCellaR_BM/VEGFAC_SingleCellaR.rds")









#VEGFAC GSEA Pre-Ranked Genes
pre_rankedGenes_for_GSEA<-identifyGSEAPrerankedGenes_for_all_clusters(VEGFAC,
                                                                      cluster.type = "louvain")
#My own curated list
fgsea_Results_1<-Run_fGSEA_for_multiple_comparisons(GSEAPrerankedGenes_list = pre_rankedGenes_for_GSEA,minSize = 10,
                                                    gmt.file = "/Users/khanas/Desktop/scRNAseq_2021/BM_scRNA_count/SingleCellaR_supp_data/GSEA_Set_1.gmt", eps = 0)
tiff("VEGFAC_GSEA_list_AB.tiff", units="in", width=10, height=10, res=300)
plot_heatmap_for_fGSEA_all_clusters(fgsea_Results_1,isApplyCutoff = TRUE,
                                    use_pvalues_for_clustering=T,
                                    show_NES_score = T,fontsize_row =7,
                                    adjusted_pval = 0.05,
                                    show_only_NES_positive_score = T,format.digits = 3,
                                    clustering_method = "ward.D",
                                    clustering_distance_rows = "euclidean",
                                    clustering_distance_cols = "euclidean",show_text_for_ns = F)
dev.off()




