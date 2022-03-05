
#This is a script for analysis of a single sample BEFORE integration
library(SingCellaR)

data_matrices_dir<-"/Users/khanas/Desktop/scRNAseq_2021/BM_scRNA_count/VEGFA/filtered_feature_bc_matrix/"
VEGFA<-new("SingCellaR")
VEGFA@dir_path_10x_matrix<-data_matrices_dir
VEGFA@sample_uniq_id<-"VEGFA"

load_matrices_from_cellranger(VEGFA,cellranger.version = 6)
#Sparse matrix created run: VEGFA to check should return 'An object of class SingleCellaR...
#Annotate MT genes
process_cells_annotation(VEGFA,mitochondiral_genes_start_with="MT-")

#Produce a QC Tiff
tiff("VEGFA_QC.tiff", units="in", width=10, height=10, res=300)
plot_cells_annotation(VEGFA,type="histogram")
dev.off()


tiff("VEGFA_QC_Boxplot.tiff", units="in", width=10, height=10, res=300)
plot_cells_annotation(VEGFA,type="boxplot")
dev.off()

tiff("VEGFA_QC_UMIvsDet.tiff", units="in", width=10, height=10, res=300)
plot_UMIs_vs_Detected_genes(VEGFA)
dev.off()

filter_cells_and_genes(VEGFA,
                       min_UMIs=500,
                       max_UMIs=50000,
                       min_detected_genes=300,
                       max_detected_genes=6000,
                       max_percent_mito=12,
                       isRemovedDoublets = FALSE,
                       genes_with_expressing_cells = 10)

normalize_UMIs(VEGFA,use.scaled.factor = T)

remove_unwanted_confounders(VEGFA,residualModelFormulaStr="~UMI_count+percent_mito")


get_variable_genes_by_fitting_GLM_model(VEGFA,mean_expr_cutoff = 0.05,disp_zscore_cutoff = 0.05)

tiff("VEGFA_Variable_Gene.tiff", units="in", width=10, height=10, res=300)
plot_variable_genes(VEGFA)
dev.off()

runPCA(VEGFA,use.components=50,use.regressout.data = T)

tiff("VEGFA_Elbow_Plot", units="in", width=10, height=10, res=300)
plot_PCA_Elbowplot(VEGFA)
dev.off()

#UMAP

runUMAP(VEGFA,dim_reduction_method = "pca",n.dims.use = 50,n.neighbors = 30,
        uwot.metric = "euclidean")
plot_umap_label_by_a_feature_of_interest(VEGFA,feature = "UMI_count",point.size = 0.1)
identifyClusters(VEGFA,n.dims.use = 50,n.neighbors = 30,knn.metric = "euclidean")

tiff("VEGFA_UMAP.tiff", units="in", width=10, height=10, res=300)
plot_umap_label_by_clusters(VEGFA,show_method = "louvain",point.size = 0.80)
dev.off()

runFA2_ForceDirectedGraph(VEGFA,n.dims.use = 50,
                          n.neighbors = 5,n.seed = 1,fa2_n_iter = 1000)

tiff("VEGFA_Force_directed.tiff", units="in", width=10, height=10, res=300)
plot_forceDirectedGraph_label_by_clusters(VEGFA,show_method = "louvain",vertex.size = 0.85,
                                          background.color = "white")
dev.off()

tiff("VEGFA_Feature_Test.tiff", units="in", width=10, height=10, res=300)
plot_umap_label_by_genes(VEGFA,gene_list = c("GATA2","CD14","TFRC","PPBP","CDH5","CXCL12","GNLY"))
dev.off()

tiff("VEGFA_UMAP_unlabelled.tiff", units="in", width=10, height=10, res=300)
plot_umap_label_by_clusters(VEGFA,show_method = "louvain",point.size = 0.80, mark.clusters = FALSE)
dev.off()

tiff("VEGFA_Force_directed_unlabelled.tiff", units="in", width=10, height=10, res=300)
plot_forceDirectedGraph_label_by_clusters(VEGFA,show_method = "louvain",vertex.size = 0.85, mark.clusters = FALSE,
                                          background.color = "white")
dev.off()


findMarkerGenes(VEGFA,cluster.type = "louvain")

tiff("VEGFA_HeatMap.tiff", units="in", width=10, height=10, res=300)
plot_heatmap_for_marker_genes(VEGFA,cluster.type = "louvain",n.TopGenes = 10,rowFont.size = 5)
dev.off()

#### Save VEFGFA object ####
save(VEGFA,file="/Users/khanas/Desktop/scRNAseq_2021/BM_scRNA_count/Count_RProject/SingleCellaR/SingleCellaR_BM/Integration_Files/VEGFA_integration.rds")



#VEGFAC GSEA Pre-Ranked Genes
pre_rankedGenes_for_GSEA_VEGFA<-identifyGSEAPrerankedGenes_for_all_clusters(VEGFA,
                                                                      cluster.type = "louvain")
#My own curated list
fgsea_Results_1_VEGFA<-Run_fGSEA_for_multiple_comparisons(GSEAPrerankedGenes_list = pre_rankedGenes_for_GSEA_VEGFA,minSize = 10,
                                                    gmt.file = "/Users/khanas/Desktop/scRNAseq_2021/BM_scRNA_count/SingleCellaR_supp_data/GSEA_Set_1.gmt", eps = 0)
tiff("VEGFA_GSEA_list_AB.tiff", units="in", width=10, height=10, res=300)
plot_heatmap_for_fGSEA_all_clusters(fgsea_Results_1_VEGFA,isApplyCutoff = TRUE,
                                    use_pvalues_for_clustering=T,
                                    show_NES_score = T,fontsize_row =7,
                                    adjusted_pval = 0.05,
                                    show_only_NES_positive_score = T,format.digits = 3,
                                    clustering_method = "ward.D",
                                    clustering_distance_rows = "euclidean",
                                    clustering_distance_cols = "euclidean",show_text_for_ns = F)
dev.off()


