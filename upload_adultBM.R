rm(list = ls())

library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(harmony)
library(symphony)
library(Seurat)
library(Matrix)
library(SingCellaR)
ht_opt$message = FALSE



#############################################################################################

setwd("~/Documents/Symphony/")
suppressPackageStartupMessages({
  source('libs.R') 
  source('utils.R')
})


############### matrix for the reference data ##########

load("~/Documents/SingCellaR_objects/R_objects/ABM_1.SingCellaR.rdata")
load(file = "~/Documents/SingCellaR_objects/R_objects/human_life_v0.2.Supervised-Harmony.SingCellaR.rdata")

plot_umap_label_by_clusters(human_life,show_method = "louvain")

humanlife.umap.results <- human_life@umap.result
humanlife.sc.clusters <- human_life@sc.clusters
humanlife.umap.results <- merge(humanlife.umap.results,humanlife.sc.clusters, by = "Cell")
dim(humanlife.umap.results)

ABM.umap.results <- humanlife.umap.results[humanlife.umap.results$sampleID %in% "1_Adult_BM",]
ref.meta.data <- ABM.umap.results

exprs_norm <- human_life@assays$data$normalized.umi
exprs_norm_ref <- exprs_norm[,ref.meta.data$Cell]

dim(exprs_norm_ref)
dim(ref.meta.data)

identical(as.character(ref.meta.data$Cell),colnames(exprs_norm_ref))
ref.meta.data %>% head(10)


########### build the reference ########

set.seed(0)
reference = symphony::buildReference(
  exprs_norm_ref,
  ref.meta.data,
  K = 50,                   # number of Harmony clusters
  verbose = TRUE,            # verbose output
  do_umap = TRUE,            # can set to FALSE if want to run umap separately later
  do_normalize = FALSE,      # set to TRUE if input counts are not normalized yet
  vargenes_method = 'vst',   # method for variable gene selection ('vst' or 'mvp')
  topn = 2000,               # number of variable genes to choose per group
  d = 40,                    # number of PCs
  save_uwot_path = './ABM_uwot_model_1'
)

reference$normalization_method = 'SingCellaR_Normalization' # optionally save normalization method in custom slot

# Save reference (modify with your desired output path)
saveRDS(reference, './ABM_ref.rds')

ref.meta.data$UMAP1 <- NULL
ref.meta.data$UMAP2 <- NULL

umap_labels = cbind(ref.meta.data, reference$umap$embedding)

head(umap_labels)

umap_labels$annotations[umap_labels$louvain_cluster == "cl1"] <- "HSC/MPP"
umap_labels$annotations[umap_labels$louvain_cluster == "cl2"] <- "Myeloid 1 (neutrophil)"
umap_labels$annotations[umap_labels$louvain_cluster == "cl3"] <- "Erythroid"
umap_labels$annotations[umap_labels$louvain_cluster == "cl4"] <- "Megakaryocyote/Erythroid progenitor"
umap_labels$annotations[umap_labels$louvain_cluster == "cl5"] <- "early myeloid"
umap_labels$annotations[umap_labels$louvain_cluster == "cl6"] <- "early lymphoid"
umap_labels$annotations[umap_labels$louvain_cluster == "cl7"] <- "MPP (myeloid)"
umap_labels$annotations[umap_labels$louvain_cluster == "cl8"] <- "B lymphoid"
umap_labels$annotations[umap_labels$louvain_cluster == "cl9"] <- "Erythroid (cycling)"
umap_labels$annotations[umap_labels$louvain_cluster == "cl10"] <- "B lymphoid (cycling)"
umap_labels$annotations[umap_labels$louvain_cluster == "cl11"] <- "Monocyte/Dendritic precursor (cycling)"
umap_labels$annotations[umap_labels$louvain_cluster == "cl12"] <- "MPP (lymphoid)"
umap_labels$annotations[umap_labels$louvain_cluster == "cl13"] <- "Megakaryocyte/Erythroid (cycling)"
umap_labels$annotations[umap_labels$louvain_cluster == "cl14"] <- "Myeloid 2 (cycling)"
umap_labels$annotations[umap_labels$louvain_cluster == "cl15"] <- "MPP"
umap_labels$annotations[umap_labels$louvain_cluster == "cl16"] <- "Myeloid 2 (monocytic)"
umap_labels$annotations[umap_labels$louvain_cluster == "cl17"] <- "Eosinophil/Basophil/Mast"
umap_labels$annotations[umap_labels$louvain_cluster == "cl18"] <- "Plasmacytoid dendritic precursor"
umap_labels$annotations[umap_labels$louvain_cluster == "cl19"] <- "Endothelial cell"
umap_labels$annotations[umap_labels$louvain_cluster == "cl20"] <- "Plasmacytoid dendritic precursor (cycling)"
umap_labels$annotations[umap_labels$louvain_cluster == "cl21"] <- "Late Erythroid"


# table(umap.results$louvain_cluster,umap.results$louvain_cluster_color)


group.colors <- c("HSC/MPP" = "gray",
                  "B lymphoid (cycling)" = "orange",
                  "Monocyte/Dendritic precursor (cycling)" = "cyan",
                  "MPP (lymphoid)" = "gray",
                  "Megakaryocyte/Erythroid (cycling)" = "purple",
                  "Myeloid 2 (cycling)" = "cyan",
                  "MPP" = "gray",
                  "Myeloid 2 (monocytic)" = "cyan",
                  "Eosinophil/Basophil/Mast" = "deeppink",
                  "Plasmacytoid dendritic precursor" = "orange",
                  "Endothelial cell" = "green",
                  "Myeloid 1 (neutrophil)" = "cyan",
                  "Plasmacytoid dendritic precursor (cycling)" = "orange",
                  "Late Erythroid" = "red",
                  "Erythroid" = "red",
                  "Megakaryocyote/Erythroid progenitor" = "purple",
                  "early myeloid" = "cyan",
                  "early lymphoid" = "orange",
                  "MPP (myeloid)" = "gray",
                  "B lymphoid"= "orange",
                  "Erythroid (cycling)" = "red")

tiff(file = "~/Documents/Abs/VEGFAC/Symphony-mapping_ABM_full.tiff",units = "in",width =6, height = 6, res=300,compression = "lzw")
plotBasic(umap_labels, color.by = 'annotations',show.labels = F,legend.position = "none",color.mapping = group.colors)
dev.off()



################################ build query dataset ######################################


load(file = "~/Documents/Abs/OneDrive_1_20-11-2021/17_11_2021_VEGFAC_Haem_SingleCellaR.rds" )

meta.data <- VEGFAC_haem@sc.clusters
expr <- get_umi_count(VEGFAC_haem)
expr <- expr[,meta.data$Cell]
dim(expr)

identical(colnames(expr),meta.data$Cell)
meta.data$state <- "organoids"

# head(ref.meta.data)

################################ run symphony to map query dataset ######################################
set.seed(2021)
query = mapQuery(expr, 
                 meta.data, 
                 ref_obj = reference, 
                 do_normalize = TRUE
)

set.seed(2021)
query = knnPredict(query, 
                   reference, 
                   reference$meta_data$louvain_cluster, 
                   k = 20,
                   confidence = TRUE,
                   seed = 1
)

save(query,file = "~/Documents/Abs/VEGFAC/query_vegfac_hemat_ABM.rdata")

r_metadata = reference$meta_data
q_metadata = query$meta_data
r_metadata$ref_query = 'reference'
q_metadata$ref_query = 'query'

head(r_metadata)
head(q_metadata)


r_metadata$cell_type_pred_knn <- r_metadata$louvain_cluster
r_metadata$cell_type_pred_knn_prob <- 1

head(r_metadata)
head(q_metadata)

r_metadata <- r_metadata[,c("Cell",'louvain_cluster',"cell_type_pred_knn","cell_type_pred_knn_prob","ref_query")]
q_metadata <- q_metadata[,c("Cell",'annotations',"cell_type_pred_knn","cell_type_pred_knn_prob","ref_query")]

head(r_metadata)
head(q_metadata)

colnames(r_metadata) <- c("Cell",'CellType',"cell_type_pred_knn","cell_type_pred_knn_prob","ref_query")
colnames(q_metadata) <- c("Cell",'CellType',"cell_type_pred_knn","cell_type_pred_knn_prob","ref_query")

meta_data_combined = rbind(q_metadata,r_metadata)
umap_combined = rbind(query$umap,reference$umap$embedding)

umap_combined_labels = cbind(meta_data_combined, umap_combined) %>% 
  mutate(cell_type_pred_knn = fct_relevel(cell_type_pred_knn, group.ordering))



umap_combined_labels %>% head(4)
umap_combined_labels$dataset[umap_combined_labels$ref_query == "reference"] <- "Roy et al"
umap_combined_labels$dataset[umap_combined_labels$ref_query == "query"] <- "Organoid-VEGFAC-hemat"

setwd("~/Documents/Abs/VEGFAC/")
tiff(file = "Symphony-mapping_ABM.tiff",units = "in",width =8, height = 6, res=300,compression = "lzw")
plotBasic(umap_combined_labels,
          title = 'Datasets comparison', 
          color.by = 'dataset',
          color.mapping = c("orange","gray"),
          show.labels = F)
dev.off()

tiff(file = "Symphony-mapping_by_dataset_hemat_202202.tiff",units = "in",width =8, height = 6, res=300,compression = "lzw")
p <- umap_combined_labels %>%
  dplyr::sample_frac(1L) %>% 
  ggplot(aes(x = UMAP1, y = UMAP2, color = dataset)) +
  geom_point(size = 0.3) +
  scale_color_manual(values = c("#69C8ECFF","gray")) +
  theme_void() +
  guides(colour = guide_legend(override.aes = list(size = 5))) 

p
dev.off()

umap_combined_query <- umap_combined_labels[umap_combined_labels$ref_query == "query",]
head(umap_combined_query)
dim(umap_combined_query)
table(umap_combined_query$cell_type_pred_knn)

umap_combined_query$annotations[umap_combined_query$cell_type_pred_knn == "cl1"] <- "HSC/MPP"
umap_combined_query$annotations[umap_combined_query$cell_type_pred_knn == "cl2"] <- "Myeloid 1 (neutrophil)"
umap_combined_query$annotations[umap_combined_query$cell_type_pred_knn == "cl3"] <- "Erythroid"
umap_combined_query$annotations[umap_combined_query$cell_type_pred_knn == "cl4"] <- "Megakaryocyote/Erythroid progenitor"
umap_combined_query$annotations[umap_combined_query$cell_type_pred_knn == "cl5"] <- "early myeloid"
umap_combined_query$annotations[umap_combined_query$cell_type_pred_knn == "cl6"] <- "early lymphoid"
umap_combined_query$annotations[umap_combined_query$cell_type_pred_knn == "cl7"] <- "MPP (myeloid)"
umap_combined_query$annotations[umap_combined_query$cell_type_pred_knn == "cl8"] <- "B lymphoid"
umap_combined_query$annotations[umap_combined_query$cell_type_pred_knn == "cl9"] <- "Erythroid (cycling)"
umap_combined_query$annotations[umap_combined_query$cell_type_pred_knn == "cl10"] <- "B lymphoid (cycling)"
umap_combined_query$annotations[umap_combined_query$cell_type_pred_knn == "cl11"] <- "Monocyte/Dendritic precursor (cycling)"
umap_combined_query$annotations[umap_combined_query$cell_type_pred_knn == "cl12"] <- "MPP (lymphoid)"
umap_combined_query$annotations[umap_combined_query$cell_type_pred_knn == "cl13"] <- "Megakaryocyte/Erythroid (cycling)"
umap_combined_query$annotations[umap_combined_query$cell_type_pred_knn == "cl14"] <- "Myeloid 2 (cycling)"
umap_combined_query$annotations[umap_combined_query$cell_type_pred_knn == "cl15"] <- "MPP"
umap_combined_query$annotations[umap_combined_query$cell_type_pred_knn == "cl16"] <- "Myeloid 2 (monocytic)"
umap_combined_query$annotations[umap_combined_query$cell_type_pred_knn == "cl17"] <- "Eosinophil/Basophil/Mast"
umap_combined_query$annotations[umap_combined_query$cell_type_pred_knn == "cl18"] <- "Plasmacytoid dendritic precursor"
umap_combined_query$annotations[umap_combined_query$cell_type_pred_knn == "cl19"] <- "Endothelial cell"
umap_combined_query$annotations[umap_combined_query$cell_type_pred_knn == "cl20"] <- "Plasmacytoid dendritic precursor (cycling)"
umap_combined_query$annotations[umap_combined_query$cell_type_pred_knn == "cl21"] <- "Late Erythroid"


#########################
tiff(file = "Symphony-mapping_predicted_celltype_hemat_202202.tiff",units = "in",width =8, height = 6, res=300,compression = "lzw")

p <- umap_combined_query %>%
  dplyr::sample_frac(1L) %>% 
  ggplot(aes(x = UMAP1, y = UMAP2, color = annotations)) +
  geom_point(size = 0.3) +
  scale_color_manual(values = group.colors) +
  theme_void() +
  theme(legend.position= "none") +
  guides(colour = guide_legend(override.aes = list(size = 5)))

labels.cent = umap_combined_query %>%
  dplyr::group_by_at("annotations") %>%  # group_by_at takes variable column name
  dplyr::select(UMAP1, UMAP2) %>%
  dplyr::summarize_all(median)

p = p + ggrepel::geom_text_repel(data = labels.cent, aes(x= UMAP1, y = UMAP2, label = get("annotations")),
                                 segment.alpha = 0.5, segment.size = 0.2, box.padding = 0.01, color = 'black')

p
dev.off()


#########################
tiff(file = "Symphony-mapping_original_stromal_v1_NI.tiff",units = "in",width =10, height = 6, res=300,compression = "lzw")
plotBasic(umap_combined_query,
          title = 'Original annotation',
          color.by = 'CellType',
          legend.position = 'right')
dev.off()



VEGFAC_haem@sc.clusters$annotations[VEGFAC_haem@sc.clusters$annotations == "Myeloid Progenitor 1"] <- "Monocyte/Neutrophil Progenitors"
VEGFAC_haem@sc.clusters$annotations[VEGFAC_haem@sc.clusters$annotations == "Myeloid Progenitor 2"] <- "Eosinophil/Basophil/Mast Cell Progenitors"

#########################
group.colors.org = c(   
  'Early Erythroid 1'='#83e3f0',
  'Early Erythroid 2'='#1d6d1f',
  'Endothelium'='#4f8c9d',
  'Fibroblast 1'='#eb1fcb',
  'Fibroblast 2'='#2f5bb1',
  'Fibroblast 3'='#f6932e',
  'HSPC'='#9698dc',
  'Late Erythroid'='#ffb2be',
  'Late Erythroid 1'='#e36f6f',
  'Megakaryocyte 1'='#5ebf72',
  'Megakaryocyte 2'='#fd5917',
  'Mid Erythroid' = '#2f5bb1',
  'Monocyte'='#d59e9a',
  'MSC 1'='#d6061a',
  'MSC 1'='#fd5917',
  'Monocyte/Neutrophil Progenitors'='#caf243',
  'Eosinophil/Basophil/Mast Cell Progenitors'='#6533ed'
  
)

#########################
umap_combined_query$CellType[umap_combined_query$CellType == "Myeloid Progenitor 1"] <- "Monocyte/Neutrophil Progenitors"
umap_combined_query$CellType[umap_combined_query$CellType == "Myeloid Progenitor 2"] <- "Eosinophil/Basophil/Mast Cell Progenitors"


#########################
tiff(file = "Symphony-mapping_original_hemat_202202.tiff",units = "in",width =8, height = 6, res=300,compression = "lzw")
p <- umap_combined_query %>%
  dplyr::sample_frac(1L) %>% 
  ggplot(aes(x = UMAP1, y = UMAP2, color = CellType)) +
  geom_point(size = 0.3) +
  scale_color_manual(values = group.colors.org) +
  theme_void() +
  theme(legend.position= "none") +
  guides(colour = guide_legend(override.aes = list(size = 5)))

labels.cent = umap_combined_query %>%
  dplyr::group_by_at("CellType") %>%  # group_by_at takes variable column name
  dplyr::select(UMAP1, UMAP2) %>%
  dplyr::summarize_all(median)

p = p + ggrepel::geom_text_repel(data = labels.cent, aes(x= UMAP1, y = UMAP2, label = get("CellType")),
                                 segment.alpha = 0.5, segment.size = 0.2, box.padding = 0.01, color = 'black')

p
dev.off()


#########################

tiff(file = "Symphony-mapping_predicted_score_hemat_202202.tiff",units = "in",width =8, height = 6, res=300,compression = "lzw")
ggplot(umap_combined_query,aes(x = UMAP1, y =UMAP2)) +
  geom_point(aes(color = cell_type_pred_knn_prob)) +
  scale_color_viridis_b() +
  theme_void() +
  theme(legend.position = "bottom")

dev.off()

write.table(umap_combined_query,file = "Symphony_mapping_hemat_202002.txt",quote = F,sep = "\t",col.names = T,row.names = F)


#########################
ggplot(umap_combined_query,aes(x =cell_type_pred_knn, y = cell_type_pred_knn_prob)) +
  geom_boxplot(aes(color = cell_type_pred_knn)) +
  # scale_color_viridis_b() +
  theme_bw()


#########################
table(umap_combined_query$CellType,umap_combined_query$cell_type_pred_knn)


cell.type <- table(umap_combined_query$CellType,umap_combined_query$cell_type_pred_knn)
rownames(cell.type)
cell.type
write.table(cell.type,file = "~/Documents/Abs/VEGFAC/original_vs_predicted_stromal_ni.txt",col.names = T,row.names = T,quote = F,sep = "\t")


