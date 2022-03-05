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

iMSC <- readRDS(file = "~/Documents/Datasets/iMSC_Nature immunology/nonhematopoietic_dejong2021.rds")

ref.meta.data <- iMSC@meta.data
table(ref.meta.data$source)
ref.meta.data <- ref.meta.data[ref.meta.data$source == "control",]
dim(ref.meta.data)
head(ref.meta.data)
table(ref.meta.data$state)
table(ref.meta.data$source)


exprs_norm <- iMSC@assays$RNA@data
exprs_norm_ref <- exprs_norm[,rownames(ref.meta.data)]

dim(exprs_norm_ref)
dim(ref.meta.data)


ref.meta.data %>% head(10)


########### build the reference ########

set.seed(0)
reference = symphony::buildReference(
  exprs_norm_ref,
  ref.meta.data,
  vars = c('state'),
  theta = c(1), # variables to integrate over
  K = 100,                   # number of Harmony clusters
  verbose = TRUE,            # verbose output
  do_umap = TRUE,            # can set to FALSE if want to run umap separately later
  do_normalize = FALSE,      # set to TRUE if input counts are not normalized yet
  vargenes_method = 'vst',   # method for variable gene selection ('vst' or 'mvp')
  vargenes_groups = 'state', # metadata column specifying groups for variable gene selection 
  topn = 1000,               # number of variable genes to choose per group
  d = 20,                    # number of PCs
  save_uwot_path = './iMSC_control_uwot_model_1'
)

reference$normalization_method = 'SingCellaR_Normalization' # optionally save normalization method in custom slot

# Save reference (modify with your desired output path)
saveRDS(reference, './iMSC_control_ref.rds')

umap_labels = cbind(ref.meta.data, reference$umap$embedding)
# reference <- readRDS(file = "~/Documents/Supat_Ontogeny/Symphony/mouse10x_GFP_ref.rds")

head(umap_labels)

tiff(file = "~/Documents/Abs/VEGFAC/Symphony-mapping_iMSC_full.tiff",units = "in",width =6, height = 6, res=300,compression = "lzw")
plotBasic(umap_labels, title = 'full', color.by = 'clusternames',show.labels = F)
dev.off()



################################ build query dataset ######################################


load(file = "~/Documents/Abs/OneDrive_1_20-11-2021/17_11_2021_VEGFAC_Stroma_SingleCellaR.rds" )

meta.data <- VEGFAC_stroma@sc.clusters
expr <- get_umi_count(VEGFAC_stroma)
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
                 vars = c('state'),
                 do_normalize = TRUE
)

set.seed(2021)
query = knnPredict(query, 
                   reference, 
                   reference$meta_data$clusternames, 
                   k = 20,
                   confidence = TRUE,
                   seed = 1
)

save(query,file = "~/Documents/Abs/VEGFAC/query_vegfac_stromal_iMSC.rdata")


############################
r_metadata = reference$meta_data
q_metadata = query$meta_data
r_metadata$ref_query = 'reference'
q_metadata$ref_query = 'query'

head(r_metadata)
head(q_metadata)

r_metadata$Cell <- rownames(r_metadata)
r_metadata$cell_type_pred_knn <- r_metadata$clusternames
r_metadata$cell_type_pred_knn_prob <- 1

head(r_metadata)
head(q_metadata)

r_metadata <- r_metadata[,c("Cell",'clusternames','state',"cell_type_pred_knn","cell_type_pred_knn_prob","ref_query")]
q_metadata <- q_metadata[,c("Cell",'annotations','state',"cell_type_pred_knn","cell_type_pred_knn_prob","ref_query")]

head(r_metadata)
head(q_metadata)

colnames(r_metadata) <- c("Cell",'CellType','state',"cell_type_pred_knn","cell_type_pred_knn_prob","ref_query")
colnames(q_metadata) <- c("Cell",'CellType','state',"cell_type_pred_knn","cell_type_pred_knn_prob","ref_query")

meta_data_combined = rbind(q_metadata,r_metadata)
umap_combined = rbind(query$umap,reference$umap$embedding)

umap_combined_labels = cbind(meta_data_combined, umap_combined) %>% 
  mutate(cell_type_pred_knn = fct_relevel(cell_type_pred_knn, group.ordering))

umap_combined_labels %>% head(4)
umap_combined_labels$dataset[umap_combined_labels$ref_query == "reference"] <- "de Jong et al"
umap_combined_labels$dataset[umap_combined_labels$ref_query == "query"] <- "Organoid-VEGFAC-stroma"


############################
setwd("~/Documents/Abs/VEGFAC/")
tiff(file = "Symphony-mapping_stromal_NI.tiff",units = "in",width =8, height = 6, res=300,compression = "lzw")
plotBasic(umap_combined_labels,
          title = 'Datasets comparison', 
          color.by = 'dataset',
          color.mapping = c("gray","orange"),
          show.labels = F)
dev.off()


############################
tiff(file = "Symphony-mapping_stromal_NI_202202.tiff",units = "in",width =8, height = 6, res=300,compression = "lzw")

p <- umap_combined_labels %>%
  dplyr::sample_frac(1L) %>% 
  ggplot(aes(x = UMAP1, y = UMAP2, color = dataset)) +
  geom_point(size = 0.3) +
  scale_color_manual(values = c("gray","orange")) +
  theme_void() +
  guides(colour = guide_legend(override.aes = list(size = 5))) 

p
dev.off()

############################
tiff(file = "Symphony-mapping_by_dataset_stromal_NI.tiff",units = "in",width =15, height = 6, res=300,compression = "lzw")
plotBasic(umap_combined_labels,
          title = 'Datasets comparison', 
          color.by = 'dataset',
          facet.by = 'dataset',
          color.mapping = c("gray","orange"),
          show.labels = F)
dev.off()

umap_combined_query <- umap_combined_labels[umap_combined_labels$ref_query == "query",]
head(umap_combined_query)
dim(umap_combined_query)



tiff(file = "Symphony-mapping_predicted_celltype_stromal_NI_202202.tiff",units = "in",width =8, height = 6, res=300,compression = "lzw")

p <- umap_combined_query %>%
  dplyr::sample_frac(1L) %>% 
  ggplot(aes(x = UMAP1, y = UMAP2, color = cell_type_pred_knn)) +
  geom_point(size = 0.3) +
  # scale_color_manual(values = group.colors) +
  theme_void() +
  theme(legend.position= "none") +
  guides(colour = guide_legend(override.aes = list(size = 5)))

labels.cent = umap_combined_query %>%
  dplyr::group_by_at("cell_type_pred_knn") %>%  # group_by_at takes variable column name
  dplyr::select(UMAP1, UMAP2) %>%
  dplyr::summarize_all(median)

p = p + ggrepel::geom_text_repel(data = labels.cent, aes(x= UMAP1, y = UMAP2, label = get("cell_type_pred_knn")),
                                 segment.alpha = 0.5, segment.size = 0.2, box.padding = 0.01, color = 'black')

p
dev.off()

############################
table(VEGFAC_stroma@sc.clusters$annotations)



############################

tiff(file = "Symphony-mapping_original_stromal_v1_NI_202202.tiff",units = "in",width =8, height = 6, res=300,compression = "lzw")
 <- umap_combined_query %>%
  dplyr::sample_frac(1L) %>% 
  ggplot(aes(x = UMAP1, y = UMAP2, color = CellType)) +
  geom_point(size = 0.3) +
  # scale_color_manual(values = group.colors.org) +
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


############################
tiff(file = "Symphony-mapping_original_stromal_v2_NI.tiff",units = "in",width =10, height = 6, res=300,compression = "lzw")
plotBasic(umap_combined_query,
          title = 'Original annotation',
          color.by = 'CellType',show.labels = F,
          legend.position = 'right')
dev.off()

tiff(file = "Symphony-mapping_predicted_score_stromal_NI_202202.tiff",units = "in",width =8, height = 6, res=300,compression = "lzw")
ggplot(umap_combined_query,aes(x = UMAP1, y =UMAP2)) +
  geom_point(aes(color = cell_type_pred_knn_prob)) +
  scale_color_viridis_b() +
  theme_void() +
  theme(legend.position = "bottom")
dev.off()


write.table(umap_combined_query,file = "Symphony_mapping_stromal_NI_202002.txt",quote = F,sep = "\t",col.names = T,row.names = F)


############################
ggplot(umap_combined_query,aes(x =cell_type_pred_knn, y = cell_type_pred_knn_prob)) +
  geom_boxplot(aes(color = cell_type_pred_knn)) +
  # scale_color_viridis_b() +
  theme_bw()



############################
table(umap_combined_query$CellType,umap_combined_query$cell_type_pred_knn)


cell.type <- table(umap_combined_query$CellType,umap_combined_query$cell_type_pred_knn)
rownames(cell.type)
cell.type
write.table(cell.type,file = "~/Documents/Abs/VEGFAC/original_vs_predicted_stromal_ni.txt",col.names = T,row.names = T,quote = F,sep = "\t")


