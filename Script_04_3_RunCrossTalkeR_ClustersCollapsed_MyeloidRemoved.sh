# Load packages
library(CrossTalkeR)
library(igraph) # Required for V function wihtin plot_cci function
library(ggraph) # Required to reproduce vignette figures

# Define input files
    # Group 1
    path.ctr <- "/Users/seanwen/Documents/Abs/CrossTalkeR/Input/"
    file.ctr <- "VEGFA_ClustersCollapsed_MyeloidRemoved.csv"
    
    # Group 2
    path.exp <- "/Users/seanwen/Documents/Abs/CrossTalkeR/Input/"
    file.exp <- "VEGFA-C_ClustersCollapsed_MyeloidRemoved.csv"
    
    # Merge
    paths <- c("VEGFA"=paste(path.ctr, file.ctr, sep=""),
               "VEGFAC"=paste(path.exp, file.exp, sep="")
                )

# Define example gene (ligands)
genes <- c("TGFB1", "CXCL12", "CD44", "ICAM", "VCAM")

# Check if gene list exist in BOTH groups
    # Group 1
    path.ctr <- "/Users/seanwen/Documents/Abs/CrossTalkeR/Input/"
    file.ctr <- "VEGFA_ClustersCollapsed_MyeloidRemoved.csv"
    df <- read.table(paste(path.ctr, file.ctr, sep=""), sep=",", header=TRUE, stringsAsFactors=FALSE)
    
    ligand.list.1 <- unique(df$Ligand)
    
    # Group 2
    path.ctr <- "/Users/seanwen/Documents/Abs/CrossTalkeR/Input/"
    file.ctr <- "VEGFA-C_ClustersCollapsed_MyeloidRemoved.csv"
    df <- read.table(paste(path.ctr, file.ctr, sep=""), sep=",", header=TRUE, stringsAsFactors=FALSE)
    
    ligand.list.2 <- unique(df$Ligand)
    
    # Check
    overlap.g1 <- intersect(genes, ligand.list.1)
    overlap.g2 <- intersect(genes, ligand.list.2)
    overlap.g1.g2 <- intersect(overlap.g1, overlap.g2)
    
    print(paste(paste(overlap.g1, collapse=", "), " ligands found in Group 1", sep=""))
    print(paste(paste(overlap.g2, collapse=", "), " ligands found in Group 1", sep=""))
    print(paste(paste(overlap.g1.g2, collapse=", "), " ligands found BOTH Groups and retained", sep=""))
    
    # Retain ligands found in BOTH groups
    genes <- overlap.g1.g2
    genes
    
# Define output folder
#output <- "/Users/seanwen/Documents/Abs/CrossTalkeR/Output/Example/"
output <- "/Users/seanwen/Documents/Abs/CrossTalkeR/Output/VEGFA-C vs VEGFA/Clusters Collapsed_Myeloid Removed/"

#########################################################################

# Generate report
data <- generate_report(paths,
                        genes,
                        out_path=paste0(output,'/'),
                        threshold=0,
                        out_file = 'vignettes_example.html',
                        output_fmt = "html_document",
                        report = FALSE
                        )
    
#########################################################################

# TGFB1
gene <- "TGFB1"

path <- "/Users/seanwen/Documents/Abs/CrossTalkeR/Output/VEGFA-C vs VEGFA/Clusters Collapsed_Myeloid Removed/Sankey/"
file <- paste(gene, ".pdf", sep="")
pdf(paste(path, file, sep=""), width=10, height=5)

plot_sankey(lrobj_tbl = data@tables$VEGFAC_x_VEGFA,
            target = c(gene),
            ligand_cluster = NULL,
            receptor_cluster = NULL,
            plt_name = gene)

dev.off()

# CXCL12
gene <- "CXCL12"

path <- "/Users/seanwen/Documents/Abs/CrossTalkeR/Output/VEGFA-C vs VEGFA/Clusters Collapsed_Myeloid Removed/Sankey/"
file <- paste(gene, ".pdf", sep="")
pdf(paste(path, file, sep=""), width=10, height=5)

plot_sankey(lrobj_tbl = data@tables$VEGFAC_x_VEGFA,
            target = c(gene),
            ligand_cluster = NULL,
            receptor_cluster = NULL,
            plt_name = gene)

dev.off()

# CD44
gene <- "CD44"

path <- "/Users/seanwen/Documents/Abs/CrossTalkeR/Output/VEGFA-C vs VEGFA/Clusters Collapsed_Myeloid Removed/Sankey/"
file <- paste(gene, ".pdf", sep="")
pdf(paste(path, file, sep=""), width=10, height=5)

plot_sankey(lrobj_tbl = data@tables$VEGFAC_x_VEGFA,
            target = c(gene),
            ligand_cluster = NULL,
            receptor_cluster = NULL,
            plt_name = gene)

dev.off()
