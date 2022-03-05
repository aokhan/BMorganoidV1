# Load packages
library(CrossTalkeR)
library(igraph) # Required for V function wihtin plot_cci function
library(ggraph) # Required to reproduce vignette figures

# Define input files
path.exp <- "/Users/seanwen/Documents/Abs/CrossTalkeR/Input/"
file.exp <- "VEGFA-C_Old_ClustersCollapsed.csv"
#file.exp <- "VEGFA-C_Old.csv"

# Define path
paths <- c("VEGFAC"=paste(path.exp, file.exp, sep=""))

# Define example gene (ligands)
genes <- c("TGFB1", "CXCL12", "CD44", "ICAM", "VCAM")

# Check if gene list exist in BOTH groups
    # Group 2
    path.ctr <- "/Users/seanwen/Documents/Abs/CrossTalkeR/Input/"
    file.ctr <- "VEGFA-C_Old_ClustersCollapsed.csv"
    df <- read.table(paste(path.ctr, file.ctr, sep=""), sep=",", header=TRUE, stringsAsFactors=FALSE)
    
    ligand.list <- unique(df$Ligand)
    
    # Check
    overlap <- intersect(genes, ligand.list)

    print(paste(paste(overlap, collapse=", "), " ligands found in Group 1", sep=""))
    
    # Retain ligands found in BOTH groups
    genes <- overlap
    
# Define output folder
output <- "/Users/seanwen/Documents/Abs/CrossTalkeR/Output/VEGFA-C_Old/Clusters Collapsed/"

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
############### CCI INTERACTION NETWORK ANALYSIS: ALL LRs ###############
#########################################################################

# Complete plot
path <- "/Users/seanwen/Documents/Abs/CrossTalkeR/Output/VEGFA-C_Old/Clusters Collapsed/CCI/"
file <- "Cell-Cell Interaction Network VEGFA-C.pdf"
pdf(paste(path, file, sep=""), width=7.5, height=7.5)

plot_cci(graph = data@graphs$VEGFAC,
        colors = data@colors,
        plt_name = "VEGFAC",
        coords = data@coords[V(data@graphs$VEGFAC)$name,],
        emax = NULL,
        leg = FALSE,
        low = 0,
        high = 0,
        ignore_alpha = FALSE,
        log = FALSE,
        efactor = 0.02,
        vfactor = 12,
        node_label_size=0.8,
        pct.interaction=FALSE,
        original.input.file=df,
        result.to.return="network.plot"
        )

dev.off()

# Avoid overlap between labels and nodes: censor labels
path <- "/Users/seanwen/Documents/Abs/CrossTalkeR/Output/VEGFA-C_Old/Clusters Collapsed/CCI/"
file <- "Cell-Cell Interaction Network VEGFA-C_No Labels.pdf"
pdf(paste(path, file, sep=""), width=7.5, height=7.5)

plot_cci(graph = data@graphs$VEGFAC,
        colors = data@colors,
        plt_name = "VEGFAC",
        coords = data@coords[V(data@graphs$VEGFAC)$name,],
        emax = NULL,
        leg = FALSE,
        low = 0,
        high = 0,
        ignore_alpha = FALSE,
        log = FALSE,
        efactor = 0.02,
        vfactor = 12,
        node_label_size=0.01,
        pct.interaction=FALSE,
        original.input.file=df,
        result.to.return="network.plot"
        )

dev.off()

# Return n interaction table (for customising legend)
freq <- plot_cci(graph = data@graphs$VEGFAC,
                 colors = data@colors,
                 plt_name = "VEGFAC",
                 coords = data@coords[V(data@graphs$VEGFAC)$name,],
                 emax = NULL,
                 leg = FALSE,
                 low = 0,
                 high = 0,
                 ignore_alpha = FALSE,
                 log = FALSE,
                 efactor = 0.02,
                 vfactor = 12,
                 node_label_size=0.01,
                 pct.interaction=FALSE,
                 original.input.file=df,
                 result.to.return="n.interaction.table"
                 )

path <- "/Users/seanwen/Documents/Abs/CrossTalkeR/Output/VEGFA-C_Old/Clusters Collapsed/CCI/"
file <- "Number of Interactions_For Creating Legends.txt"
write.table(freq, paste(path, file, sep=""), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
