# Define gene names to plot at line 24

# Load packages
library(ggplot2)
library(ggrepel)

# Read R object
path <- "/Users/seanwen/Documents/SingCellaR/Custom Functions/In/"
file <- "DE_cl11_VEGFA_vs_VEGFAC.txt"
df <- read.table(paste(path, file, sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Recode Inf,-Inf values with arbitrary values
df$log2FC[which(df$FoldChange=="Inf")] <- log2(df$ExpA[which(df$FoldChange=="Inf")] + 1)
df$log2FC[which(df$FoldChange==0)] <- log2(df$ExpB[which(df$FoldChange==0)] + 1) * -1

# Indicate sig genes + direction
df$sig.genes.direction <- NA
df$sig.genes.direction[which(df$adjusted.pval < 0.05 & df$log2FC > 0.5)] <- "up"
df$sig.genes.direction[which(df$adjusted.pval < 0.05 & df$log2FC < -0.5)] <- "down"
df$sig.genes.direction[is.na(df$sig.genes.direction)] <- "n.s."
df$sig.genes.direction <- factor(df$sig.genes.direction, levels=c("up", "down", "n.s."))

# Genes to label
#label <- df$Gene[c(1:15)]
label <- c("ANGPT2", "CD34")
intersect(label, df$Gene)
df$label <- ifelse(df$Gene %in% label, df$Gene, "")

######################################################################
############################ VOLCANO PLOT ############################
######################################################################

# Definition
data <- df
x <- data$log2FC
y <- -log10(data$adjusted.pval)
z <- data$sig.genes.direction
#label <- data$gene_name
maintitle <- ""
xtitle <- "log2(FC)"
ytitle <- "-log10(FDR)"

xmin <- floor(min(x)) ; xmax <- ceiling(max(x)) ; xinterval <- 2

xintercept <- c(-0.5, 0.5)
yintercept <- -log10(0.05)
    
sig.up <- which(data$sig.genes.direction=="up")
sig.down <- which(data$sig.genes.direction=="down")

if(length(sig.up) != 0 & length(sig.down) != 0) {

    col.breaks <- c("red", "blue", "gray")
    
} else if(length(sig.up) != 0 & length(sig.down) == 0) {

    col.breaks <- c("red", "gray")
    
} else if(length(sig.up) == 0 & length(sig.down) != 0) {

    col.breaks <- c("blue", "gray")

} else if(length(sig.up) == 0 & length(sig.down) == 0) {

    col.breaks <- "gray"

}

# Plot
plot <- ggplot(data, aes(x=x, y=y)) +
           geom_point(aes(color=z), shape=20, alpha = 0.25) +
           geom_hline(yintercept=yintercept, linetype="dashed", color="black", size=0.25) +
           geom_vline(xintercept=xintercept, linetype="dashed", color="black", size=0.25) +
           geom_text_repel(aes(label=label), max.overlaps = Inf, box.padding = 0.5, size=2.5, max.time = 1, max.iter = 1e5, segment.alpha=0.5, segment.size=0.1) +
           scale_colour_manual(values=col.breaks) +
           scale_x_continuous(breaks=seq(xmin, xmax, by=xinterval), limits=c(xmin, xmax)) +
           #scale_y_continuous(breaks=seq(ymin, ymax, by=yinterval), limits=c(ymin, ymax)) +
           labs(title=maintitle, x=xtitle, y=ytitle) +
            theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border=element_blank(),
                plot.title=element_text(hjust = 0.5, size=15),
                plot.subtitle=element_text(hjust = 0.5, size=15),
                axis.line.y.left = element_line(color="black"),
                axis.line.x = element_line(color="black"),
                axis.title=element_text(size=12),
                axis.text=element_text(size=12),
                axis.text.x=element_text(size=8, colour="black"),
                axis.text.y=element_text(size=8, colour="black"),
                legend.position="none",
                legend.title=element_text(size=8),
                legend.text=element_text(size=8)
                )
                
# Save plot
path <- "/Users/seanwen/Documents/SingCellaR/Custom Functions/Out/"
file <- "DE_cl11_VEGFA_vs_VEGFAC_Volcano Plot_Top Genes Annotated.png"
ggsave(paste(path, file, sep=""), plot, width=3, height=3)
