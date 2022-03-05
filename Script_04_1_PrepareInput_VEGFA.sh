# Load packages
library(data.table)
library(plyr)

# Read count file
path <- "/Users/seanwen/Documents/Abs/CellPhoneDB/Output/VEGFA CellPhoneDB/out/"
file <- "significant_means.txt"
df <- as.data.frame(fread(paste(path, file, sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE))

# Extend each LR interaction by sig. cell pairs
.list <- list()

pb <- txtProgressBar(1, nrow(df), style=3)

for(i in 1:nrow(df)) {

    # Retrieve LR
    df.small <- df[i, ]
    
    # Retrieve sig cell pairs
    df.small. <- df.small[,-c(1:12)]
    df.small. <- as.data.frame(t(df.small.))
    df.small. <- na.omit(df.small.)
    names(df.small.) <- "MeanLR"
    
    if(nrow(df.small.) >= 1) {
    
        # Retrieve L,R clusters
        . <- strsplit(row.names(df.small.), split="|", fixed=TRUE)
        df.small.$Ligand.Cluster <- sapply(., function(x) {x[1]})
        df.small.$Receptor.Cluster <- sapply(., function(x) {x[2]})
        
        # Annotate LR id for tracking
        df.small.$id_cp_interaction <- df.small$id_cp_interaction
        
        # Save into list
        .list[[i]] <- df.small.
    
    }
    
    # Track progress
    setTxtProgressBar(pb, i)

}

. <- do.call(rbind.data.frame, .list)
dim(.) ; dim(df)

# Annotate metadata
df.small <- unique(df[,c("id_cp_interaction", "interacting_pair", "partner_a", "partner_b", "receptor_a", "receptor_b")])
df <- join(., df.small, by="id_cp_interaction", type="left")
nrow(.)==nrow(df)

# Create L-R column
df$Ligand <- NA
df$Receptor <- NA
df$interaction_type <- NA

# Retrieve L,R names
    # Retrieve complex names (with "_") in partner_a
    index.a <- grep("_", df$partner_a, fixed=TRUE)
    index.b <- grep("_", df$partner_b, fixed=TRUE)
    index.a.only <- setdiff(index.a, index.b)
    index.b.only <- setdiff(index.b, index.a)
    index.ab <- intersect(index.a, index.b)
    length(index.a.only) ; length(index.b.only); length(index.ab)
    
    # Fix partner_a
    df.small <- df[index.a, ]
    unique(df.small$partner_a)
    
    interacting_pair.new <- NULL
    
    for(i in 1:nrow(df.small)) {
    
        df.small. <- df.small[i,]
    
        partner_a <- df.small.$partner_a
        partner_a <- gsub("complex:", "", partner_a)
        
        partner_a.new <- gsub("_", " ", partner_a)
        
        interacting_pair <- df.small.$interacting_pair
        interacting_pair.new[i] <- gsub(partner_a, partner_a.new, interacting_pair, fixed=TRUE)
     
    }
    
    df.small$interacting_pair <- interacting_pair.new
    
    df.small.partner_a_fixed <- df.small
    
    # Fix partner_b
    df.small <- df[index.b, ]
    unique(df.small$partner_b)
    
    interacting_pair.new <- NULL
    
    for(i in 1:nrow(df.small)) {
    
        df.small. <- df.small[i,]
    
        partner_b <- df.small.$partner_b
        partner_b <- gsub("complex:", "", partner_b)
        
        partner_b.new <- gsub("_", " ", partner_b)
        
        interacting_pair <- df.small.$interacting_pair
        interacting_pair.new[i] <- gsub(partner_b, partner_b.new, interacting_pair, fixed=TRUE)
     
    }
    
    df.small$interacting_pair <- interacting_pair.new
    
    df.small.partner_b_fixed <- df.small
    
    # Merge
    df.small <- df[-c(index.a, index.b), ]
    nrow(df) == nrow(df.small) + nrow(df.small.partner_a_fixed) + nrow(df.small.partner_b_fixed)
    df <- rbind.data.frame(df.small, df.small.partner_a_fixed, df.small.partner_b_fixed)
    
    # Check that L and R have ONLT 1 each
    . <- strsplit(df$interacting_pair, split="_", fixed=TRUE)
    unique(lapply(., length)) # Should return only 2
    
# Assign LR
    # Split pair
    . <- strsplit(df$interacting_pair, split="_", fixed=TRUE)
    df$interacting_pair_receptor_a <- sapply(., function(x) {x[1]})
    df$interacting_pair_receptor_b <- sapply(., function(x) {x[2]})

    # L-R
    . <- which(is.na(df$Ligand) & is.na(df$Ligand) &
               df$receptor_a==FALSE & df$receptor_b==TRUE
               )
    df$Ligand[.] <- df$interacting_pair_receptor_a[.]
    df$Receptor[.] <- df$interacting_pair_receptor_b[.]
    df$interaction_type[.] <- "L-R"
    
    df.small <- df[which(df$interaction_type=="L-R"), c("Ligand", "Receptor", "Ligand.Cluster", "Receptor.Cluster")]
    nrow(df.small)
    nrow(unique(df.small))
    
    df.small$id <- paste(df.small[,1], df.small[,2], df.small[,3], df.small[,4], sep="_")
    . <- as.data.frame(table(df.small$id))
    . <- .[order(.$Freq, decreasing=TRUE), ]
    head(.)
    
    # L-R (inverted)
    . <- which(is.na(df$Ligand) & is.na(df$Ligand) &
               df$receptor_a==TRUE & df$receptor_b==FALSE
               )
    df$Ligand[.] <- df$interacting_pair_receptor_b[.]
    df$Receptor[.] <- df$interacting_pair_receptor_a[.]
    df$interaction_type[.] <- "L-R (inverted)"
    
    df.small <- df[which(df$interaction_type=="L-R (inverted)"), c("Ligand", "Receptor", "Ligand.Cluster", "Receptor.Cluster")]
    nrow(df.small)
    nrow(unique(df.small))
    
    df.small$id <- paste(df.small[,1], df.small[,2], df.small[,3], df.small[,4], sep="_")
    . <- as.data.frame(table(df.small$id))
    . <- .[order(.$Freq, decreasing=TRUE), ]
    head(.)
    
    # L-L
    . <- which(is.na(df$Ligand) & is.na(df$Ligand) &
               df$receptor_a==FALSE & df$receptor_b==FALSE
               )
    df$Ligand[.] <- df$interacting_pair_receptor_a[.]
    df$Receptor[.] <- df$interacting_pair_receptor_b[.]
    df$interaction_type[.] <- "L-L"
    
    df.small <- df[which(df$interaction_type=="L-L"), c("Ligand", "Receptor", "Ligand.Cluster", "Receptor.Cluster")]
    nrow(df.small)
    nrow(unique(df.small))
    
    df.small$id <- paste(df.small[,1], df.small[,2], df.small[,3], df.small[,4], sep="_")
    . <- as.data.frame(table(df.small$id))
    . <- .[order(.$Freq, decreasing=TRUE), ]
    head(.)
    
    # R-R
    . <- which(is.na(df$Ligand) & is.na(df$Ligand) &
               df$receptor_a==TRUE & df$receptor_b==TRUE
               )
    df$Ligand[.] <- df$interacting_pair_receptor_a[.]
    df$Receptor[.] <- df$interacting_pair_receptor_b[.]
    df$interaction_type[.] <- "R-R"
    
    df.small <- df[which(df$interaction_type=="R-R"), c("Ligand", "Receptor", "Ligand.Cluster", "Receptor.Cluster")]
    nrow(df.small)
    nrow(unique(df.small))
    
    df.small$id <- paste(df.small[,1], df.small[,2], df.small[,3], df.small[,4], sep="_")
    . <- as.data.frame(table(df.small$id))
    . <- .[order(.$Freq, decreasing=TRUE), ]
    head(.)
    
    # Check for missing assignments
    sum(is.na(df$Ligand))
    sum(is.na(df$Receptor))
    sum(is.na(df$interaction_type))
    table(df$interaction_type)

    # Check for unique L-R and cluster pair
    nrow(df)==nrow(unique(df[,c("Ligand", "Receptor", "Ligand.Cluster", "Receptor.Cluster")]))
    
    # Edit space in L-R naming
    #df$Ligand <- gsub(" ", "_", df$Ligand, fixed=TRUE)
    #df$Receptor <- gsub(" ", "_", df$Receptor, fixed=TRUE)
    
# Reorder columns to match input format
col.main <- c("Ligand", "Receptor", "Ligand.Cluster", "Receptor.Cluster", "MeanLR")
col.others <- names(df)[-which(names(df) %in% col.main)]
df <- df[,c(col.main, col.others)]

# Recode commas (to avoid conflict with saving table as csv)
#df$annotation_strategy <- gsub(",", ";", df$annotation_strategy, fixed=TRUE)

# Save file
#path <- "/Users/seanwen/Documents/Abs/CrossTalkeR/Input/"
#file <- "VEGFA.csv"
#write.table(df, paste(path, file, sep=""), sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)

#########################################################################
########################### COLLAPSE CLUSTERS ###########################
#########################################################################

# Subset most relevant columns
cols <- c("Ligand", "Receptor", "Ligand.Cluster", "Receptor.Cluster", "MeanLR")
df <- df[, cols]

# Collapse clusters
    # Remove number suffix
        # Ligand
        df$Ligand.Cluster <- gsub(" [0-9]", "", df$Ligand.Cluster)
        unique(df$Ligand.Cluster)
    
        # Receptor
        df$Receptor.Cluster <- gsub(" [0-9]", "", df$Receptor.Cluster)
        unique(df$Receptor.Cluster)

    # Pool ery
        # Ligand
        index <- which(df$Ligand.Cluster %in% c("Early Erythroid", "Mid Erythroid", "Late Erythroid"))
        df$Ligand.Cluster[index] <- "Erythroid"
        unique(df$Ligand.Cluster)

        # Receptor
        index <- which(df$Receptor.Cluster %in% c("Early Erythroid", "Mid Erythroid", "Late Erythroid"))
        df$Receptor.Cluster[index] <- "Erythroid"
        unique(df$Receptor.Cluster)
    
# Rename remaining clusters
    # Ligand
    table(df$Ligand.Cluster)
    df$Ligand.Cluster[which(df$Ligand.Cluster=="Endothelium")] <- "Endo"
    df$Ligand.Cluster[which(df$Ligand.Cluster=="Erythroid")] <- "Ery"
    df$Ligand.Cluster[which(df$Ligand.Cluster=="Fibroblast")] <- "Fibro"
    df$Ligand.Cluster[which(df$Ligand.Cluster=="Megakaryocyte")] <- "Mk"
    df$Ligand.Cluster[which(df$Ligand.Cluster=="Monocyte")] <- "Mono"
    df$Ligand.Cluster[which(df$Ligand.Cluster=="Myeloid Progenitor")] <- "Mye"
    table(df$Ligand.Cluster)
    
    # Receptor
    table(df$Receptor.Cluster)
    df$Receptor.Cluster[which(df$Receptor.Cluster=="Endothelium")] <- "Endo"
    df$Receptor.Cluster[which(df$Receptor.Cluster=="Erythroid")] <- "Ery"
    df$Receptor.Cluster[which(df$Receptor.Cluster=="Fibroblast")] <- "Fibro"
    df$Receptor.Cluster[which(df$Receptor.Cluster=="Megakaryocyte")] <- "Mk"
    df$Receptor.Cluster[which(df$Receptor.Cluster=="Monocyte")] <- "Mono"
    df$Receptor.Cluster[which(df$Receptor.Cluster=="Myeloid Progenitor")] <- "Mye"
    table(df$Receptor.Cluster)

# Aggregate mean
df$id <- paste(df$Ligand, df$Receptor, df$Ligand.Cluster, df$Receptor.Cluster, sep="|")
df <- df[,c("id", "MeanLR")]
df.collapsed <- aggregate(MeanLR ~ id, data=df, mean)

. <- strsplit(df.collapsed$id, split="|", fixed=TRUE)

df.collapsed <- data.frame("Ligand"=sapply(., function(x) {x[1]}),
                           "Receptor"=sapply(., function(x) {x[2]}),
                           "Ligand.Cluster"=sapply(., function(x) {x[3]}),
                           "Receptor.Cluster"=sapply(., function(x) {x[4]}),
                           "MeanLR"=df.collapsed$MeanLR,
                           stringsAsFactors=FALSE
                           )
                           
nrow(df) ; nrow(df.collapsed)

# Save file
#path <- "/Users/seanwen/Documents/Abs/CrossTalkeR/Input/"
#file <- "VEGFA_ClustersCollapsed.csv"
#write.table(df.collapsed, paste(path, file, sep=""), sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)

#########################################################################
####################### REMOVE MYELOID CLUSTERS #########################
#########################################################################

# Remove myeloid clusters
    # Ligand
    df.collapsed <- df.collapsed[which(df.collapsed$Ligand.Cluster != "Mye"), ]
    
    # Receptor
    df.collapsed <- df.collapsed[which(df.collapsed$Receptor.Cluster != "Mye"), ]

# Save file
path <- "/Users/seanwen/Documents/Abs/CrossTalkeR/Input/"
file <- "VEGFA_ClustersCollapsed_MyeloidRemoved.csv"
write.table(df.collapsed, paste(path, file, sep=""), sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)
