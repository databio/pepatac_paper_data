library(data.table)
library(ggplot2)
library(reshape2)

theme_PEPATAC <- function(base_family = "sans", ...){
    theme_classic(base_family = base_family, base_size = 14, ...) +
        theme(
            axis.line = element_line(size = 0.5),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "transparent"),
            plot.background = element_rect(fill = "transparent", color = NA),
            legend.background = element_rect(fill = "transparent", color = NA),
            legend.box.background = element_rect(fill = "transparent", color = NA),
            aspect.ratio = 1,
            legend.position = "none",
            plot.title = element_text(hjust = 0.5),
            panel.border = element_rect(colour = "black", fill=NA, size=0.5)
        )
}

inter <- fread('consensus_peak_caller_intersections.csv')

svg("peak_caller_jaccard_plot1.svg")
ggplot(data = inter, aes(x=peak_caller_a, y=peak_caller_b,
                         fill=jaccard)) + 
    geom_tile(color = "white") +
    xlab('Peak caller A') +
    ylab('Peak caller B') +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0.5, limit = c(0,1), space = "Lab", 
                         name="Jaccard\nIndex") +
    coord_fixed() +
    geom_text(aes(peak_caller_a, peak_caller_b, label = signif(jaccard,2)),
              color = "black", size = 4) +
    theme_PEPATAC()
dev.off()

cortable <- fread('consensus_peak_caller_intersections_cormat.csv')
cormat <- as.matrix(cortable[,2:7])
rownames(cormat) <- cortable$V1

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
}
upper_tri <- get_upper_tri(cormat)
upper_tri
melted_cormat <- melt(upper_tri, na.rm = TRUE)

ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    xlab('Peak caller A') +
    ylab('Peak caller B') +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0.5, limit = c(0,1), space = "Lab", 
                         name="Jaccard\nIndex") +
    coord_fixed() +
    geom_text(aes(Var2, Var1, label = signif(value,2)),
              color = "black", size = 4) +
    theme_PEPATAC()

reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
}

cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0.5, limit = c(0,1), space = "Lab", 
                         name="Jaccard\nIndex") +
    xlab('Peak caller A') +
    ylab('Peak caller B') +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1))+
    coord_fixed()  +
    geom_text(aes(Var2, Var1, label = signif(value,2)),
              color = "black", size = 4) +
    theme_PEPATAC()
# Print the heatmap
print(ggheatmap)

svg("peak_caller_jaccard_plot2.svg")
ggheatmap
dev.off()

# Look at just a single sample's peaks and compare peak callers
inter <- fread('gold1_peak_caller_intersections.csv')

svg("gold1_peak_caller_jaccard_plot1.svg")
ggplot(data = inter, aes(x=peak_caller_a, y=peak_caller_b,
                         fill=jaccard)) + 
    geom_tile(color = "white") +
    xlab('Peak caller A') +
    ylab('Peak caller B') +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0.5, limit = c(0,1), space = "Lab", 
                         name="Jaccard\nIndex") +
    coord_fixed() +
    geom_text(aes(peak_caller_a, peak_caller_b, label = signif(jaccard,2)),
              color = "black", size = 4) +
    theme_PEPATAC()
dev.off()

cortable <- fread('gold1_peak_caller_intersections_cormat.csv')
cormat <- as.matrix(cortable[,2:7])
rownames(cormat) <- cortable$V1

upper_tri <- get_upper_tri(cormat)
upper_tri

melted_cormat <- melt(upper_tri, na.rm = TRUE)

ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    xlab('Peak caller A') +
    ylab('Peak caller B') +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0.5, limit = c(0,1), space = "Lab", 
                         name="Jaccard\nIndex") +
    coord_fixed() +
    geom_text(aes(Var2, Var1, label = signif(value,2)),
              color = "black", size = 4) +
    theme_PEPATAC()

reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
}

cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0.5, limit = c(0,1), space = "Lab", 
                         name="Jaccard\nIndex") +
    xlab('Peak caller A') +
    ylab('Peak caller B') +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1))+
    coord_fixed()  +
    geom_text(aes(Var2, Var1, label = signif(value,3)),
              color = "black", size = 4) +
    theme_PEPATAC()
ggheatmap
# Print the heatmap
svg("gold1_peak_caller_jaccard_plot2.svg")
ggheatmap
dev.off()
