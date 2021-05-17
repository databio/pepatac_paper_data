library(data.table)
library(ggplot2)
library(stringr)

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

stats <- fread('gold_alt_annotation_stats.tsv')

g1 <- ggplot(stats, aes(x=sample_name, y=TSS_score, col=annotation)) +
    geom_point() +
    theme_PEPATAC() +
    theme(legend.position = "right")

svg("PEPATAC_alt_annotation_TSS-score_1.svg")
g1
dev.off()


g2 <- ggplot(stats, aes(x=annotation, y=TSS_score, col=annotation)) +
    geom_boxplot() +
    theme_PEPATAC() +
    theme(legend.position = "right")
svg("PEPATAC_alt_annotation_TSS-score_2.svg")
g2
dev.off()

g3 <- ggplot(data=stats, aes(x=annotation, y=TSS_score)) +
    geom_violin(width = 0.1, draw_quantiles = c(0.25,0.75),
                linetype="dashed", alpha=0.5) +
    geom_violin(width=0.1, fill="transparent",
                draw_quantiles = 0.5, alpha=0.5) +
    stat_summary(fun = "mean", geom = "point",
                 shape = 1, size = 2) +
    geom_boxplot(width = 0.25, alpha=0.1) +
    geom_point(aes(col=stats$sample_name)) +
    theme_PEPATAC() +
    labs(x="", y="TSS score") +
    theme(legend.position = "right")
svg("PEPATAC_alt_annotation_TSS-score_3.svg")
g3
dev.off()

