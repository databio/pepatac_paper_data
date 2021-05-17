# install.packages("devtools")
# devtools::install_github("pepkit/pepr")
# devtools::install_github("databio/GenomicDistributions")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("BSgenome")
# BiocManager::install("GenomicFeatures")
# BiocManager::install("ensembldb")
# BiocManager::install("AnnotationHub")
# BiocManager::install("ExperimentHub")
# install.packages("http://big.databio.org/GenomicDistributionsData/GenomicDistributionsData_0.0.1.tar.gz", repos=NULL)
# devtools::install_github("databio/peppro", subdir="PEPPROr")
library(PEPPROr)
library(ggpubr)
library(ggrepel)
theme_PEPPRO <- function(base_family = "sans", ...){
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
            rect = element_rect(fill = "transparent"),
            aspect.ratio = 1,
            legend.position = "none",
            plot.title = element_text(hjust = 0.5),
            panel.border = element_rect(colour = "black", fill=NA, size=0.5)
        )
}
setwd("C:/Users/srroj/Downloads/papers/PEPPRO_reviews/")

################################################################################
# DM6
dm6 <- fread('PEPPRO_dm6_stats_summary.tsv')

# PI
q <- ggplot(dm6, aes(x=Pause_index, y=sample_name)) +
    geom_bar(stat="identity", position = position_dodge()) +
    theme_PEPPRO() +
    #coord_cartesian(xlim=c(0,100), breaks=c(0,20,40,60,80)) +
    scale_x_continuous(limits = c(0, 90), breaks=c(0,20,40,60,80)) +
    labs(x="Median pause index", y='') +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
    geom_vline(xintercept=10, color='#333333ff', linetype='dashed')

basename <- "dm6"
filename <- paste0(basename, "_pause_index.svg")
svg(filename, height=7.5, width=7.5, bg = 'transparent')
q
dev.off()


################################################################################
# Annotation comparison plots
cmp <- fread("alt_annotations_comparison.tsv")
cmp$sample_name <- factor(cmp$sample_name,
                          levels = rev(c('K562_PRO-seq_100',
                                         'K562_RNA-seq_0',
                                         'K562_RNA-seq_10',
                                         'K562_RNA-seq_20',
                                         'K562_RNA-seq_30',
                                         'K562_RNA-seq_40',
                                         'K562_RNA-seq_50',
                                         'K562_RNA-seq_60',
                                         'K562_RNA-seq_70',
                                         'K562_RNA-seq_80',
                                         'K562_RNA-seq_90',
                                         'K562_RNA-seq_100',
                                         'K562_GRO-seq',
                                         'HelaS3_GRO-seq',
                                         'Jurkat_ChRO-seq_1',
                                         'Jurkat_ChRO-seq_2',
                                         'HEK_PRO-seq',
                                         'HEK_ARF_PRO-seq',
                                         'H9_PRO-seq_1',
                                         'H9_PRO-seq_2',
                                         'H9_PRO-seq_3',
                                         'H9_treated_PRO-seq_1',
                                         'H9_treated_PRO-seq_2',
                                         'H9_treated_PRO-seq_3')))

# mRNA
p <- ggplot(cmp, aes(x=mRNA_contamination, y=sample_name, fill=source)) +
    geom_bar(stat="identity", position = position_dodge()) +
    scale_fill_manual(values = c("red", "darkgreen", "blue")) +
    theme_PEPPRO() +
    coord_cartesian(xlim=c(0,15)) +
    labs(x="Median mRNA contamination", y='') +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
    annotate("rect", ymin=0, ymax=Inf, xmin=1.0, xmax=1.8, alpha=0.5) +
    annotate(geom = 'text', x = 13, y = 24, 
             label = "RefSeq", color = 'blue', hjust=0) +
    annotate(geom = 'text', x = 13, y = 23, 
             label = "GENCODE", color = 'darkgreen', hjust=0) +
    annotate(geom = 'text', x = 13, y = 22,
             label = "Ensembl", color = 'red', hjust=0)
basename <- "annotation_comparison"
filename <- paste0(basename, "_mRNA_contamination.svg")
svg(filename, height=7.5, width=7.5, bg = 'transparent')
p
dev.off()

# PI
q <- ggplot(cmp, aes(x=pause_index, y=sample_name, fill=source)) +
    geom_bar(stat="identity", position = position_dodge()) +
    scale_fill_manual(values = c("red", "darkgreen", "blue")) +
    theme_PEPPRO() +
    scale_x_continuous(limits = c(0, 90), breaks=c(0,20,40,60,80)) +
    labs(x="Median pause index", y='') +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
    geom_vline(xintercept=10, color='#333333ff', linetype='dashed') + 
    annotate(geom = 'text', x = 78, y = 24, 
             label = "RefSeq", color = 'blue', hjust=0) +
    annotate(geom = 'text', x = 78, y = 23, 
             label = "GENCODE", color = 'darkgreen', hjust=0) +
    annotate(geom = 'text', x = 78, y = 22,
             label = "Ensembl", color = 'red', hjust=0)

basename <- "annotation_comparison"
filename <- paste0(basename, "_pause_index.svg")
svg(filename, height=7.5, width=7.5, bg = 'transparent')
q
dev.off()

################################################################################

stats <- fread("PEPPRO_stats_summary_2021.tsv")

stats$sample_name <- factor(stats$sample_name,
                            levels = rev(c('K562_PRO-seq_02',
                                           'K562_PRO-seq_04',
                                           'K562_PRO-seq_06',
                                           'K562_PRO-seq_08',
                                           'K562_PRO-seq_10',
                                           'K562_PRO-seq_20',
                                           'K562_PRO-seq_30',
                                           'K562_PRO-seq_40',
                                           'K562_PRO-seq_50',
                                           'K562_PRO-seq_60',
                                           'K562_PRO-seq_70',
                                           'K562_PRO-seq_80',
                                           'K562_PRO-seq_90',
                                           'K562_PRO-seq_100',
                                           'K562_RNA-seq_0',
                                           'K562_RNA-seq_10',
                                           'K562_RNA-seq_20',
                                           'K562_RNA-seq_30',
                                           'K562_RNA-seq_40',
                                           'K562_RNA-seq_50',
                                           'K562_RNA-seq_60',
                                           'K562_RNA-seq_70',
                                           'K562_RNA-seq_80',
                                           'K562_RNA-seq_90',
                                           'K562_RNA-seq_100',
                                           'K562_GRO-seq',
                                           'HelaS3_GRO-seq',
                                           'Jurkat_ChRO-seq_1',
                                           'Jurkat_ChRO-seq_2',
                                           'HEK_PRO-seq',
                                           'HEK_ARF_PRO-seq',
                                           'H9_PRO-seq_1',
                                           'H9_PRO-seq_2',
                                           'H9_PRO-seq_3',
                                           'H9_treated_PRO-seq_1',
                                           'H9_treated_PRO-seq_2',
                                           'H9_treated_PRO-seq_3',
                                           'H9_PRO-seq_10',
                                           'H9_PRO-seq_20',
                                           'H9_PRO-seq_30',
                                           'H9_PRO-seq_40',
                                           'H9_PRO-seq_50',
                                           'H9_PRO-seq_60',
                                           'H9_PRO-seq_70',
                                           'H9_PRO-seq_80',
                                           'H9_PRO-seq_90',
                                           'H9_PRO-seq_100')))

################################################################################
# rDNA ratios
stats$rDNA_ratio_raw <- stats$Aligned_reads_human_rDNA/stats$Raw_read

# Correlate rDNA ratio to mRNA contamination

cor(stats$rDNA_ratio_raw, stats$mRNA_contamination,
    method = "pearson", use = "complete.obs")
cor.test(stats$rDNA_ratio_raw, stats$mRNA_contamination,
         method = "pearson", use = "complete.obs")
# Drop synthetic spike-in and subsetted samples
core_samples <- stats[-c(1:13,15:25,38:47),]
corplot1 <- ggscatter(core_samples,
          x = "rDNA_ratio_raw", y = "mRNA_contamination",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          cor.coef.coord = c(0.05, 3.27),
          xlab = expression(frac("rDNA aligned reads", "Total reads")),
          ylab = "Median mRNA contamination") +
    geom_text_repel(data=core_samples[3,],
                    aes(rDNA_ratio_raw,mRNA_contamination,label=sample_name),
                    nudge_y = 0.025, nudge_x = -0.025) +
    theme_PEPPRO() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=2))
filename = "Core-Samples_wHelaS3_rDNA-to-mRNA.svg"
svg(filename, height=7.5, width=7.5, bg = 'transparent')
corplot1
dev.off()    

# Drop spike-in, subsets, AND HelaS3 outlier
core_no_outlier <- core_samples[-3,]
corplot2 <- ggscatter(core_no_outlier,
                     x = "rDNA_ratio_raw", y = "mRNA_contamination",
                     add = "reg.line", conf.int = TRUE, 
                     cor.coef = TRUE, cor.method = "pearson",
                     cor.coef.coord = c(0.05, 1.45),
                     xlab = expression(frac("rDNA aligned reads", "Total reads")),
                     ylab = "Median mRNA contamination") +
    theme_PEPPRO() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=2))
filename = "Core-Samples_rDNA-to-mRNA.svg"
svg(filename, height=7.5, width=7.5, bg = 'transparent')
corplot2
dev.off()

# Add additional samples to increase power (from: GSE126919)
hek_samples <- fread('guertin_hek_stats_summary.tsv')
hek_samples$rDNA_ratio_raw <- hek_samples$Aligned_reads_human_rDNA/hek_samples$Raw_read

# remove the two samples that are already part of the original core samples
hek_samples <- hek_samples[-c(20,24),]
# label clones
hek_samples$group = c(rep("GSE126919_A",8), rep("GSE126919_B", 4),
                      rep("GSE126919_C",7), rep("GSE126919_B", 3))

# Combine the core_samples (no outlier) with the HEK samples
combined_samples <- data.table(
    sample_name = c(core_no_outlier$sample_name,
                    hek_samples$sample_name),
    rDNA_ratio_raw = c(core_no_outlier$rDNA_ratio_raw,
                       hek_samples$rDNA_ratio_raw),
    mRNA_contamination = c(core_no_outlier$mRNA_contamination,
                           hek_samples$mRNA_contamination),
    group = c(rep("core_samples", nrow(core_no_outlier)), hek_samples$group))
cor.test(combined_samples$rDNA_ratio_raw,
         combined_samples$mRNA_contamination,
         method = "pearson", use = "complete.obs")
corplot3 <- ggscatter(combined_samples, 
                      x = "rDNA_ratio_raw", y = "mRNA_contamination",
                      add = "reg.line", conf.int = TRUE, 
                      cor.coef = TRUE, cor.method = "pearson",
                      cor.coef.coord = c(0.05, 1.55),
                      xlab = expression(frac("rDNA aligned reads", "Total reads")),
                      ylab = "Median mRNA contamination") +
    geom_point(data=combined_samples, aes(x=rDNA_ratio_raw, 
                                          y=mRNA_contamination,
                                          color=group)) +
    theme_PEPPRO() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=2),
          legend.position = "right", legend.title = element_blank())
filename = "Core-Samples_HEK-Samples_rDNA-to-mRNA.svg"
svg(filename, height=7.5, width=7.5, bg = 'transparent')
corplot3
dev.off()

################################################################################
# Overall plots
ExpectedUniquePlot <- ggplot(core_samples[core_samples$umi_len > 0,],
                             aes(x=sample_name, y=Frac_exp_unique_at_10M)) +
    geom_bar(stat="identity",
             fill = c("#f8766dff", "#d89000ff", "#a3a500ee", "#39b600cf",
                      "#00bf7df7", "#00bfc4ff", "#00b0f6ff", "#9590ffff",
                      "#e76bf3ff", "#ff62bcff")) +
    theme_PEPPRO()
p <- ExpectedUniquePlot + labs(y="Fraction expected unique at 10M reads", x='') +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
filename = "Core-samples_frac-exp-unique-at-10M.svg"
svg(filename, height=7.5, width=7.5, bg = 'transparent')
p + geom_hline(yintercept=0.75, color='#333333ff')
dev.off()

# Add back the RNA-seq spike-in samples for this plot
# Drop subsetted samples
nosubset_samples <- stats[-c(1:13,38:47),]
TooShortPlot <- ggplot(nosubset_samples,
                       aes(x=Pct_uninformative_adapter_reads, y=sample_name)) +
    geom_bar(stat="identity") +
    theme_PEPPRO()
p <- TooShortPlot + labs(x="Reads too short (%)", y='') +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
filename = "H9_subset_pct_uninformative_adapter_reads.svg"
svg(filename, height=7.5, width=7.5, bg = 'transparent')
p + geom_vline(xintercept=5, color='#333333ff', linetype='dashed') +
    coord_cartesian(xlim=c(0,50)) +
    xlab('Uninformative adapter reads (%)') +
    annotate("rect", ymin=0, ymax=Inf, xmin=20, xmax=50, alpha=0.5)
dev.off()

TSSplusplot <- ggplot(nosubset_samples, aes(x=TSS_coding_score, y=sample_name)) +
    geom_bar(stat="identity", fill=rgb(0,139,69, max=255)) +
    theme_PEPPRO()
tssPlus <- TSSplusplot + labs(x="TSS enrichment score \n(+ strand)", y='') +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

TSSminusplot <- ggplot(nosubset_samples, aes(x=sample_name, y=`TSS_non-coding_score`)) +
    geom_bar(stat="identity", fill=rgb(0,0,255, max=255)) +
    theme_PEPPRO()
tssMinus <- TSSminusplot + labs(y="TSS enrichment score \n(- strand)", x='') +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

# Combined TSS plot
tssPlot <- tssPlus + geom_bar(aes(x=`TSS_non-coding_score`, y=sample_name),
                              stat="identity", fill=rgb(0,0,255, max=255)) +
    labs(x="TSS enrichment score", y='') +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
    annotate("text", x=60, y=Inf, hjust=0, vjust=2,
             size=theme_get()$text[["size"]]/3,
             label="coding", angle=0, col=rgb(0,139,69, max=255)) +
    annotate("text", x=60, y=Inf, hjust=0, vjust=4,
             size=theme_get()$text[["size"]]/3,
             label="non-coding", angle=0, col=rgb(0,0,255, max=255))
filename = "Core-Samples_RNA-spike_TSS-scores.svg"
svg(filename, height=7.5, width=7.5, bg = 'transparent')
tssPlot + 
    #geom_vline(xintercept=10, color='darkgreen') +
    #geom_vline(xintercept=5, color='darkblue') +
    coord_cartesian(xlim=c(0,70))
dev.off()

Alignment_rate_human_rDNAplot <- ggplot(nosubset_samples,
                                        aes(x=Alignment_rate_human_rDNA,
                                            y=sample_name)) +
    geom_bar(stat="identity") +
    theme_PEPPRO()
p <- Alignment_rate_human_rDNAplot +
    coord_cartesian(xlim=c(0,60)) +
    labs(x="Alignment rate rDNA", y='') +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
filename = "Core-Samples_RNA-spike-in_rDNA-alignmentRate.svg"
svg(filename, height=7.5, width=7.5, bg = 'transparent')
p + geom_vline(xintercept=20, color='#333333ff', linetype="dashed")
dev.off()

# Degradation plot
degradationPlot <- ggplot(core_samples,
                          aes(x=`Degradation_ratio`, y=sample_name)) +
    geom_bar(stat="identity") +
    theme_PEPPRO()
p <- degradationPlot +
    labs(x='Degradation ratio', y='') +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
filename = "Core-Samples_degradation-ratio.svg"
svg(filename, height=7.5, width=7.5, bg = 'transparent')
p + geom_vline(xintercept=1, color='#333333ff', linetype="dashed") +
    coord_cartesian(xlim=c(0,1.4))
dev.off()

################################################################################
# H9 subsetting

H9_subset <- rbind(stats[grepl("subset H9", stats$sample_desc),],
                   stats[sample_name == "H9_PRO-seq_100",])

H9_subset$sample_name <- factor(H9_subset$sample_name,
                                levels = rev(c('H9_PRO-seq_10',
                                               'H9_PRO-seq_20',
                                               'H9_PRO-seq_30',
                                               'H9_PRO-seq_40',
                                               'H9_PRO-seq_50',
                                               'H9_PRO-seq_60',
                                               'H9_PRO-seq_70',
                                               'H9_PRO-seq_80',
                                               'H9_PRO-seq_90',
                                               'H9_PRO-seq_100')))

mRNAplot <- ggplot(H9_subset, aes(x=mRNA_contamination, y=sample_name)) +
    geom_bar(stat="identity") +
    theme_PEPPRO()
mRNAplot <- mRNAplot +
    coord_cartesian(xlim=c(0,15)) +
    labs(x=expression(median ((over(exon[RPKM], intron[RPKM]))~X~Gene)),
         y='') +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
mRNAplot <- mRNAplot + annotate("rect", ymin=0, ymax=Inf, xmin=1, xmax=1.5, alpha=0.5)
filename = "H9_mRNA-contamination.svg"
svg(filename, height=7.5, width=7.5, bg = 'transparent')
mRNAplot + xlab('Median mRNA contamination')
dev.off()

PIplot <- ggplot(H9_subset, aes(x=Pause_index, y=sample_name)) +
    geom_bar(stat="identity") +
    theme_PEPPRO()
p <- PIplot + labs(x="median(pause index)", y='') +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
p <- p + geom_vline(xintercept=10, color='#333333ff')
filename = "H9_pause-index.svg"
svg(filename, height=7.5, width=7.5, bg = 'transparent')
p + coord_cartesian(xlim=c(0,80))
dev.off()

ExpectedUniquePlot <- ggplot(H9_subset[H9_subset$umi_len > 0,],
                             aes(x=sample_name, y=Frac_exp_unique_at_10M)) +
    geom_bar(stat="identity",
             fill = c("#f8766dff", "#d89000ff", "#a3a500ee", "#39b600cf",
                      "#00bf7df7", "#00bfc4ff", "#00b0f6ff", "#9590ffff",
                      "#e76bf3ff", "#ff62bcff")) +
    theme_PEPPRO()
p <- ExpectedUniquePlot + labs(y="Fraction expected unique at 10M reads", x='') +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
filename = "H9_frac-exp-unique-at-10M.svg"
svg(filename, height=7.5, width=7.5, bg = 'transparent')
p + geom_hline(yintercept=0.75, color='#333333ff')
dev.off()

TooShortPlot <- ggplot(H9_subset,
                       aes(x=Pct_uninformative_adapter_reads, y=sample_name)) +
    geom_bar(stat="identity") +
    theme_PEPPRO()
p <- TooShortPlot + labs(x="Reads too short (%)", y='') +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
filename = "H9_pct-uninformative-adapter-reads.svg"
svg(filename, height=7.5, width=7.5, bg = 'transparent')
p + geom_vline(xintercept=5, color='#333333ff', linetype='dashed') +
    coord_cartesian(xlim=c(0,50)) +
    xlab('Uninformative adapter reads (%)') +
    annotate("rect", ymin=0, ymax=Inf, xmin=20, xmax=50, alpha=0.5)
dev.off()

TSSplusplot <- ggplot(H9_subset, aes(x=TSS_coding_score, y=sample_name)) +
    geom_bar(stat="identity", fill=rgb(0,139,69, max=255)) +
    theme_PEPPRO()
tssPlus <- TSSplusplot + labs(x="TSS enrichment score \n(+ strand)", y='') +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

TSSminusplot <- ggplot(H9_subset, aes(x=sample_name, y=`TSS_non-coding_score`)) +
    geom_bar(stat="identity", fill=rgb(0,0,255, max=255)) +
    theme_PEPPRO()
tssMinus <- TSSminusplot + labs(y="TSS enrichment score \n(- strand)", x='') +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

# Combined TSS plot
tssPlot <- tssPlus + geom_bar(aes(x=`TSS_non-coding_score`, y=sample_name),
                              stat="identity", fill=rgb(0,0,255, max=255)) +
    labs(x="TSS enrichment score", y='') +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
    annotate("text", x=60, y=Inf, hjust=0, vjust=2,
             size=theme_get()$text[["size"]]/3,
             label="coding", angle=0, col=rgb(0,139,69, max=255)) +
    annotate("text", x=60, y=Inf, hjust=0, vjust=4,
             size=theme_get()$text[["size"]]/3,
             label="non-coding", angle=0, col=rgb(0,0,255, max=255))
filename = "H9_TSS-scores.svg"
svg(filename, height=7.5, width=7.5, bg = 'transparent')
tssPlot + coord_cartesian(xlim=c(0,70))
dev.off()

Alignment_rate_human_rDNAplot <- ggplot(H9_subset,
                                        aes(x=Alignment_rate_human_rDNA,
                                            y=sample_name)) +
    geom_bar(stat="identity") +
    theme_PEPPRO()
p <- Alignment_rate_human_rDNAplot +
    coord_cartesian(xlim=c(0,60)) +
    labs(x="Alignment rate rDNA", y='') +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
filename = "H9_rDNA-alignmentRate.svg"
svg(filename, height=7.5, width=7.5, bg = 'transparent')
p + geom_vline(xintercept=20, color='#333333ff')
dev.off()

# Degradation plot
degradationPlot <- ggplot(H9_subset,
                          aes(x=`Degradation_ratio`, y=sample_name)) +
    geom_bar(stat="identity") +
    theme_PEPPRO()
p <- degradationPlot +
    labs(x='Degradation ratio', y='') +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
filename = "H9_degradation-ratio.svg"
svg(filename, height=7.5, width=7.5, bg = 'transparent')
p + geom_vline(xintercept=1, color='#333333ff', linetype="dashed") +
    coord_cartesian(xlim=c(0,1.4))
dev.off()

## FRiF plot
h9_frif <- fread('H9_subset_FRiF_values.csv')
g <- ggplot(h9_frif, aes(x=pct, y=value, color=feature)) +
    geom_line(size=1.5) +
    geom_point(shape=0, size=3)
g <- g + scale_color_manual(values=c("#fe3b2fff",
                                     "#ff9900ff",
                                     "#00b300ff",
                                     "#005f94ff",
                                     "#00b5ebff",
                                     "#fe8affff",
                                     "#6c01e4ff"))
g <- g + scale_x_continuous(breaks=c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100))

g <- g + theme(
    axis.line = element_line(size = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA),
    legend.key=element_rect(fill = "transparent", color = NA),
    legend.title = element_blank(),
    aspect.ratio = 0.3,
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5)
)

svg("H9_FRiF.svg")
g + labs(x="H9 PRO-seq subsample %", y="FRiF") +
    coord_cartesian(ylim=c(0,0.6)) +
    scale_y_continuous(breaks=c(0.2, 0.4, 0.6))
dev.off()

## plot original pause index calculation (prior to restricting by expression quartile)
## For H9, load a file with just the H9 samples' pause index value when using
## the previous approach (all genes included).
pi <- fread('PEPPRO_original_pause_index_calculation.csv')
pi$sample_name <- factor(pi$sample_name,
                         levels = rev(c('H9_PRO-seq_10',
                                        'H9_PRO-seq_20',
                                        'H9_PRO-seq_30',
                                        'H9_PRO-seq_40',
                                        'H9_PRO-seq_50',
                                        'H9_PRO-seq_60',
                                        'H9_PRO-seq_70',
                                        'H9_PRO-seq_80',
                                        'H9_PRO-seq_90',
                                        'H9_PRO-seq_100')))
PIplot <- ggplot(pi, aes(x=Pause_index, y=sample_name)) +
    geom_bar(stat="identity") +
    theme_PEPPRO()
p <- PIplot + labs(x="median(pause index)", y='') +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
p <- p + geom_vline(xintercept=10, color='#333333ff', linetype='dashed')
filename = "H9_pause-index_original-calculation.svg"
svg(filename, height=7.5, width=7.5, bg = 'transparent')
p + coord_cartesian(xlim=c(0,80))
dev.off()

################################################################################
# K562 subsetting

K562_subset <- stats[grep('K562_PRO-seq_*', stats$sample_name),]

mRNAplot <- ggplot(K562_subset, aes(x=mRNA_contamination, y=sample_name)) +
    geom_bar(stat="identity") +
    theme_PEPPRO()
mRNAplot <- mRNAplot +
    coord_cartesian(xlim=c(0,15)) +
    labs(x=expression(median ((over(exon[RPKM], intron[RPKM]))~X~Gene)),
         y='') +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
mRNAplot <- mRNAplot + annotate("rect", ymin=0, ymax=Inf, xmin=1, xmax=1.5, alpha=0.5)
filename = "K562_mRNA-contamination.svg"
svg(filename, height=7.5, width=7.5, bg = 'transparent')
mRNAplot + xlab('Median mRNA contamination')
dev.off()

PIplot <- ggplot(K562_subset, aes(x=Pause_index, y=sample_name)) +
    geom_bar(stat="identity") +
    theme_PEPPRO()
p <- PIplot + labs(x="median(pause index)", y='') +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
p <- p + geom_vline(xintercept=10, color='#333333ff')
filename = "K562_pause-index.svg"
svg(filename, height=7.5, width=7.5, bg = 'transparent')
p + coord_cartesian(xlim=c(0,80))
dev.off()

TooShortPlot <- ggplot(K562_subset,
                       aes(x=Pct_uninformative_adapter_reads, y=sample_name)) +
    geom_bar(stat="identity") +
    theme_PEPPRO()
p <- TooShortPlot + labs(x="Reads too short (%)", y='') +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
filename = "K562_pct-uninformative-adapter-reads.svg"
svg(filename, height=7.5, width=7.5, bg = 'transparent')
p + geom_vline(xintercept=5, color='#333333ff', linetype='dashed') +
    coord_cartesian(xlim=c(0,50)) +
    xlab('Uninformative adapter reads (%)') +
    annotate("rect", ymin=0, ymax=Inf, xmin=20, xmax=50, alpha=0.5)
dev.off()

TSSplusplot <- ggplot(K562_subset, aes(x=TSS_coding_score, y=sample_name)) +
    geom_bar(stat="identity", fill=rgb(0,139,69, max=255)) +
    theme_PEPPRO()
tssPlus <- TSSplusplot + labs(x="TSS enrichment score \n(+ strand)", y='') +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

TSSminusplot <- ggplot(K562_subset, aes(x=sample_name, y=`TSS_non-coding_score`)) +
    geom_bar(stat="identity", fill=rgb(0,0,255, max=255)) +
    theme_PEPPRO()
tssMinus <- TSSminusplot + labs(y="TSS enrichment score \n(- strand)", x='') +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

# Combined TSS plot
tssPlot <- tssPlus + geom_bar(aes(x=`TSS_non-coding_score`, y=sample_name),
                              stat="identity", fill=rgb(0,0,255, max=255)) +
    labs(x="TSS enrichment score", y='') +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
    annotate("text", x=60, y=Inf, hjust=0, vjust=2,
             size=theme_get()$text[["size"]]/3,
             label="coding", angle=0, col=rgb(0,139,69, max=255)) +
    annotate("text", x=60, y=Inf, hjust=0, vjust=4,
             size=theme_get()$text[["size"]]/3,
             label="non-coding", angle=0, col=rgb(0,0,255, max=255))
filename = "K562_TSS-scores.svg"
svg(filename, height=7.5, width=7.5, bg = 'transparent')
tssPlot + coord_cartesian(xlim=c(0,70))
dev.off()

Alignment_rate_human_rDNAplot <- ggplot(K562_subset,
                                        aes(x=Alignment_rate_human_rDNA,
                                            y=sample_name)) +
    geom_bar(stat="identity") +
    theme_PEPPRO()
p <- Alignment_rate_human_rDNAplot +
    coord_cartesian(xlim=c(0,60)) +
    labs(x="Alignment rate rDNA", y='') +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
filename = "K562_rDNA-alignmentRate.svg"
svg(filename, height=7.5, width=7.5, bg = 'transparent')
p + geom_vline(xintercept=20, color='#333333ff', linetype="dashed")
dev.off()

# Degradation plot
degradationPlot <- ggplot(K562_subset,
                          aes(x=`Degradation_ratio`, y=sample_name)) +
    geom_bar(stat="identity") +
    theme_PEPPRO()
p <- degradationPlot +
    labs(x='Degradation ratio', y='') +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
filename = "K562_degradation-ratio.svg"
svg(filename, height=7.5, width=7.5, bg = 'transparent')
p + geom_vline(xintercept=1, color='#333333ff', linetype="dashed") +
    coord_cartesian(xlim=c(0,1.4))
dev.off()

## FRiF plot
k562_frif <- fread('K562_PRO-seq_subset_FRiF_values.csv')
g <- ggplot(k562_frif, aes(x=pct, y=value, color=feature)) +
    geom_line(size=1.5) +
    geom_point(shape=0, size=3)
g <- g + scale_color_manual(values=c("#fe3b2fff",
                                     "#ff9900ff",
                                     "#00b300ff",
                                     "#005f94ff",
                                     "#00b5ebff",
                                     "#fe8affff",
                                     "#6c01e4ff"))
g <- g + scale_x_continuous(breaks=c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100))

g <- g + theme(
        axis.line = element_line(size = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        legend.key=element_rect(fill = "transparent", color = NA),
        legend.title = element_blank(),
        aspect.ratio = 0.3,
        legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)
    )

svg("K562_FRiF.svg")
g + labs(x="K562 PRO-seq subsample %", y="FRiF")
dev.off()

###############################################################################
# H9 low complexity samples

lc <- fread('PEPPRO_low_complexity_stats_summary.tsv')
lc$sample_name <- factor(lc$sample_name,
                         levels = rev(c('H9_PRO-seq_50pct',
                                        'H9_PRO-seq_60pct',
                                        'H9_PRO-seq_70pct',
                                        'H9_PRO-seq_80pct',
                                        'H9_PRO-seq_90pct',
                                        'H9_PRO-seq_92pct',
                                        'H9_PRO-seq_94pct',
                                        'H9_PRO-seq_96pct',
                                        'H9_PRO-seq_98pct',
                                        'H9_PRO-seq_100pct')))
lc$rDNA_ratio_raw <- lc$Aligned_reads_human_rDNA/lc$Raw_read

mRNAplot <- ggplot(lc, aes(x=mRNA_contamination, y=sample_name)) +
    geom_bar(stat="identity") +
    theme_PEPPRO()
mRNAplot <- mRNAplot +
    coord_cartesian(xlim=c(0,15)) +
    labs(x=expression(median ((over(exon[RPKM], intron[RPKM]))~X~Gene)),
         y='') +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
mRNAplot <- mRNAplot + annotate("rect", ymin=0, ymax=Inf, xmin=1, xmax=1.5, alpha=0.5)
filename = "Low-Complexity_mRNA-contamination.svg"
svg(filename, height=7.5, width=7.5, bg = 'transparent')
mRNAplot + xlab('Median mRNA contamination')
dev.off()

PIplot <- ggplot(lc, aes(x=Pause_index, y=sample_name)) +
    geom_bar(stat="identity") +
    theme_PEPPRO()
p <- PIplot + labs(x="median(pause index)", y='') +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
p <- p + geom_vline(xintercept=10, color='#333333ff', linetype="dashed")
filename = "Low-Complexity_pause-index.svg"
svg(filename, height=7.5, width=7.5, bg = 'transparent')
p + coord_cartesian(xlim=c(0,80))
dev.off()

TooShortPlot <- ggplot(lc,
                       aes(x=Pct_uninformative_adapter_reads, y=sample_name)) +
    geom_bar(stat="identity") +
    theme_PEPPRO()
p <- TooShortPlot + labs(x="Reads too short (%)", y='') +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
filename = "Low-Complexity_pct-uninformative-adapter-reads.svg"
svg(filename, height=7.5, width=7.5, bg = 'transparent')
p + geom_vline(xintercept=5, color='#333333ff', linetype='dashed') +
    coord_cartesian(xlim=c(0,50)) +
    xlab('Uninformative adapter reads (%)') +
    annotate("rect", ymin=0, ymax=Inf, xmin=20, xmax=50, alpha=0.5)
dev.off()

TSSplusplot <- ggplot(lc, aes(x=TSS_coding_score, y=sample_name)) +
    geom_bar(stat="identity", fill=rgb(0,139,69, max=255)) +
    theme_PEPPRO()
tssPlus <- TSSplusplot + labs(x="TSS enrichment score \n(+ strand)", y='') +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

TSSminusplot <- ggplot(lc, aes(x=sample_name, y=`TSS_non-coding_score`)) +
    geom_bar(stat="identity", fill=rgb(0,0,255, max=255)) +
    theme_PEPPRO()
tssMinus <- TSSminusplot + labs(y="TSS enrichment score \n(- strand)", x='') +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

# Combined TSS plot
tssPlot <- tssPlus + geom_bar(aes(x=`TSS_non-coding_score`, y=sample_name),
                              stat="identity", fill=rgb(0,0,255, max=255)) +
    labs(x="TSS enrichment score", y='') +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
    annotate("text", x=60, y=Inf, hjust=0, vjust=2,
             size=theme_get()$text[["size"]]/3,
             label="coding", angle=0, col=rgb(0,139,69, max=255)) +
    annotate("text", x=60, y=Inf, hjust=0, vjust=4,
             size=theme_get()$text[["size"]]/3,
             label="non-coding", angle=0, col=rgb(0,0,255, max=255))
filename = "Low-Complexity_TSS-scores.svg"
svg(filename, height=7.5, width=7.5, bg = 'transparent')
tssPlot + coord_cartesian(xlim=c(0,70))
dev.off()

Alignment_rate_rDNAplot <- ggplot(lc,
                                  aes(x=Alignment_rate_human_rDNA,
                                      y=sample_name)) +
    geom_bar(stat="identity") +
    theme_PEPPRO()
p <- Alignment_rate_rDNAplot +
    coord_cartesian(xlim=c(0,60)) +
    labs(x="Alignment rate rDNA", y='') +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
filename = "Low-Complexity_rDNA-alignmentRate.svg"
svg(filename, height=7.5, width=7.5, bg = 'transparent')
p + geom_vline(xintercept=20, color='#333333ff', linetype="dashed")
dev.off()

# Degradation plot
degradationPlot <- ggplot(lc,
                          aes(x=`Degradation_ratio`, y=sample_name)) +
    geom_bar(stat="identity") +
    theme_PEPPRO()
p <- degradationPlot +
    labs(x='Degradation ratio', y='') +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
filename = "Low-Complexity_degradation-ratio.svg"
svg(filename, height=7.5, width=7.5, bg = 'transparent')
p + geom_vline(xintercept=1, color='#333333ff', linetype="dashed") +
    coord_cartesian(xlim=c(0,1.4))
dev.off()

ExpectedUniquePlot <- ggplot(lc[lc$umi_len > 0,],
                             aes(x=sample_name, y=Frac_exp_unique_at_10M)) +
    geom_bar(stat="identity",
             fill = c("#f8766dff", "#b79f00ff", "#00ba38ff", "#00bfc4ff",
                      "#619cffff", "#f564e3d9")) +
    theme_PEPPRO()
p <- ExpectedUniquePlot + labs(y="Fraction expected unique at 10M reads", x='') +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
filename = "Low-Complexity_frac-exp-unique-at-10M.svg"
svg(filename, height=7.5, width=7.5, bg = 'transparent')
p +
    geom_hline(yintercept=0.75, color='#333333ff', linetype="dashed") +
    annotate("text", x=6.9, y=0, label="*", size=8, col="#333333ff") +
    annotate("text", x=7.9, y=0, label="*", size=8, col="#333333ff") +
    annotate("text", x=8.9, y=0, label="*", size=8, col="#333333ff") +
    annotate("text", x=9.9, y=0, label="*", size=8, col="#333333ff")
dev.off()

## FRiF plot
lc_frif <- fread('low_complexity_FRiF_values.csv')
g <- ggplot(lc_frif, aes(x=pct, y=value, color=feature)) +
    geom_line(size=1.5) +
    geom_point(shape=0, size=3)
g <- g + scale_color_manual(values=c("#fe3b2fff",
                                     "#ff9900ff",
                                     "#00b300ff",
                                     "#005f94ff",
                                     "#00b5ebff",
                                     "#fe8affff",
                                     "#6c01e4ff"))
g <- g + scale_x_continuous(breaks=c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100))

g <- g + theme(
        axis.line = element_line(size = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        legend.key=element_rect(fill = "transparent", color = NA),
        legend.title = element_blank(),
        aspect.ratio = 0.3,
        legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)
    ) +
    labs(x="H9 PRO-seq complexity %", y="FRiF") +
    coord_cartesian(ylim=c(0,0.6)) +
    scale_y_continuous(breaks=c(0.2, 0.4, 0.6))

svg("Low-Complexity_FRiF.svg")
g
dev.off()

################################################################################
# RNA-seq spike-in

## FRiF
rna_frif <- fread('K562_RNA-seq_spike-in_FRiF_values.csv')
g <- ggplot(rna_frif, aes(x=pct, y=value, color=feature)) +
    geom_line(size=1.5) +
    geom_point(shape=0, size=3)
g <- g + scale_color_manual(values=c("#fe3b2fff",
                                     "#ff9900ff",
                                     "#00b300ff",
                                     "#005f94ff",
                                     "#00b5ebff",
                                     "#fe8affff",
                                     "#6c01e4ff"))
g <- g + scale_x_continuous(breaks=c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100))

g <- g + theme(
    axis.line = element_line(size = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA),
    legend.key=element_rect(fill = "transparent", color = NA),
    legend.title = element_blank(),
    aspect.ratio = 0.3,
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5)
) 

svg("K562_RNA-seq_FRiF.svg")
g + labs(x="K562 RNA-seq %", y="FRiF")
dev.off()
