# Prep

Make sure the sample_name method matches the protocol_type method. There are fast-ATAC samples that were labeled Omni and some standard-ATAC labeled as Omni.

```
R
library(data.table)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(gridExtra)

pre             <- fread('prealignments_stats_summary.tsv')
no_pre          <- fread('no_prealignments_stats_summary.tsv')

colnames(pre) <- c("sample_name", "srr", "protocol", "read_type", "read1", "read2", "organism", "protocol_type", "gsm_accession", "cell_line_or_source", "cell_type", "Raw_reads", "Fastq_reads", "Trimmed_reads", "Trim_loss_rate", "Mapped_reads", "QC_filtered_reads", "Aligned_reads", "Alignment_rate", "Total_efficiency", "Mitochondrial_reads", "NRF", "PBC1", "PBC2", "Unmapped_reads", "Duplicate_reads", "Dedup_aligned_reads", "Dedup_alignment_rate", "Dedup_total_efficiency", "Read_length", "Genome_size", "Frac_exp_unique_at_10M", "TSS_score", "Peak_count", "FRiP", "File_mb", "Read_type", "Genome", "Aligned_reads_human_repeats", "Alignment_rate_human_repeats", "Aligned_reads_rCRSd", "Alignment_rate_rCRSd", "FRiP_ref", "Time", "Success", "Aligned_reads_mm_mtDNAd", "Alignment_rate_mm_mtDNAd", "Aligned_reads_at_mtDNAd", "Alignment_rate_at_mtDNAd")
pre$method <- "pre"

colnames(no_pre) <- c("sample_name", "srr", "protocol", "read_type", "read1", "read2", "organism", "protocol_type", "gsm_accession", "cell_line_or_source", "cell_type", "Raw_reads", "Fastq_reads", "Trimmed_reads", "Trim_loss_rate", "Mapped_reads", "QC_filtered_reads", "Aligned_reads", "Alignment_rate", "Total_efficiency", "Mitochondrial_reads", "NRF", "PBC1", "PBC2", "Unmapped_reads", "Duplicate_reads", "Dedup_aligned_reads", "Dedup_alignment_rate", "Dedup_total_efficiency", "Read_length", "Genome_size", "Frac_exp_unique_at_10M", "TSS_score", "Peak_count", "FRiP", "File_mb", "Read_type", "Genome", "FRiP_ref", "Time", "Success")
no_pre$method <- "no-pre"

stats <- plyr::rbind.fill(pre, no_pre)

stats$protocol_type <- ordered(stats$protocol_type,
                               levels = c("Standard-ATAC", "Fast-ATAC", "Omni-ATAC"))

stats$method <- ordered(stats$method, levels = c("pre", "no-pre"))
                           
t <- data.table(sample_name=rep(stats$sample_name, times=14),
                method=rep(stats$method, times=14),
                organism=rep(stats$organism, times=14),
                library_protocol=rep(stats$protocol_type, times=14),
                type=rep(c("Raw_reads", "Fastq_reads", "Trimmed_reads",
                           "Mapped_reads", "QC_filtered_reads", "Aligned_reads",
                           "Unmapped_reads", "Duplicate_reads",
                           "Dedup_aligned_reads", "Aligned_reads_human_repeats",
                           "Aligned_reads_rCRSd", "Aligned_reads_mm_mtDNAd",
                           "Aligned_reads_at_mtDNAd", "Mitochondrial_reads"),
                         each=nrow(stats)))
t$reads <- c(stats$Raw_reads, stats$Fastq_reads, stats$Trimmed_reads,
             stats$Mapped_reads, stats$QC_filtered_reads, stats$Aligned_reads,
             stats$Unmapped_reads, stats$Duplicate_reads,
             stats$Dedup_aligned_reads, stats$Aligned_reads_human_repeats,
             stats$Aligned_reads_rCRSd, stats$Aligned_reads_mm_mtDNAd,
             stats$Aligned_reads_at_mtDNAd, stats$Mitochondrial_reads)

# Replace NA with 0
t[is.na(t)] <- 0

t$pct_aligned <- 100*(t$reads/t$reads[t$type == "Raw_reads"])

total <- data.table(sample_name=stats$sample_name,
                    method=stats$method,
                    organism=stats$organism,
                    library_protocol=stats$protocol_type,
                    type="Total_aligned_reads",
                    reads=(t[t$type == "Aligned_reads", reads] +
                           t[t$type == "Aligned_reads_human_repeats", reads] +
                           t[t$type == "Aligned_reads_rCRSd", reads] +
                           t[t$type == "Aligned_reads_mm_mtDNAd", reads] +
                           t[t$type == "Aligned_reads_at_mtDNAd", reads]),
                    pct_aligned=100*((t[t$type == "Aligned_reads", reads] +
                                      t[t$type == "Aligned_reads_human_repeats", reads] +
                                      t[t$type == "Aligned_reads_rCRSd", reads] +
                                      t[t$type == "Aligned_reads_mm_mtDNAd", reads] +
                                      t[t$type == "Aligned_reads_at_mtDNAd", reads])/
                                      t[t$type == "Raw_reads", reads]))

t2 <- rbind(t, total)
t2$method <- ordered(t2$method, levels = c("pre", "no-pre"))

# add a total mitochondrial aligned reads just like this
mtDNA <- data.table(sample_name=stats$sample_name,
                    method=stats$method,
                    organism=stats$organism,
                    library_protocol=stats$protocol_type,
                    type="Total_mtDNA_reads",
                    reads=(t[t$type == "Aligned_reads_rCRSd", reads] +
                           t[t$type == "Aligned_reads_mm_mtDNAd", reads] +
                           t[t$type == "Aligned_reads_at_mtDNAd", reads] +
                           t[t$type == "Mitochondrial_reads", reads]),
                    pct_aligned=100*((t[t$type == "Aligned_reads_rCRSd", reads] +
                                      t[t$type == "Aligned_reads_mm_mtDNAd", reads] +
                                      t[t$type == "Aligned_reads_at_mtDNAd", reads] +
                                      t[t$type == "Mitochondrial_reads", reads])/
                                      t[t$type == "Raw_reads", reads]))

t3 <- rbind(t2, mtDNA)
t3$method <- ordered(t3$method, levels = c("pre", "no-pre""))
# add all reads together to show 100% total (should match raw_reads total)
allReads <- data.table(sample_name=stats$sample_name,
                       method=stats$method,
                       organism=stats$organism,
                       library_protocol=stats$protocol_type,
                       type="Summed_reads",
                       reads=(t[t$type == "Aligned_reads_rCRSd", reads] +
                              t[t$type == "Aligned_reads_mm_mtDNAd", reads] +
                              t[t$type == "Aligned_reads_at_mtDNAd", reads] +
                              t[t$type == "Aligned_reads_human_repeats", reads] +
                              #t[t$type == "Mitochondrial_reads", reads] +
                              t[t$type == "Aligned_reads", reads] +
                              t[t$type == "QC_filtered_reads", reads] +
                              t[t$type == "Unmapped_reads", reads] +
                              (t[t$type == "Raw_reads", reads] - t[t$type == "Trimmed_reads", reads]) ),
                       pct_aligned=100*((t[t$type == "Aligned_reads_rCRSd", reads] +
                                         t[t$type == "Aligned_reads_mm_mtDNAd", reads] +
                                         t[t$type == "Aligned_reads_at_mtDNAd", reads] +
                                         t[t$type == "Aligned_reads_human_repeats", reads] +
                                         #t[t$type == "Mitochondrial_reads", reads] +
                                         t[t$type == "Aligned_reads", reads] +
                                         t[t$type == "QC_filtered_reads", reads] +
                                         t[t$type == "Unmapped_reads", reads] +
                                         (t[t$type == "Raw_reads", reads] - t[t$type == "Trimmed_reads", reads]) )/
                                         t[t$type == "Raw_reads", reads]))
t4 <- rbind(t3, allReads)
t4$method <- ordered(t4$method, levels = c("pre", "no-pre"))

t6 <- t5[t5$type == "Total_mtDNA_reads"]
t6_all <- t4[t4$type == "Total_mtDNA_reads"]

mtDNA_pre_no_pre_human <- data.table(
    sample_name = t6[method == "pre"]$sample_name,
    method = "pre",
    organism = t6[method == "pre"]$organism,
    library_protocol = t6[method == "pre"]$library_protocol,
    type = "Total_mtDNA_reads",
    reads = t6[method == "pre"]$reads,
    pct_aligned = t6[method == "pre"]$pct_aligned,
    diff = t6[method == "pre"]$pct_aligned - t6[method == "no-pre"]$pct_aligned
)

mtDNA_pre_no_pre_human$library_protocol <- ordered(mtDNA_pre_no_pre_human$library_protocol,
    levels = c("Standard-ATAC", "Fast-ATAC", "Omni-ATAC"))

stat.test <- mtDNA_pre_no_pre_human %>%
    group_by(library_protocol) %>%
    t_test(data=., diff ~ 1, mu=0, p.adjust.method = "BH") %>%
    add_xy_position(x = "library_protocol")

mtDNA_diff_plot_human <- ggplot(mtDNA_pre_no_pre_human, aes(x=library_protocol, y=diff)) +
    geom_violin(alpha=.6, width=.5, aes(fill=library_protocol)) + 
    geom_boxplot(outlier.shape=NA, width=.5, alpha=0, aes(fill=library_protocol)) +
    geom_point(aes(color=library_protocol), position=position_jitterdodge(jitter.width=0.75), alpha=.3) +
    ylab('Additional mtDNA Aligned Reads with Prealignments (%)') +
    xlab('ATAC-seq Library Protocol') +
    stat_pvalue_manual(stat.test, hide.ns = TRUE, alpha=0.75) +
    ylim(0,10) +
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
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)
    )
mtDNA_diff_plot_human

# Let's confirm prealignments always improves mtDNA alignment rates
mtDNA_pre_no_pre_human[diff < 0,] # Empty as we would expect (prealignment always improves it)

svg("human_mtDNA_differences-between-prealignments-or-none.svg")
mtDNA_diff_plot_human
dev.off()

# Now just look at all combined and see

mtDNA_pre_no_pre <- data.table(
    sample_name = t6_all[method == "pre"]$sample_name,
    method = "pre",
    organism = t6_all[method == "pre"]$organism,
    library_protocol = t6_all[method == "pre"]$library_protocol,
    type = "Total_mtDNA_reads",
    reads = t6_all[method == "pre"]$reads,
    pct_aligned = t6_all[method == "pre"]$pct_aligned,
    diff = t6_all[method == "pre"]$pct_aligned - t6_all[method == "no-pre"]$pct_aligned
)

mtDNA_pre_no_pre$library_protocol <- ordered(mtDNA_pre_no_pre$library_protocol,
    levels = c("Standard-ATAC", "Fast-ATAC", "Omni-ATAC"))

# Drop samples without a known "library_protocol" (i.e. the plant samples)
mtDNA_pre_no_pre <- mtDNA_pre_no_pre[complete.cases(mtDNA_pre_no_pre[, 4]), ]

stat.test <- mtDNA_pre_no_pre %>%
    group_by(library_protocol) %>%
    t_test(data=., diff ~ 1, mu=0, p.adjust.method = "BH") %>%
    add_xy_position(x = "library_protocol")

mtDNA_diff_plot <- ggplot(mtDNA_pre_no_pre, aes(x=library_protocol, y=diff)) +
    geom_violin(alpha=.6, width=.5, aes(fill=library_protocol)) + 
    geom_boxplot(outlier.shape=NA, width=.5, alpha=0, aes(fill=library_protocol)) +
    geom_point(aes(color=library_protocol), position=position_jitterdodge(jitter.width=0.75), alpha=.3) +
    ylab('Additional mtDNA Aligned Reads with Prealignments (%)') +
    xlab('ATAC-seq Library Protocol') +
    stat_pvalue_manual(stat.test, hide.ns = TRUE, alpha=0.75) +
    #ylim(0,10) +
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
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)
    )

mtDNA_pre_no_pre[diff < 0,] # Two samples got worse, but overall trend is the same

svg("all_mtDNA_differences-between-prealignments-or-none.svg")
mtDNA_diff_plot
dev.off()

t7 <- t6[method == "pre" | method == "no-pre"]
t7$library_protocol <- ordered(t7$library_protocol,
                               levels = c("Standard-ATAC", "Fast-ATAC", "Omni-ATAC"))
                               
ba <- data.table(
    sample_name=t7[method == "no-pre"]$sample_name,
    library_protocol=t7[method == "no-pre"]$library_protocol,
    before_reads=t7[method == "no-pre"]$reads,
    before_pct_aligned=t7[method == "no-pre"]$pct_aligned,
    after_reads=t7[method == "pre"]$reads,
    after_pct_aligned=t7[method == "pre"]$pct_aligned,
    diff=(t7[method == "pre"]$pct_aligned - t7[method == "no-pre"]$pct_aligned))

standard <- ba[library_protocol == "Standard-ATAC",]
fast <- ba[library_protocol == "Fast-ATAC",]
omni <- ba[library_protocol == "Omni-ATAC",]

# Drop outliers
standard <- standard[!(before_reads %in% boxplot.stats(standard$before_reads)$out),]

stat.test <- standard %>%
    t_test(data=., diff ~ 1, mu=0, p.adjust.method = "BH") %>%
    add_xy_position(x = "library_protocol")

stat.test$y.position <- 80

standard_mtDNA_beforeAfter <- ggplot(standard) + 
    geom_segment(aes(x = 1, xend = 2, y = before_reads/1000000,
                     yend = after_reads/1000000)) + 
    geom_segment(data = standard, aes(x=1, xend=2, y = mean(before_reads/1000000),
                 yend = mean(after_reads/1000000)), col='blue', size=2) +
    stat_pvalue_manual(stat.test, hide.ns = TRUE, alpha=0.75) +
    theme(
        axis.line = element_line(size = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)
    ) +
    scale_x_discrete(
        breaks = c("1", "2"),
        labels = c("None", "Prealignments"),
        limits = c(1, 2)
    ) + 
    labs(x = "Standard-ATAC", y = "mtDNA Reads (M)")

svg("mtDNA_standard_human_beforeAfter.svg")
standard_mtDNA_beforeAfter
dev.off()

# Drop outliers
fast <- fast[!(before_reads %in% boxplot.stats(fast$before_reads)$out),]

stat.test <- fast %>%
    t_test(data=., diff ~ 1, mu=0, p.adjust.method = "BH") %>%
    add_xy_position(x = "library_protocol")
    
stat.test$y.position <- 7

fast_mtDNA_beforeAfter <- ggplot(fast) + 
    geom_segment(aes(x = 1, xend = 2, y = before_reads/1000000,
                     yend = after_reads/1000000)) + 
    geom_segment(data = fast, aes(x=1, xend=2, y = mean(before_reads/1000000),
                 yend = mean(after_reads/1000000)), col='blue', size=2) +
    stat_pvalue_manual(stat.test, hide.ns = TRUE, alpha=0.75) +
    theme(
        axis.line = element_line(size = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)
    ) +
    scale_x_discrete(
        breaks = c("1", "2"),
        labels = c("None", "Prealignments"),
        limits = c(1, 2)
    ) + 
    labs(x = "Fast-ATAC", y = "mtDNA Reads (M)")

svg("mtDNA_fast_human_beforeAfter_noOutliers.svg")
fast_mtDNA_beforeAfter
dev.off()

# Drop outliers
omni <- omni[!(before_reads %in% boxplot.stats(omni$before_reads)$out),]

stat.test <- omni %>%
    t_test(data=., diff ~ 1, mu=0, p.adjust.method = "BH") %>%
    add_xy_position(x = "library_protocol")

stat.test$y.position <- 1.5

omni_mtDNA_beforeAfter <- ggplot(omni) + 
    geom_segment(aes(x = 1, xend = 2, y = before_reads/1000000,
                     yend = after_reads/1000000)) + 
    geom_segment(data = omni, aes(x=1, xend=2, y = mean(before_reads/1000000),
                 yend = mean(after_reads/1000000)), col='blue', size=2) +
    stat_pvalue_manual(stat.test, hide.ns = TRUE, alpha=0.75) +
    theme(
        axis.line = element_line(size = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)
    ) +
    scale_x_discrete(
        breaks = c("1", "2"),
        labels = c("None", "Prealignments"),
        limits = c(1, 2)
    ) + 
    labs(x = "Omni-ATAC", y = "mtDNA Reads (M)")

svg("mtDNA_omni_human_beforeAfter.svg")
omni_mtDNA_beforeAfter
dev.off()

################################################################################
### Default FRiP
FRiP_diff <- data.table(
    sample_name = stats[stats$method == "pre", ]$sample_name,
    organism = stats[stats$method == "pre", ]$organism,
    library_protocol = stats[stats$method == "pre", ]$protocol_type,
    type = "FRiP",
    none_FRiP = stats[stats$method == "no-pre", ]$FRiP,
    pre_FRiP = stats[stats$method == "pre", ]$FRiP,
    diff = stats[stats$method == "pre", ]$FRiP - stats[stats$method == "no-pre", ]$FRiP
)

FRiP_diff$library_protocol <- ordered(FRiP_diff$library_protocol,
                                      levels = c("Standard-ATAC", "Fast-ATAC", "Omni-ATAC"))

FRiP_diff[diff < 0,] # 2 samples in mouse and 2 in human get worse

# Drop samples without library protocol
FRiP_diff <- FRiP_diff[complete.cases(FRiP_diff), ]

# Drop outliers
FRiP_diff <- FRiP_diff[!(none_FRiP %in% boxplot.stats(FRiP_diff$none_FRiP)$out),]

stat.test <- FRiP_diff %>%
    group_by(library_protocol) %>%
    t_test(data=., diff ~ 1, mu=0, p.adjust.method = "BH") %>%
    add_xy_position(x = "library_protocol")
    
# This uses FRiP where the denominator is "Aligned_reads"
FRiP_diff_plot <- ggplot(FRiP_diff, aes(x=library_protocol, y=diff)) +
    geom_violin(alpha=.6, width=.5, aes(fill=library_protocol)) + 
    geom_boxplot(outlier.shape=NA, width=.5, alpha=0, aes(fill=library_protocol)) +
    geom_point(aes(color=library_protocol, shape=organism), position=position_jitterdodge(jitter.width=0.75), alpha=.3) +
    ylab('FRiP') +
    xlab('ATAC-seq Library Protocol') +
    stat_pvalue_manual(stat.test, hide.ns = TRUE, alpha=0.75) +
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
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)
    )
 
svg("FRiP_aligned_reads_diff.svg")
FRiP_diff_plot
dev.off()

### Default FRiP before and after plot (standard, fast, omni)
standard <- FRiP_diff[library_protocol == "Standard-ATAC",]
fast     <- FRiP_diff[library_protocol == "Fast-ATAC",]
omni     <- FRiP_diff[library_protocol == "Omni-ATAC",]

# Drop outliers
standard <- standard[!(diff %in% boxplot.stats(standard$diff)$out),]

stat.test <- standard %>%
    t_test(data=., diff ~ 1, mu=0, p.adjust.method = "BH") %>%
    add_xy_position(x = "library_protocol")

stat.test$y.position <- 1

standard_FRiP_beforeAfter <- ggplot(standard) + 
    geom_segment(aes(x = 1, xend = 2, y = none_FRiP,
                     yend = pre_FRiP, col = organism)) + 
    geom_segment(data = standard, aes(x=1, xend=2, y = mean(none_FRiP),
                 yend = mean(pre_FRiP)), col='blue', size=2) +
    stat_pvalue_manual(stat.test, hide.ns = TRUE, alpha=0.75) +
    theme(
        axis.line = element_line(size = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)
    ) +
    scale_x_discrete(
        breaks = c("1", "2"),
        labels = c("None", "Prealignments"),
        limits = c(1, 2)
    ) + 
    labs(x = "Standard-ATAC", y = "FRiP")

svg("FRiP_standard_aligned_reads_diff_beforeAfter.svg")
standard_FRiP_beforeAfter
dev.off()

# Drop outliers
fast <- fast[!(diff %in% boxplot.stats(fast$diff)$out),]

stat.test <- fast %>%
    t_test(data=., diff ~ 1, mu=0, p.adjust.method = "BH") %>%
    add_xy_position(x = "library_protocol")

stat.test$y.position <- 1

fast_FRiP_beforeAfter <- ggplot(fast) + 
    geom_segment(aes(x = 1, xend = 2, y = none_FRiP,
                     yend = pre_FRiP, col = organism)) + 
    geom_segment(data = fast, aes(x=1, xend=2, y = mean(none_FRiP),
                 yend = mean(pre_FRiP)), col='blue', size=2) +
    stat_pvalue_manual(stat.test, hide.ns = TRUE, alpha=0.75) +
    theme(
        axis.line = element_line(size = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)
    ) +
    scale_x_discrete(
        breaks = c("1", "2"),
        labels = c("None", "Prealignments"),
        limits = c(1, 2)
    ) + 
    labs(x = "Fast-ATAC", y = "FRiP")

svg("FRiP_fast_aligned_reads_diff_beforeAfter.svg")
fast_FRiP_beforeAfter
dev.off()

# Drop outliers
omni <- omni[!(diff %in% boxplot.stats(omni$diff)$out),]

stat.test <- omni %>%
    t_test(data=., diff ~ 1, mu=0, p.adjust.method = "BH") %>%
    add_xy_position(x = "library_protocol")

stat.test$y.position <- 1

omni_FRiP_beforeAfter <- ggplot(omni) + 
    geom_segment(aes(x = 1, xend = 2, y = none_FRiP,
                     yend = pre_FRiP, col = organism)) + 
    geom_segment(data = omni, aes(x=1, xend=2, y = mean(none_FRiP),
                 yend = mean(pre_FRiP)), col='blue', size=2) +
    stat_pvalue_manual(stat.test, hide.ns = TRUE, alpha=0.75) +
    theme(
        axis.line = element_line(size = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)
    ) +
    scale_x_discrete(
        breaks = c("1", "2"),
        labels = c("None", "Prealignments"),
        limits = c(1, 2)
    ) + 
    labs(x = "Omni-ATAC", y = "FRiP")

svg("FRiP_omni_aligned_reads_diff_beforeAfter.svg")
omni_FRiP_beforeAfter
dev.off()

################################################################################
### Total aligned reads corrected FRiP
stats[is.na(stats)] <- 0

corrected_FRiP_diff <- data.table(
    sample_name = stats[stats$method == "pre", ]$sample_name,
    organism = stats[stats$method == "pre", ]$organism,
    library_protocol = stats[stats$method == "pre", ]$protocol_type,
    type = "FRiP",
    none_FRiP = stats[stats$method == "no-pre", ]$FRiP,
    pre_FRiP = ((stats[stats$method == "pre", ]$FRiP * stats[stats$method == "pre", ]$Aligned_reads)/ (stats[stats$method == "pre", ]$Aligned_reads + stats[stats$method == "pre", ]$Aligned_reads_human_repeats + stats[stats$method == "pre", ]$Aligned_reads_rCRSd + stats[stats$method == "pre", ]$Aligned_reads_mm_mtDNAd + stats[stats$method == "pre", ]$Aligned_reads_at_mtDNAd))
)

corrected_FRiP_diff$diff <- (corrected_FRiP_diff$pre_FRiP - corrected_FRiP_diff$none_FRiP)

# Drop samples without library protocol
corrected_FRiP_diff <- corrected_FRiP_diff[complete.cases(corrected_FRiP_diff), ]

# Drop outliers
corrected_FRiP_diff <- corrected_FRiP_diff[!(none_FRiP %in% boxplot.stats(corrected_FRiP_diff$none_FRiP)$out),]

stat.test <- corrected_FRiP_diff %>%
    group_by(library_protocol) %>%
    t_test(data=., diff ~ 1, mu=0, p.adjust.method = "BH") %>%
    add_xy_position(x = "library_protocol")
    
# This uses FRiP where the denominator is "Aligned_reads"
FRiP_diff_plot <- ggplot(corrected_FRiP_diff, aes(x=library_protocol, y=diff)) +
    geom_violin(alpha=.6, width=.5, aes(fill=library_protocol)) + 
    geom_boxplot(outlier.shape=NA, width=.5, alpha=0, aes(fill=library_protocol)) +
    geom_point(aes(color=library_protocol, shape=organism), position=position_jitterdodge(jitter.width=0.75), alpha=.3) +
    ylab('FRiP') +
    xlab('ATAC-seq Library Protocol') +
    stat_pvalue_manual(stat.test, hide.ns = TRUE, alpha=0.75) +
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
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)
    )
 
svg("FRiP_total_aligned_reads_diff.svg")
FRiP_diff_plot
dev.off()

### Default FRiP before and after plot (standard, fast, omni)
standard <- corrected_FRiP_diff[library_protocol == "Standard-ATAC",]
fast     <- corrected_FRiP_diff[library_protocol == "Fast-ATAC",]
omni     <- corrected_FRiP_diff[library_protocol == "Omni-ATAC",]

# Drop outliers
standard <- standard[!(diff %in% boxplot.stats(standard$diff)$out),]

stat.test <- standard %>%
    t_test(data=., diff ~ 1, mu=0, p.adjust.method = "BH") %>%
    add_xy_position(x = "library_protocol")

stat.test$y.position <- 1

standard_FRiP_beforeAfter <- ggplot(standard) + 
    geom_segment(aes(x = 1, xend = 2, y = none_FRiP,
                     yend = pre_FRiP, col = organism)) + 
    geom_segment(data = standard, aes(x=1, xend=2, y = mean(none_FRiP),
                 yend = mean(pre_FRiP)), col='blue', size=2) +
    stat_pvalue_manual(stat.test, hide.ns = TRUE, alpha=0.75) +
    theme(
        axis.line = element_line(size = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)
    ) +
    scale_x_discrete(
        breaks = c("1", "2"),
        labels = c("None", "Prealignments"),
        limits = c(1, 2)
    ) + 
    labs(x = "Standard-ATAC", y = "FRiP")

svg("FRiP_standard_total_aligned_reads_diff_beforeAfter.svg")
standard_FRiP_beforeAfter
dev.off()

# Drop outliers
fast <- fast[!(diff %in% boxplot.stats(fast$diff)$out),]

stat.test <- fast %>%
    t_test(data=., diff ~ 1, mu=0, p.adjust.method = "BH") %>%
    add_xy_position(x = "library_protocol")

stat.test$y.position <- 1

fast_FRiP_beforeAfter <- ggplot(fast) + 
    geom_segment(aes(x = 1, xend = 2, y = none_FRiP,
                     yend = pre_FRiP, col = organism)) + 
    geom_segment(data = fast, aes(x=1, xend=2, y = mean(none_FRiP),
                 yend = mean(pre_FRiP)), col='blue', size=2) +
    stat_pvalue_manual(stat.test, hide.ns = TRUE, alpha=0.75) +
    theme(
        axis.line = element_line(size = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)
    ) +
    scale_x_discrete(
        breaks = c("1", "2"),
        labels = c("None", "Prealignments"),
        limits = c(1, 2)
    ) + 
    labs(x = "Fast-ATAC", y = "FRiP")

svg("FRiP_fast_total_aligned_reads_diff_beforeAfter.svg")
fast_FRiP_beforeAfter
dev.off()

# Drop outliers
omni <- omni[!(diff %in% boxplot.stats(omni$diff)$out),]

stat.test <- omni %>%
    t_test(data=., diff ~ 1, mu=0, p.adjust.method = "BH") %>%
    add_xy_position(x = "library_protocol")

stat.test$y.position <- 1

omni_FRiP_beforeAfter <- ggplot(omni) + 
    geom_segment(aes(x = 1, xend = 2, y = none_FRiP,
                     yend = pre_FRiP, col = organism)) + 
    geom_segment(data = omni, aes(x=1, xend=2, y = mean(none_FRiP),
                 yend = mean(pre_FRiP)), col='blue', size=2) +
    stat_pvalue_manual(stat.test, hide.ns = TRUE, alpha=0.75) +
    theme(
        axis.line = element_line(size = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)
    ) +
    scale_x_discrete(
        breaks = c("1", "2"),
        labels = c("None", "Prealignments"),
        limits = c(1, 2)
    ) + 
    labs(x = "Omni-ATAC", y = "FRiP")

svg("FRiP_omni_total_aligned_reads_diff_beforeAfter.svg")
omni_FRiP_beforeAfter
dev.off()

################################################################################
### Raw reads FRiP
stats[is.na(stats)] <- 0

raw_FRiP_diff <- data.table(
    sample_name = stats[stats$method == "pre", ]$sample_name,
    organism = stats[stats$method == "pre", ]$organism,
    library_protocol = stats[stats$method == "pre", ]$protocol_type,
    type = "FRiP",
    none_FRiP = stats[stats$method == "no-pre", ]$FRiP,
    pre_FRiP = ((stats[stats$method == "pre", ]$FRiP * stats[stats$method == "pre", ]$Aligned_reads)/ (stats[stats$method == "pre", ]$Raw_reads))
)

raw_FRiP_diff$diff <- (raw_FRiP_diff$pre_FRiP - raw_FRiP_diff$none_FRiP)

# Drop samples without library protocol
raw_FRiP_diff <- raw_FRiP_diff[complete.cases(raw_FRiP_diff), ]

# Drop outliers
raw_FRiP_diff <- raw_FRiP_diff[!(none_FRiP %in% boxplot.stats(raw_FRiP_diff$none_FRiP)$out),]

stat.test <- raw_FRiP_diff %>%
    group_by(library_protocol) %>%
    t_test(data=., diff ~ 1, mu=0, p.adjust.method = "BH") %>%
    add_xy_position(x = "library_protocol")
    
# This uses FRiP where the denominator is "Aligned_reads"
FRiP_diff_plot <- ggplot(raw_FRiP_diff, aes(x=library_protocol, y=diff)) +
    geom_violin(alpha=.6, width=.5, aes(fill=library_protocol)) + 
    geom_boxplot(outlier.shape=NA, width=.5, alpha=0, aes(fill=library_protocol)) +
    geom_point(aes(color=library_protocol, shape=organism), position=position_jitterdodge(jitter.width=0.75), alpha=.3) +
    ylab('FRiP') +
    xlab('ATAC-seq Library Protocol') +
    stat_pvalue_manual(stat.test, hide.ns = TRUE, alpha=0.75) +
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
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)
    )
 
svg("FRiP_raw_aligned_reads_diff.svg")
FRiP_diff_plot
dev.off()

### Default FRiP before and after plot (standard, fast, omni)
standard <- raw_FRiP_diff[library_protocol == "Standard-ATAC",]
fast     <- raw_FRiP_diff[library_protocol == "Fast-ATAC",]
omni     <- raw_FRiP_diff[library_protocol == "Omni-ATAC",]

# Drop outliers
standard <- standard[!(diff %in% boxplot.stats(standard$diff)$out),]

stat.test <- standard %>%
    t_test(data=., diff ~ 1, mu=0, p.adjust.method = "BH") %>%
    add_xy_position(x = "library_protocol")

stat.test$y.position <- 1

standard_FRiP_beforeAfter <- ggplot(standard) + 
    geom_segment(aes(x = 1, xend = 2, y = none_FRiP,
                     yend = pre_FRiP, col = organism)) + 
    geom_segment(data = standard, aes(x=1, xend=2, y = mean(none_FRiP),
                 yend = mean(pre_FRiP)), col='blue', size=2) +
    stat_pvalue_manual(stat.test, hide.ns = TRUE, alpha=0.75) +
    theme(
        axis.line = element_line(size = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)
    ) +
    scale_x_discrete(
        breaks = c("1", "2"),
        labels = c("None", "Prealignments"),
        limits = c(1, 2)
    ) + 
    labs(x = "Standard-ATAC", y = "FRiP")

svg("FRiP_standard_raw_aligned_reads_diff_beforeAfter.svg")
standard_FRiP_beforeAfter
dev.off()

# Drop outliers
fast <- fast[!(diff %in% boxplot.stats(fast$diff)$out),]

stat.test <- fast %>%
    t_test(data=., diff ~ 1, mu=0, p.adjust.method = "BH") %>%
    add_xy_position(x = "library_protocol")

stat.test$y.position <- 1

fast_FRiP_beforeAfter <- ggplot(fast) + 
    geom_segment(aes(x = 1, xend = 2, y = none_FRiP,
                     yend = pre_FRiP, col = organism)) + 
    geom_segment(data = fast, aes(x=1, xend=2, y = mean(none_FRiP),
                 yend = mean(pre_FRiP)), col='blue', size=2) +
    stat_pvalue_manual(stat.test, hide.ns = TRUE, alpha=0.75) +
    theme(
        axis.line = element_line(size = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)
    ) +
    scale_x_discrete(
        breaks = c("1", "2"),
        labels = c("None", "Prealignments"),
        limits = c(1, 2)
    ) + 
    labs(x = "Fast-ATAC", y = "FRiP")

svg("FRiP_fast_raw_aligned_reads_diff_beforeAfter.svg")
fast_FRiP_beforeAfter
dev.off()

# Drop outliers
omni <- omni[!(diff %in% boxplot.stats(omni$diff)$out),]

stat.test <- omni %>%
    t_test(data=., diff ~ 1, mu=0, p.adjust.method = "BH") %>%
    add_xy_position(x = "library_protocol")

stat.test$y.position <- 1

omni_FRiP_beforeAfter <- ggplot(omni) + 
    geom_segment(aes(x = 1, xend = 2, y = none_FRiP,
                     yend = pre_FRiP, col = organism)) + 
    geom_segment(data = omni, aes(x=1, xend=2, y = mean(none_FRiP),
                 yend = mean(pre_FRiP)), col='blue', size=2) +
    stat_pvalue_manual(stat.test, hide.ns = TRUE, alpha=0.75) +
    theme(
        axis.line = element_line(size = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)
    ) +
    scale_x_discrete(
        breaks = c("1", "2"),
        labels = c("None", "Prealignments"),
        limits = c(1, 2)
    ) + 
    labs(x = "Omni-ATAC", y = "FRiP")

svg("FRiP_omni_raw_aligned_reads_diff_beforeAfter.svg")
omni_FRiP_beforeAfter
dev.off()

```


