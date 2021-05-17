library(data.table)
library(lubridate)
library(ggplot2)

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

filelist <- c("gold1_profile.tsv",
              "gold2_profile.tsv",
              "gold3_profile.tsv",
              "gold4_profile.tsv",
              "gold5_profile.tsv")

profiles <- lapply(filelist, fread)
names(profiles) <- c("gold1", "gold2", "gold3", "gold4", "gold5")

totals <- data.table(run=as.character(),
                     total_runtime=as.numeric(),
                     trd=period(),
                     dhms=as.character(),
                     name=as.character())

for (sample in names(profiles)) {
    tmp  <- unique(profiles[[sample]], by="hash", fromLast=T)
    dt <- tmp[, sum(runtime_seconds), by = run]
    colnames(dt) <- c("run","total_runtime")
    dt[,trd:=seconds_to_period(total_runtime)]
    dt[,dhms:=sprintf('%02g %02g:%02g:%02g', day(trd), trd@hour, minute(trd), second(trd))]
    dt[,name:=sample]
    totals <- rbind(totals, dt)
}

g <- ggplot(data=totals, aes(x=name, y=total_runtime, group=name)) +
    geom_boxplot(width = 0.25) +
    geom_violin(width = 1, draw_quantiles = c(0.25,0.75),
                linetype="dashed") +
    geom_violin(width=1, fill="transparent", draw_quantiles = 0.5) +
    stat_summary(fun = "mean", geom = "point",
                 shape = 1, size = 2) +
    theme_PEPATAC() +
    labs(x="", y="Total runtime (s)")
svg("gold_runtime_plot.svg")
g
dev.off()


