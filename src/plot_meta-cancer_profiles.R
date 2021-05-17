library(data.table)
library(lubridate)
library(ggplot2)
library(ggpubr)
library(dplyr)
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

# Take all the meta-cancer samples and 1) bin by total reads 2) plot average total runtime by those bins

filelist <- list.files(pattern="*.tsv")
profiles <- lapply(filelist, fread, fill=TRUE)
names(profiles) <- str_remove(filelist, "_profile.tsv")

sizes <- data.table(name=as.character(),
                    file_size=as.numeric(),
                    raw_reads=as.numeric())

for (sample in profiles) {
    s <- data.table(name=sample[1, sample],
                    file_size=sample[1, file_mb],
                    raw_reads=sample[1, raw_reads])
    sizes <- rbind(sizes, s)
}

# order by input file size
sizes  <- sizes[order(sizes$file_size),]
breaks <- c(0,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,
            5500,6000,6500,7000,7500,8000,8500,9000,9500,10000,20000)
bins   <- c("[0-500)", "[500-1000)", "[1000-1500)", "[1500-2000)",
            "[2000-2500)", "[2500-3000)","[3000-3500)", "[3500-4000)",
            "[4000-4500)", "[4500-5000)", "[5000-5500)", "[5500-6000)",
            "[6000-6500)", "[6500-7000)", "[7000-7500)", "[7500-8000)",
            "[8000-8500)", "[8500-9000)", "[9000-9500)", "[9500-10000)",
            "[10000-20000)")
group_tags <- cut(sizes$file_size, breaks=breaks, labels=bins,
                  include.lowest=TRUE, right=FALSE)
summary(group_tags)

bsizes <- sizes %>% 
    mutate(bin = case_when(
        file_size < 500 ~ bins[1],
        file_size >= 500 & file_size < 1000 ~ bins[2],
        file_size >= 1000 & file_size < 1500 ~ bins[3],
        file_size >= 1500 & file_size < 2000 ~ bins[4],
        file_size >= 2000 & file_size < 2500 ~ bins[5],
        file_size >= 2500 & file_size < 3000 ~ bins[6],
        file_size >= 3000 & file_size < 3500 ~ bins[7],
        file_size >= 3500 & file_size < 4000 ~ bins[8],
        file_size >= 4000 & file_size < 4500 ~ bins[9],
        file_size >= 4500 & file_size < 5000 ~ bins[10],
        file_size >= 5000 & file_size < 5500 ~ bins[11],
        file_size >= 5500 & file_size < 6000 ~ bins[12],
        file_size >= 6000 & file_size < 6500 ~ bins[13],
        file_size >= 6500 & file_size < 7000 ~ bins[14],
        file_size >= 7000 & file_size < 7500 ~ bins[15],
        file_size >= 7500 & file_size < 8000 ~ bins[16],
        file_size >= 8000 & file_size < 8500 ~ bins[17],
        file_size >= 8500 & file_size < 9000 ~ bins[18],
        file_size >= 9000 & file_size < 9500 ~ bins[19],
        file_size >= 9500 & file_size < 10000 ~ bins[20],
        file_size >= 10000 & file_size < 20000 ~ bins[21]
    ))

runtimes <- data.table(name=as.character(),
                       total_runtime=as.numeric())

for (sample in profiles) {
    tmp  <- unique(sample, by="hash", fromLast=T)
    times <- data.table(name=tmp[1, sample],
                        total_runtime=tmp[, sum(runtime_seconds)])
    runtimes <- rbind(runtimes, times)
}


dt <- merge(bsizes, runtimes, by="name")
dt <- dt[order(dt$file_size),]
dt$bin <- factor(dt$bin, levels=unique(dt$bin))

give.n <- function(x){
    return(c(y = median(x)*1.05, label = length(x))) 
}

g <- ggplot(dt, aes(x=bin, y=total_runtime)) +
    geom_boxplot() +
    stat_summary(fun.data = give.n, geom = "text", fun = median, col="gray40",
                 position = position_dodge(width = 0.75)) +
    geom_point(alpha=0.05) + 
    theme_PEPATAC() +
    labs(x="Input file size bins (mb)", y="Total runtime (s)")

g + geom_smooth(method = "lm", se=FALSE, formula = y ~ x,
                color="gray50", aes(group=1)) + stat_cor()

lm_eqn <- function(df, y, x){
    formula = as.formula(sprintf('%s ~ %s', y, x))
    m <- lm(formula, data=df);
    # formating the values into a summary string to print out
    # ~ give some space, but equal size and comma need to be quoted
    eq <- substitute(italic(target) == a + b %.% italic(input)*","~~italic(r)^2~"="~r2*","~~p~"="~italic(pvalue), 
                     list(target = y,
                          input = x,
                          a = format(as.vector(coef(m)[1]), digits = 2), 
                          b = format(as.vector(coef(m)[2]), digits = 2), 
                          r2 = format(summary(m)$r.squared, digits = 3),
                          # getting the pvalue is painful
                          pvalue = format(summary(m)$coefficients[2,'Pr(>|t|)'], digits=1)
                     )
    )
    as.character(as.expression(eq));                 
}

svg('meta-cancer_runtime_plot.svg')
g + geom_text(x=8.5,y=62000,label=lm_eqn(dt, 'total_runtime','file_size'),
              color='gray30',parse=T,check_overlap = TRUE)
    #geom_smooth(method = "lm", se=FALSE, formula = y ~ x,
    #            alpha=0.5, color='gray50', aes(group=1), size=0.1)
dev.off()

# Plot it in hours

dt[,total_runtime_hr:=total_runtime/3600]
svg('meta-cancer_runtime_hr_plot.svg')
ggplot(dt, aes(x=bin, y=total_runtime_hr)) +
    geom_boxplot() +
    stat_summary(fun.data = give.n, geom = "text", fun = median, col="gray40",
                 position = position_dodge(width = 0.75)) +
    geom_point(alpha=0.05) + 
    theme_PEPATAC() +
    labs(x="Input file size bins (mb)", y="Total runtime (h)") +
    geom_text(x=8.5,y=62000,label=lm_eqn(dt, 'total_runtime_hr','file_size'),
              color='gray30',parse=T,check_overlap = TRUE)
dev.off()


# Plot memory use
memory <- data.table(name=as.character(),
                     file_size=as.numeric(),
                     peak_memory=as.numeric())

for (sample in profiles) {
    tmp  <- unique(sample, by="hash", fromLast=T)
    times <- data.table(name=tmp[1, sample],
                        file_size=tmp[1, file_mb],
                        peak_memory=tmp[, max(as.numeric(mem), na.rm=TRUE)])
    memory <- rbind(memory, times)
}

memory <- memory %>% 
    mutate(bin = case_when(
        file_size < 500 ~ bins[1],
        file_size >= 500 & file_size < 1000 ~ bins[2],
        file_size >= 1000 & file_size < 1500 ~ bins[3],
        file_size >= 1500 & file_size < 2000 ~ bins[4],
        file_size >= 2000 & file_size < 2500 ~ bins[5],
        file_size >= 2500 & file_size < 3000 ~ bins[6],
        file_size >= 3000 & file_size < 3500 ~ bins[7],
        file_size >= 3500 & file_size < 4000 ~ bins[8],
        file_size >= 4000 & file_size < 4500 ~ bins[9],
        file_size >= 4500 & file_size < 5000 ~ bins[10],
        file_size >= 5000 & file_size < 5500 ~ bins[11],
        file_size >= 5500 & file_size < 6000 ~ bins[12],
        file_size >= 6000 & file_size < 6500 ~ bins[13],
        file_size >= 6500 & file_size < 7000 ~ bins[14],
        file_size >= 7000 & file_size < 7500 ~ bins[15],
        file_size >= 7500 & file_size < 8000 ~ bins[16],
        file_size >= 8000 & file_size < 8500 ~ bins[17],
        file_size >= 8500 & file_size < 9000 ~ bins[18],
        file_size >= 9000 & file_size < 9500 ~ bins[19],
        file_size >= 9500 & file_size < 10000 ~ bins[20],
        file_size >= 10000 & file_size < 20000 ~ bins[21]
    ))

memory <- memory[order(memory$file_size),]
memory$bin <- factor(memory$bin, levels=unique(memory$bin))

mem_plot <- ggplot(memory, aes(x=bin, y=peak_memory)) +
    geom_boxplot() +
    stat_summary(fun.data = give.n, geom = "text", fun = median, col="gray40",
                 position = position_dodge(width = 0.75)) +
    geom_point(alpha=0.05) + 
    theme_PEPATAC() +
    labs(x="Input file size bins (mb)", y="Peak memory (GB)") +
    geom_text(x=8.75,y=10.5,label=lm_eqn(memory, 'peak_memory','file_size'),
              color='gray30',parse=T,check_overlap = TRUE)
svg('meta-cancer_peak_memory_plot.svg')
mem_plot
dev.off()

# Interesting consequence of very small input files with <1000 aligned reads is that MACS2 peaks at 6GB of RAM use
# What function is typically the peak memory user

memory <- data.table(name=as.character(),
                     file_size=as.numeric(),
                     cmd=as.character(),
                     peak_memory=as.numeric())

for (sample in profiles) {
    tmp  <- unique(sample, by="hash", fromLast=T)
    times <- data.table(name=tmp[1, sample],
                        file_size=tmp[1, file_mb],
                        cmd=tmp[mem==max(as.numeric(mem), na.rm=TRUE),]$cmd,
                        peak_memory=tmp[, max(as.numeric(mem), na.rm=TRUE)])
    memory <- rbind(memory, times)
}

#peak_mem_cmd <- list()
peak_mem_cmd <- data.table(bin=as.character(),
                           cmd=as.character())
for (b in unique(memory$bin)) {
    #peak_mem_cmd[[b]] <- unique(memory[bin==b,]$cmd)
    tmp <- data.table(bin=b,
                      cmd=memory[bin==b & peak_memory==max(memory[bin==b,]$peak_memory),]$cmd)
    peak_mem_cmd <- rbind(peak_mem_cmd, tmp)
}

# Remove the path from the matching strings (will need to adjust for users paths)
peak_mem_cmd$cmd <- str_remove(peak_mem_cmd$cmd, "/scratch/jps3dp/tools/databio//pepatac/tools/")

mem_plot + stat_summary(geom = 'text', label = peak_mem_cmd$cmd,
                        fun=median, vjust = 0, hjust=-0.6, angle=90,
                        size=3, col='gray40')

library(gridExtra)
mytheme <- gridExtra::ttheme_default(
    core = list(padding = unit(c(2.5, 2.5), "mm")))
tbl <- tableGrob(peak_mem_cmd, theme = mytheme, rows = NULL)

grid.arrange(mem_plot,
             tbl,
             ncol = 2,
             as.table = TRUE,
             heights = c(4, 1))

fwrite(peak_mem_cmd, 'meta-cancer_peak_memory_cmd.csv', sep=',')

memory[peak_memory==max(memory$peak_memory, na.rm=T),]  # The outlier is bamQC.py at 10.77 GB
