#! /usr/bin/env Rscript



##### LOAD ARGUMENTPARSER #####
loadLibrary <- tryCatch (
    {
        suppressWarnings(suppressPackageStartupMessages(library(argparser)))
    },
    error=function(e) {
        message("Error: Install the \"argparser\"",
                " library before proceeding.")
        return(NULL)
    },
    warning=function(e) {
        message(e)
        return(TRUE)
    }
)
if (length(loadLibrary)!=0) {
    suppressWarnings(library(argparser))
} else {
    quit()
}

# Create a parser
p <- arg_parser("Produce ATAC-seq Prealignment Comparison Plots")

# Add command line arguments
p <- add_argument(p, arg="with", short="-y",
                  help="PEPATAC project_config.yaml with prealignments")
p <- add_argument(p, arg="without", short="-n",
                  help="PEPATAC project_config.yaml without prealignments")
p <- add_argument(p, arg="resultsw", short="-r",
                  help="Project results output subdirectory path with prealignments")
p <- add_argument(p, arg="resultswo", short="-s",
                  help="Project results output subdirectory path without prealignments")
p <- add_argument(p, arg="output", short="-o",
                  help="Project parent output directory path")
# Parse the command line arguments
argv <- parse_args(p)

###############################################################################
##### LOAD DEPENDENCIES #####
required_libraries <- c("pepr", "data.table", "ggplot2", "lubridate", "tidyr")
for (i in required_libraries) {
    loadLibrary <- tryCatch (
        {
            suppressPackageStartupMessages(
                suppressWarnings(library(i, character.only=TRUE)))
        },
        error=function(e) {
            message("Error: Install the \"", i,
                    "\" library before proceeding.")
            return(NULL)
        },
        warning=function(e) {
            message(e)
            return(1)
        }
    )
    if (length(loadLibrary)!=0) {
        suppressWarnings(library(i, character.only=TRUE))
    } else {
        quit()
    }
}

################################################################################

aggregateProfiles <- function(project, results_subdir, output_dir) {
    project_name <- pepr::config(project)$name

    # Create assets_summary file
    project_samples <- pepr::sampleTable(project)$sample_name
    missing_files   <- 0
    profile <- data.table(sample_name=character(),
                          pid=numeric(),
                          hash=character(),
                          cid=character(),
                          runtime=character(),
                          mem=numeric(),
                          cmd=character(),
                          lock=character())
    write(paste0("Creating profile summary..."), stdout())

    for (sample in project_samples) {
        sample_output_folder <- file.path(results_subdir, sample)
        sample_profile_file  <- file.path(sample_output_folder,
                                          "PEPATAC_profile.tsv")

        if (!file.exists(sample_profile_file)) {
            missing_files <- missing_files + 1
            next
        }

        t <- fread(sample_profile_file, header=TRUE, skip=1)
        if(ncol(t) == 7) {
            colnames(t) <- c('pid', 'hash', 'cid', 'runtime',
                             'mem', 'cmd', 'lock')
        } else {
            missing_files <- missing_files + 1
            next
        }
       
        t <- t[!duplicated(t[, c('pid', 'hash', 'cid', 'runtime',
                                 'mem', 'cmd', 'lock')],
               fromLast=TRUE),]
        t[,sample_name:=sample]
        profile = rbind(profile, t)
    }
    if (missing_files > 0) {
        warning(sprintf("Profile files missing for %s samples.", missing_files))
    }

    return(profile)
}

# Set the output directory
summary_dir <- suppressMessages(file.path(argv$output, "summary"))
# Produce output directory (if needed)
dir.create(summary_dir, showWarnings = FALSE)

project <- Project(file = argv$with)
profile <- aggregateProfiles(project, argv$resultsw, summary_dir)

# Plot runtime
total_yes   <- profile[, sum(period_to_seconds(hms(runtime))), by = sample_name]
totals2_yes <- separate(total_yes, sample_name, into=c("pct", "total"), sep="_")
totals2_yes$pct   <- ordered(totals2_yes$pct,
                             levels = c("0-10", "1-9", "2-8", "3-7", "4-6",
                                        "5-5", "6-4", "7-3", "8-2", "9-1",
                                        "10-0"))
totals2_yes$total <- ordered(totals2_yes$total,
                             levels = c("10M", "20M", "40M", "60M", "80M",
                                        "100M", "120M", "140M", "160M",
                                        "180M", "200M"))
p <- ggplot(totals2_yes, aes(x=pct, y=V1/1000, group=total)) +
 geom_line(aes(col=total)) +
 ylab('Seconds (K)') +
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
  )
pdf(paste0(summary_dir, "/prealignments_total_runtimes.pdf"))
print(p)
dev.off()

# Add data from run without prealignments
project   <- Project(file = argv$without)
profile   <- aggregateProfiles(project, argv$resultswo, output_dir)
totals_no   <- profile[, sum(period_to_seconds(hms(runtime))), by = sample_name]
totals2_no  <- separate(totals_no, sample_name, into=c("pct", "total"), sep="_")
totals2_no$pct   <- ordered(totals2_no$pct,
                            levels = c("0-10", "1-9", "2-8", "3-7", "4-6",
                                       "5-5", "6-4", "7-3", "8-2", "9-1",
                                       "10-0"))
totals2_no$total <- ordered(totals2_no$total,
                            levels = c("10M", "20M", "40M", "60M", "80M",
                                       "100M", "120M", "140M", "160M",
                                       "180M", "200M"))
p <- ggplot(totals2_no, aes(x=pct, y=V1/1000, group=total)) +
 geom_line(aes(col=total)) +
 ylab('Seconds (K)') +
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
  )
pdf(paste0(summary_dir, "/no_prealignments_total_runtimes.pdf"))
print(p)
dev.off()

# Combine prealignments with no-prealignments
totals2_yes[ ,prealignment := "yes"]
totals2_no[ ,prealignment := "no"]
total <- rbind(totals2_yes, totals2_no)
total$prealignment <- ordered(total$prealignment, levels = c("yes", "no"))

p <- ggplot(total, aes(x=pct, y=V1/1000, group=interaction(total, prealignment))) +
 geom_line(aes(col=total, linetype=prealignment)) +
 ylab('Seconds (K)') +
 geom_point(aes(col=total, shape=prealignment)) +
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
  )
pdf(paste0(summary_dir, "/total_runtimes.pdf"))
print(p)
dev.off()

# Log comparison
total_log <- data.table(
    pct=total[prealignment == "yes",]$pct,
    total=total[prealignment == "yes",]$total,
    log_diff=log(total[prealignment == "yes",]$V1/
                 total[prealignment == "no",]$V1)
)

p <- ggplot(total_log, aes(x=pct, y=log_diff)) +
 geom_boxplot(outlier.shape=NA, width=.5, alpha=0, aes(x=pct, y=log_diff)) +
 geom_violin(alpha=.25, col='gray', fill='#A4A4A4') + 
 geom_point(aes(color=total),
            position=position_jitterdodge(jitter.width=0.1),
            alpha=.75) +
 geom_hline(yintercept=0, col='red') + 
 ylab(expression(log(over("time w/", "time w/o")))) +
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
  )
pdf(paste0(summary_dir, "/log_total_runtimes.pdf"))
print(p)
dev.off()
