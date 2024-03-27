#!/usr/bin/env Rscript

# Suppress warnings
options(warn=-1)

args <- commandArgs(trailingOnly = TRUE)

# Check if sufficient arguments are provided
if(length(args) < 2) {
  stop("Usage: Rscript script_name.R <folder_path> <file_prefix>", call. = FALSE)
}
# Assign arguments to variables
folder <- args[1]
file_pref <- args[2]
# Load packages
library(ggplot2)
library(ggseqlogo)
library(gridExtra)



# plot for genome mapped reads' size distro
genome_sd <- read.table(file.path(folder, "tables",
                                  paste0(file_pref,
                                          "_genome_sd.txt")))
mapperc <- read.table(file.path(folder, "tables", paste0(
  file_pref, "_map_stats.log"
)))
p1 <- ggplot(genome_sd, aes(x = V2, y = V1/sum(V1)*100))+
  geom_bar(stat = 'identity', fill = 'orange')+
  ylab("% reads")+
  xlab(NULL)+
  scale_x_continuous(breaks = genome_sd$V2)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  ggtitle(paste0("Size distribution for genome mapped reads\n",
                 mapperc[1,2], " of reads mapped"))

# plot for TE mapped reads' size distro
te_sd_p <- read.table(file.path(folder, "tables",
                                paste0(file_pref, "_te_sense_sd.txt")))
te_sd_m <- read.table(file.path(folder, "tables", paste0(
  file_pref, "_te_asense_sd.txt"
)))
te_sd <- merge(te_sd_p, te_sd_m, by = 'V2', all = T)
colnames(te_sd) <- c("len", "sense", "asense")
p2 <- ggplot(te_sd, aes(x = len, y = sense/sum(c(sense,asense))*100))+
  geom_bar(stat = 'identity', fill = 'darkred')+
  geom_bar(stat = 'identity', aes(x = len, y = -asense/sum(c(sense,asense))*100),
           fill = 'darkblue')+
  ylab("% reads")+
  xlab(NULL)+
  scale_x_continuous(breaks = te_sd$len)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  ggtitle(paste0("Size distribution for TE mapped reads\n",
                 mapperc[2,2], " of reads mapped"))

# Nucleotide composition for +
pfm_sense <- as.matrix(read.table(file.path(folder, "tables", paste0(
  file_pref,  "_te_map_sense.pfm"
)),
                                  row.names = 1))
rownames(pfm_sense) <- c("A", "C", "G", "U")
p3 <- ggplot()+
  geom_logo(pfm_sense, seq_type = "RNA")+
  theme_logo()+
  ggtitle("Weblogo for TE sense mapped reads")

# Nucleotide composition for -
pfm_asense <- as.matrix(read.table(file.path(folder, "tables", paste0(
  file_pref, "_te_map_asense.pfm"
)),
                                  row.names = 1))
rownames(pfm_asense) <- c("A", "C", "G", "U")
p4 <- ggplot()+
  geom_logo(pfm_asense, seq_type = "RNA")+
  theme_logo()+
  ggtitle("Weblogo for TE antisense mapped reads")

# Ping pong signature profile
te_pp <- read.table(file.path(folder, "tables", paste0(
  file_pref, "_te_mapped.pp"
)))
z10 <- (te_pp[, 2][te_pp$V1 == 10] - mean(te_pp$V2))/sd(te_pp$V2)
p5 <- ggplot(te_pp, aes(x = V1, y = V2/sum(V2)*100))+
  geom_bar(stat = 'identity', fill = '#5D00CD')+
  ylab("% overlapping reads")+
  xlab(NULL)+
  scale_x_continuous(breaks = te_pp$V1)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  ggtitle(paste0("Ping-pong signature for TE mapped reads\n",
                 "Z10 = ", format(z10, digits = 3)))

pdf(file.path(folder,
              "plots",
              paste0(file_pref, "_QC_plot.pdf")), width = 10, height = 16)
grid.arrange(
  p1,
  p2,
  p3,
  p4,
  p5,
  ncol = 2
)
bye_little_message <- dev.off()
