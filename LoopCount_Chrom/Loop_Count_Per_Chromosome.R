#!/usr/bin/env Rscript

#============================
# this script counts loops per chromosome
# and plots 
# parameters:
# 1) input loop file (assuming first 6 columns denote the interacting bins)
# 2) output plot file (absolute path)
# 3) if specified, a boolean variable indicating whether the input file has header or not. Default = 1
#============================

library(ggpubr)

args <- commandArgs(TRUE)

inploopfile <- args[1]
outplotfile <- args[2]
if (length(args) > 2) {
	headerinp <- as.integer(args[3])
} else {
	headerinp <- 1
}

outdir <- dirname(outplotfile)
system(paste("mkdir -p", outdir))

tempfile <- paste0(outdir, '/temp_dump.bed')

if (headerinp == 1) {
	system(paste("awk \'{if (NR>1) {print $1}}\' ", inploopfile, " | sort -k1,1 | uniq -c | awk \'{print $2\"\t\"$1}\' - > ", tempfile))
} else {
	system(paste("cut -f1", inploopfile, " | sort -k1,1 | uniq -c | awk \'{print $2\"\t\"$1}\' - > ", tempfile))
}

tempdata <- read.table(tempfile, header=F, sep="\t", stringsAsFactors=F)
colnames(tempdata) <- c('chr', 'loopcnt')

tot_loop_cnt <- sum(tempdata$loopcnt)
fracloopcnt_vec <- (tempdata$loopcnt * 1.0) / tot_loop_cnt

tempdata <- cbind.data.frame(tempdata, fracloopcnt_vec)
colnames(tempdata) <- c('chr', 'loopcnt', 'fracloopcnt')

plotfile <- paste0(outdir, '/loopcount_per_chr.pdf')
currplot <- ggpubr::ggbarplot(tempdata, x="chr", y="loopcnt", color="chr", fill="chr", add.params=list(size = 0.01), title="chromsome specific loop count", xlab="chr", ylab="loopcnt", legend = "right") + rotate_x_text(45)
ggarrange(currplot, ncol = 1) %>% ggexport(filename = plotfile)

plotfile <- paste0(outdir, '/loopcount_fraction_per_chr.pdf')
currplot <- ggpubr::ggbarplot(tempdata, x="chr", y="fracloopcnt", color="chr", fill="chr", add.params=list(size = 0.01), title="chromsome specific loop count (fraction)", xlab="chr", ylab="fraction_loopcnt", legend = "right") + rotate_x_text(45)
ggarrange(currplot, ncol = 1) %>% ggexport(filename = plotfile)

system(paste("rm", tempfile))


