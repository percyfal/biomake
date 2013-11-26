#! /usr/bin/Rscript --vanilla
# File: plotThetas.R
# Created: Tue Oct 29 17:12:00 2013
# $Id: $
#
# Copyright (C) 2013 by Per Unneberg
#
# Author: Per Unneberg
#
# Description:
#

library(utils)
library(lattice)
library(RColorBrewer)
cpal <- colorRampPalette(brewer.pal(9,"Paired"))(1000)

args <- commandArgs(TRUE)
if (length(args) != 2) {
    message("Usage: plotThetas.R infile outfile!")
    quit("yes")
}

infile <- args[1]

dir <- basename(dirname(normalizePath(infile)))
if (grep("w[0-9]+_[0-9]+", dir)) {
    message("We have a window match")
    window <- gsub("w([0-9]+)_.*", "\\1", dir)
    step <- gsub(".*_([0-9]+)", "\\1", dir)
} else {
    window <- 50000
    step <- 50000
}
pop = gsub("(^[A-Z]+)_.*", "\\1", infile)

d <- read.table(infile)
names(d) <- c("range", "chr", "pos", "tW", "tP", "tF", "tH", "tL", "tajD", "fulif", "fuliD", "fayH", "zengsE", "numSites")
d$pop = pop
d.stack <- cbind(stack(d[,c("tajD", "fuliD", "fulif", "fayH", "zengsE", "tW", "tP")]), pop=d$pop, pos=d$pos, chr=d$chr)
outfile <- args[2]

pdf(outfile)
xyplot(values ~ pos/1e6 | ind, data=d.stack, type="l", layout=c(1,7), xlab="pos (Mb)", scales=list(y=(relation="free")), strip=FALSE, strip.left=TRUE, main=paste(basename(infile), ", w:", window, ", s:", step, sep=""))
dev.off()
