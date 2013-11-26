#! /usr/bin/Rscript --vanilla
# File: plotNgstk.R
# Created: Tue Oct 29 16:26:25 2013
# $Id: $
#
# Copyright (C) 2013 by Per Unneberg
#
# Author: Per Unneberg
#
# Description:
#

library(utils)
library(fields)
library(RColorBrewer)
cpal <- colorRampPalette(brewer.pal(9,"Paired"))(1000)

args <- commandArgs(TRUE)
if (length(args) != 2) {
    message("Usage: plotNgstk.R infile outfile!")
    quit("yes")
}

ngstkfile <- args[1]
ngstkdata <- read.table(ngstkfile)
outfile <- args[2]

coverage = gsub(".*C([0-9]+).*", "\\1", ngstkfile)
fraction = gsub(".*F([0-9]+\\.[0-9]+).*", "\\1", ngstkfile)
popa = gsub("(^[A-Z]+)_.*", "\\1", ngstkfile)
popb = gsub("^[A-Z]+_([A-Z])_.*", "\\1", ngstkfile)

chr = gsub(".ngstk.txt", "", ngstkfile)

pdf(outfile)
tmp = as.matrix( log10( ngstkdata ) )
image.plot(tmp, col=cpal, zlim=c(0,max(tmp)), xlab=paste("allele frequency, pop", popa), ylab=paste("allele frequency, pop", popb),
           main=paste(chr, paste("coverage>=", coverage, sep=""), paste("fraction>=", fraction, sep=""), sep=", "))

dev.off()
