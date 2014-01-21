#! /usr/bin/Rscript --vanilla
# File: plotNgstk.R
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
cpal.out <- colorRampPalette(brewer.pal(9,"Paired"))(10)

# Function definitions
dup <- function(y)  {
    x <- y[is.na(y$FLOWCELL),]
    x <- x[order(x$SAMPLE_ID),]
    print(
        stripplot(100 * PERCENT_DUPLICATION ~ SAMPLE_ID, data=x,
                  scales=(list(x=list(rot=45))), par.settings=simpleTheme(pch=19),
                  ylim=c(-1,100), xlab="Sample", ylab="Percent duplication")
        )
}

insert <- function(y) {

}

hs <- function(y) {
    x <- y[is.na(y$FLOWCELL),]
    x <- x[order(x$SAMPLE_ID),]
    i <- c("ZERO_CVG_TARGETS_PCT", "PCT_TARGET_BASES_2X", "PCT_TARGET_BASES_10X", "PCT_TARGET_BASES_20X", "PCT_TARGET_BASES_30X")
    ## First make outlier plot
    x.st <- stack(x[,i])
    x.st <- cbind(x.st, SAMPLE_ID=x$SAMPLE_ID)
    x.st$ind <- factor(x.st$ind, levels=levels(x.st$ind)[c(5,3,1,2,4)])
    levels(x.st$ind) <- c("0X", "2X", "10X", "20X", "30X")
    g <- (100*x.st$values) %in% boxplot(100*values ~ ind, data=x.st, plot=FALSE)$out
    qc <- rep("Not outlier", length(x.st$SAMPLE_ID))
    qc[g] <- as.character(x.st$SAMPLE_ID[g])

    print(
        bwplot(100*values ~ ind , data=x.st,
               par.settings=simpleTheme(pch=19, col=cpal.out),
               xlab="Coverage", ylab="Percent on target", main="Percentage bases with a given coverage.",
               auto.key=list(space="right"),
               scales=(list(x=list(rot=45))),
               groups = qc)
        )

    ## Then plot per sample
    PLOTS_PER_PANEL = 9
    groups <- cut(1:dim(x)[1], ceiling(dim(x)[1]/PLOTS_PER_PANEL), labels=FALSE)
    print(
        by(x, groups, function(z){
            x.st <- stack(z[,i])
            x.st <- cbind(x.st, SAMPLE_ID=z$SAMPLE_ID)
            x.st$ind <- factor(x.st$ind, levels=levels(x.st$ind)[c(5,3,1,2,4)])
            levels(x.st$ind) <- c("0X", "2X", "10X", "20X", "30X")
            stripplot(100*values ~ ind | SAMPLE_ID, data=x.st,
                      scales=list(x=list(rot=45)),
                      main="Percentage bases with a given coverage.",
                      xlab="Coverage", ylab="Percentage bases",
                      ylim=c(-1,100), par.settings=simpleTheme(pch=19)
                      )
        })
        )
}

align <- function(y) {
    # Sample-based plots
    x <- y[is.na(y$FLOWCELL) & y$CATEGORY=="PAIR",]
    x <- x[order(x$SAMPLE_ID),]
    g <- log10(x$TOTAL_READS) %in% boxplot(log10(x$TOTAL_READS), plot=FALSE)$out
    qc <- rep("Not outlier", length(x$SAMPLE_ID))
    qc[g] <- as.character(x$SAMPLE_ID[g])

    print(
        stripplot(log10(TOTAL_READS) ~ SAMPLE_ID, data=x,
                  scales=list(x=list(rot=45, draw=FALSE)), ylab="Reads (log10 millions)", xlab="Sample", main="Reads per sample. Outliers indicate low read count.",
                  par.settings=simpleTheme(pch=19, col=cpal.out),
                  auto.key=list(space="right"),
                  groups = qc,
                  )
        )

    g1 <- log10(x$TOTAL_READS) %in% boxplot(log10(x$TOTAL_READS), plot=FALSE)$out
    g2 <- (100*x$PCT_PF_READS_ALIGNED) %in% boxplot(100*(x$PCT_PF_READS_ALIGNED), plot=FALSE)$out
    qc <- rep("Not outlier", length(x$SAMPLE_ID))
    qc[g1 | g2] <- as.character(x$SAMPLE_ID[g1 | g2])

    print(
        xyplot(100*PCT_PF_READS_ALIGNED ~ log10(TOTAL_READS), data=x,
               scales=list(y=list(rot=45)), xlab="Reads (log millions)", ylab="% reads aligned", main="Aligned reads. Outlier indicate low read count and alignment rate.",
               par.settings=simpleTheme(pch=19, col=cpal.out),
               auto.key=list(space="right"),
               groups = qc,
               )
        )

    ## Based on all categories
    x <- y[is.na(y$FLOWCELL),]
    x <- x[order(x$SAMPLE_ID),]

    print(
        stripplot(log10(TOTAL_READS) + log10(PF_READS_ALIGNED)  ~ CATEGORY | SAMPLE_ID, data=x,
                  auto.key=list(text=c("Reads", "Aligned")),
                  scales=list(x=list(rot=45)),
                  ylab="Reads (log10 millions)", xlab="Category", main="Mapping statistics, reads per sample",
                  par.settings=simpleTheme(col=cpal[1:2], pch=21)
                  )
        )
        
    g <- (100*x$PCT_PF_READS_ALIGNED) %in% boxplot(100*(x$PCT_PF_READS_ALIGNED), plot=FALSE)$out
    qc <- rep("Not outlier", length(x$SAMPLE_ID))
    qc[g] <- as.character(x$SAMPLE_ID[g])

    print(
        stripplot(100 * PCT_PF_READS_ALIGNED ~ SAMPLE_ID | CATEGORY, groups=qc,
                  data=x, auto.key=list(space="right"), scales=list(x=list(rot=45, draw=FALSE)),
                  ylab="Aligned reads (percent).", xlab="Sample", main="Mapping statistics, percentage. Outliers indicate low alignment rate.",
                  par.settings=simpleTheme(pch=19, col=cpal.out)
                  )
        )

    # Flowcell-based plots
    x <- y[!is.na(y$FLOWCELL),]
    x <- x[order(x$SAMPLE_ID),]

    if(nrow(x) > 0) {
      print(
        stripplot(100 * PCT_PF_READS_ALIGNED ~ SAMPLE_ID | CATEGORY,
                  groups=FLOWCELL, data=x, auto.key=list(space="right"),
                  scales=list(x=list(rot=45, draw=FALSE)),
                  ylab="Aligned reads (percent)", xlab="Sample", main="Mapping statistics, percentage, grouped by flowcell",
                  par.settings=simpleTheme(pch=19, col=cpal)
                  )
        )
    }
}

args <- commandArgs(TRUE)
if (length(args) != 3) {
    message("Usage: plotMetrics.R infile outfile type!")
    quit("yes")
}

infile <- args[1]
outfile <- args[2]
type <- args[3]

x <- read.table(infile, header=TRUE, sep="\t")

pdf(outfile)
type <- match.arg(type, choices = c("align", "insert", "dup", "hs"))
switch(type,
       align = align(x),
       insert = insert(x),
       hs = hs(x),
       dup = dup(x)
       )
dev.off()

