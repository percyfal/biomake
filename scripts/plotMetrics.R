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
cpal.out <- c("grey", rep(brewer.pal(10, "Paired"), 10))

# Function definitions

label.outliers <- function(x, x.labels, f=rep(1,nrow(x))) {
  ## Make sure x is an array (we might get a vector)
  x <- cbind(x) 
  ## Check arguments
  stopifnot(nrow(x) == length(x.labels))
  stopifnot(nrow(x) == length(f))
  ## Identify outliers for each column in x and factor level in f
  outlier <- apply(x, 2, function(x) {
    min.bound <- max.bound <- rep(NA, length(x)) 
    for(i in split(1:length(x), f)) {
      box.stats <- boxplot.stats(x[i])$stats
      min.bound[i] <- box.stats[1]
      max.bound[i] <- box.stats[5]
    }
    return(x < min.bound | x > max.bound)
  })
  outlier <- apply(outlier, 1, any)
  ## Set labels accordingly
  x.labels <- as.character(x.labels)
  x.labels[ !outlier ] <- "Not outlier"
  return(factor(x.labels, levels=unique(x.labels)))
}   

dup <- function(y)  {
    x <- y[is.na(y$FLOWCELL),]
    x <- x[order(x$SAMPLE_ID),]

    print(
      stripplot(100 * PERCENT_DUPLICATION ~ SAMPLE_ID, data=x,
                scales=list(x=list(draw=FALSE)), ylab="Percent duplication", xlab="Sample", main="",
                par.settings=simpleTheme(pch=19, col=cpal.out),
                auto.key=list(space="right"),
                groups = label.outliers(x$PERCENT_DUPLICATION, x$SAMPLE_ID)
                )
      )
   
}

insert <- function(y) {

}

hs <- function(y) {
    x <- y[is.na(y$FLOWCELL),]
    x <- x[order(x$SAMPLE_ID),]

    ## Percent missed targets
    print(
        stripplot(100*ZERO_CVG_TARGETS_PCT ~ SAMPLE_ID, data=x,
                  scales=list(x=list(draw=FALSE)),
                  ylab="Missed target regions (%)", xlab="Sample", main="Percent missed targets\n(missed = coverage < 2 at each base)",
                  par.settings=simpleTheme(pch=19, col=cpal.out),
                  auto.key=list(space="right"),
                  groups = label.outliers(x$ZERO_CVG_TARGETS_PCT, x$SAMPLE_ID)
                  )
        )   

    ## Mean target coverage
    print(
        stripplot(MEAN_TARGET_COVERAGE ~ SAMPLE_ID, data=x,
                  scales=list(x=list(draw=FALSE)),
                  ylab="Mean target coverage", xlab="Sample", main="Mean coverage (excluding missed targets)",
                  par.settings=simpleTheme(pch=19, col=cpal.out),
                  auto.key=list(space="right"),
                  groups = label.outliers(x$MEAN_TARGET_COVERAGE, x$SAMPLE_ID)
                  )
        )   

    ## Fold enrichment
    print(
        stripplot(FOLD_ENRICHMENT ~ SAMPLE_ID, data=x,
                  scales=list(x=list(draw=FALSE)),
                  ylab="Fold enrichment", xlab="Sample", main="Fold enrichment for baited regions over genomic background",
                  par.settings=simpleTheme(pch=19, col=cpal.out),
                  auto.key=list(space="right"),
                  groups = label.outliers(x$FOLD_ENRICHMENT, x$SAMPLE_ID)
                  )
        )   

    ## Percent usable bases on target (out of all sequenced bases)
    print(
        stripplot(100 * PCT_USABLE_BASES_ON_TARGET ~ SAMPLE_ID, data=x,
                  scales=list(x=list(draw=FALSE)),
                  ylab="Percent of sequenced bases", xlab="Sample", main="Percent usable bases\n(usable = matching target and, if applicable, de-duped)",
                  par.settings=simpleTheme(pch=19, col=cpal.out),
                  auto.key=list(space="right"),
                  groups = label.outliers(x$PCT_USABLE_BASES_ON_TARGET, x$SAMPLE_ID)
                  )
        )   

    ## Percent usable bases on target (out of all sequenced bases)
    print(
        stripplot(100 * ON_TARGET_BASES / PF_UQ_BASES_ALIGNED ~ SAMPLE_ID, data=x,
                  scales=list(x=list(draw=FALSE)),
                  ylab="Percent of mapped bases", xlab="Sample", main="Percent of mapped bases on target",
                  par.settings=simpleTheme(pch=19, col=cpal.out),
                  auto.key=list(space="right"),
                  groups = label.outliers(x$ON_TARGET_BASES / x$PF_UQ_BASES_ALIGNED, x$SAMPLE_ID)
                  )
        )   

    ## Target coverage as a function of read count
    print(
        xyplot(100*PCT_TARGET_BASES_10X ~ TOTAL_READS/1e6, data=x,
               ylab="Percent of targeted bases with 10X coverage", xlab="Sequenced reads (millions)",
               main="Total sequenced reads versus target coverage",
               par.settings=simpleTheme(pch=19, col=cpal.out),
               auto.key=list(space="right"),
               groups = label.outliers(cbind(PCT_TARGET_BASES_10X, x$TOTAL_READS), x$SAMPLE_ID)
               )
        )
    
    ## Arrange data for basewise coverage plots
    cov.levels <- c("2X", "10X", "20X", "30X", "40X", "50X", "100X")
    i <- paste("PCT_TARGET_BASES", cov.levels, sep="_")
    x.st <- stack(x[,i])
    x.st <- cbind(x.st, SAMPLE_ID=x$SAMPLE_ID)
    x.st$ind <- factor(x.st$ind, levels=i)
    levels(x.st$ind) <- cov.levels
    ylim <- c(100 * min(x.st$values) - 5, 100)
    
    ## Coverage boxplot over all samples, showing several coverage levels (2X, 10X etc)
    ## Note: the "groups" argument was previously used in the bwplot call below, but it seems bwplot does not properly support this.
    print(
        bwplot(100*values ~ ind , data=x.st,
               par.settings=simpleTheme(pch=19, col=cpal.out),
               xlab="Coverage", ylab="Percent of targeted bases", main="Coverage of targeted bases",
               scales=(list(x=list(rot=45))), ylim=ylim,
               )
        )

    ## Detailed coverage with outliers named. One plot per level (2X, 10X etc).
    for(cov.level in c("2X", "10X", "20X", "30X")) {
      x.st.subset <- x.st[x.st$ind == cov.level, ]
      print(
        stripplot(100 * values ~ SAMPLE_ID, data=x.st.subset,
                  scales=list(x=list(draw=FALSE)), ylim=ylim,
                  xlab="Sample", ylab="Percent of targeted bases",
                  main=paste("Percentage bases with", cov.level, "coverage"),
                  par.settings=simpleTheme(pch=19, col=cpal.out),
                  auto.key=list(space="right"),
                  groups = label.outliers(x.st.subset$values, x.st.subset$SAMPLE_ID))
      )
    }   
    
    ## Coverage, one panel per sample
    print(
      stripplot(100*values ~ ind | SAMPLE_ID, data=x.st,
                scales=list(x=list(rot=90, cex=0.7)), ylim=ylim,
                main="Percentage bases with a given coverage",
                xlab="Coverage", ylab="Percent of targeted bases",
                par.settings=simpleTheme(pch=19), as.table=TRUE,
                layout=c(6,3), par.strip.text=list(cex=0.5)
                )
      )

}

align <- function(y) {

    ## Fix order of category factor levels
    stopifnot(all(levels(y$CATEGORY) == c("FIRST_OF_PAIR", "PAIR", "SECOND_OF_PAIR")))
    levels(y$CATEGORY) <- c("Read 1", "Pair", "Read 2") # Rename levels
    y$CATEGORY <- factor(y$CATEGORY, levels=c("Read 1", "Read 2", "Pair")) # Reorder levels
  
    ## Sample-based plots
    x <- y[is.na(y$FLOWCELL) & y$CATEGORY=="Pair",]
    x <- x[order(x$SAMPLE_ID),]

    print(
        stripplot(TOTAL_READS/1e6 ~ SAMPLE_ID, data=x,
                  scales=list(x=list(rot=45, draw=FALSE)),
                  ylab="Sequenced reads (millions)", xlab="Sample", main="Reads per sample",
                  par.settings=simpleTheme(pch=19, col=cpal.out),
                  auto.key=list(space="right"),
                  groups = label.outliers(x$TOTAL_READS, x$SAMPLE_ID)
                  )
        )
    
    print(
        stripplot(TOTAL_READS/1e6 ~ SAMPLE_ID, data=x,
                  scales=list(x=list(rot=45, draw=FALSE), y=list(log=TRUE, equispaced.log=FALSE)),
                  ylab="Sequenced reads (millions)", xlab="Sample", main="Reads per sample (log scale)",
                  par.settings=simpleTheme(pch=19, col=cpal.out),
                  auto.key=list(space="right"),
                  groups = label.outliers(log10(x$TOTAL_READS), x$SAMPLE_ID)
                  )
        )
   
    print(
        xyplot(100*PCT_PF_READS_ALIGNED ~ TOTAL_READS/1e6, data=x,
               scales=list(x=list(log=TRUE, equispaced.log=FALSE), y=list(rot=45)),
               xlab="Sequenced reads (millions)", ylab="Aligned reads (%)",
               main="Total sequenced reads versus percent aligned",
               par.settings=simpleTheme(pch=19, col=cpal.out),
               auto.key=list(space="right"),
               groups = label.outliers(cbind(PCT_PF_READS_ALIGNED, log10(x$TOTAL_READS)), x$SAMPLE_ID)
               )
        )
    
    ## Based on all categories
    x <- y[is.na(y$FLOWCELL),]
    x <- x[order(x$SAMPLE_ID),]

    print(
        stripplot(100 * PCT_PF_READS_ALIGNED ~ SAMPLE_ID | CATEGORY, data=x,
                  auto.key=list(space="bottom", columns=3), scales=list(x=list(rot=45, draw=FALSE)),
                  ylab="Aligned reads (%)", xlab="", main="Percent aligned reads per sample",
                  par.settings=simpleTheme(pch=19, col=cpal.out), as.table=TRUE, layout=c(3,1),
                  groups=label.outliers(x$PCT_PF_READS_ALIGNED, x$SAMPLE_ID, f=x$CATEGORY),
                  )
        )

    print(
        stripplot(TOTAL_READS/1e6 + PF_READS_ALIGNED/1e6 ~ CATEGORY | SAMPLE_ID, data=x,
                  auto.key=list(text=c("Sequenced", "Aligned")),
                  scales=list(x=list(rot=45), y=list(log=TRUE, equispaced.log=FALSE)),
                  ylab="Reads (millions)", xlab="Category", main="Mapping statistics, reads per sample",
                  par.settings=simpleTheme(col=cpal[1:2], pch=21), as.table=TRUE,
                  layout=c(6,3), par.strip.text=list(cex=0.5)
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
                  ylab="Aligned reads (percent)", xlab="Sample",
                  main="Mapping statistics, percentage, grouped by flowcell",
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

