#!/usr/bin/env Rscript

library(kinship2)
library(jsonlite)

f <- file("stdin")
data <- read.table(f, header=FALSE)
ped <- pedigree(id=data[,2], dadid=data[,3], momid=data[,4], sex=data[,5],
                famid=data[,1])

ped1 <- ped['1']
outData <- list()
outData$layout <- align.pedigree(ped1)
pedigree <- list()
pedigree$id <- ped1$id
pedigree$famid <- ped1$famid
pedigree$sex <- ped1$sex
pedigree$findex <- ped1$findex
pedigree$mindex <- ped1$mindex
outData$pedigree <- pedigree

jsonOut <- (toJSON(outData, pretty=TRUE, digits=NA))
cat(jsonOut)
