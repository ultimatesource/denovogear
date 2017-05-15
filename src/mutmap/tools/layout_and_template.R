#!@RSCRIPT_COMMAND@ --vanilla

computeLayout <- function(pedFilePath) {

    library(kinship2, quietly=TRUE)
    library(jsonlite, quietly=TRUE)

    data <- read.table(pedFilePath, header=FALSE)

    ped <- pedigree(id=data[,2], dadid=data[,3], momid=data[,4], sex=data[,5],
                    famid=data[,1])

    #familyId <- data[1,1]
    #ped1 <- ped[familyId]
    ped1 <- ped['1463']
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
    return(jsonOut)
}

readFile <- function(filePath) {
    return(readChar(filePath, file.info(filePath)$size))
}

escapeAsciiForJavaScript <- function(string) {
    return(gsub('\n', '\\\\x0A', gsub('"', '\\\\x22', string)))
}

build <- function(sourceDirectory, pedFilePath, dngOutputFilePath,
                  outputFilePath) {

    pedFileText <- escapeAsciiForJavaScript(readFile(pedFilePath))
    bundleFilePath <- file.path(sourceDirectory, 'bundle.js')
    bundleFileText <- readFile(bundleFilePath)

    templateFilePath <- file.path(sourceDirectory, 'index.html')
    templateFileText <- readFile(templateFilePath)
    outputText <- sub("<!--BUNDLE_JS_PLACEHOLDER-->", bundleFileText,
                      templateFileText, fixed=TRUE)

    outputText <- sub("/*PEDIGREE_FILE_TEXT_PLACEHOLDER*/",
                      paste("var pedigreeFileText = \"", pedFileText, "\";",
                            sep=""),
                      outputText, fixed=TRUE)

    layoutData <- computeLayout(pedFilePath)
    outputText <- sub("/*LAYOUT_DATA_PLACEHOLDER*/",
                      paste("var layoutData = ", layoutData, ";", sep=""),
                      outputText, fixed=TRUE)

    dngOutputFileText <- escapeAsciiForJavaScript(readFile(dngOutputFilePath))
    outputText <- sub("/*DNG_VCF_DATA_PLACEHOLDER*/",
                      paste("var dngOutputFileText = \"", dngOutputFileText,
                            "\";", sep=""),
                      outputText, fixed=TRUE)

    write(outputText, outputFilePath)
}
