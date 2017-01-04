#!@RSCRIPT_COMMAND@ --vanilla

argErrorAndExit <- function(errorMessage) {
    message <- c(errorMessage, ". Use --help to get usage information")
    stop(message)
}

parseArguments <- function() {

    library(argparser, quietly=TRUE)

    argParser <- arg_parser("Visualization for DeNovoGear")
    argParser <- add_argument(argParser, "--pedFilePath", short='-p',
                              help="Path to pedigree file", type="character")
    argParser <- add_argument(argParser, "--dngOutputFilePath", short='-d',
                              help="DeNovoGear output VCF", type="character")
    argParser <- add_argument(argParser, "--outputFilePath", short='-o',
                              help="Where to write the output",
                              type="character")
    argv <- parse_args(argParser)

    if (is.na(argv$pedFilePath)) {
        argErrorAndExit("Pedigree file path not provided")
    }
    else if (is.na(argv$dngOutputFilePath)) {
        argErrorAndExit("DeNovoGear output VCF file path not provided")
    }
    else if (is.na(argv$outputFilePath)) {
        argErrorAndExit("Output file path not provided")
    }
    return(argv)
}

computeLayout <- function(pedFilePath) {

    library(kinship2, quietly=TRUE)
    library(jsonlite, quietly=TRUE)

    data <- read.table(pedFilePath, header=FALSE)

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

# designed to be run from within the build/src/mutmap directory
test <- function() {
    sourceDirectory <- file.path('.')
    pedFilePath <- file.path('..', '..', 'testdata', 'relationship_graph',
                             'relationship_graph.ped')
    dngOutputFilePath <- file.path('..', '..', 'dng_test_output.vcf')
    outputFilePath <- file.path('mutmap.html')

    build(sourceDirectory, pedFilePath, dngOutputFilePath, outputFilePath)
}

testServer <- function() {

    # use a socket as a simple IPC signalling mechanism. When a connection is
    # initiated run the test function
    while(TRUE){
        con <- socketConnection(host="localhost", port = 6011, blocking=TRUE,
                                server=TRUE, open="r+", timeout=86400)
        data <- readLines(con, 1)
        print("Signal received")

        test()

        response <- "Test build completed successfully"
        writeLines(response, con) 

        close(con)
    }
}

main <- function() {
    argv <- parseArguments()
    sourceDirectory <- file.path('@MUTMAP_DIR@')
    pedFilePath <- argv$pedFilePath
    dngOutputFilePath <- argv$dngOutputFilePath
    outputFilePath <- argv$outputFilePath

    build(sourceDirectory, pedFilePath, dngOutputFilePath, outputFilePath)
}

if (!interactive()) {
    main()
}

