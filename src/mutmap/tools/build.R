#!@RSCRIPT_COMMAND@ --vanilla

source(file.path('@MUTMAP_DIR@', 'layout_and_template.R'))

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

