#!@RSCRIPT_COMMAND@ --vanilla

source(file.path('layout_and_template.R'))

# designed to be run from within the build/src/mutmap directory
test <- function() {
    sourceDirectory <- file.path('.')
    pedFilePath <- file.path('..', '..', 'testdata', 'human_trio',
                             'trio.ped')
    #pedFilePath <- file.path('..', '..', 'testdata', 'relationship_graph',
    #                         'relationship_graph.ped')
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

if (!interactive()) {
    testServer()
}
