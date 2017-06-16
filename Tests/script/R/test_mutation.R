source("mutation_function.R")
options(digits=10)

## Test F81
mu6 <- 1e-6
mu3 <- 1e-3
equal_freq <- rep(0.25,4)
unequal_freq <- seq(0.1, 0.4, by=0.1)
crazy_freq <- c(0.01, 0.1, 0.19, 0.7)


f81_core( mu6, equal_freq)
f81_core( mu6, unequal_freq)
f81_core( mu3, crazy_freq)



catCppMatrix<-function(matrix, file=""){
    s<- paste("{",
        apply(matrix,1,function(y){
            paste(y, sep="", collapse=", ") 
        })
        ,"}", collapse=",\n", sep="")
    cat("{", paste(s, collapse="}, \n{"), "};\n", sep="", file=file )
}

catEigenMatrix<-function(matrix, file=""){
    s<- paste(
            apply(matrix,1,function(y){
                paste(y, sep="", collapse=", ") 
            }), collapse=",\n", sep="")
    cat("<< ", paste(s, collapse=", \n{"), ";\n", sep="", file=file )
}

## Test transition_matrix
F81<- f81_full(1e-6, unequal_freq)

catEigenMatrix(mitosis_haploid_matrix(F81, -1))
catEigenMatrix(mitosis_haploid_matrix(F81, 0))
catEigenMatrix(mitosis_haploid_matrix(F81, 1))


catEigenMatrix(mitosis_diploid_matrix(F81, -1))
catEigenMatrix(mitosis_diploid_matrix(F81, 0))
catEigenMatrix(mitosis_diploid_matrix(F81, 1))
catEigenMatrix(mitosis_diploid_matrix(F81, 2))


catEigenMatrix(meiosis_haploid_matrix(F81, -1))
catEigenMatrix(meiosis_haploid_matrix(F81, 0))
catEigenMatrix(meiosis_haploid_matrix(F81, 1))

catEigenMatrix(meiosis_diploid_matrix(F81, F81, -1), file="temp.txt")
catEigenMatrix(meiosis_diploid_matrix(F81, F81, 0), file="temp.txt")
catEigenMatrix(meiosis_diploid_matrix(F81, F81, 1), file="temp.txt")
catEigenMatrix(meiosis_diploid_matrix(F81, F81, 2), file="temp.txt")
