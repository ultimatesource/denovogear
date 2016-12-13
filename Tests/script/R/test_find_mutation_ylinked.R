source("mutation_function.R")
options(digits=15)


## test_prior
theta<- 0.001
freq<- c(0.3, 0.2, 0.2, 0.3)
ref_weight<- 1

freq_sum<- sum(freq)
normalised_theta<- theta/freq_sum

genotype_prior<- vector(length=5,mode="list")

for(i in 1:5){
    alpha<- freq * normalised_theta
    if(i<5){
        alpha[i] <- alpha[i] + ref_weight
    }
    alpha_sum<- sum(alpha)
    genotype_prior[[i]] <- alpha/alpha_sum
}    


for(i in 1:5){
    cat("{", paste(genotype_prior[[i]], collapse=", "), "},", sep="", fill=T)
}



## y-linked M12
rd<-vector(length=3, mode="list")
genotype<-vector(length=3, mode="list")
rd[[1]]<-c(0, 1, 25, 29) #878
genotype[[1]]<-c(1.9681812586082493465e-19, 5.8848619632386626597e-17, 9.3761538576728770258e-11, 9.3067830049776481198e-11, 1.1750042113891158116e-16, 2.571990828847844964e-08, 2.7827281184883153803e-08, 1.817757848706868286e-07, 1, 8.3192969181890452131e-07)

rd[[2]]<-c(0, 0, 57, 0) #891
genotype[[2]]<-c(2.9503055407498715637e-18, 2.9503055407498715637e-18, 2.2779180406756631894e-07, 2.9503055407498715637e-18, 2.9503055407498715637e-18, 2.2779180406756631894e-07, 2.9503055407498715637e-18, 1, 2.2779180406756631894e-07, 2.9503055407498715637e-18)

rd[[3]]<-c(0, 0, 76, 1) #892
genotype[[3]]<-c(1.6296875623801300125e-19, 1.6296875623801300125e-19, 7.848024187308119499e-08, 4.8727658115165870234e-17, 1.6296875623801300125e-19, 7.848024187308119499e-08, 4.8727658115165870234e-17, 1, 2.152806634867116697e-05, 9.7292347474093002355e-17)


genotype4<- sapply(genotype, function(x){
    temp<-log(x[c(1,5,8,10)])
    temp<- temp-max(temp)
    return(exp(temp))
}, simplify=F)
upper<- vector(length=3, mode="list")
upper[[1]]<- c(0.0002997002997003, 0.0001998001998002, 0.999200799200799, 0.0002997002997003)
upper[[2]]<- upper[[1]]
upper[[3]]<- upper[[1]]

lower<- vector(length=10, mode="list")
lower[[6]]<- genotype4[[1]]  #L1   Node1
lower[[7]]<- genotype4[[3]]  #L3   Node7
lower[[8]]<- genotype4[[1]]  #L4   Node3
lower[[9]]<- genotype4[[1]]  #L10  Node6
lower[[10]]<- genotype4[[2]] #L11  Node9


catCppGenotype<- function(geno){
    l<- length(geno[[1]])
    s<- sapply(geno, function(x){
        paste(x, sep="", collapse=", ")
    })
    cat("std::vector<std::array<double, ", l ,">> expected_genotype {\n\t{", paste(s, collapse="}, \n\t{"), "}\n};\n", sep="" )
}



##
numMutation<- 0
numMutation<- -1
peel_forward_m12 <- function(numMutation){
    F81_mu1<- f81_full(1e-8,freq)
    F81_mu2<- f81_full(2e-8,freq)

    mut_somatic<- mitosis_haploid_matrix(F81_mu1, numMutation)
    mut_library<- mitosis_haploid_matrix(F81_mu2, numMutation)

    lower[[5]] <- mut_library %*% lower[[10]]
    lower[[3]] <- mut_somatic %*% lower[[5]] 
    lower[[3]] <- lower[[3]] * mut_library %*% lower[[9]]

    lower[[2]] <- mut_library %*% lower[[8]]


    lower[[4]] <- mut_library %*% lower[[7]]
    lower[[1]] <- mut_somatic %*% lower[[4]] 
    lower[[1]] <- lower[[1]] * mut_library %*% lower[[6]]

    total <- 0
    for(i in 1:3){
        total <- total + log(sum(lower[[i]] * upper[[i]]))
    }
    return(total)
}

log_nomut<- peel_forward_m12(0)
log_fullmut<- peel_forward_m12(-1)
ret<- -expm1(log_nomut - log_fullmut)


cat("double expected_log_nomut = ", log_nomut, ";\n", sep="")
cat("double expected_log_fullmut = ", log_fullmut, ";\n", sep="")
cat("double expected_mup = ", ret, ";\n", sep="")


catCppGenotype(lower)

