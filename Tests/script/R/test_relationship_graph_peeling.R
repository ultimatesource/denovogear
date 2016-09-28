source("mutation_function.R")
options(digits=10)

rd<-vector(length=3, mode="list")
genotype<-vector(length=3, mode="list")
rd[[1]]<-c(0, 1, 25, 29) #878
genotype[[1]]<-c(1.9681812586082493465e-19, 5.8848619632386626597e-17, 9.3761538576728770258e-11, 9.3067830049776481198e-11, 1.1750042113891158116e-16, 2.571990828847844964e-08, 2.7827281184883153803e-08, 1.817757848706868286e-07, 1, 8.3192969181890452131e-07)

rd[[2]]<-c(0, 0, 57, 0) #891
genotype[[2]]<-c(2.9503055407498715637e-18, 2.9503055407498715637e-18, 2.2779180406756631894e-07, 2.9503055407498715637e-18, 2.9503055407498715637e-18, 2.2779180406756631894e-07, 2.9503055407498715637e-18, 1, 2.2779180406756631894e-07, 2.9503055407498715637e-18)

rd[[3]]<-c(0, 0, 76, 1) #892
genotype[[3]]<-c(1.6296875623801300125e-19, 1.6296875623801300125e-19, 7.848024187308119499e-08, 4.8727658115165870234e-17, 1.6296875623801300125e-19, 7.848024187308119499e-08, 4.8727658115165870234e-17, 1, 2.152806634867116697e-05, 9.7292347474093002355e-17)

upper<-vector(length=2, mode="list")  #1 dad, 2 mom
upper[[1]]<- c(0.00014982019479770606772, 5.9910104887616127319e-08, 0.00029961043454296829892, 8.9865157331424197595e-08, 9.9870144847656117711e-05, 0.00019974028969531223542, 5.9910104887616127319e-08, 0.99880131862140886234, 0.00029961043454296829892, 0.00014982019479770606772)
upper[[2]]<- upper[[1]]

lower<-vector(length=2, mode="list") #1 dad, 2 mom

freq<-c(0.3, 0.2, 0.2, 0.3)
F81_mu3<- f81_full(3e-8,freq)
F81_mu2<- f81_full(2e-8,freq)

numMutation<- 0; #0, 1, -1, 2
mut_germ<- meiosis_diploid_matrix(F81_mu3, F81_mu3, numMutation)
mut_somatic<- mitosis_diploid_matrix(F81_mu2, numMutation)

#up_fast (mom)
lower[[2]]<- mut_somatic %*% genotype[[3]] 

#to_father_fast
paired_buffer<- mut_germ %*% genotype[[1]]
paired_buffer_matirx<- matrix(paired_buffer, nrow=10, ncol=10, byrow=T)
lower[[1]]<- paired_buffer_matirx %*% (upper[[2]]*lower[[2]])

#up (dad)
lower_root<- lower[[1]] * (mut_somatic %*% genotype[[2]]) 

# roots
root_product<- lower_root*upper[[1]]
ret<- log(sum(root_product))

cat("expected_lower_mom << ", paste(lower[[2]], collapse=", "), ";\n", sep="")
cat("expected_lower_dad << ", paste(lower[[1]], collapse=", "), ";\n", sep="")
cat("expected_lower_root << ", paste(lower_root, collapse=", "), ";\n", sep="")
cat("expected_root_product << ", paste(root_product, collapse=", "), ";\n", sep="")
cat("expected_ret = ", ret, ";\n", sep="");
