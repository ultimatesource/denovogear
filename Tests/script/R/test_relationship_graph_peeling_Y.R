source("mutation_function.R")
options(digits=15)

f10to4<- c(1,5,8,10)

rd<-vector(length=1, mode="list")
genotype<-vector(length=1, mode="list")

rd[[1]]<-c(0, 0, 57, 0) #891
genotype[[1]]<-c(2.9503055407498715637e-18, 2.9503055407498715637e-18, 2.2779180406756631894e-07, 2.9503055407498715637e-18, 2.9503055407498715637e-18, 2.2779180406756631894e-07, 2.9503055407498715637e-18, 1, 2.2779180406756631894e-07, 2.9503055407498715637e-18)[f10to4]


upper<-vector(length=2, mode="list")  #1 
upper[[1]]<- c(0.0002997002997003, 0.0001998001998002, 0.999200799200799, 0.0002997002997003)

##
F81_mu2<- f81_full(2e-8,freq)

##
numMutation<- 0
mut_somatic<- mitosis_haploid_matrix(F81_mu2, numMutation)

# upFast
lower_root<-  (mut_somatic %*% genotype[[1]]) 
# roots
root_product<- lower_root*upper[[1]]
nomut_ret<- log(sum(root_product))

## 
numMutation<- -1
mut_somatic<- mitosis_haploid_matrix(F81_mu2, numMutation)

# upFast
lower_root<-  (mut_somatic %*% genotype[[1]]) 
# roots
root_product<- lower_root*upper[[1]]
fullmut_ret<- log(sum(root_product))

mup = -(expm1(nomut_ret - fullmut_ret) )

# cat("expected_lower_dad << ", paste(lower[[1]], collapse=", "), ";\n", sep="")

cat("{", paste(lower_root, collapse=", "), "}\n", sep="")
# cat("{", paste(root_product, collapse=", "), "}\n", sep="")
cat("double expected_nomut_ret = ", nomut_ret, ";\n", sep="");
cat("double expected_full_ret = ", fullmut_ret, ";\n", sep="");
cat("double expected_mup = ", mup, ";\n", sep="");


















F81_mu2<- f81_full(2e-1,freq)

numMutation<- -1; #0, 1, -1, 2
mut_germ<- meiosis_diploid_matrix(F81_mu3, F81_mu3, numMutation)
mut_somatic<- 

mitosis_diploid_matrix(F81_mu2, numMutation)

mitosis_haploid_matrix(F81_mu2, numMutation)
f10to4<- c(1,5,8,10)


f1<- c(freq[1]+1,freq[-1])
f1/sum(f1)
theta<- 0.001
theta2<- theta/sum(freq)
f1<-theta2*freq; 
f1[1]<-f1[1]+1

alpha<-f1
alpha_sum<- sum(alpha)

g<- c(
alpha[1]*(1.0 + alpha[1]) / alpha_sum / (1.0 + alpha_sum),
2.0 * alpha[1]*(alpha[2]) / alpha_sum / (1.0 + alpha_sum), 
2.0 * alpha[1]*(alpha[3]) / alpha_sum / (1.0 + alpha_sum), 
2.0 * alpha[1]*(alpha[4]) / alpha_sum / (1.0 + alpha_sum), 
alpha[2]*(1.0 + alpha[2]) / alpha_sum / (1.0 + alpha_sum), 
2.0 * alpha[2]*(alpha[3]) / alpha_sum / (1.0 + alpha_sum), 
2.0 * alpha[2]*(alpha[4]) / alpha_sum / (1.0 + alpha_sum), 
alpha[3]*(1.0 + alpha[3]) / alpha_sum / (1.0 + alpha_sum), 
2.0 * alpha[3]*(alpha[4]) / alpha_sum / (1.0 + alpha_sum), 
alpha[4]*(1.0 + alpha[4]) / alpha_sum / (1.0 + alpha_sum))


g
sum(g)
g4<-g[f10to4]
g4/sum(g4)
alpha/alpha_sum


