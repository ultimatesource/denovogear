options(digits=10)

## Test F81

mu6 <- 1e-6
mu3 <- 1e-3
equal_freq <- rep(0.25,4)
unequal_freq <- seq(0.1, 0.4, by=0.1)
crazy_freq <- c(0.01, 0.1, 0.19, 0.7)

f81<- function(mu, freq){
    beta <- 1/(1- sum(freq^2) )
    exp_beta <- exp(-beta*mu)
    diag<- exp_beta+ freq*(1- exp_beta )
    off_diag<- freq*(1-exp_beta)
    return( list(diag=diag, off_diag=off_diag))
}


f81( mu6, equal_freq)
f81( mu6, unequal_freq)
f81( mu3, crazy_freq)

### Result
# > f81( mu6, equal_freq)
# $diag
# [1] 0.999999 0.999999 0.999999 0.999999
# 
# $off_diag
# [1] 3.333331111e-07 3.333331111e-07 3.333331111e-07 3.333331111e-07
# 
# > f81( mu6, unequal_freq)
# $diag
# [1] 0.9999987143 0.9999988571 0.9999990000 0.9999991429
# 
# $off_diag
# [1] 1.428570408e-07 2.857140816e-07 4.285711225e-07 5.714281633e-07
# 
# > f81( mu3, crazy_freq)
# $diag
# [1] 0.9978677587 0.9980615989 0.9982554390 0.9993538663
# 
# $off_diag
# [1] 0.0000215377905 0.0002153779050 0.0004092180195 0.0015076453352
