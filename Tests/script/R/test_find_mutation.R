## Test prior
options(digits=15)

prior0<- c(0,0,0,0)
freq<- c(0.3,0.2,0.2,0.3)
theta <- 0.001
tempComb<- combn(4,2)
allComb<-cbind(c(1,1), tempComb[,1:3], c(2,2), tempComb[,4:5], c(3,3), tempComb[,6], c(4,4))
result<- vector(mode="list", length=5)
for(i in 1:5){
    prior<- prior0
    if(i<5){
        prior[i]<- prior[i]+1
    }
    alpha <- freq*theta + prior
    alpha_sum <- sum(alpha)
    scale<- alpha_sum * (1+alpha_sum)
    result[[i]]<- apply(allComb, 2, function(x){
        if(x[1]==x[2]){
            return(alpha[x[1]]*(1+alpha[x[1]])/scale)
        }
        else{
            return(2*alpha[x[1]]*alpha[x[2]]/scale)
        }
    })

}

s<- sapply(result, function(x){
    paste(x, sep="", collapse=", ")
})
cat("{", paste(s, collapse="}, \n{"), "}\n" )

