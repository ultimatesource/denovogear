
f81_core<- function(mu, freq){
    beta <- 1/(1- sum(freq^2) )
    exp_beta <- exp(-beta*mu)
    diag<- exp_beta+ freq*(1- exp_beta )
    off_diag<- freq*(1-exp_beta)
    return( list(diag=diag, off_diag=off_diag))
}

f81_full<- function(mu, freq){
    temp<- f81_core(mu, freq)
    m<- matrix(temp$off_diag, nrow=4, ncol=4, byrow=T)
    diag(m)<- temp$diag
    return(m)
}

mitotic_haploid_mutation_counts<- rbind(
    c(0, 1, 1, 1), c(1, 0, 1, 1), c(1, 1, 0, 1), c(1, 1, 1, 0))


mitosis_haploid_matrix<- function(mut, numMutation){
    mitosis<- matrix(0, nrow=4, ncol=4)
    if(numMutation<0){
        mitosis<- mut
    }
    else if(numMutation==0){
        diag(mitosis)<- diag(mut)
    }
    else{
        index<- numMutation==mitotic_haploid_mutation_counts
        for(i in 1:4){
            mitosis[i,index[i,]]<- mut[i,index[i,]]
        }
    }
    return(mitosis)

}

mitotic_diploid_mutation_counts<- rbind(
    c(0, 1, 1, 1, 1, 2, 2, 2, 1, 2, 2, 2, 1, 2, 2, 2),
    c(1, 0, 1, 1, 2, 1, 2, 2, 2, 1, 2, 2, 2, 1, 2, 2),
    c(1, 1, 0, 1, 2, 2, 1, 2, 2, 2, 1, 2, 2, 2, 1, 2),
    c(1, 1, 1, 0, 2, 2, 2, 1, 2, 2, 2, 1, 2, 2, 2, 1),
    c(1, 2, 2, 2, 0, 1, 1, 1, 1, 2, 2, 2, 1, 2, 2, 2),
    c(2, 1, 2, 2, 1, 0, 1, 1, 2, 1, 2, 2, 2, 1, 2, 2),
    c(2, 2, 1, 2, 1, 1, 0, 1, 2, 2, 1, 2, 2, 2, 1, 2),
    c(2, 2, 2, 1, 1, 1, 1, 0, 2, 2, 2, 1, 2, 2, 2, 1),
    c(1, 2, 2, 2, 1, 2, 2, 2, 0, 1, 1, 1, 1, 2, 2, 2),
    c(2, 1, 2, 2, 2, 1, 2, 2, 1, 0, 1, 1, 2, 1, 2, 2),
    c(2, 2, 1, 2, 2, 2, 1, 2, 1, 1, 0, 1, 2, 2, 1, 2),
    c(2, 2, 2, 1, 2, 2, 2, 1, 1, 1, 1, 0, 2, 2, 2, 1),
    c(1, 2, 2, 2, 1, 2, 2, 2, 1, 2, 2, 2, 0, 1, 1, 1),
    c(2, 1, 2, 2, 2, 1, 2, 2, 2, 1, 2, 2, 1, 0, 1, 1),
    c(2, 2, 1, 2, 2, 2, 1, 2, 2, 2, 1, 2, 1, 1, 0, 1),
    c(2, 2, 2, 1, 2, 2, 2, 1, 2, 2, 2, 1, 1, 1, 1, 0) )

# Lower Triangle
unfold_index = c()
for(a in 0:3) {
    for(b in 0:a) {
        unfold_index = c(unfold_index, 1+a*4+b)
    } 
}

fold_index = rep(-1,16)
for(i in seq_along(unfold_index)) {
    h = unfold_index[i]
    fold_index[h] = i
    a = (h-1)%/%4
    b = (h-1)%%4
    fold_index[1+b*4+a] = i
}

mitosis_diploid_matrix<- function(mut, numMutation){

    mitosis<- matrix(0, nrow=10, ncol=10)
    kroneckerMut<- kronecker(mut,mut)
    if(numMutation==0){
        diag(mitosis)<- diag(kroneckerMut)[unfold_index]
    }
    else {
        for(p in 1:10){
            i<- unfold_index[p]
            for(j in 1:16){
                q<- fold_index[j]
                if( numMutation<0 | numMutation==mitotic_diploid_mutation_counts[i,j]){
                    mitosis[p,q] = mitosis[p,q] + kroneckerMut[i,j]
                }
            
            }
        }
    }
    return(mitosis)

}

meiosis_haploid_matrix<- function(mut, numMutation=0){
    mm<- mitosis_haploid_matrix(mut, numMutation)
    meiosis<- matrix(0,nrow=10, ncol=4)
    for(i in 1:10) {
        h = unfold_index[i]
        a = (h-1)%/%4
        b = (h-1)%%4
        for(j in 1:4) {
            meiosis[i,j] = 0.5*(mm[a+1,j]+mm[b+1,j])
        }
    }
    return(meiosis)
}

meiosis_diploid_matrix<- function(dad, mom, numMutation){
    if(numMutation <=0){
        meiosis_dad<- meiosis_haploid_matrix(dad, numMutation)
        meiosis_mom<- meiosis_haploid_matrix(mom, numMutation)
        temp<- kronecker(meiosis_dad, meiosis_mom)
    }
    else{
        temp<- matrix(0, nrow=100, ncol=16)
        for(i in 0:(numMutation)){
            meiosis_dad<- meiosis_haploid_matrix(dad, i)
            meiosis_mom<- meiosis_haploid_matrix(mom, numMutation-i)
            temp<- temp + kronecker(meiosis_dad, meiosis_mom)
        }
    }


#     meiosis_diploid<- matrix(0,nrow=100, ncol=10)
#     for(i in 1:16){
#         meiosis_diploid[,fold_index[i]] <- meiosis_diploid[,fold_index[i]]+temp[,i]
#     }
    meiosis_diploid<- t(rowsum(t(temp), fold_index, reorder=TRUE))
    return(meiosis_diploid)
}
