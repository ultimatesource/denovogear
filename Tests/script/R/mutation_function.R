
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
unfold_index<-c(1, 2, 3, 4, 6, 7, 8, 11, 12, 16)
fold_index<-c(1, 2, 3, 4, 2, 5, 6, 7, 3, 6, 8, 9, 4, 7, 9, 10)


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
#     meiosis<- matrix(0,nrow=10, ncol=4)
#     index<- 1
#     for(i in 1:4){
#         for(j in i:4){
#             meiosis[index,]<-0.5*(mm[i,]+mm[j,])
#             index<- index+1
#         }
#     }
    meiosis<- 0.5*(mm[c(1,1,1,1,2,2,2,3,3,4),]+mm[c(1:4,2:4,3:4,4),])
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
    meiosis_diploid<- t(rowsum(t(temp), fold_index, reorder=F))
    return(meiosis_diploid)
}
