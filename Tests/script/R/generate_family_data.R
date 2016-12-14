options(digits=15)

set.seed(0)
dip_count<- 10

# family order: dad, mom, child1, child2
lower<- sapply(1:4, function(x){runif(dip_count)}, simplify=F)
upper<- sapply(1:4, function(x){runif(dip_count)}, simplify=F)

transition<- sapply(1:2, function(x){ matrix(runif(dip_count^2), nrow=dip_count, byrow=T) }, simplify=F)
transition[[3]]<- matrix(runif(dip_count^3), nrow=dip_count*dip_count, byrow=T)
transition[[4]]<- matrix(runif(dip_count^3), nrow=dip_count*dip_count, byrow=T)

## up
up_fast <- transition[[2]] %*% lower[[2]]
up<- lower[[1]] * up_fast

## to father
buffer <- transition[[3]] %*% lower[[3]]
buffer_matrix <- matrix(buffer,nrow=dip_count, byrow=T)
mom_lower_upper <- lower[[2]] * upper[[2]]
father_fast <- buffer_matrix %*% mom_lower_upper
father <- lower[[1]] * father_fast

## to mother
buffer <- transition[[3]] %*% lower[[3]]
buffer_matrix <- matrix(buffer,nrow=dip_count, byrow=F)##by column here
dad_lower_upper <- lower[[1]] * upper[[1]]
mother_fast <- buffer_matrix %*% dad_lower_upper
mother <- lower[[2]] * mother_fast

## to child
buffer<- (lower[[1]]*upper[[1]]) %x% ((lower[[2]]*upper[[2]]))
other_child <- transition[[4]] %*% lower[[4]]
updated_buffer <- buffer * as.vector(other_child)

child_fast<- t(transition[[3]]) %*% buffer
child<- t(transition[[3]]) %*% updated_buffer

## output to c++
formatList<-function(name,v){
    paste(name, "[", 0:(length(v)-1), "] << ", 
        sapply(v,function(x){
            paste(x,collapse=", ")
        }), 
    ";", sep="")
}

formatVector<-function(name,v){
    paste(name, " << ", 
        paste(v,collapse=", "),
        ";", sep="")
}

outputToCPP<- function(filename=""){
    cat("", fill=T, file=filename)
    cat(formatList("lower_array", lower), fill=T, file=filename, append=T)
    cat(formatList("upper_array", upper), fill=T, file=filename, append=T)
    
    cat(paste("trans_matrix[",0:(length(transition)-1),"] << ", 
        sapply(transition,function(x){
            paste(t(x),collapse=", ") ## eigen << init take it by row
        }), 
        ";", sep=""),
    fill=T, file=filename, append=T)
    
    cat(formatVector("expected_up_fast", up_fast), fill=T, file=filename, append=T)
    cat(formatVector("expected_up", up), fill=T, file=filename, append=T)

    cat(formatVector("expected_father_fast", father_fast), fill=T, file=filename, append=T)
    cat(formatVector("expected_father", father), fill=T, file=filename, append=T)
    
    cat(formatVector("expected_mother_fast", mother_fast), fill=T, file=filename, append=T)
    cat(formatVector("expected_mother", mother), fill=T, file=filename, append=T)
    
    cat(formatVector("expected_child_fast", child_fast), fill=T, file=filename, append=T)
    cat(formatVector("expected_child", child), fill=T, file=filename, append=T)
}

outputToCPP("generateFamilyData")

