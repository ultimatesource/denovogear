#! Rscript --vanilla

library(doParallel)

registerDoParallel()

DNGLL_BIN = "dng"

fixed_pars = c(
    "theta" = 0.001
)

init_pars = c(
        "lib-bias" = 1.009997e+00,
        "lib-error" =  3.200240e-03,
        "lib-overdisp" = 2.526826e-02,
        "mu-somatic" = 3.276369e-04,
        "ref-weight" = 1
)

upper_pars = c(
    "lib-bias" = 1.01,
    "lib-error" = 0.1,
    "lib-overdisp" = 0.5,
    "mu-somatic" = 0.1,
    "theta" = 10,
    "ref-weight" = 10000
)

lower_pars = c(
    "lib-bias" = 0.5,
    "lib-error" = 1e-10,
    "lib-overdisp" = 1e-10,
    "mu-somatic" = 1e-15,
    "theta" = 1e-10,
    "ref-weight" = 1e-10
)

# Run dng loglike and extract the log like
run_once = function(pars,input) {
    pars = c(pars,fixed_pars)
    args = paste("--", names(pars), "=", pars,sep="")
    if(grepl("dng$",DNGLL_BIN)) {
        args = c("loglike",args)
    }
    out = system2(DNGLL_BIN, args=c(args,input),stdout=TRUE,stderr=FALSE)
    if(is.numeric(out)) {
        return(NA)
    }
    score = strsplit(tail(out,1),"\t")[[1]][1]
    -as.numeric(score)
}

loglike = function(pars,peds,inputs) {
    n = names(pars)
    if(any(pars < lower_pars[n] | pars > upper_pars[n])) {
        return(NA)
    }
    #print(pars)
    results = foreach(i=seq_along(inputs),.combine=c) %dopar% run_once(c(pars, ped=peds[i]),inputs[i])
    total = sum(results)
    total
}

main = function(peds,inputs) {
    pars = init_pars
    o = optim(pars,loglike,peds=peds,inputs=inputs,method="Nelder-Mead",control=list(
         trace=6,REPORT=1,reltol=1e-8
    ))
    o
}

if(!interactive()) {
    ARGS = commandArgs(trailingOnly=TRUE)
    data = read.csv(ARGS[1],header=F,stringsAsFactors=FALSE)
    o = main(data$V1,data$V2)
    print(o)
}
