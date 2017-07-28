#! Rscript --vanilla

DNGLL_BIN = "dng loglike"
DNGLL_BIN = "./src/dng-loglike"

# Run dng loglike and extract the log like
run_once = function(pars,input) {
    args = paste("--", names(pars), "=", pars,sep="")
    out = system2(DNGLL_BIN, args=c(args,input),stdout=TRUE,stderr=FALSE)
    if(is.numeric(out)) {
        return(NA)
    }
    score = strsplit(tail(out,1),"\t")[[1]][1]
    -as.numeric(score)
}

loglike = function(pars,peds,inputs) {
    total = c()
    for(i in seq_along(inputs)) {
        p = c(pars, ped=peds[i])
        score = run_once(p,inputs[i])
        total = c(total,score)
    }
    print(total)
    total = sum(total)
    print(c(total,pars))

    total
}

main = function(peds,inputs) {
    pars = make_pars(c(
        "lib-bias" = 1.01,
        "lib-error" = 0.0036,
        "lib-overdisp" = 0.038,
        "mu-somatic" = 0.00055,
        "theta" = 0.054,
        "ref-weight" = 1.01
    ))
    o = optim(pars,loglike,peds=peds,inputs=inputs,method="Nelder-Mead",control=list(
        trace=6,REPORT=1,reltol=1e-4
    ))
    o
}

if(!interactive()) {
    ARGS = commandArgs(trailingOnly=TRUE)
    data = read.csv(ARGS[1],header=F,stringsAsFactors=FALSE)
    main(data$V1,data$V2)
}
