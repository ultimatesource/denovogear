#! Rscript --vanilla

options(stringsAsFactors=FALSE)

library(doParallel)

registerDoParallel()

fixed_pars = c(
)

init_pars = c(
    "theta" = 0.247461702,
    "lib-error" =  0.002311513,
    "lib-overdisp" = 0.184123736,
    "lib-bias" = 2.473174701,
    "ref-weight" = 4.929903491
)

link_func = list(
    logit = list(
        do = function(x) {
            log(x)-log1p(-x)
        },
        undo = function(x) {
            y = exp(x)
            y/(y+1)
        },
    ),
    log = list(do = log, undo = exp)
)

link = c(
    "lib-bias" = "log",
    "lib-error" = "logit",
    "lib-overdisp" = "logit",
    "mu-library" = "logit",
    "mu-somatic" = "logit",
    "mu-germline" = "logit",
    "theta" = "log",
    "ref-weight" = "log"
)

makepars = function(x) {
    for(n in names(x)) {
        x[n] = link_func[[link[n]]]$do(x[n])
    }
    x
}

unmakepars = function(x) {
    for(n in names(x)) {
        x[n] = link_func[[link[n]]]$undo(x[n])
    }
    x
}

# Run dng loglike and extract the log like
run_once = function(cmd, pars) {
    scanned_args = scan(text=cmd,what=character())
    prog = scanned_args[1]
    scanned_args = scanned_args[-1]
    pars = c(fixed_pars, pars)
    args = paste("--", names(pars), "=", pars,sep="")
    o = which(scanned_args == "--")[1]
    if(is.na(o)) {
        args = c(scanned_args,args)
    } else {
        args = append(scanned_args,args,after=o-1)
    }
    out = system2(prog, args=args,stdout=TRUE,stderr=FALSE)
    if(is.numeric(out)) {
        return(NA)
    }
    score = strsplit(tail(out,3),"\t")[[1]][2]
    -as.numeric(score)
}

loglike = function(pars,cmds) {
    pars = unmakepars(pars)
    results = foreach(cmd=cmds,.combine=c) %dopar% run_once(cmd,pars)
    total = sum(results)
    total
}

main = function(cmds) {
    pars = makepars(init_pars)
    o = optim(pars,loglike,cmds=cmds,method="BFGS",control=list(
         trace=6,REPORT=1,reltol=1e-8,maxit=1000
    ))
    o$par = unmakepars(o$par)
    o
}

if(!interactive()) {
    ARGS = commandArgs(trailingOnly=TRUE)
    cmds = readLines(ARGS[1])
    o = main(cmds)
    print(o)
}
