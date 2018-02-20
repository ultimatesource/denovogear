#! Rscript --vanilla
#
# Copyright (c) 2017-2018 Reed A. Cartwright
# Copyright (c) 2018      Adam J. Orr
#
# Authors:  Reed A. Cartwright <reed@cartwrig.ht>
#           Adam J. Orr <adamjorr@gmail.com>
#
# This file is part of DeNovoGear.
#
# DeNovoGear is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/>.
#


options(stringsAsFactors=FALSE)

library(doParallel)

registerDoParallel()

fixed_pars = c(
    "lib-overdisp-hom" = 1e-6,
    "lib-overdisp-het" = 1e-6,
    "lib-error" =  0,
    "mu" = 1e-8,
    "lib-bias" = 1
)

init_pars = c(
    "theta" = 0.01,
    "asc-hom" = 0.1,
    "asc-het" = 0.1
)

link_func = list(
    logit = list(
        "do" = function(x) {
            log(x)-log1p(-x)
        },
        "undo" = function(x) {
            y = exp(x)
            y/(y+1)
        }
    ),
    log = list("do" = log, "undo" = exp)
)

link = c(
    "mu" = "logit",
    "mu-library" = "logit",
    "mu-somatic" = "logit",
    "mu-entropy" = "log",
    "theta" = "log",
    "asc-hom" = "logit",
    "asc-het" = "logit",
    "asc-hap" = "logit",

    "lib-bias" = "log",
    "lib-error" = "logit",
    "lib-error-entropy" = "log",
    "lib-overdisp-hom" = "logit",
    "lib-overdisp-het" = "logit"
)

parscale = c(
    "mu" = 1e-5,
    "mu-library" = 1e-5,
    "mu-somatic" = 1e-5,
    "mu-entropy" = 1e-5,
    "theta" = 1e-5,
    "asc-hom" = 0.01,
    "asc-het" = 0.01,
    "asc-hap" = 0.01,
    "lib-bias" = 1e-3,
    "lib-error" = 1e-5,
    "lib-error-entropy" = 1e-5,
    "lib-overdisp-hom" = 1e-5,
    "lib-overdisp-het" = 1e-5
)
parscale_multiplier = 1 #increase this value for additional sites

makeparscale = function(x){
    for(n in names(x)) {
        x[n] = parscale_multiplier * parscale[n]
    }
    x
}

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
    scanned_args = scan(text=cmd,what=character(),quiet=TRUE)
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
    pscale = makeparscale(pars)
    o = optim(pars,loglike,cmds=cmds,method="BFGS",control=list(
         trace=6,REPORT=1,reltol=1e-8,maxit=1000,parscale=pscale
    ))
    o$par = unmakepars(o$par)
    o$parscale = pscale
    o
}

if(!interactive()) {
    ARGS = commandArgs(trailingOnly=TRUE)
    cmds = readLines(ARGS[1])
    o = main(cmds)
    print(o)
}
