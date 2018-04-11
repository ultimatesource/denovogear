#!/usr/bin/Rscript --vanilla
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
library(dfoptim)

registerDoParallel()

fixed_pars = c(
    "mu" = 0,
    "lib-bias" = 1
)

init_pars = c(
    "theta" = 0.021359190,
    "lib-error" = 0.001256881,
    "ref-bias-hom" = 0.938929627,
    "ref-bias-het" = 0.757601819,
    "lib-overdisp-hom" = 0.001,
    "lib-overdisp-het" = 0.001
)

lower_pars = c(
    "mu" = 0,
    "mu-library" = 0,
    "mu-somatic" = 0,   
    "theta" = 0,
    "ref-bias-hom" = -1,
    "ref-bias-het" = -1,
    "ref-bias-hap" = -1,
    "lib-bias" = 0,
    "lib-error" = 0,
    "lib-overdisp-hom" = 0,
    "lib-overdisp-het" = 0  
)

upper_pars = c(
    "mu" = 1,
    "mu-library" = 1,
    "mu-somatic" = 1,
    "theta" = 1,
    "ref-bias-hom" = 1,
    "ref-bias-het" = 1,
    "ref-bias-hap" = 1,
    "lib-bias" = 2,
    "lib-error" = 1,
    "lib-overdisp-hom" = 1,
    "lib-overdisp-het" = 1
)

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

loglike = function(pars,cmds,names) {
    names(pars) = names
    results = foreach(cmd=cmds,.combine=c) %dopar% run_once(cmd,pars)
    total = sum(results)
    total
}

main = function(cmds) {
    pars = init_pars
    par_names=names(pars)
    o = nmkb(pars,loglike,lower=lower_pars[par_names],upper=upper_pars[par_names],
         cmds=cmds,names=par_names,control=list(
         trace=TRUE,tol=1e-6
    ))
    names(o$par) = par_names
    o
}

if(!interactive()) {
    ARGS = commandArgs(trailingOnly=TRUE)
    cmds = readLines(ARGS[1])
    o = main(cmds)
    print(o)
}

