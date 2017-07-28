#! Rscript --vanilla

DNGLL_BIN = "dng loglike"

# Run dng loglike and extract the log like
run_once = function(x) {

}

main = function(args) {

}

# GTEX-13OVK.20170717.vcf
# GTEX-WOFM.20170717.vcf


if(!interactive()) {
    ARGS = commandArgs(trailingOnly=TRUE)
    data = readLInes(ARGS[1])
    main(data)
}
