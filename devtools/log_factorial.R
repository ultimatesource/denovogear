g = lgamma(1:1024)
s = sprintf("%0.17g", g)

if(all(as.numeric(s) == g)) {
	m = matrix(s,ncol=4,byrow=TRUE)
	ss = apply(m, 1, paste, collapse=", ")
	ss = paste(ss, collapse=",\n\t")
	cat("constexpr double log_factorial[128] = {\n\t")
	cat(ss)
	cat("\n};\n")
}

