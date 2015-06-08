n = c("A","C","G","T")
h = outer(n,n,function(x,y) { sprintf("%s>%s",x,y) })
d = kronecker(h,h,function(x,y) { sprintf("%s,%s",x,y) })
s = strsplit(d,",|>")
z = sapply(s, function(x) {
	if(x[1] != x[2] && x[3] != x[4]) {
		return(sprintf("%s%s>%s%s", x[1],x[3],x[2],x[4]))
	} else if(x[1] != x[2]) {
		return(sprintf("%s>%s", x[1],x[2]))		
	} else if(x[3] != x[4]) {
		return(sprintf("%s>%s", x[3],x[4]))	
	}
	return("")
})
z = matrix(z,16,16)

g = outer(n,n,function(x,y) { sprintf("%s%s",x,y) })
f = t(g)[lower.tri(g,TRUE)]
m = outer(f,f,function(x,y) { sprintf("%s>%s",x,y) })

p = kronecker(f,f,function(x,y) { sprintf("%sx%s",x,y) })
q = outer(p,f,function(x,y) { sprintf("%s>%s",x,y)})

a = apply(m,1,function(x) {paste(sprintf("\"%s\"",x), collapse=", ")})
a = sprintf("    {%s}",a)
a = paste(a,collapse=",\n")
a = sprintf("{\n%s\n};\n",a)

cat("const char dng::mitotic_diploid_mutation_labels[10][10][6] = ")
cat(a)

a = apply(q,1,function(x) {paste(sprintf("\"%s\"",x), collapse=", ")})
a = sprintf("    {%s}",a)
a = paste(a,collapse=",\n")
a = sprintf("{\n%s\n};\n",a)

cat("const char dng::meiotic_diploid_mutation_labels[100][10][9] = ")
cat(a)
