CC=gcc -g -Wall -O2
LIB=

all: bgzf.o dict.o array.o kstring.o bcf.o knetfile.o vcf.o bcfutils.o

bgzf.o: bgzf.c bgzf.h
	$(CC) -c $<

array.o: array.c array.h
	$(CC) -c $<

dict.o: dict.c dict.h
	$(CC) -c $<

readseq.o: readseq.c readseq.h
	$(CC) -c $<

kstring.o: kstring.c 
	$(CC) -c $< 

bcf.o: bcf.c 
	$(CC) -c $< 

knetfile.o: knetfile.c 
	$(CC) -c $<

vcf.o: vcf.c
	$(CC) -c $< 

bcfutils.o: bcfutils.c
	$(CC) -c $<


