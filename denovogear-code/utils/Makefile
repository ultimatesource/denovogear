CC=gcc -g -Wall -O2
LIB=

all: bgzf.o dict.o array.o kstring.o bcf.o knetfile.o vcf.o bcfutils.o sam_view.o sam.o bam.o bedidx.o bam_aux.o sam_header.o bam_import.o bam_pileup.o bam_index.o faidx.o razf.o

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

sam_view.o: sam_view.c 
	$(CC) -c $< 

sam.o: sam.c 
	$(CC) -c $< 

bam.o: bam.c 
	$(CC) -c $< 

bedidx.o: bedidx.c 
	$(CC) -c $< 

bam_aux.o: bam_aux.c
	$(CC) -c $< 

sam_header.o: sam_header.c
	$(CC) -c $< 

bam_import.o: bam_import.c
	$(CC) -c $< 

bam_pileup.o: bam_pileup.c
	$(CC) -c $< 

bam_index.o: bam_index.c
	$(CC) -c $< 

faidx.o: faidx.c
	$(CC) -c $< 

razf.o: razf.c
	$(CC) -c $< 
