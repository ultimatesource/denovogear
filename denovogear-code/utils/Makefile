CC=gcc
LIB=

all: bgzf.o dict.o array.o readseq.o

bgzf.o: bgzf.c bgzf.h
	$(CC) -c bgzf.c

array.o: array.c array.h
	$(CC) -c array.c

dict.o: dict.c dict.h
	$(CC) -c dict.c

readseq.o: readseq.c readseq.h
	$(CC) -c readseq.c
