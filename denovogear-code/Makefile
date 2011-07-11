CC = g++
LIB=-lm -lz -lnewmat
NEWMAT=newmat
DIRS=utils


all: denovogear

denovogear: denovogear.cc glf.o utils
	$(CC) -Iutils -I$(NEWMAT) -o denovogear denovogear.cc glf.o utils/bgzf.o utils/array.o utils/dict.o $(LIB)  -L$(NEWMAT)

glf.o: glf.c glf.h
	$(CC) -Iutils -c glf.c

utils: force_update
	set -e; $(MAKE) -C utils


NM := $(wildcard $(NEWMAT)*.o)

libnewmat : $(NM)
	ar rc libnewmat.a $(NM)

clean:
	rm *.o
	rm utils/*.o

force_update:
	true
