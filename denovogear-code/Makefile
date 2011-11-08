CC=gcc
CPP=g++ -g -Wall -O2
LIB=-lm -lz -lnewmat
NEWMAT=newmat
DIRS=utils

all: denovogear

denovogear: denovogear.cc denovogear.h bcf2Qcall.o bcf2Paired.o pedParser.o snpLike.o indelLike.o  pairLike.o makeLookup.o utils
	$(CPP) -Iutils -I$(NEWMAT) -o denovogear denovogear.cc snpLike.o indelLike.o pairLike.o bcf2Qcall.o bcf2Paired.o pedParser.o makeLookup.o utils/bcf.o utils/vcf.o utils/bcfutils.o utils/knetfile.o utils/kstring.o utils/bgzf.o utils/array.o utils/dict.o $(LIB) -L. -L$(NEWMAT)

bcf2Qcall.o: bcf2Qcall.c pedParser.o
	$(CPP) -c -Iutils -o bcf2Qcall.o bcf2Qcall.c 

bcf2Paired.o: bcf2Paired.c pedParser.o
	$(CPP) -c -Iutils -o bcf2Paired.o bcf2Paired.c 

snpLike.o: snpLike.cc
	$(CPP) -c -I$(NEWMAT) -Iutils -o snpLike.o snpLike.cc

indelLike.o: indelLike.cc
	$(CPP) -c -I$(NEWMAT) -Iutils -o indelLike.o indelLike.cc

pairLike.o: pairLike.cc
	$(CPP) -c -I$(NEWMAT) -Iutils -o pairLike.o pairLike.cc

pedParser.o: pedParser.c pedParser.h
	$(CPP) -c -Iutils -o pedParser.o pedParser.c

makeLookup.o: makeLookup.cc makeLookup.h
	$(CPP) -c -I$(NEWMAT) -o makeLookup.o makeLookup.cc

utils: force_update
	set -e; $(MAKE) -C utils

NM := $(wildcard $(NEWMAT)*.o)

libnewmat : $(NM)
	ar rc libnewmat.a $(NM)

clean:
	rm -f denovogear
	rm -f *.o
	rm -f utils/*.o
	rm -f $(NEWMAT)/*.[ao]

force_update:
	true