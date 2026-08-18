CC=gcc

FFTW_HOME=/afs/cern.ch/sw/lcg/external/fftw3/3.1.2/slc4_amd64_gcc34

CCOPT=-O3 -ffast-math -fno-keep-inline-functions -fomit-frame-pointer -I$(FFTW_HOME)/include -fPIC
LLOPT=-L$(FFTW_HOME)/lib -lfftw3 -lm

guinea : guinea_pig.c background.c grv.c file.d memory.o switches.d rndm.d histogram.o physconst.h spline.h hit.c version.h store_particle.c fourtrans3.c
	$(CC) $(CCOPT) -D USE_FFTW -o guinea guinea_pig.c histogram.o memory.o $(LLOPT)

guinea_nofftw : guinea_pig.c background.c grv.c file.d memory.o switches.d rndm.d histogram.o physconst.h spline.h hit.c version.h store_particle.c
	$(CC) $(CCOPT) -o guinea_nofftw guinea_pig.c histogram.o memory.o -lm

guinea.so : guinea_pig.o memory.o histogram.o
	gcc -O3 -shared -o guinea.so guinea_pig.o memory.o histogram.o -L$(FFTW_HOME) -lfftw -lm

guinea.o : guinea_pig.c background.c grv.c file.d memory.o switches.d rndm.d histogram.o physconst.h spline.h hit.c version.h store_particle.c
	$(CC) $(CCOPT) -c guinea_pig.c

memory.o : memory.c memory.h
	$(CC) $(CCOPT) -O3 -c memory.c

histogram.o : histogram.c memory.h
	$(CC) $(CCOPT) -c -O3 histogram.c -lm

gpv : view.c histogram.o histogram.h memory.o memory.h file.d
	$(CC) $(CCOPT) -o gpv view.c histogram.o memory.o -lm

scatter : convert.c histogram.h histogram.o rndm.d memory.h memory.o
	$(CC) $(CCOPT) -pipe -O3 -o scatter  convert.c histogram.o memory.o -lm

gp.tar.gz : guinea_pig.c background.c grv.c file.d memory.c memory.h switches.d rndm.d histogram.c histogram.h physconst.h spline.h hit.c makefile view.c hist.d rndm.d gpv.tcl lumi_extr.tcl extract.tcl acc.dat version.h store_particle.c coherent.c bessel.c dump2.c convert.c fourtrans3.c
	tar -cf gp.tar guinea_pig.c background.c grv.c file.d memory.c memory.h switches.d rndm.d histogram.c histogram.h physconst.h spline.h hit.c makefile view.c hist.d rndm.d gpv.tcl lumi_extr.tcl extract.tcl acc.dat version.h store_particle.c coherent.c bessel.c dump2.c convert.c fourtrans3.c
	gzip -9 gp.tar

clean:
	rm -f guinea memory.o histogram.o
