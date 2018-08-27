#ARMAlib=/home/zlebcr/libs/arma/armalib/lib64
#ARMAinc=/home/zlebcr/libs/arma/armalib/include

ARMAlib=/home/radek/Dropbox/Krakow/bkNew/bkEvol/arma/install/lib
ARMAinc=/home/radek/Dropbox/Krakow/bkNew/bkEvol/arma/install/include

CFLAGS = $(shell root-config --cflags)
LIBS := $(shell root-config --libs)


fitting: fitting.cpp
	g++ -fopenmp -std=c++11 -g -O3 -I$(ARMAinc) $(CFLAGS)  $<  -Wl,-R$(ARMAlib)  -L$(ARMAlib) $(LIBS)  -larmadillo -o $@

plotKernel: plotKernel.C
	g++ -fopenmp -std=c++11 -g -O3 -I$(ARMAinc) $(CFLAGS)  $<  -Wl,-R$(ARMAlib)  -L$(ARMAlib) $(LIBS)  -larmadillo -o $@
