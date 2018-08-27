#ROOT libraries
rootCFLAGS := $(shell root-config --cflags)
rootLibs := $(shell root-config --libs)


#ARMADILLO path (install by ./installARMA.sh)
armaDir :=$(shell  echo \${ARMADIR})
ifeq ($(strip $(armaDir)),)
	    armaDir := $(shell  echo \${PWD}/../arma)
	endif
armaLib := $(shell echo \${armaDir}/install/*lib*)  #Currently only include used
armaInc := $(shell echo \${armaDir}/install/include)


OBJ=obj/cheb.o obj/integration.o obj/nodes.o obj/Spline.o obj/alphaSpline.o
cheb: ${OBJ} obj/main.o
	g++ -g -O2 $^ -o $@  $(rootLibs)

solve: ${OBJ} obj/solve.o
	g++ -g -O2 $^ -o $@  $(rootLibs)  -lopenblas


testmain: obj/testmain.o obj/tests.o obj/cheb.o obj/nodes.o  obj/Spline.o obj/alphaSpline.o
	g++ -g  -O2 $^ -o $@

obj/%.o: src/%.cpp
	g++ -g -O2 -Iinc  $(rootCFLAGS) -DARMA_DONT_USE_WRAPPER  -I$(armaInc)  -c $^ -o $@  
