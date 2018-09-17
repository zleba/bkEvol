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


OBJ=obj/cheb.o obj/integration.o obj/nodes.o obj/Spline.o obj/alphaSpline.o obj/kernels.o 
cheb: ${OBJ} obj/main.o 
	g++ -g -O2 $^ -o $@  $(rootLibs) -lopenblas

one: ${OBJ} obj/one.o 
	g++ -g -O2 $^ -o $@  $(rootLibs) -lopenblas  -fopenmp 

solve: ${OBJ} obj/solve.o
	g++ -g -O2 $^ -o $@  $(rootLibs)  -lopenblas

solveOne: ${OBJ} obj/solveOne.o
	g++ -g -O2 $^ -o $@  $(rootLibs)  -lopenblas


testmain: obj/testmain.o obj/tests.o  $(OBJ)
	g++ -g  -O2 $^ -o $@  -lopenblas

obj/%.o: src/%.cpp 
	g++ -g -O2 -fopenmp -Iinc  $(rootCFLAGS) -DARMA_DONT_USE_WRAPPER  -I$(armaInc)  -c $^    -o $@  
