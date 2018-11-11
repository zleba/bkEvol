#Run with nVidia GPU
hasGPU=false

#ROOT libraries
rootCFLAGS := $(shell root-config --cflags)
rootLibs := $(shell root-config --libs)
$(info ${rootCFLAGS} ${rootLibs})

#ARMADILLO path (install by ./installARMA.sh)
armaDir :=$(shell  echo \${ARMADIR})
ifeq ($(strip $(armaDir)),)
    armaDir := $(shell  echo \${PWD}/arma)
endif
armaLib := $(shell echo \${armaDir}/install/*lib*)  #Currently only include used
armaInc := $(shell echo \${armaDir}/install/include)
$(info $$armaLib is [${armaLib}])


#HDF5 library (to store data tables to files)
hdf5Lib :=-lhdf5_serial
ifeq ($(shell ldconfig -p | grep hdf5_serial.so),)
	hdf5Lib:=-lhdf5
    $(info Inside if)
endif



#CC=mpic++ #MPI compiler
CC=g++ #MPI compiler
myFlags=-std=c++11 -O2 -g  -fopenmp  -I${CMAKE_PREFIX_PATH}/include


ifeq ($(strip $(hasGPU)),true)
	CUDAbase=/usr/local/cuda-9.1
	CUDAcompiler=$(CUDAbase)/bin/nvcc
	#CUDAcompiler=/usr/bin/nvcc
	CUDAlib=-L$(CUDAbase)/targets/x86_64-linux/lib/ -lcudart  -lcublas  
	CUDAINC=-I$(CUDAbase)/targets/x86_64-linux/include/ 
endif

GSL=-lgsl  -lgslcblas

#iter: iterate.cpp
	#g++ -g -O3 $(CFLAGS)  $^ $(LIBS) -lgsl -lgslcblas -fopenmp -o $@ 


SRCS = src/Solver.cpp  src/iterate.cpp  src/main.cpp src/integration.cpp src/Fitter.cpp src/alphaSpline.cpp src/Spline.cpp src/newKernels.cpp 
OBJS = obj/Solver.o obj/iterate.o   obj/integration.o  obj/Fitter.o obj/alphaSpline.o obj/Spline.o  obj/newKernels.o
ifeq ($(strip $(hasGPU)),true)
	OBJS += obj/gpuBooster.o
endif

INCS = inc/Solver.h inc/integration.h inc/gpuBooster.h inc/Fitter.h inc/Settings.h inc/kernels.h


bkEvol: $(OBJS) obj/main.o
	$(CC)  $(myFlags)  $^ $(rootLibs) -lMinuit   $(hdf5Lib)  -lopenblas  $(CUDAlib) $(GSL) -lm   -L${CMAKE_PREFIX_PATH}/lib -lmpi  -lmpi_cxx    -o $@  #-llapack 

bkevol.so:   $(OBJS) obj/pyBind.o
	$(CC) $(myFlags) -shared -std=c++11   $^ -o $@  $(rootLibs) -lMinuit    -lopenblas  $(CUDAlib) $(GSL) -lm  -L${CMAKE_PREFIX_PATH}/lib $(hdf5Lib)    -fPIC

#`python -m pybind11 --includes` 
obj/pyBind.o: src/pyBind.cpp $(INCS)
	$(CC) $(myFlags)   -I$(armaInc)  `python-config --includes` -Ipybind  -DpwdDir=$(shell pwd) -DhasGPU=$(hasGPU) -I./inc   -DARMA_DONT_USE_WRAPPER   $(CUDAINC)  -c   -o $@ $< $(rootCFLAGS) -fPIC

obj/%.o: src/%.cpp $(INCS)
	$(CC) $(myFlags)   -I$(armaInc)   -DpwdDir=$(shell pwd) -DhasGPU=$(hasGPU) -I./inc   -DARMA_DONT_USE_WRAPPER   $(CUDAINC)  -c   -o $@ $< $(rootCFLAGS) -fPIC

obj/gpuBooster.o: src/gpuBooster.cxx
	 $(CUDAcompiler) -arch=sm_52 -ccbin g++ -I./inc -I$(CUDAbase)/samples/common/inc/ -I$(armaInc) -DARMA_DONT_USE_WRAPPER  -m64  -Xcompiler  -fopenmp    -std=c++11  -o obj/gpuBooster.o -c src/gpuBooster.cxx

clean:
	rm obj/*.o bkEvol bkevol.so
