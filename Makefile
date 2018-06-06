<<<<<<< HEAD
CFLAGS = $(shell root-config --cflags)
LIBS := $(shell root-config --libs)


armaLib=/home/zlebcr/libs/arma/armalib/lib64
armaInc=/home/zlebcr/libs/arma/armalib/include

CC=g++
CC=mpic++

CUDA=/usr/local/cuda-8.0/bin/nvcc
#CUDA=/usr/bin/nvcc
=======
#Run with nVidia GPU
hasGPU=false

#ROOT libraries
rootCFLAGS := $(shell root-config --cflags)
rootLibs := $(shell root-config --libs)

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



CC=mpic++ #MPI compiler
myFlags=-std=c++11 -O3 -g -fopenmp 


ifeq ($(strip $(hasGPU)),true)
	CUDAbase=/usr/local/cuda-9.1
	CUDAcompiler=$(CUDAbase)/bin/nvcc
	#CUDAcompiler=/usr/bin/nvcc
	CUDAlib=-L$(CUDAbase)/targets/x86_64-linux/lib/ -lcudart  -lcublas  
	CUDAINC=-I$(CUDAbase)/targets/x86_64-linux/include/ 
endif

GSL=-lgsl  -lgslcblas
>>>>>>> 3218ac9a0ae3dc07b896fb991fac340ae8181511

#iter: iterate.cpp
	#g++ -g -O3 $(CFLAGS)  $^ $(LIBS) -lgsl -lgslcblas -fopenmp -o $@ 


SRCS = src/Solver.cpp  src/iterate.cpp src/kernels.cpp src/main.cpp src/integration.cpp src/Fitter.cpp src/alphaSpline.cpp src/Spline.cpp  
<<<<<<< HEAD
OBJS = obj/Solver.o obj/iterate.o obj/kernels.o obj/main.o  obj/integration.o  obj/Fitter.o obj/alphaSpline.o obj/Spline.o  #obj/gpuBooster.o
=======
OBJS = obj/Solver.o obj/iterate.o obj/kernels.o obj/main.o  obj/integration.o  obj/Fitter.o obj/alphaSpline.o obj/Spline.o  
ifeq ($(strip $(hasGPU)),true)
	OBJS += obj/gpuBooster.o
endif
>>>>>>> 3218ac9a0ae3dc07b896fb991fac340ae8181511



#obj/%.o: src/%.cpp inc/Solver.h inc/integration.h
	#$(CC) -g  -fopenmp  -I$(armaInc) -DARMA_DONT_USE_WRAPPER   -I./inc -c   -o $@ $< $(CFLAGS)
#iter: $(OBJS) 
	#$(CC) -g  $^ $(LIBS)  -fopenmp  /opt/anaconda/3/pkgs/mkl-2017.0.3-0/lib/libmkl_rt.so -o $@

<<<<<<< HEAD

iter: $(OBJS) 
	$(CC) -g  $^ $(LIBS) -lMinuit  -fopenmp   -Wl,-R$(armaLib) -L$(armaLib) -larmadillo -L/usr/local/cuda-8.0/targets/x86_64-linux/lib/ -lcudart  -lcublas  -lgsl -lgslcblas -lm   -o $@


obj/%.o: src/%.cpp inc/Solver.h inc/integration.h inc/gpuBooster.h inc/Fitter.h inc/Settings.h
	$(CC) -g  -fopenmp  -I$(armaInc)  -DpwdDir=$(shell pwd)  -I./inc  -I/usr/local/cuda-8.0/targets/x86_64-linux/include/   -c   -o $@ $< $(CFLAGS)

obj/gpuBooster.o: src/gpuBooster.cxx
	 $(CUDA) -arch=sm_52 -ccbin g++ -I./inc -I/usr/local/cuda-8.0/samples/common/inc/ -I/home/zlebcr/libs/arma/armalib/include -DARMA_DONT_USE_WRAPPER  -m64  -Xcompiler  -fopenmp    -std=c++11  -o obj/gpuBooster.o -c src/gpuBooster.cxx


=======
#-Wl,-R$(armaLib) -L$(armaLib)

iter: $(OBJS) 
	$(CC)  $(myFlags)  $^ $(rootLibs) -lMinuit   $(hdf5Lib)  -lopenblas  $(CUDAlib) $(GSL) -lm   -o $@  #-llapack 


obj/%.o: src/%.cpp inc/Solver.h inc/integration.h inc/gpuBooster.h inc/Fitter.h inc/Settings.h
	$(CC) $(myFlags)   -I$(armaInc)  -DpwdDir=$(shell pwd) -DhasGPU=$(hasGPU) -I./inc   -DARMA_DONT_USE_WRAPPER   $(CUDAINC)  -c   -o $@ $< $(rootCFLAGS)

obj/gpuBooster.o: src/gpuBooster.cxx
	 $(CUDAcompiler) -arch=sm_52 -ccbin g++ -I./inc -I$(CUDAbase)/samples/common/inc/ -I$(armaInc) -DARMA_DONT_USE_WRAPPER  -m64  -Xcompiler  -fopenmp    -std=c++11  -o obj/gpuBooster.o -c src/gpuBooster.cxx

clean:
	rm obj/*.o iter
>>>>>>> 3218ac9a0ae3dc07b896fb991fac340ae8181511
