# bkEvol
## The Calculation Procedure
The program is using several ways to boost the computation speed.

When fitting the HERA data there are 3 elementary steps which can be computationally intense:
1. Calculation of the evolution kernel (initialPDF -> PDF)
2. Calculation of the convolution kernel (PDF -> F2, FL, sigmaR)
3. Applying convolution and evolution kernel (called several hundreds of times when fitting)

The parallelisation of the step 1) and 2) is trivial as the kernel cubes (3D matrices) can be calculated independently for each z (evolution kernel) and each Q2 (convolution kernel).
These kernels can be consequently stored into the files (currently HDF5 is used).
Fitting then represents only loading the kernels form hdf5 files and applying them multiple times to the inputPDF to find the minimum.

### Parallelisation applied in each step:
1. Evolution kernel
Currently openMP is coded as default to use all CPU cores.
In addition the openMPI parallelisation is also implemented for systems with shared memory (~100 cores).  In case that only one CPU is available (e.g. running on laptop) this parallelisation has no effect.
The batch calculation of the evolution kernel can be easily implemented but as there is not large time issue, it's not used yet.

2. Convolution kernel
The calculation requires multidimensional integration which is quite time consuming.
To boost it, each q2 node is calculated as one job at batch system (~50jobs).
Each such job employs openMP to make use of all cores of the machine.

3. Application of the kernels
This last step is purely about linear algebra but cannot be parallelised trivially.
The PDF always depends only on PDF at higher x-values.
Therefore, the next, lower x-value depends on the previous, higher, one.
To boost such calculation the special library for linear algebra (armadillo using BLAS procedures) is used.
The BLAS procedures internally use openMP and SIMD parallelisation. 
Furthermore, the evaluation of PDF for single x-value is parallelised by openMPI.
As an alternative, also the GPU-based kernel application is implemented.
In this case, the BLAS procedures are evaluated at GPU within CUDA.
It can speed up the calculation substantially but nVidia GPU is required.

## Installation:
### Required libraries
The code is written in such a way to be able to be executed on usual laptop.
In such case the running speed is sufficient for single evolution but too slow for fitting.
Following libraries are used
-  installation of ROOT, as Minuit is used for minimisation.
-  GSL for alternative minimisation
-  HDF5 library for writing kernels to the files
-  openMPI libraries
-  BLAS + Armadillo for linear algebra

### Four Lines to Get Result from any OS
To simplify the installation Docker image including all libraries above was created, so the simplest way to get results is:
```
git clone https://github.com/zleba/bkEvol.git
cd bkEvol
./rd make
./rd ./iter \< steer/config.ini
```
Where the first command download the source code (bkEvol directory), then one runs the compilation using compiler in the Docker image and finally runs the program itself with the given steering file.
The `rd` script will first automatically download the image which can take a while.
Note that Docker must be installed. Instructions for various operation systems are given in 

https://docs.docker.com/install/

For example in Ubuntu one can get Docker and git by this command
```
sudo apt-get install docker.io git
```
Notice, that running via Docker is as quick as the native running. 

### Compiling using Classical Way
Without using Docker all the packages above must be installed which can be tedious, especially the installation of ROOT. The Armadillo can be installed using the attached script `installARMA.sh`.
For other libraries I recommend as an inspiration the content of the `Dockerfile` and of file containing the libraries required by ROOT

https://github.com/root-project/rootspi/blob/master/docker/ubuntu16-base/packages


## Configuration
The interface using the steering file is under development.
Currently one can set few parameters, as the alphaS of the fit and the parameters of the grids in kT and x.

When fitting the HERA data, the nice feature is that the fitted formula for initial PDF parametrisation as well as the starting values of the parameters and their ranges can be entered in the steering.
