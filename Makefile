COMPILER_OPTIONS= -Wall -fopenmp -O3 -std=c++11 -g
CC              = g++-6 $(COMPILER_OPTIONS)
# Change this to:
# -lafcpu    - For CPU backend
# -lafcuda   - For CUDA backend
# -lafopencl - For OpenCL backend
# -laf       - For unified backend
LIBS            = -laf -lhdf5
LIB_PATHS       =  -L $(HDF5_DIR)/lib -L $(AF_PATH)/lib
INCLUDES        = -I $(AF_PATH)/include -I $(EIGEN_PATH) -I $(HIGHFIVE_PATH)/include -I ./header -I $(HDF5_DIR)/include
SOURCES         = ./test/testFMM/testFMM2DAdaptive.cpp
EXECUTABLE      = ./exec/testFMM2D

all:
	$(CC) $(LIB_PATHS) $(SOURCES) -o $(EXECUTABLE) $(INCLUDES) $(LIBS) 

clean:
	rm -rf ./exec/*
