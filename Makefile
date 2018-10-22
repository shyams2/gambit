COMPILER_OPTIONS= -Wall -fopenmp -O3 -std=c++11 -g
# Using h5c++ to make possible file writing using HDF5:
CC              = g++ $(COMPILER_OPTIONS)
# Change this to:
# -lafcpu    - For CPU backend
# -lafcuda   - For CUDA backend
# -lafopencl - For OpenCL backend
# -laf       - For unified backend
LIBS            = -laf -lhdf5
LIB_PATHS       =  -L $(HDF5_DIR)/lib -L $(AF_PATH)/lib
INCLUDES        = -I $(AF_PATH)/include -I $(EIGEN_PATH) -I $(HIGHFIVE_PATH)/include -I ./header -I $(HDF5_DIR)/include
SOURCES         = ./test/testMatrixData.cpp
EXECUTABLE      = ./exec/testMatrixData

all:
	$(CC) $(LIB_PATHS) $(SOURCES) -o $(EXECUTABLE) $(INCLUDES) $(LIBS) 

clean:
	rm -rf ./exec/*
