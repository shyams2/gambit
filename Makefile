COMPILER_OPTIONS= -Wall -fopenmp -O3 -std=c++11 -g
# Using h5c++ to make possible file writing using HDF5:
CC              = h5c++ $(COMPILER_OPTIONS)
# Change this to:
# -lafcpu    - For CPU backend
# -lafcuda   - For CUDA backend
# -lafopencl - For OpenCL backend
# -laf       - For unified backend
LIBS            = -laf
LIB_PATHS       = -L $(AF_PATH)/lib
INCLUDES        = -I $(AF_PATH)/include -I $(EIGEN_PATH) -I $(HIGHFIVE_PATH)/include -I ./header
SOURCES         = ./test/testTree/testTree2D.cpp
EXECUTABLE      = ./exec/testTree

all:
	$(CC) $(SOURCES) -o $(EXECUTABLE) $(INCLUDES) $(LIBS) $(LIB_PATHS)

clean:
	rm -rf ./exec/*
