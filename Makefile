COMPILER_OPTIONS= -Wall -O3 -std=c++11 -g
CC              = g++ $(COMPILER_OPTIONS)
# Change this to:
# -lafcpu    - For CPU backend
# -lafcuda   - For CUDA backend
# -lafopencl - For OpenCL backend
# -laf       - For unified backend
LIBS            = -laf
LIB_PATHS       = -L $(AF_PATH)/lib
INCLUDES        = -I $(AF_PATH)/include -I $(EIGEN_PATH) -I ./header
SOURCES         = ./test/testID.cpp
EXECUTABLE      = ./exec/testID

all:
	$(CC) $(SOURCES) -o $(EXECUTABLE) $(INCLUDES) $(LIBS) $(LIB_PATHS)

clean:
	rm -rf ./exec/*
