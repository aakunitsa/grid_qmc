CFLAGS = -march=native -ffast-math -O2 -std=c++17 -fopenmp -I$(HOME)/local/include  -DMT64 -DARMA_DONT_USE_WRAPPER -DNORMAL_POLYMER1
FFLAGS = -O2 -cpp
LFLAGS =   -lgsl -lgslcblas -lm -llapack -lblas -lgfortran -lgomp
LIBS = -L$(HOME)/local/lib

