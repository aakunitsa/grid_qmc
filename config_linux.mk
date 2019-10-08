CFLAGS = -march=native -ffast-math -O2 -std=c++17 -fopenmp -I$(HOME)/local/include -I$(HOME)/arma_latest/include -DMT64 -DARMA_DONT_USE_WRAPPER -DNORMAL_POLYMER1
#CFLAGS = -O0 -std=c++17 -I$(HOME)/local/include -I$(HOME)/arma_latest/include -DMT64 -DARMA_DONT_USE_WRAPPER -DAUXBAS
FFLAGS = -O2 -cpp
#FFLAGS = -O0 
#LFLAGS = $(HOME)/local/lib/libopenblas.a 
#LFLAGS = $(HOME)/local/lib/libgsl.a $(HOME)/local/lib/libgslcblas.a -larmadillo -L$(HOME)/local/lib
#LFLAGS =   -lgsl -lgslcblas -lm -larmadillo -llapack -lblas -lgfortran
LFLAGS =   -lgsl -lgslcblas -lm -llapack -lblas -lgfortran -lgomp
#LFLAGS =   -lgsl -lgslcblas -lm -llapack -lblas # Would not work 
LIBS = -L$(HOME)/local/lib

