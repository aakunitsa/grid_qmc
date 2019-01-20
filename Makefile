MF = Makefile
CC = g++
FC = gfortran
#CFLAGS = -march=native -ffast-math -O3 -std=c++14 -I$(HOME)/local/include -DMT64 
CFLAGS = -O3 -std=c++17 -I$(HOME)/local/include -DMT64 
FFLAGS = -O3
#LFLAGS = $(HOME)/local/lib/libopenblas.a 
#LFLAGS = $(HOME)/local/lib/libgsl.a $(HOME)/local/lib/libgslcblas.a -larmadillo -L$(HOME)/local/lib
#LFLAGS =   -lgsl -lgslcblas -larmadillo -lm
LFLAGS =   -lgsl -lgslcblas -lm -larmadillo -llapack -lblas -lgfortran
LIBS = -L$(HOME)/local/lib 

EXE = vmc.x

SRC= $(wildcard *.cpp)
SRCEXT= lebedev_grid/sphere_lebedev_rule.cpp

FORTRAN_SRC= $(wildcard *.f90)

#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .cpp .o

OBJ = $(SRC:.cpp=.o)  $(SRCEXT:.cpp=.o) $(FORTRAN_SRC:.f90=.o)

.f90.o:
	$(FC) $(FFLAGS) -c $< -o $@

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@

all:	$(EXE)

$(EXE):	$(OBJ)
	$(CC) $(LIBS) -o $@ $(OBJ) $(LFLAGS) 

$(OBJ):	$(MF)

clean:
	rm -f $(OBJ) $(EXE) 
