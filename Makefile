MF = Makefile
CC = g++
FC = gfortran
CFLAGS = -march=native -ffast-math -O2 -std=c++17 -I$(HOME)/local/include -I$(HOME)/arma_latest/include -DMT64 -DARMA_DONT_USE_WRAPPER -DAUXBAS
FFLAGS = -O2 
#LFLAGS = $(HOME)/local/lib/libopenblas.a 
#LFLAGS = $(HOME)/local/lib/libgsl.a $(HOME)/local/lib/libgslcblas.a -larmadillo -L$(HOME)/local/lib
#LFLAGS =   -lgsl -lgslcblas -lm -larmadillo -llapack -lblas -lgfortran
LFLAGS =   -lgsl -lgslcblas -lm -llapack -lblas -lgfortran
LIBS = -L$(HOME)/local/lib 

EXE = vmc.x

SOLVER = ./poisson_solver

SRC= $(wildcard *.cpp)
SRCEXT= lebedev_grid/sphere_lebedev_rule.cpp
FORTRAN_SRC_TOP=$(wildcard *.f90)

FORTRAN_SRC= $(wildcard $(SOLVER)/*.f90)  $(wildcard $(SOLVER)/*.f)
FORTRAN_OBJ_=$(FORTRAN_SRC:.f90=.o)
FORTRAN_OBJ=$(FORTRAN_OBJ_:.f=.o)


#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .f90 .cpp .o

OBJ = $(SRC:.cpp=.o)  $(SRCEXT:.cpp=.o) $(FORTRAN_SRC_TOP:.f90=.o)

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@

.f90.o:
	$(FC) $(FFLAGS) -c $< -o $@

all:	$(EXE)

test:
	@echo $(OBJ)

$(EXE):	$(OBJ)
	cd $(SOLVER); make
	$(CC) $(LIBS) -o $@ $(OBJ) $(FORTRAN_OBJ) $(LFLAGS) 

$(OBJ):	$(MF)

clean:
	cd $(SOLVER); make clean
	rm -fv $(OBJ) $(EXE) *~
