MF = Makefile
CC = g++
#CFLAGS = -march=native -ffast-math -O3 -std=c++14 -I$(HOME)/local/include -DMT64 
CFLAGS = -O3 -std=c++17 -I$(HOME)/local/include -DMT64 
#LFLAGS = $(HOME)/local/lib/libopenblas.a 
#LFLAGS = $(HOME)/local/lib/libgsl.a $(HOME)/local/lib/libgslcblas.a -larmadillo -L$(HOME)/local/lib
#LFLAGS =   -lgsl -lgslcblas -larmadillo -lm
LFLAGS =   -lgsl -lgslcblas -lm -larmadillo -llapack -lblas 
LIBS = -L$(HOME)/local/lib

EXE = vmc.x

SRC= $(wildcard *.cpp)
SRCEXT= lebedev_grid/sphere_lebedev_rule.cpp

#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .cpp .o

OBJ = $(SRC:.cpp=.o)  $(SRCEXT:.cpp=.o)

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@

all:	$(EXE)

$(EXE):	$(OBJ)
	$(CC) $(LIBS) -o $@ $(OBJ) $(LFLAGS) 

$(OBJ):	$(MF)

clean:
	rm -f $(OBJ) $(EXE) 
