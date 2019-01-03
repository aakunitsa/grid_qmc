MF = Makefile
CC = g++
#CFLAGS = -march=native -ffast-math -O3 -std=c++14 -I$(HOME)/local/include -DMT64 
CFLAGS = -ffast-math -O3 -std=c++14 -I$(HOME)/local/include -DMT64 
#CFLAGS = -O0 -std=c++14 -I$(HOME)/local/include -DMT64 
#LFLAGS = -L$(HOME)/local/lib -lopenblas -lpthread -lgfortran
#LFLAGS = $(HOME)/local/lib/libopenblas.a 
LFLAGS =  $(HOME)/local/lib/libgslcblas.a $(HOME)/local/lib/libgsl.a

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
	$(CC) $(LFLAGS) -o $@ $(OBJ) 

$(OBJ):	$(MF)

clean:
	rm -f $(OBJ) $(EXE) 
