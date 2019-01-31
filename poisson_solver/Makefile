MF = Makefile
FC = gfortran
#FFLAGS = -O0 -fcheck=all
FFLAGS = -O2

SRC= pmodules_stripped.f90 ppoisson.f90 plebedev.f

#
# No need to edit below this line
#

.SUFFIXES: 
.SUFFIXES: .f .f90 .o

OBJF90  = $(SRC:.f90=.o) 
OBJ = $(OBJF90:.f=.o)

.f90.o:
	$(FC) $(FFLAGS) -c $< -o $@


all:	$(OBJ)

test:
	@echo $(OBJ)

clean:
	rm -fv *.o *.mod *~


$(OBJ):	$(MF)

