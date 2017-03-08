# Use this Makefile with make

# Executable name
CMD = align.exe
INC=-I$(MKLROOT)/include -I/export/zimmerman/mjafari/ZStruct-S/source-code/align-structures -I.
#INC=-I$(MKLROOT)/include

# -------- description of DFLAGS ---------------


# -------- Define environmental variable C_COMPILER -----------
# Make sure it is defined
#          ifeq ($(strip$(FORTRAN_COMPILER)),)
# Otherwise you can define it here also by uncommenting next line
#FC = icpc -gcc-name=gcc-4.5 -std=c++0x -openmp -I$(MKLROOT)/include
#FC = g++ -gdwarf-3 -std=c++11 -I$(MKLROOT)/include
FC = g++ -Wall -gdwarf-3 -std=c++11 $(INC)
# FC = g++ -fopenmp -I$(MKLROOT)/include
DFLAGS =  #-Define the cpp flags to be used
#DFLAGS =  #-Define the cpp flags to be used
OFLAGS =  # optimization

#Intel parallel openmp (only w/icpc compiler)
#LINKERFLAGS =  -L$(MKLROOT)/lib/em64t -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lm
#LINKERFLAGS =  -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread 
LINKERFLAGS =  -Wl,--no-as-needed -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm
# MAC OS linkers
#LINKERFLAGS = -lm -framework Accelerate



#
# Implicit rules for handling src files
#  ( uses static pattern rules, see info for make )
.c.o:
	$(FC) -c -g $(DFLAGS) -Wimplicit $<
.cpp.o:
	$(FC) -c -g $(DFLAGS) $<

OBJECTS = main.o icoord.o pTable.o stringtools.o utils.o  mm_grad.o align.o  Surface.o BindingSite.o Atom.o Molecule.o Coordinates.o

$(CMD) : $(OBJECTS)
	$(FC) $(DEBUG_FLAGS) $(OFLAGS) $(LINKERFLAGS) $(OBJECTS)  -o ./$(CMD)

clean:
	/bin/rm -f *.o *.i *.mod *.exe output-* a.out make.log
	
all: clean build

build: $(OBJECTS)
	$(FC) $(DEBUG_FLAGS) $(OFLAGS) $(LINKERFLAGS) $(OBJECTS)  -o ./$(CMD)

run:
	./$(CMD)

cleano:
	rm -f *.o *.i

depend :
	g++ -MM *.cpp *.c >> Makefile 

# DO NOT DELETE created with g++ -MM *.cpp *.c
main.o: main.cpp icoord.h align.h
align.o: align.cpp align.h icoord.h utils.h
icoord.o: icoord.cpp icoord.h 
mm_grad.o: mm_grad.cpp icoord.h
pTable.o: pTable.cpp pTable.h
stringtools.o: stringtools.cpp stringtools.h
utils.o: utils.cpp utils.h
Surface.o: Surface.cpp Surface.h 
BindingSite.o: BindingSite.cpp BindingSite.h
Atom.o: Atom.cpp Atom.h
Molecule.o: Molecule.cpp Molecule.h
Coordinates.o: Coordinates.cpp Coordinates.h
