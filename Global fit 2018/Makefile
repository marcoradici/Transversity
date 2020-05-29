CORE = Evolution.cc
COREo = Evolution.o

PDF = u_valence.f90 d_valence.f90
PDFo = u_valence.o d_valence.o

MODULES = params.f90 evo_grid.f90
MODULESo = params.o evo_grid.o

SOURCE = h1grid.f90
SOURCEo = h1grid.o

EXE = h1grid.exe

F90 = gfortran
CXX = g++
FLAGS = -O3

compile:
	$(F90) $(FLAGS) -c $(PDF) $(MODULES) $(SOURCE)
	$(CXX) $(FLAGS) -fPIC -std=c++11 `apfelxx-config --cppflags` -c $(CORE)

link:
	$(F90) -o $(EXE) $(COREo) $(PDFo) $(MODULESo) $(SOURCEo) -lc++ -lstdc++ `apfelxx-config --ldflags`

cleanmod:
	rm *.mod
	
cleano:
	rm *.o
	
cleanall:
	rm *.o *.mod h1grid.exe
