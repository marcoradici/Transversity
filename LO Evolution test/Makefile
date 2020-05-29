# order of modules is important for dependencies
CORE = Evolution.cc
COREo = Evolution.o

PDF = u_valence.f90 d_valence.f90
PDFo = u_valence.o d_valence.o

SOURCE = evotest.f90
SOURCEo = evotest.o

EXE = evotest.exe

F90 = gfortran
CXX = g++
FLAGS = -O3

compile:
	$(F90) $(FLAGS) -c $(SOURCE) $(PDF)
	$(CXX) $(FLAGS) -fPIC -std=c++11 `apfelxx-config --cppflags` -c $(CORE)
	
link:
	$(F90) -o $(EXE) $(COREo) $(PDFo) $(SOURCEo) -lc++ -lstdc++ `apfelxx-config --ldflags`

cleanobj:
	rm *.o 

allclean:
	rm *.o *.txt
