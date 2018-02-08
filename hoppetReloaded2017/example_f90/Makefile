# Makefile generated automatically with
# ../scripts/makef95makefile "tabulation_example" "-L../src" "-lhoppet_v1" "-I../src" "FC=gfortran" "FFLAGS=-fno-whole-file"
# default program to compile
PROG =	tabulation_example

ALLPROG = 	a sumrules tabulation_example tabulation_example_2 \
	tabulation_example_streamlined

# This will be used one day...
ALLPROGSRC =	a.f90 sumrules.f90 tabulation_example.f90 tabulation_example_2.f90 \
	tabulation_example_streamlined.f90

ALLPROGOBJ =	a.o sumrules.o tabulation_example.o tabulation_example_2.o \
	tabulation_example_streamlined.o

SRCS =	

POSTSRCS =	

OBJS =	

POSTOBJS =	
POSTLIB = 

LIBS = -L../src -lhoppet_v1  

CC = cc
CFLAGS = -O
FC = gfortran
FFLAGS = -fno-whole-file -I../src 
F90 = gfortran
F90FLAGS = -fno-whole-file -I../src 
LDFLAGS =  

all: $(ALLPROG)

ALL:  $(ALLPROG)

a: a.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o a a.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

sumrules: sumrules.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o sumrules sumrules.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

tabulation_example: tabulation_example.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o tabulation_example tabulation_example.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

tabulation_example_2: tabulation_example_2.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o tabulation_example_2 tabulation_example_2.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

tabulation_example_streamlined: tabulation_example_streamlined.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o tabulation_example_streamlined tabulation_example_streamlined.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

libclean:
	rm -f  $(ALLPROGOBJ) $(OBJS) $(POSTOBJS) 

clean:
	rm -f  $(ALLPROGOBJ) $(OBJS) $(POSTOBJS) *.mod *.d

distclean: realclean
	rm -f Makefile
realclean:
	rm -f  $(ALLPROG) $(ALLPROGOBJ) $(OBJS) $(POSTOBJS) *.mod *.d

make:
	../scripts/makef95makefile  "tabulation_example" "-L../src" "-lhoppet_v1" "-I../src" "FC=gfortran" "FFLAGS=-fno-whole-file"

install:

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

a.o: 
sumrules.o: 
tabulation_example.o: 
tabulation_example_2.o: 
tabulation_example_streamlined.o: 
