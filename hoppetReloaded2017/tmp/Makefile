# Makefile generated automatically with
# ../scripts/makef95makefile coeff_fn -L../src -lhoppet_v1 -I../src
# default program to compile
PROG =	coeff_fn

ALLPROG = 	coeff_fn delta_all delta_all_orig

# This will be used one day...
ALLPROGSRC =	coeff_fn.f90 delta_all.f90 delta_all_orig.f90

ALLPROGOBJ =	coeff_fn.o delta_all.o delta_all_orig.o

SRCS =	

POSTSRCS =	

OBJS =	

POSTOBJS =	
POSTLIB = 

LIBS = -L../src -lhoppet_v1  

CC = cc
CFLAGS = -O
FC = lf95
FFLAGS = -O -I../src 
F90 = lf95
F90FLAGS = -O -I../src 
LDFLAGS =  

all: $(ALLPROG)

ALL:  $(ALLPROG)

coeff_fn: coeff_fn.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o coeff_fn coeff_fn.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

delta_all: delta_all.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o delta_all delta_all.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

delta_all_orig: delta_all_orig.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o delta_all_orig delta_all_orig.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

libclean:
	rm -f  $(ALLPROGOBJ) $(OBJS) $(POSTOBJS) 

clean:
	rm -f  $(ALLPROGOBJ) $(OBJS) $(POSTOBJS) *.mod *.d

realclean:
	rm -f  $(ALLPROG) $(ALLPROGOBJ) $(OBJS) $(POSTOBJS) *.mod *.d

make:
	../scripts/makef95makefile  coeff_fn -L../src -lhoppet_v1 -I../src

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90 
	$(F90) $(F90FLAGS) -c $<

coeff_fn.o: 
delta_all.o: 
delta_all_orig.o: 
