.SUFFIXES: .f90 .mod

LIBDIR = $(GRASP)/lib
F90 = $(FC)
F90_FLAGS = $(FC_FLAGS)
LIBA = liblapack90.a

install:
	make -f Makefileu.inc install
	make -f Makefiled.inc install
	make -f Makefileb.inc install
	@echo " Building " $(LIBA)
	ar -crvs $(LIBA) *.o 
	mv $(LIBA) $(LIBDIR)

clean:
	make -f Makefileu.inc clean
	make -f Makefiled.inc clean
	make -f Makefileb.inc clean


