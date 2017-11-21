.SUFFIXES: .f90 .mod  
EXE = rcsfsymexpandfull
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
GRASPBLOCKD = ${GRASPLIB}/graspBD.o
BINFILE = $(BINDIR)/$(EXE)
FC = gfortran
FC_FLAGS = -fno-automatic

APP_OBJ= \
       symexpand_mod.o rcsfsymexpandfull.o testmodule.o

$(EXE): $(APP_OBJ)
	$(FC) -o $(BINFILE) $(APP_OBJ) -L$(GRASPLIB) $(GRASPLIBS) 

.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@
.f90.o:
	$(FC) -c $(FC_FLAGS) $< -o $@
clean:
	-rm -f *.o core 
