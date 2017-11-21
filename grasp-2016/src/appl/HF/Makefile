.SUFFIXES: .f90 .mod

BIN = $(GRASP)/bin
UTIL = HF

install: EXE 

EXE : $(BIN)/HF  \

OBJ : HF.o

$(BIN)/HF: HF.o
	$(FC) -o $(BIN)/HF HF.o 

.f90.o:
	$(FC) -c $(FC_FLAGS) -I $(GRASP_INCLUDES) $< -o $@

clean:
	-rm -f *.o core *.mod

