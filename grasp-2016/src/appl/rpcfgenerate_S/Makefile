EXE = rpcfgenerate_S
BINDIR = ${GRASP}/bin
BINFILE = $(BINDIR)/$(EXE)

APP_OBJ= \
        rpcfgenerate_S.o rcsfblock.o 

$(EXE): $(APP_OBJ)
	$(FC) -o $(BINFILE) $(FC_LD) $(APP_OBJ) 

.f.o:
	$(FC) -c $(FC_FLAGS) -I $(GRASP_INCLUDES) $< -o $@

clean:
	-rm -f *.o core 


