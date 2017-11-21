.SUFFIXES: .f90 .mod 
EXE = rci_symb
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
GRASPBLOCKD = ${GRASPLIB}/graspBD.o
BINFILE = $(BINDIR)/$(EXE)

APP_LIBS = ${GRASPBLOCKD} -L${GRASPLIB} ${NEWGRASPLIBS} 

APP_OBJ= \
        symexpand_mod.o  testmodule.o rcsfsymexpand.o \
        auxblk.o bessel.o breid.o brint1.o brint2.o brint3.o brint4.o \
        brint5.o brint6.o brintf.o brra.o cxk.o engout.o evcout.o funk.o \
        funl.o fzalf.o genintrk.o genmat.o hmout.o hovlap.o iabint.o indtpi.o \
        keint.o klamaq.o mohr.o ncharg.o qed.o rkint.o rkintc.o setdbg.o \
        setsum.o skint.o snrc.o talk.o triangrk.o vac2.o vac4.o vacpol.o \
        vint.o vinti.o vpint.o vpintf.o wghtd5.o zkf.o \
        dnicmv.o genmat2.o getcid.o iniestdm.o iniestsd.o lodmix.o lodres.o \
        maneig.o matrix.o rci92.o  setcsl.o setham_gg.o setmix.o setres.o \
        spodmv.o strsum.o shield.o qed_slfen.o lodstatetune.o genintbreit1.o \
        genintbreit2.o triangbreit1.o triangbreit2.o findtype.o onescalar11.o \
		  rkco_gg11.o transfer.o matrixblock12.o matrixblock2.o matrixblock22.o matrixblock13.o \
		  matrixblock23.o matrixblock24.o matrixblock3.o matrixblock14.o matrixblock15.o matrixblock25.o


$(EXE): $(APP_OBJ)
	$(FC) -o $(BINFILE) $(FC_LD) $(APP_OBJ) $(APP_LIBS) \
                 -L ${LAPACK_DIR} ${LAPACK_LIB}

.f.o:
	$(FC) -c $(FC_FLAGS) -I $(GRASP_INCLUDES) $< -o $@

clean:
	-rm -f *.o core  *.mod

