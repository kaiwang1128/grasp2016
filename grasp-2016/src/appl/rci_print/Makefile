EXE = rci_print
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
GRASPBLOCKD = ${GRASPLIB}/graspBD.o
BINFILE = $(BINDIR)/$(EXE)

APP_LIBS = ${GRASPBLOCKD} -L${GRASPLIB} ${NEWGRASPLIBS} 

APP_OBJ= \
        auxblk.o bessel.o breid.o brint1.o brint2.o brint3.o brint4.o \
        brint5.o brint6.o brintf.o brra.o cxk.o engout.o evcout.o funk.o \
        funl.o fzalf.o genintrk.o genmat.o hmout.o hovlap.o iabint.o indtpi.o \
        keint.o klamaq.o mohr.o ncharg.o qed.o rkint.o rkintc.o setdbg.o \
        setsum.o skint.o snrc.o talk.o triangrk.o vac2.o vac4.o vacpol.o \
        vint.o vinti.o vpint.o vpintf.o wghtd5.o zkf.o \
        dnicmv.o genmat2.o getcid.o iniestdm.o iniestsd.o lodmix.o lodres.o \
        maneig.o matrix.o rci92.o  setcsl.o setham_gg.o setmix.o setres.o \
        spodmv.o strsum.o shield.o qed_slfen.o lodstatetune.o genintbreit1.o \
        genintbreit2.o triangbreit1.o triangbreit2.o


$(EXE): $(APP_OBJ)
	$(FC) -o $(BINFILE) $(FC_LD) $(APP_OBJ) $(APP_LIBS) \
                 -L ${LAPACK_DIR} ${LAPACK_LIB}

.f.o:
	$(FC) -c $(FC_FLAGS) -I $(GRASP_INCLUDES) $< -o $@

clean:
	-rm -f *.o core 

