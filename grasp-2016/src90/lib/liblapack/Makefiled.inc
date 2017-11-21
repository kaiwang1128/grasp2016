.SUFFIXES: .f90 .mod

MODDIR = ../libcommons

RM= /bin/rm -f

VASTO = ../libmod/vast_kind_param_M.o
MODOBJ = \

OBJ = \
dlae2.o   dlaebz.o  dlaev2.o  dlagtf.o  dlagts.o  dlansp.o  dlanst.o  dlapy2.o\
dlarf.o   dlarfg.o  dlarnv.o  dlartg.o  dlaruv.o  dlascl.o  dlaset.o  dlasr.o\
dlasrt.o  dlassq.o  dopgtr.o  dopmtr.o  dorg2l.o  dorg2r.o  dspevx.o  dsptrd.o\
dstebz.o  dstein.o  dsteqr.o  dsterf.o

install:$(OBJ)
	@echo " Building install lapack90"
#	ar -crvs $(LIB) $? # $(LIBOBJ)

$(OBJ) : ${MODOBJ}
${MODOBJ} : ${VASTO}

.f90.o:
	$(F90) -c $(F90_FLAGS) -I ${MODDIR} $< # -o $@

clean:
	-$(RM) $(OBJ)

