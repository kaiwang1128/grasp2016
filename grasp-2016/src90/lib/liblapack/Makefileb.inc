.SUFFIXES: .f90 .mod

MODDIR = ../libcommons
VASTO = $(MODDIR)/vast_kind_param_M.o

RM= /bin/rm -f

MODOBJ = \


OBJ = \
   dasum.o    daxpy.o    dcopy.o    ddot.o     dgbmv.o    dgemm.o    dgemv.o \
   dger.o     dmach.o    dnrm2.o    drot.o     drotg.o    dsbmv.o    dscal.o \
   dspmv.o    dspr.o     dspr2.o    dswap.o    dsymm.o    dsymv.o    dsyr.o \
   dsyr2.o    dsyr2k.o   dsyrk.o    dtbmv.o    dtbsv.o    dtpmv.o    dtpsv.o \
   dtrmm.o    dtrmv.o    dtrsm.o    dtrsv.o    gather.o   idamax.o   iinit.o \
   scatter.o  iinit.o    dgthr.o    droti.o    dsctr.o

install: $(OBJ)
	@echo " Building " $(LIB)
#	ar -crvs $(LIB) $? # $(LIBOBJ)

$(OBJ) : ${MODOBJ}
${MODOBJ} : ${VASTO}

.f90.o:
	$(F90) -c $(F90_FLAGS) -I ${MODDIR} $< # -o $@

clean:
	-$(RM) -f $(OBJ) 

