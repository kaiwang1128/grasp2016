.SUFFIXES: .f90 .mod

MODDIR = ../libcommons
VASTO = $(MODDIR)/vast_kind_param_M.o

MODOBJ = \

F90_FLAGSU= $(F90_FLAGS) -O0

OBJ = \
xerbla.o slamch.o dlamc1.o  dlamc2.o  dlamc3.o  dlamc4.o  dlamc5.o  dlamch.o  ilaenv.o  lsame.o

install:$(OBJ)
	@echo " Building install lapack90"
#	ar -crvs $(LIB) $?

$(OBJ) : ${MODOBJ}
${MODOBJ} : ${VASTO}

.f90.o:
	$(F90) -c $(F90_FLAGSU) -I ${MODDIR} $< # -o $@

clean:
	-rm -f $(OBJ) $(LIBDIR)/liblapack90.a core 

