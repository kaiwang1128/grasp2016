.SUFFIXES: .f90 .mod

EXE = rcsfzerofirst
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
BINFILE = $(BINDIR)/$(EXE)
SRCLIBDIR = ../../lib
LAPACK90_DIR = $(GRASPLIB)
MODDIR = ${SRCLIBDIR}/libcommons
MODL9290 = ${SRCLIBDIR}/lib9290
F90 = $(FC)
F90_FLAGS = $(FC_FLAGS)
#GRASPLIBS =-l9290 -lcommons90 -lrang90_inter
GRASPLIBS =-l9290 -lcommons90
LAPACK90_LIB = -llapack90

APP_LIBS = -L${GRASPLIB} ${GRASPLIBS}

#   Define data types
#VASTO = ${MODDIR}/vast_kind_param_M.o

APP_OBJ= \
	chkplt_I.o \
	lodcsl_Zero_I.o lodcsl_Part_I.o \
	set_CSF_ZFlist_I.o set_CSF_number_I.o \
\
        RCSFzerofirst.o chkplt.o \
        lodcsl_Zero.o lodcsl_Part.o \
	set_CSF_ZFlist.o set_CSF_number.o

$(EXE): $(APP_OBJ)
	$(F90) -o $(BINFILE) $(APP_OBJ) $(APP_LIBS) \
              -L ${LAPACK90_DIR} ${LAPACK90_LIB} $(APP_LIBS)

.f90.o:
	$(F90) -c $(F90_FLAGS) $< -I $(MODDIR) -I ${MODL9290} -I $(MODDIR) -o $@

.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

clean:
	-rm -f *.o core  *.mod
