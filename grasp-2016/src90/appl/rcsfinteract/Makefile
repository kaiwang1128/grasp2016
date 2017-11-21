.SUFFIXES: .f90 .mod

EXE = rcsfinteract
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
BINFILE = $(BINDIR)/$(EXE)
SRCLIBDIR = ../../lib
LAPACK90_DIR = $(GRASPLIB)
MODDIR = ${SRCLIBDIR}/libcommons
MODL9290 = ${SRCLIBDIR}/lib9290
MODLRANG90 = ${SRCLIBDIR}/librang90_inter
F90 = $(FC)
F90_FLAGS = $(FC_FLAGS)
# MODLNJG = ${SRCLIBDIR}/libnjgraf
#GRASPLIBS =-l9290 -lcommons90 -lnjgraf90 -llapack90
#GRASPLIBS =-l9290 -lcommons90 -lrang90_inter -lmcp90
GRASPLIBS =-l9290 -lcommons90 -lrang90_inter
LAPACK90_LIB = -llapack90

APP_LIBS = -L${GRASPLIB} ${GRASPLIBS}

#   Define data types
VASTO = ${MODDIR}/vast_kind_param_M.o

APP_OBJ= \
	chkplt_I.o getinf_I.o lodres_I.o  Interact_MR_I.o \
	lodcsl_MR_I.o lodcsl_CSF_I.o \
	set_CSF_list_I.o set_CSF_number_I.o \
\
        RCSFinteract.o chkplt.o getinf.o lodres.o Interact_MR.o \
        lodcsl_MR.o lodcsl_CSF.o \
	set_CSF_list.o set_CSF_number.o

$(EXE): $(APP_OBJ)
	$(F90) -o $(BINFILE) $(APP_OBJ) $(APP_LIBS) \
              -L ${LAPACK90_DIR} ${LAPACK90_LIB} $(APP_LIBS)

.f90.o:
	$(F90) -c $(F90_FLAGS) $< -I $(MODDIR) -I ${MODL9290} -I $(MODLRANG90) -I $(MODDIR) -o $@

.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

clean:
	-rm -f *.o core  *.mod
