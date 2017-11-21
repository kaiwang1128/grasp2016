.SUFFIXES: .f90 .mod

LIBDIR = $(GRASP)/lib
LIBA = $(LIBDIR)/librang90_inter.a

LIBMOD = -lcommons90
LIBOBJ = $(OBJ)
MODDIR = ../libcommons
MODL92 = ../lib9290
#MODLMCP90 = ../libmcp90
VASTO = $(MODDIR)/vast_kind_param_M.o
F90 = $(FC)
F90_FLAGS = $(FC_FLAGS)

RM= /bin/rm -f

MODOBJ = \


OBJ = \
   nmtejj_I.o perko1_I.o perko2_I.o \
   itrexg_I.o eile_I.o \
   dracah_I.o ittk_I.o nine0_I.o sixj1_I.o sixj35_I.o sixj4_I.o sixj_I.o \
   ixjtik_I.o nine_I.o sixj2_I.o sixj3_I.o sixj5_I.o \
   onescalarINT_I.o onescalar1INT_I.o onescalar2INT_I.o \
   el1INT_I.o el2INT_I.o el3INT_I.o el31INT_I.o el32INT_I.o el33INT_I.o el4INT_I.o el41INT_I.o \
   el5INT_I.o el51INT_I.o el52INT_I.o el53INT_I.o Interact_csf_I.o \
   coulomb_I.o diaga1_I.o diaga2_I.o diaga3_I.o diaga4_I.o \
   reco_I.o rec3_I.o reco3_I.o \
   reco2_I.o reco4_I.o recoonescalar_I.o snrc_I.o \
   \
   nmtejj.o perko1.o perko2.o \
   itrexg.o eile.o \
   dracah.o ittk.o nine0.o sixj1.o sixj35.o sixj4.o sixj.o \
   ixjtik.o nine.o sixj2.o sixj3.o sixj5.o \
   onescalarINT.o onescalar1INT.o onescalar2INT.o \
   el1INT.o el2INT.o el3INT.o el31INT.o el32INT.o el33INT.o el4INT.o el41INT.o \
   el5INT.o el51INT.o el52INT.o el53INT.o Interact_csf.o \
   coulomb.o diaga1.o diaga2.o diaga3.o diaga4.o \
   reco.o rec3.o reco3.o \
   reco2.o reco4.o recoonescalar.o snrc.o

install : $(LIBA)
$(LIBA) : $(LIBOBJ)
	@echo " Building " $(LIBA)
	ar -crvs $(LIBA)  $?

.f90.o:
	$(F90) -c $(F90_FLAGS) $< -I $(MODDIR) -I $(MODL92) -I $(MODDIR) -I .  -o $@

clean:
	-rm -f *.o *.mod core
