.SUFFIXES: .f90 .mod

LIBDIR = $(GRASP)/lib
LIBA = $(LIBDIR)/lib9290.a

LIBMOD = -lcommons90
LIBOBJ = $(OBJ)
MODDIR = ../libcommons
VASTO = $(MODDIR)/vast_kind_param_M.o
F90 = $(FC)
F90_FLAGS = $(FC_FLAGS)

RM = /bin/rm -f

MODOBJ = \


OBJ = \
   alcbuf_I.o clrx_I.o convrt_double_I.o convrt_I.o cord_I.o cre_I.o factt_I.o \
   getyn_I.o ichkq1_I.o ichkq2_I.o ichop_I.o iq_I.o  \
   ispar_I.o itjpo_I.o itrig_I.o jcup_I.o jqs_I.o  \
   lodcsh_I.o lodcsl_I.o lodplt_I.o lodiso_I.o openfl_I.o pack_I.o \
   parsjl_I.o prsrcn_I.o prsrsl_I.o speak_I.o\
   setiso_I.o setcsla_I.o setqna_I.o starttime_I.o stoptime_I.o \
   \
   alcbuf.o clrx.o convrt_double.o convrt.o cord_I.o cre.o factt.o \
   getyn.o ichkq1.o ichkq2.o ichop.o iq.o  \
   ispar.o itjpo.o itrig.o jcup.o jqs.o  \
   lodcsh.o lodcsl.o lodplt.o lodiso.o openfl.o pack.o \
   parsjl.o prsrcn.o prsrsl.o speak.o\
   setiso.o setcsla.o setqna.o starttime.o stoptime.o

install : $(LIBA)
$(LIBA) : $(LIBOBJ) 
	@echo " Building " $(LIBA)
	ar -crvs $(LIBA)  $?

clean:
	-rm -f *.o *.mod core

.f90.o:
	$(F90) -c $(F90_FLAGS) $< -I  $(MODDIR) -I .  -o $@
