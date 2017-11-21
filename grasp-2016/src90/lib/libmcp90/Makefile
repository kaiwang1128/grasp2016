.SUFFIXES: .f90 .mod

LIBDIR = $(GRASP)/lib
#LIBA = $(LIBDIR)/libmcp90.a
LIBA = libmcp90.a

LIBMOD = -lcommons90
LIBOBJ = $(OBJ)
MODDIR = ../libcommons
MODL92 = ../lib9290
VASTO = $(MODDIR)/vast_kind_param_M.o
F90 = $(FC)
F90_FLAGS = $(FC_FLAGS)

RM = /bin/rm -f

MODOBJ = \


OBJ = \
   cxk_I.o snrc_I.o talk_I.o \
   \
   cxk.o snrc.o talk.o

install : $(LIBA)
$(LIBA) : $(LIBOBJ)
	@echo " Building " $(LIBA)
	ar -crvs $(LIBA)
	mv $(LIBA) $(LIBDIR)

clean:
	-rm -f *.o *.mod core

.f90.o:
	$(F90) -c $(F90_FLAGS) $< -I  $(MODDIR) -I $(MODL92) -I .  -o $@
