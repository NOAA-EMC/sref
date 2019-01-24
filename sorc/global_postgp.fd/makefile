SHELL	=/bin/sh
SRCM	=postgp.f90
INCM	=
SRCSP	=postgp_module.f90 funcphys.f90 mptran.f90 physcons_v.f90 machine.f90
OBJSP	=postgp_module.o funcphys.o mptran.o physcons_v.o machine.o
SRCSS	=postgp_module.f90 funcphys.f90 mptran.std.f90 physcons_v.f90 machine.f90
OBJSS	=postgp_module.o funcphys.o mptran.std.o physcons_v.o machine.o
FCP	=ifort
FCS	=ifort
OPTS	= -O -fp-model strict
FREE	= -FR
TRAPS	= 
DEBUG	= -g -traceback
FFLAGS	= $(OPTS) $(FREE) $(TRAPS) $(DEBUG)
LDFLAGS	= -qopenmp
LIBS =  ${SIGIO_LIB4} ${BACIO_LIB4} ${W3NCO_LIB4} ${IP_LIB4} ${SP_LIB4}
CMDP	=global_postgp
CMDS	=global_postgs
both:	$(CMDP) $(CMDS)
$(CMDP):	$(SRCM) $(INCM) $(OBJSP)
	$(FCP) -I $(SIGIO_INC4) $(FFLAGS) $(LDFLAGS) $(SRCM) $(OBJSP) $(LIBS) -o $(CMDP)
$(CMDS):	$(SRCM) $(INCM) $(OBJSS)
	$(FCS) -I $(SIGIO_INC4) $(FFLAGS) $(LDFLAGS) $(SRCM) $(OBJSS) $(LIBS) -o $(CMDS)
postgp_module.o:	postgp_module.f90
	$(FCS) $(FFLAGS) -c postgp_module.f90
funcphys.o:	funcphys.f90 physcons_v.o machine.o
	$(FCS) $(FFLAGS) -c funcphys.f90
physcons_v.o:	physcons_v.f90 machine.o
	$(FCS) $(FFLAGS) -c physcons_v.f90
machine.o:	machine.f90
	$(FCS) $(FFLAGS) -c machine.f90
mptran.o:	mptran.f90 machine.o
	$(FCP) $(FFLAGS) -c mptran.f90
mptran.std.o:	mptran.std.f90 machine.o
	$(FCS) $(FFLAGS) -c mptran.std.f90
clean:
	-rm -f $(OBJSS) $(OBJSP) *.mod
clobber:        clean
	-rm -f $(CMDS) $(CMDP)
void:   clobber
	-rm -f $(SRCS) $(SRCM) $(INCM) makefile
