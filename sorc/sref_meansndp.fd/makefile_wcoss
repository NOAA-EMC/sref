SHELL=/bin/sh
#
SRCS=	BFRHDR.f BFRIZE.f CALHEL.f CALWXT.f SMIX.f SNDPST.f VAP.f \
        WETBLB.f CALVIS.f CALWXT_RAMER.f CALWXT_EXPLICIT.f \
        CALWXT_REVISED.f CALWXT_BOURG.f CALWXT_DOMINANT.f \
        CALWXT_DOMINANT2.f

OBJS=	BFRHDR.o BFRIZE.o CALHEL.o CALWXT.o SMIX.o SNDPST.o VAP.o \
        WETBLB.o CALVIS.o CALWXT_RAMER.o CALWXT_EXPLICIT.o \
        CALWXT_REVISED.o CALWXT_BOURG.o CALWXT_DOMINANT.o \
        CALWXT_DOMINANT2.o

# Tunable parameters
#
# FC		Name of the fortran compiling system to use
# LDFLAGS	Flags to the loader
# LIBS		List of libraries
# CMD		Name of the executable
# PROFLIB	Library needed for profiling
#
FC =		ifort -convert big_endian
LDFLAGS = 
CMD =	 sref_meansndp	
PROFLIB =

LIBS =           ${BUFR_LIB4} 

FFLAGS = -O

# Lines from here on down should not need to be changed.  They are the
# actual rules which make uses to build a.out.
#
all:		$(CMD)

$(CMD):		$(OBJS)
	$(FC) $(LDFLAGS) -o $(@) $(OBJS) $(LIBS)

# Make the profiled version of the command and call it a.out.prof
#
$(CMD).prof:	$(OBJS)
	$(FC) $(LDFLAGS) -o $(@) $(OBJS) $(PROFLIB) $(LIBS)

clean:
	-rm -f $(OBJS) *.lst

clobber:	clean
	-rm -f $(CMD) $(CMD).prof

void:	clobber
	-rm -f $(SRCS) makefile
