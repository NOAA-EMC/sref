SHELL=/bin/sh
#
SRCS=	WGTMKR.f POLATWGT.f makgds.f w3fi63_g218_12km.f w3fi71_g218_12km.f 

OBJS=	WGTMKR.o POLATWGT.o makgds.o w3fi63_g218_12km.o w3fi71_g218_12km.o 

# Tunable parameters
#
# CF		Name of the fortran compiling system to use
# LDFLAGS	Flags to the loader
# LIBS		List of libraries
# CMD		Name of the executable
# PROFLIB	Library needed for profiling
#

FC =		ifort -convert big_endian -g -traceback
LDFLAGS =
LIBS = ${IP_LIB4} ${W3NCO_LIB4} ${SP_LIB4}

CMD =		wgtmkr_nmmb.x
PROFLIB =	

# To perform the default compilation, use the first line
# To compile with flowtracing turned on, use the second line
# To compile giving profile additonal information, use the third line
# WARNING:  SIMULTANEOUSLY PROFILING AND FLOWTRACING IS NOT RECOMMENDED 
#FFLAGS =  -qrealsize=4 -qmaxmem=-1 -qsource -qnosave -qspillsize=548
FFLAGS =  -O2 -I${IP_INC4}

# Lines from here on down should not need to be changed.  They are the
# actual rules which make uses to build awprgn.xc.
#
all:		$(CMD)

$(CMD):		$(OBJS)
	$(FC) $(LDFLAGS) -o $(@) $(OBJS) $(LIBS)

# Make the profiled version of the command and call it awprgn.xc.prof
#
$(CMD).prof:	$(OBJS)
	$(FC) $(LDFLAGS) -o $(@) $(OBJS) $(PROFLIB) $(LIBS)

clean:
	-rm -f $(OBJS) *.lst $(CMD)

clobber:	clean
	-rm -f $(CMD) $(CMD).prof

void:	clobber
	-rm -f $(SRCS) makefile
