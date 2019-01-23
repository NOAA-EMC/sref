SHELL=	/bin/sh
LIBS=  ${W3NCO_LIB4} ${BACIO_LIB4}
SRCS=	sref_dvrtma_bias.f90
CMD=	sref_dvrtma_bias 
FC=	ifort
FFLAGS=	-O3
#FFLAGS=	-O3 -qsuffix=f=f90 -qlist -qsource -bnoquiet 

$(CMD):	$(SRCS)
	$(FC) $(SRCS) $(FFLAGS) -o $(@) $(OBJS) $(LIBS)

clean:
	rm -f *.o $(CMD)
