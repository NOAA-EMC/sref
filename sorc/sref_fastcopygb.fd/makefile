################################################################################
# This makefile was produced by fmjif at 17:21:11 on 97/02/24.
SHELL=	/bin/sh
CMD=	fastcopygb
FOPTIP= -O3
FOPTW3=	-O3 ${W3EMC_INCd}
FOPTS=	-O3 
# For old sp lib:
#LOPTS=	
# For new sp lib (multi thread)
LOPTS= -openmp -auto
# For new sp lib (single thread)
#LOPTS= -mkl
INCS=	
OBJS=	copygb.o sv.polates0.o sv.polatev0.o sv.polates1.o sv.polatev1.o sv.polates2.o sv.polatev2.o putgben.o
LIBS=  ${W3NCO_LIBd}  ${W3EMC_LIBd} ${BACIO_LIB4} ${IP_LIBd} ${SP_LIBd}
################################################################################
$(CMD):	$(OBJS)
	ifort $(LOPTS) $(OBJS) $(LIBS) -o $(CMD)
$(OBJS):	$(INCS)
sv.polates0.o:	sv.polates0.f
	ifort -c $(FOPTIP) $*.f
sv.polatev0.o:	sv.polatev0.f
	ifort -c $(FOPTIP) $*.f
sv.polates1.o:	sv.polates1.f
	ifort -c $(FOPTIP) $*.f
sv.polatev1.o:	sv.polatev1.f
	ifort -c $(FOPTIP) $*.f
sv.polates2.o:	sv.polates2.f
	ifort -c $(FOPTIP) $*.f
sv.polatev2.o:	sv.polatev2.f
	ifort -c $(FOPTIP) $*.f
putgben.o:	putgben.f
	ifort -c $(FOPTW3) $*.f
copygb.o:	copygb.f
	ifort -c $(FOPTS) $*.f
clean:	
	-rm -f $(OBJS)
	ifort -c \$(FFLAGS) \$<
