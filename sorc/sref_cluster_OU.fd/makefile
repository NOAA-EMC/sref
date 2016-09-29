#!/bin/sh
FCMP    = ifort -O3
#FCMP    = ifort -O3 -mcmodel=medium
CMD 	= sref_cluster_OU
SRCS 	= main.f90 tree.f90 cophenetic.f90 resemblance.f90 CAlib.f90 Global.f90 agglomerative.f90 
FOPTS   = -convert big_endian -g -traceback
#FOPTS   = 
LIBS    =   ${BACIO_LIB4} ${W3NCO_LIB4} ${W3EMC_LIB4}
$(CMD):	 ${SRCS}
	$(FCMP) $(FOPTS) $(SRCS) $(LIBS) -o $(CMD)

#clean:
#	rm -f *.o

