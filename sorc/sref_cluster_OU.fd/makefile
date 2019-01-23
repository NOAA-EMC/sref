#!/bin/sh
FCMP    = ifort -O3
CMD 	= sref_cluster_OU
SRCS 	= Global.f90 main.f90 tree.f90 cophenetic.f90 resemblance.f90 CAlib.f90 agglomerative.f90 
FOPTS   = -convert big_endian -g -traceback -I.
LIBS    =   ${BACIO_LIB4} ${W3NCO_LIB4}
$(CMD):	 ${SRCS}
	$(FCMP) $(FOPTS) $(SRCS) $(LIBS) -o $(CMD)
clean:
	rm -f *.o *.mod $(CMD)
