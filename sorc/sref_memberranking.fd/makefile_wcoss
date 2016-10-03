SHELL=/bin/sh
SRCS 	= ranking.f
OBJS	= ranking.o

#FC     = ifort -O3 -convert big_endian
#FC     = ifort -O3 -stack_size:0x2400000000
#FC     = ifort -O3 -mcmodel=medium
FC 	= ifort
LDFLAGS =
LIBS 	= ${W3NCO_LIB4} ${W3EMC_LIB4} ${BACIO_LIB4}
CMD 	= sref_ranking

all:	$(CMD)

$(CMD):	 ${OBJS}
	$(FC) $(LDFLAGS) -o $(@) $(OBJS) $(LIBS)

clean:
	rm -f *.o
