SHELL=/bin/sh
SRCS 	= estimate_bias.f
OBJS	= estimate_bias.o

FC 	= ifort
LDFLAGS =
LIBS 	= ${W3NCO_LIB4} ${BACIO_LIB4}
CMD 	= sref_estimate_bias

all:	$(CMD)

$(CMD):	 ${OBJS}
	$(FC) $(LDFLAGS) -o $(@) $(OBJS) $(LIBS)

clean:
	-rm -f *.o
