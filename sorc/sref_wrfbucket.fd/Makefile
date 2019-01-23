SHELL=/bin/sh

SRCS = WRFBUCKET.F

FC =   ifort
LIBS = ${W3NCO_LIB4} ${BACIO_LIB4}

TARGETS =  sref_pcpbucket_g212 sref_pcpbucket_g216 sref_pcpbucket_g221 sref_pcpbucket_g243 sref_pcpbucket_g132

all: $(TARGETS)


sref_pcpbucket_g212: $(SRCS)
	cp parmeta_g212 parmeta
	$(FC) -o $(@) $(SRCS) $(LIBS)
	rm parmeta

sref_pcpbucket_g216: $(SRCS) 
	cp parmeta_g216 parmeta
	$(FC) -o $(@) $(SRCS) $(LIBS)
	rm parmeta

sref_pcpbucket_g221: $(SRCS) 
	cp parmeta_g221 parmeta
	$(FC) -o $(@) $(SRCS) $(LIBS)
	rm parmeta

sref_pcpbucket_g243: $(SRCS) 
	cp parmeta_g243 parmeta
	$(FC) -o $(@) $(SRCS) $(LIBS)
	rm parmeta

sref_pcpbucket_g132: $(SRCS) 
	cp parmeta_g132 parmeta
	$(FC) -o $(@) $(SRCS) $(LIBS)
	rm parmeta

clean:
	-rm -f *.o *.lst $(TARGETS)
