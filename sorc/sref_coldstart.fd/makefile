SHELL=/bin/sh

#----------------------------------------------------------------
#  Makefile for emcsfc_coldstart program.
#
#  Must be invoked from the 'make.sh' driver script.  See the
#  README.build file for details.
#
#  This program works best with 8 byte floats!!
#----------------------------------------------------------------

EXE = sref_coldstart_wrf

OBJS = program_setup.o \
       coldstart_wrf.o \
       read_data.o \
       interp_data.o \
       write_data.o \
       grib_utils.o \
       surface_chgres.o \
       machine.o

LIBS= $(LANDSFCUTIL_LIBd) $(G2_LIBd) $(NEMSIO_LIB) $(W3NCO_LIBd) \
      $(BACIO_LIB4) $(IP_LIBd) $(SFCIO_LIB4) $(SP_LIBd) $(NETCDF_LDFLAGS_F) \
      $(JASPER_LIB) $(PNG_LIB) $(Z_LIB)

$(EXE): ${OBJS}
	$(FCOMP) $(FFLAGS) $(LFLAGS) -o $(EXE) $(OBJS) $(LIBS)

coldstart_wrf.o: program_setup.o read_data.o interp_data.o write_data.o coldstart_wrf.f
	${FCOMP} ${FFLAGS} -c coldstart_wrf.f

program_setup.o: program_setup.f
	${FCOMP} ${FFLAGS} -I$(LANDSFCUTIL_INCd) -I$(G2_INCd) -I$(NEMSIO_INC) $(NETCDF_INCLUDE) -c program_setup.f

grib_utils.o: grib_utils.f
	${FCOMP} ${FFLAGS} -I$(G2_INCd) -c grib_utils.f

machine.o: machine.f
	${FCOMP} ${FFLAGS} -c machine.f

surface_chgres.o: machine.o surface_chgres.F
	${FCOMP} ${FFLAGS} ${FPPFLAGS} -I$(IP_INCd) -I$(LANDSFCUTIL_INCd) -I$(G2_INCd) -c surface_chgres.F

read_data.o: program_setup.o surface_chgres.o read_data.f
	${FCOMP} ${FFLAGS} -I$(NEMSIO_INC) -I$(LANDSFCUTIL_INCd) -I$(SFCIO_INC4) -c read_data.f

interp_data.o: program_setup.o read_data.o surface_chgres.o interp_data.f
	${FCOMP} ${FFLAGS} -I$(NEMSIO_INC) -c interp_data.f

write_data.o: interp_data.o program_setup.o write_data.f
	${FCOMP} ${FFLAGS} -I$(NEMSIO_INC) -I$(LANDSFCUTIL_INCd) $(NETCDF_INCLUDE) -c write_data.f

clean:
	rm -f *.o *.mod surface_chgres.i ${EXE} ../../exec/${EXE}
install:
	/usr/bin/install -c -d ../../exec
	/usr/bin/install -c $(EXE) ../../exec
check_prereqs:
	/nwprod/spa_util/check_libs.bash $(LIBS)
#	/nwprod/spa_util/check_incs.bash $(LANDSFCUTIL_INCd) $(G2_INCd) $(NEMSIO_INC) $(NETCDF_INCLUDE) $(SFCIO_INC4)
