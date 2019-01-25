#
# Define the name of the executable
#
TARGET  = sref_bufr 
#
# CPP, Compiler, and Linker Options
#
FC       = $(COMP_MPI)
CPP      = /lib/cpp -P
ARCH     =
CPPFLAGS = 
#OPTS     = -O3 -assume byterecl -convert big_endian
OPTS     = -O0 -convert big_endian
LIST     = 
FREE     = -FR
#TRAPS    = -g -traceback -debug all -check all -O0 -ftrapuv -fpe0 -fp-stack-check -fstack-protector
TRAPS    = -g -traceback
PROFILE  = 
DEBUG    = 
MEM      =
MAP      = 
ESSL     =
MASS     =

#W3FLAGS = -I/nwprod/lib

NEMSIOFLAG = -I$(NEMSIO_INC)
NETCDF_LIB = ${NETCDF}/lib/libnetcdf.a

WRFFLAGS = -I$(WRFPATH)/frame
WRFLIB   = $(WRFPATH)/frame/pack_utils.o  \
           $(WRFPATH)/frame/module_internal_header_util.o \
           $(WRFPATH)/main/libwrflib.a \
           $(WRFPATH)/external/io_grib1/libio_grib1.a \
           $(WRFPATH)/external/io_grib_share/libio_grib_share.a \
           $(WRFPATH)/external/io_int/libwrfio_int.a \
           $(WRFPATH)/external/io_netcdf/libwrfio_nf.a \
           $(WRFPATH)/external/esmf_time_f90/libesmf_time.a \
           $(WRFPATH)/external/RSL_LITE/librsl_lite.a

SEARCH   =
#
## Assemble Options
#
FFLAGS   = $(OPTS) $(LIST) $(TRAPS) $(PROFILE) $(DEBUG) $(WRFFLAGS) $(W3FLAGS) $(NEMSIOFLAG)
FFLAGST  = $(OPTS) $(LIST) $(FREE) $(TRAPS) $(PROFILE) $(DEBUG) $(WRFFLAGS) $(W3FLAGS) $(NEMSIOFLAG)
LDFLAGS  = $(MEM) $(MAP) $(SMP) $(PROFILE)
LIBS     = $(ESSL) $(MASS) $(SEARCH) $(WRFLIB) \
           ${NEMSIO_LIB} ${W3EMC_LIB4} ${W3NCO_LIB4} \
           ${BACIO_LIB4} $(NETCDF_LDFLAGS_F) $(NETCDF_LDFLAGS_C)
  
#
# free-format object files
#
OBJST=	getIVariableB.o getIVariable.o getVariableB.o getVariable.o \
        kinds_mod.o \
        count_recs_wrf_binary_file.o inventory_wrf_binary_file.o \
        next_buf.o retrieve_index.o GETNEMSANDPLACE.o
#
# fixed-format object files
#
OBJS=	VRBLS2D_mod.o VRBLS3D_mod.o MASKS_mod.o SOIL_mod.o \
 WRFBUFR.o PROF_EM.o PROF_NMM.o PROF_NMM_NET.o PROF_EM_NET.o PROF_NMMB.o \
iw3jdn.o w3fs26.o w3movdat.o w3reddat.o gdswiz.o gdswiz01.o \
gdswiz03.o gdswiz05.o
#
# Includes
#
INCLUDES= parmsoil
#
# Common Blocks
#
COMMS=

DEPS= $(COMMS) $(INCLUDES)

.SUFFIXES:	.F .f .o

.F.f:
	$(CPP) $(CPPFLAGS) $< > $*.f

$(TARGET):	$(OBJST) $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJST) $(OBJS) $(LIBS) 

$(OBJS):	$(DEPS)
	$(FC) $(FFLAGS) -c $(@:.o=.f)

$(OBJST):	$(DEPS)
	$(FC) $(FFLAGST) -c $(@:.o=.f)

clean:	
	/bin/rm -f  $(TARGET) *.lst *.o *.mod
#
