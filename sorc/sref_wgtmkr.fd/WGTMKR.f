      PROGRAM WGTMKR
C$$$  MAIN PROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C MAIN PROGRAM: WGTMKR       MAKE INTERP WGTS FOR PRODUCT GENERATOR 
C   PRGMMR: BALDWIN/BRILL    ORG: NP22        DATE: 97-12-01  
C
C ABSTRACT: AWPRGN PRODUCES FILES THAT HAVE BEEN INTERPOLATED (USING
C   IPLIB) TO VARIOUS OUTPUT GRIDS WITH OPTIONAL WMO HEADERS.  AWPRGN
C   READS THROUGH A MASTER INPUT GRIB FILE, DETERMINES WHAT GRIDS TO
C   INTERPOLATE TO, PERFORMS PRE- AND POST-INTERPOLATION SMOOTHING,
C   PACKS THE DATA INTO GRIB, ADDS A WMO HEADER, THEN WRITES THE
C   PACKED DATA TO AN OUTPUT FILE.  AN INPUT CONTROL FILE DETERMINES 
C   THE OUTPUT GRID NUMBER, WMO HEADER TYPE, OUTPUT FILE NAME, PACKING 
C   PRECISION, AND NUMBER OF PRE- AND POST-INTERPOLATION SMOOTHING 
C   PASSES FOR EACH GRIB FIELD THAT IS DESIRED FOR POSTING.  THE
C   MASTER INPUT GRID SHOULD BE LARGE ENOUGH TO ENCOMPASS ALL OF THE
C   REQUESTED OUTPUT GRIDS, AND ALSO SHOULD CONTAIN ALL OF THE
C   GRIB PARAMETERS REQUIRED, SINCE AWPRGN DOES NOT PROVIDE
C   AND DIAGNOSTIC CALCULATIONS.
C
C PROGRAM HISTORY LOG:
C   97-12-01  BALDWIN,    ORIGINATOR
C             BRILL
C
C USAGE:  MAIN PROGRAM
C
C   INPUT FILES:
C
C   OUTPUT FILES:
C
C           UNIT 51 - INTERPOLATION WEIGHTS
C
C   SUBPROGRAMS CALLED:
C     UNIQUE:
C       POLATWGT - COMPUTE INTERPOLATION WEIGHTS
C     LIBRARY:
C       W3LIB: W3FI72
C       IPLIB: MAKGDS, IPXETAS
C       SPLIB: (FOR INTERPOLATION)
C
C   EXIT STATES:
C     COND =   1 - NORMAL EXIT
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN
C   MACHINE : CRAY J-916
C
C$$$
C
C   INTERP WGTS GENERATOR
C
C
C   IMAXIN, JJMAXIN ARE MAX DIMENSIONS OF INPUT GRID
C   IMAXOT, JJMAXOT ARE MAX DIMENSIONS OF OUTPUT GRID
C    

      PARAMETER(IMAXIN=1000,JJMAXIN=1201,JMAXIN=IMAXIN*JJMAXIN)
      PARAMETER(IMAXOT=1480,JJMAXOT=1030,JMAXOT=IMAXOT*JJMAXOT)

      INTEGER KGDSIN(200)
      INTEGER KGDSIN2(200),KGDSOUT(200)

      INTEGER N11(JMAXOT),N21(JMAXOT),
     &              N12(JMAXOT),N22(JMAXOT),
     &              NPP(JMAXOT,25)
      INTEGER NV11(JMAXOT),NV21(JMAXOT),
     &              NV12(JMAXOT),NV22(JMAXOT)

      REAL RLAT(JMAXOT),RLON(JMAXOT)
      REAL CROT(JMAXOT),SROT(JMAXOT)
      REAL W11(JMAXOT),W21(JMAXOT),
     &     W12(JMAXOT),W22(JMAXOT)
      REAL WV11(JMAXOT),WV21(JMAXOT),
     &     WV12(JMAXOT),WV22(JMAXOT)
      REAL C11(JMAXOT),C21(JMAXOT),
     &     C12(JMAXOT),C22(JMAXOT)
      REAL S11(JMAXOT),S21(JMAXOT),
     &     S12(JMAXOT),S22(JMAXOT)

      CHARACTER GDSO(400)

C
      CALL W3TAGB('INTWGTGN',0097,0355,0068,'NP22   ')

       LUNOUT=51
C
C     READ INPUT/OUTPUT GRID NUMBERS
C     
      READ(5,*) KGRIDIN
      READ(5,*) KGRIDOT
C
C  MAKE GDS FOR INPUT/OUTPUT GRIDS
C
       IF (KGRIDIN.EQ.255) THEN
C
         IUGDF=15
         READ(IUGDF)KGDSIN(1)
         READ(IUGDF)KGDSIN(2)
         READ(IUGDF)KGDSIN(3)
         READ(IUGDF)KGDSIN(4)
         READ(IUGDF)KGDSIN(5)
         READ(IUGDF)KGDSIN(6)
         READ(IUGDF)KGDSIN(7)
         READ(IUGDF)KGDSIN(8)
         READ(IUGDF)KGDSIN(9)
         READ(IUGDF)KGDSIN(10)
         READ(IUGDF)KGDSIN(11)
         READ(IUGDF)KGDSIN(12)
         READ(IUGDF)KGDSIN(13)
         READ(IUGDF)KGDSIN(20)
!cdusan         KGDSIN(8) = KGDSIN(8) + 360000
       ENDIF

       DO ii=1,20
        print*,'KGDSIN(',ii,')= ',KGDSIN(ii)
       end do
       CALL MAKGDS(KGRIDIN,KGDSIN,GDSO,LENGDS,IRET)
        write(6,*) ' kgdsin ',kgdsin
	write(6,*) 'KGRIDOT: ', KGRIDOT
       CALL MAKGDS(KGRIDOT,KGDSOUT,GDSO,LENGDS,IRET)
        write(6,*) ' kgdsout ',kgdsout
C
C  COMPUTE INTERPOLATION WEIGHTS
C
         NOUT=KGDSOUT(2)*KGDSOUT(3)
         NO=NOUT
         CALL POLATWGT(KGDSIN,KGDSOUT,
     &   JMAXIN,JMAXOT,KGRIDOT,
     &   N11,N21,N12,N22,NPP,
     &   NV11,NV21,NV12,NV22,
     &   W11,W21,W12,W22,
     &   WV11,WV21,WV12,WV22,
     &   C11,C21,C12,C22,
     &   S11,S21,S12,S22,
     &   NO,RLAT,RLON,CROT,SROT,IRET2)

C
         WRITE(LUNOUT) KGRIDOT,NOUT
         WRITE(LUNOUT) (KGDSOUT(I),I=1,22)
         WRITE(LUNOUT) (N11(I),I=1,NOUT)
         WRITE(LUNOUT) (N12(I),I=1,NOUT)
         WRITE(LUNOUT) (N21(I),I=1,NOUT)
         WRITE(LUNOUT) (N22(I),I=1,NOUT)
         WRITE(LUNOUT) (NV11(I),I=1,NOUT)
         WRITE(LUNOUT) (NV12(I),I=1,NOUT)
         WRITE(LUNOUT) (NV21(I),I=1,NOUT)
         WRITE(LUNOUT) (NV22(I),I=1,NOUT)
         WRITE(LUNOUT) (C11(I),I=1,NOUT)
         WRITE(LUNOUT) (C12(I),I=1,NOUT)
         WRITE(LUNOUT) (C21(I),I=1,NOUT)
         WRITE(LUNOUT) (C22(I),I=1,NOUT)
         WRITE(LUNOUT) (S11(I),I=1,NOUT)
         WRITE(LUNOUT) (S12(I),I=1,NOUT)
         WRITE(LUNOUT) (S21(I),I=1,NOUT)
         WRITE(LUNOUT) (S22(I),I=1,NOUT)
         WRITE(LUNOUT) (W11(I),I=1,NOUT)
         WRITE(LUNOUT) (W12(I),I=1,NOUT)
         WRITE(LUNOUT) (W21(I),I=1,NOUT)
         WRITE(LUNOUT) (W22(I),I=1,NOUT)
         WRITE(LUNOUT) (WV11(I),I=1,NOUT)
         WRITE(LUNOUT) (WV12(I),I=1,NOUT)
         WRITE(LUNOUT) (WV21(I),I=1,NOUT)
         WRITE(LUNOUT) (WV22(I),I=1,NOUT)
         WRITE(LUNOUT) (RLAT(I),I=1,NOUT)
         WRITE(LUNOUT) (RLON(I),I=1,NOUT)
         WRITE(LUNOUT) (SROT(I),I=1,NOUT)
         WRITE(LUNOUT) (CROT(I),I=1,NOUT)
         WRITE(LUNOUT) ((NPP(I,J),I=1,NOUT),J=1,25)
C
      CALL W3TAGE('INTWGTGN')
      STOP
      END


