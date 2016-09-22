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

      PARAMETER(IMAXIN=1500,JJMAXIN=1500,JMAXIN=IMAXIN*JJMAXIN)
      PARAMETER(IMAXOT=2500,JJMAXOT=1500,JMAXOT=IMAXOT*JJMAXOT)

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
c     DATA  GRD94 / 0, 255,203,345,569,  -3441, -148799, 136,   50000,
c    & -111000,    154,141,64, 0, 0, 0, 0, 0/

       IF (KGRIDIN.EQ.255) THEN
CHC        KGDSIN=0
CHC        KGDSIN(1)=3
CHC        KGDSIN(2)=259
CHC        KGDSIN(3)=163
CHC        KGDSIN(4)=19881
CHC        KGDSIN(5)=235911
CHC        KGDSIN(6)=8
CHC        KGDSIN(7)=262000
CHC        KGDSIN(8)=22000
CHC        KGDSIN(9)=22000
CHC        KGDSIN(10)=0
CHC        KGDSIN(11)=64
CHC        KGDSIN(12)=60000
CHC        KGDSIN(13)=30000
CHC        KGDSIN(20)=255
C
!         IUGDF=15
!         READ(IUGDF)KGDSIN(1)
!         READ(IUGDF)KGDSIN(2)
!         READ(IUGDF)KGDSIN(3)
!         READ(IUGDF)KGDSIN(4)
!         READ(IUGDF)KGDSIN(5)
!         READ(IUGDF)KGDSIN(6)
!         READ(IUGDF)KGDSIN(7)
!         READ(IUGDF)KGDSIN(8)
!         READ(IUGDF)KGDSIN(9)
!         READ(IUGDF)KGDSIN(10)
!         READ(IUGDF)KGDSIN(11)
!         READ(IUGDF)KGDSIN(12)
!         READ(IUGDF)KGDSIN(13)
!         READ(IUGDF)KGDSIN(20)
       ENDIF

!-- na12aq domain:
!        kgdsin(1)=205
!        kgdsin(2)=706       !- Number of x grid points
!        kgdsin(3)=564       !- Number of y grid points
!        kgdsin(4)=9222      !- 1000*Latitude of SW corner; see (A) below
!        kgdsin(5)=211657    !- 1000*(360+Longitude of SW corner); see (A) below 
!        kgdsin(6)=136       !- unchanged
!        kgdsin(7)=50000     !- 1000*(Latitude of center point)
!        kgdsin(8)=249000    !- 1000*(360+Longitude of center point)
!        kgdsin(9)=124       !- 1000*DX in longitude; see (B) below
!        kgdsin(10)=106      !- 1000*DY in latitude; see (B) below
!
!-- conus4 domain:
!        kgdsin(2)=1300      !- Number of x grid points
!        kgdsin(3)=875       !- Number of y grid points
!        kgdsin(4)=19984     !- 1000*Latitude of SW corner; see (A) below
!        kgdsin(5)=237514    !- 1000*(360+Longitude of SW corner); see (A) below
!        kgdsin(6)=136       !- unchanged ?
!        kgdsin(7)=39000     !- 1000*(Latitude of center point)
!        kgdsin(8)=262000    !- 1000*(360+Longitude of center point)
!        kgdsin(9)=37        !- 1000*DX in longitude; see (B) below
!        kgdsin(10)=36       !- 1000*DY in latitude; see (B) below
!
!-- cent4 domain:
!       kgdsin(2)=1050      !- Number of x grid points
!       kgdsin(3)=660       !- Number of y grid points
!       kgdsin(4)=23934     !- 1000*Latitude of SW corner; see (A) below
!       kgdsin(5)=247296    !- 1000*(360+Longitude of SW corner); see (A) below
!       kgdsin(6)=136       !- unchanged ?
!       kgdsin(7)=38000     !- 1000*(Latitude of center point)
!       kgdsin(8)=268000    !- 1000*(360+Longitude of center point)
!       kgdsin(9)=37        !- 1000*DX in longitude; see (B) below
!       kgdsin(10)=36       !- 1000*DY in latitude; see (B) below
!
!-- new expanded NAM domain:
!       kgdsin(1)=205
!       kgdsin(2)=961       !- Number of x grid points
!       kgdsin(3)=901       !- Number of y grid points
!       kgdsin(4)=-7450     !- 1000*Latitude of SW corner; see (A) below
!       kgdsin(5)=215860    !- 1000*(360+Longitude of SW corner); see (A) below
!       kgdsin(6)=136       !- unchanged ?
!       kgdsin(7)=54000     !- 1000*(Latitude of center point)
!       kgdsin(8)=254000    !- 1000*(360+Longitude of center point)
!       kgdsin(9)=125       !- 1000*DX in longitude; see (B) below
!       kgdsin(10)=100      !- 1000*DY in latitude; see (B) below
!
!-- new expanded NAM domain:
!-- new iplib now needs lat/lon of NE corner point
        kgdsin(1)=205
        kgdsin(2)=954       !- Number of x grid points
        kgdsin(3)=835       !- Number of y grid points
        kgdsin(4)=-7491     !- 1000*Latitude of SW corner; see (A) below
        kgdsin(5)=215866    !- 1000*(360+Longitude of SW corner); see (A) below 
        kgdsin(6)=136       !- unchanged ?
        kgdsin(7)=54000     !- 1000*(Latitude of center point)
        kgdsin(8)=254000    !- 1000*(360+Longitude of center point)
        kgdsin(9)=126       !- 1000*DX in longitude; see (B) below 
        kgdsin(10)=108      !- 1000*DY in latitude; see (B) below
        kgdsin(12)=44539    !- 1000*Latitude of NE corner
        kgdsin(13)=14802    !- 1000*Longitude of NE Corner
!
!-- new 4 km CONUS NEST
        kgdsin(1)=205
        kgdsin(2)=1371      !- Number of x grid points
        kgdsin(3)=1100       !- Number of y grid points
        kgdsin(4)=15947     !- 1000*Latitude of SW corner; see (A) below
        kgdsin(5)=234532    !- 1000*(360+Longitude of SW corner); see (A) below
        kgdsin(6)=136       !- unchanged ?
        kgdsin(7)=54000     !- 1000*(Latitude of center point)
        kgdsin(8)=254000    !- 1000*(360+Longitude of center point)
        kgdsin(9)=042       !- 1000*DX in longitude; see (B) below
        kgdsin(10)=036      !- 1000*DY in latitude; see (B) below
        kgdsin(12)=45407    !- 1000*Latitude of NE corner
        kgdsin(13)=307610   !- 1000*Longitude of NE Corner
!-- new 4 km Alaska NEST
        kgdsin(1)=205
        kgdsin(2)=595      !- Number of x grid points
        kgdsin(3)=625       !- Number of y grid points
        kgdsin(4)=34921     !- 1000*Latitude of SW corner; see (A) below
        kgdsin(5)=198337    !- 1000*(360+Longitude of SW corner); see (A) below
        kgdsin(6)=136       !- unchanged ?
        kgdsin(7)=54000     !- 1000*(Latitude of center point)
        kgdsin(8)=254000    !- 1000*(360+Longitude of center point)
        kgdsin(9)=063       !- 1000*DX in longitude; see (B) below
        kgdsin(10)=054      !- 1000*DY in latitude; see (B) below
        kgdsin(12)=83771    !- 1000*Latitude of NE corner
        kgdsin(13)=208279  !- 1000*Longitude of NE Corner
!-- new 3 km Hawaii NEST
        kgdsin(1)=205
        kgdsin(2)=373      !- Number of x grid points
        kgdsin(3)=561       !- Number of y grid points
        kgdsin(4)=11625     !- 1000*Latitude of SW corner; see (A) below
        kgdsin(5)=203661    !- 1000*(360+Longitude of SW corner); see (A) below
        kgdsin(6)=136       !- unchanged ?
        kgdsin(7)=54000     !- 1000*(Latitude of center point)
        kgdsin(8)=254000    !- 1000*(360+Longitude of center point)
        kgdsin(9)=032       !- 1000*DX in longitude; see (B) below
        kgdsin(10)=027      !- 1000*DY in latitude; see (B) below
        kgdsin(12)=30429    !- 1000*Latitude of NE corner
        kgdsin(13)=202173   !- 1000*Longitude of NE Corner
!-- new 3 km Puerto Rico NEST
        kgdsin(1)=205
        kgdsin(2)=241      !- Number of x grid points
        kgdsin(3)=241       !- Number of y grid points
        kgdsin(4)=17500     !- 1000*Latitude of SW corner; see (A) below
        kgdsin(5)=288950    !- 1000*(360+Longitude of SW corner); see (A) below
        kgdsin(6)=136       !- unchanged ?
        kgdsin(7)=54000     !- 1000*(Latitude of center point)
        kgdsin(8)=254000    !- 1000*(360+Longitude of center point)
        kgdsin(9)=032       !- 1000*DX in longitude; see (B) below
        kgdsin(10)=027      !- 1000*DY in latitude; see (B) below
        kgdsin(12)=18926    !- 1000*Latitude of NE corner
        kgdsin(13)=298808   !- 1000*Longitude of NE Corner
!-- enlarged 3 km Puerto Rico NEST
        kgdsin(1)=205
        kgdsin(2)=401      !- Number of x grid points
        kgdsin(3)=325       !- Number of y grid points
        kgdsin(4)=17609     !- 1000*Latitude of SW corner; see (A) below
        kgdsin(5)=283672    !- 1000*(360+Longitude of SW corner); see (A) below
        kgdsin(6)=136       !- unchanged ?
        kgdsin(7)=54000     !- 1000*(Latitude of center point)
        kgdsin(8)=254000    !- 1000*(360+Longitude of center point)
        kgdsin(9)=032       !- 1000*DX in longitude; see (B) below
        kgdsin(10)=027      !- 1000*DY in latitude; see (B) below
        kgdsin(12)=18840    !- 1000*Latitude of NE corner
        kgdsin(13)=298739   !- 1000*Longitude of NE Corner
!  Rapid Refresh grid
        kgdsin(1)=205
        kgdsin(2)=758
        kgdsin(3)=567
        kgdsin(4)=2228
        kgdsin(5)=219519
        kgdsin(6)=136
        kgdsin(7)=47500
        kgdsin(8)=256000
        kgdsin(9)=121
        kgdsin(10)=121
        kgdsin(12)=53492
        kgdsin(13)=349016
!-- new CONUS DGEX domain:
!-- new iplib now needs lat/lon of NE corner point
        kgdsin(1)=205
        kgdsin(2)=421       !- Number of x grid points
        kgdsin(3)=302       !- Number of y grid points
        kgdsin(4)=19757     !- 1000*Latitude of SW corner; see (A) below
        kgdsin(5)=234966    !- 1000*(360+Longitude of SW corner); see (A) below
        kgdsin(6)=136       !- unchanged ?
        kgdsin(7)=40000     !- 1000*(Latitude of center point)
        kgdsin(8)=262000    !- 1000*(360+Longitude of center point)
        kgdsin(9)=126       !- 1000*DX in longitude; see (B) below
        kgdsin(10)=108      !- 1000*DY in latitude; see (B) below
        kgdsin(12)=50073    !- 1000*Latitude of NE corner
        kgdsin(13)=303797   !- 1000*Longitude of NE Corner
!
!-- new Alaska DGEX domain:
!-- new iplib now needs lat/lon of NE corner point
!       kgdsin(1)=205
!       kgdsin(2)=325       !- Number of x grid points
!       kgdsin(3)=231       !- Number of y grid points
!       kgdsin(4)=44127     !- 1000*Latitude of SW corner; see (A) below
!       kgdsin(5)=174672    !- 1000*(360+Longitude of SW corner); see (A) below
!       kgdsin(6)=136       !- unchanged ?
!       kgdsin(7)=61000     !- 1000*(Latitude of center point)
!       kgdsin(8)=203000    !- 1000*(360+Longitude of center point)
!       kgdsin(9)=126       !- 1000*DX in longitude; see (B) below
!       kgdsin(10)=108      !- 1000*DY in latitude; see (B) below
!       kgdsin(12)=64795    !- 1000*Latitude of NE corner
!       kgdsin(13)=256112   !- 1000*Longitude of NE Corner
!
!-- new 4 km CONUS NEST
        kgdsin(1)=205
        kgdsin(2)=1371      !- Number of x grid points
        kgdsin(3)=1100       !- Number of y grid points
        kgdsin(4)=15947     !- 1000*Latitude of SW corner; see (A) below
        kgdsin(5)=234532    !- 1000*(360+Longitude of SW corner); see (A) below
        kgdsin(6)=136       !- unchanged ?
        kgdsin(7)=54000     !- 1000*(Latitude of center point)
        kgdsin(8)=254000    !- 1000*(360+Longitude of center point)
        kgdsin(9)=042       !- 1000*DX in longitude; see (B) below
        kgdsin(10)=036      !- 1000*DY in latitude; see (B) below
        kgdsin(12)=45407    !- 1000*Latitude of NE corner
        kgdsin(13)=307610   !- 1000*Longitude of NE Corner
!
!-- new 4 km Alaska NEST
        kgdsin(1)=205
        kgdsin(2)=595      !- Number of x grid points
        kgdsin(3)=625       !- Number of y grid points
        kgdsin(4)=34921     !- 1000*Latitude of SW corner; see (A) below
        kgdsin(5)=198337    !- 1000*(360+Longitude of SW corner); see (A) below
        kgdsin(6)=136       !- unchanged ?
        kgdsin(7)=54000     !- 1000*(Latitude of center point)
        kgdsin(8)=254000    !- 1000*(360+Longitude of center point)
        kgdsin(9)=063       !- 1000*DX in longitude; see (B) below
        kgdsin(10)=054      !- 1000*DY in latitude; see (B) below
        kgdsin(12)=83771    !- 1000*Latitude of NE corner
        kgdsin(13)=208279  !- 1000*Longitude of NE Corner
!
! SREF 16km NMMB domain
!
        kgdsin(1)=205
        kgdsin(2)=716       !- Number of x grid points
        kgdsin(3)=627       !- Number of y grid points
        kgdsin(4)=-7527     !- 1000*Latitude of SW corner; see (A) below
        kgdsin(5)=215881    !- 1000*(360+Longitude of SW corner); see (A) below
        kgdsin(6)=136       !- unchanged 
        kgdsin(7)=54000     !- 1000*(Latitude of center point)
        kgdsin(8)=254000    !- 1000*(360+Longitude of center point)
        kgdsin(9)=168       !- 1000*DX in longitude; see (B) below
        kgdsin(10)=144      !- 1000*DY in latitude; see (B) below
        kgdsin(12)=44532    !- 1000*Latitude of NE corner
        kgdsin(13)=14855    !- 1000*Longitude of NE Corner
!
! SREF 20km NMMB domain
!
!       kgdsin(1)=205
!       kgdsin(2)=572       !- Number of x grid points
!       kgdsin(3)=501       !- Number of y grid points
!       kgdsin(4)=-7424     !- 1000*Latitude of SW corner; see (A) below
!       kgdsin(5)=215883    !- 1000*(360+Longitude of SW corner); see (A) below
!       kgdsin(6)=136       !- unchanged
!       kgdsin(7)=54000     !- 1000*(Latitude of center point)
!       kgdsin(8)=254000    !- 1000*(360+Longitude of center point)
!       kgdsin(9)=210       !- 1000*DX in longitude; see (B) below
!       kgdsin(10)=180      !- 1000*DY in latitude; see (B) below
!       kgdsin(12)=44591    !- 1000*Latitude of NE corner
!       kgdsin(13)=14737    !- 1000*Longitude of NE Corner
!
! SREF 25km NMMB domain
!
!       kgdsin(1)=205
!       kgdsin(2)=458       !- Number of x grid points
!       kgdsin(3)=401       !- Number of y grid points
!       kgdsin(4)=-7437     !- 1000*Latitude of SW corner; see (A) below
!       kgdsin(5)=215870    !- 1000*(360+Longitude of SW corner); see (A) below
!       kgdsin(6)=136       !- unchanged 
!       kgdsin(7)=54000     !- 1000*(Latitude of center point)
!       kgdsin(8)=254000    !- 1000*(360+Longitude of center point)
!       kgdsin(9)=262       !- 1000*DX in longitude; see (B) below
!       kgdsin(10)=225      !- 1000*DY in latitude; see (B) below
!       kgdsin(12)=44573    !- 1000*Latitude of NE corner
!       kgdsin(13)=14741    !- 1000*Longitude of NE Corner

!-- Get the following information from the NPS directory:
!
!-- (A) Get the lat/lon for the SW corner by searching for the string 
!       "SW corner" from geogrid.log.
!
!-- (B) Get the dx (delta longitude) and dy (delta latitude) from namelist.nps.
!
        kgdsin(11)=64
        kgdsin(20)=255

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
      STOP1
      END


