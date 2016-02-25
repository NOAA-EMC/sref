C-----------------------------------------------------------------------
      SUBROUTINE GETGBEXM(LUGB,LUGI,JF,J,JPDS,JGDS,JENS,
     &                    MBUF,CBUF,NLEN,NNUM,MNUM,
     &                    KF,K,KPDS,KGDS,KENS,KPROB,XPROB,KCLUST,KMEMBR,
     &                    LB,F,IRET)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM: GETGBEXM       FINDS AND UNPACKS A GRIB MESSAGE
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 94-04-01
C
C ABSTRACT: FIND AND UNPACK A GRIB MESSAGE.
C   READ A GRIB INDEX FILE (OR OPTIONALLY THE GRIB FILE ITSELF)
C   TO GET THE INDEX BUFFER (I.E. TABLE OF CONTENTS) FOR THE GRIB FILE.
C   FIND IN THE INDEX BUFFER A REFERENCE TO THE GRIB MESSAGE REQUESTED.
C   THE GRIB MESSAGE REQUEST SPECIFIES THE NUMBER OF MESSAGES TO SKIP
C   AND THE UNPACKED PDS AND GDS PARAMETERS.  (A REQUESTED PARAMETER
C   OF -1 MEANS TO ALLOW ANY VALUE OF THIS PARAMETER TO BE FOUND.)
C   IF THE REQUESTED GRIB MESSAGE IS FOUND, THEN IT IS READ FROM THE
C   GRIB FILE AND UNPACKED.  ITS MESSAGE NUMBER IS RETURNED ALONG WITH
C   THE UNPACKED PDS AND GDS PARAMETERS, THE UNPACKED BITMAP (IF ANY),
C   AND THE UNPACKED DATA.  IF THE GRIB MESSAGE IS NOT FOUND, THEN THE
C   RETURN CODE WILL BE NONZERO.
C
C PROGRAM HISTORY LOG:
C   94-04-01  IREDELL
C   95-10-31  IREDELL     MODULARIZED PORTIONS OF CODE INTO SUBPROGRAMS
C                         AND ALLOWED FOR UNSPECIFIED INDEX FILE
C   97-02-11  Y.ZHU       INCLUDED PROBABILITY AND CLUSTER ARGUMENTS
C
C USAGE:    CALL GETGBEXM(LUGB,LUGI,JF,J,JPDS,JGDS,JENS,
C    &                    MBUF,CBUF,NLEN,NNUM,MNUM,
C    &                    KF,K,KPDS,KGDS,KENS,KPROB,XPROB,KCLUST,KMEMBR,
C    &                    LB,F,IRET)
C   INPUT ARGUMENTS:
C     LUGB         INTEGER UNIT OF THE UNBLOCKED GRIB DATA FILE
C     LUGI         INTEGER UNIT OF THE UNBLOCKED GRIB INDEX FILE
C                  (=0 TO GET INDEX BUFFER FROM THE GRIB FILE)
C     JF           INTEGER MAXIMUM NUMBER OF DATA POINTS TO UNPACK
C     J            INTEGER NUMBER OF MESSAGES TO SKIP
C                  (=0 TO SEARCH FROM BEGINNING)
C                  (<0 TO READ INDEX BUFFER AND SKIP -1-J MESSAGES)
C     JPDS         INTEGER (200) PDS PARAMETERS FOR WHICH TO SEARCH
C                  (=-1 FOR WILDCARD)
C          (1)   - ID OF CENTER
C          (2)   - GENERATING PROCESS ID NUMBER
C          (3)   - GRID DEFINITION
C          (4)   - GDS/BMS FLAG (RIGHT ADJ COPY OF OCTET 8)
C          (5)   - INDICATOR OF PARAMETER
C          (6)   - TYPE OF LEVEL
C          (7)   - HEIGHT/PRESSURE , ETC OF LEVEL
C          (8)   - YEAR INCLUDING (CENTURY-1)
C          (9)   - MONTH OF YEAR
C          (10)  - DAY OF MONTH
C          (11)  - HOUR OF DAY
C          (12)  - MINUTE OF HOUR
C          (13)  - INDICATOR OF FORECAST TIME UNIT
C          (14)  - TIME RANGE 1
C          (15)  - TIME RANGE 2
C          (16)  - TIME RANGE FLAG
C          (17)  - NUMBER INCLUDED IN AVERAGE
C          (18)  - VERSION NR OF GRIB SPECIFICATION
C          (19)  - VERSION NR OF PARAMETER TABLE
C          (20)  - NR MISSING FROM AVERAGE/ACCUMULATION
C          (21)  - CENTURY OF REFERENCE TIME OF DATA
C          (22)  - UNITS DECIMAL SCALE FACTOR
C          (23)  - SUBCENTER NUMBER
C          (24)  - PDS BYTE 29, FOR NMC ENSEMBLE PRODUCTS
C                  128 IF FORECAST FIELD ERROR
C                   64 IF BIAS CORRECTED FCST FIELD
C                   32 IF SMOOTHED FIELD
C                  WARNING: CAN BE COMBINATION OF MORE THAN 1
C          (25)  - PDS BYTE 30, NOT USED
C     JGDS         INTEGER (200) GDS PARAMETERS FOR WHICH TO SEARCH
C                  (ONLY SEARCHED IF JPDS(3)=255)
C                  (=-1 FOR WILDCARD)
C          (1)   - DATA REPRESENTATION TYPE
C          (19)  - NUMBER OF VERTICAL COORDINATE PARAMETERS
C          (20)  - OCTET NUMBER OF THE LIST OF VERTICAL COORDINATE
C                  PARAMETERS
C                  OR
C                  OCTET NUMBER OF THE LIST OF NUMBERS OF POINTS
C                  IN EACH ROW
C                  OR
C                  255 IF NEITHER ARE PRESENT
C          (21)  - FOR GRIDS WITH PL, NUMBER OF POINTS IN GRID
C          (22)  - NUMBER OF WORDS IN EACH ROW
C       LATITUDE/LONGITUDE GRIDS
C          (2)   - N(I) NR POINTS ON LATITUDE CIRCLE
C          (3)   - N(J) NR POINTS ON LONGITUDE MERIDIAN
C          (4)   - LA(1) LATITUDE OF ORIGIN
C          (5)   - LO(1) LONGITUDE OF ORIGIN
C          (6)   - RESOLUTION FLAG (RIGHT ADJ COPY OF OCTET 17)
C          (7)   - LA(2) LATITUDE OF EXTREME POINT
C          (8)   - LO(2) LONGITUDE OF EXTREME POINT
C          (9)   - DI LONGITUDINAL DIRECTION OF INCREMENT
C          (10)  - DJ LATITUDINAL DIRECTION INCREMENT
C          (11)  - SCANNING MODE FLAG (RIGHT ADJ COPY OF OCTET 28)
C       GAUSSIAN  GRIDS
C          (2)   - N(I) NR POINTS ON LATITUDE CIRCLE
C          (3)   - N(J) NR POINTS ON LONGITUDE MERIDIAN
C          (4)   - LA(1) LATITUDE OF ORIGIN
C          (5)   - LO(1) LONGITUDE OF ORIGIN
C          (6)   - RESOLUTION FLAG  (RIGHT ADJ COPY OF OCTET 17)
C          (7)   - LA(2) LATITUDE OF EXTREME POINT
C          (8)   - LO(2) LONGITUDE OF EXTREME POINT
C          (9)   - DI LONGITUDINAL DIRECTION OF INCREMENT
C          (10)  - N - NR OF CIRCLES POLE TO EQUATOR
C          (11)  - SCANNING MODE FLAG (RIGHT ADJ COPY OF OCTET 28)
C          (12)  - NV - NR OF VERT COORD PARAMETERS
C          (13)  - PV - OCTET NR OF LIST OF VERT COORD PARAMETERS
C                             OR
C                  PL - LOCATION OF THE LIST OF NUMBERS OF POINTS IN
C                       EACH ROW (IF NO VERT COORD PARAMETERS
C                       ARE PRESENT
C                             OR
C                  255 IF NEITHER ARE PRESENT
C       POLAR STEREOGRAPHIC GRIDS
C          (2)   - N(I) NR POINTS ALONG LAT CIRCLE
C          (3)   - N(J) NR POINTS ALONG LON CIRCLE
C          (4)   - LA(1) LATITUDE OF ORIGIN
C          (5)   - LO(1) LONGITUDE OF ORIGIN
C          (6)   - RESOLUTION FLAG  (RIGHT ADJ COPY OF OCTET 17)
C          (7)   - LOV GRID ORIENTATION
C          (8)   - DX - X DIRECTION INCREMENT
C          (9)   - DY - Y DIRECTION INCREMENT
C          (10)  - PROJECTION CENTER FLAG
C          (11)  - SCANNING MODE (RIGHT ADJ COPY OF OCTET 28)
C       SPHERICAL HARMONIC COEFFICIENTS
C          (2)   - J PENTAGONAL RESOLUTION PARAMETER
C          (3)   - K      "          "         "
C          (4)   - M      "          "         "
C          (5)   - REPRESENTATION TYPE
C          (6)   - COEFFICIENT STORAGE MODE
C       MERCATOR GRIDS
C          (2)   - N(I) NR POINTS ON LATITUDE CIRCLE
C          (3)   - N(J) NR POINTS ON LONGITUDE MERIDIAN
C          (4)   - LA(1) LATITUDE OF ORIGIN
C          (5)   - LO(1) LONGITUDE OF ORIGIN
C          (6)   - RESOLUTION FLAG (RIGHT ADJ COPY OF OCTET 17)
C          (7)   - LA(2) LATITUDE OF LAST GRID POINT
C          (8)   - LO(2) LONGITUDE OF LAST GRID POINT
C          (9)   - LATIT - LATITUDE OF PROJECTION INTERSECTION
C          (10)  - RESERVED
C          (11)  - SCANNING MODE FLAG (RIGHT ADJ COPY OF OCTET 28)
C          (12)  - LONGITUDINAL DIR GRID LENGTH
C          (13)  - LATITUDINAL DIR GRID LENGTH
C       LAMBERT CONFORMAL GRIDS
C          (2)   - NX NR POINTS ALONG X-AXIS
C          (3)   - NY NR POINTS ALONG Y-AXIS
C          (4)   - LA1 LAT OF ORIGIN (LOWER LEFT)
C          (5)   - LO1 LON OF ORIGIN (LOWER LEFT)
C          (6)   - RESOLUTION (RIGHT ADJ COPY OF OCTET 17)
C          (7)   - LOV - ORIENTATION OF GRID
C          (8)   - DX - X-DIR INCREMENT
C          (9)   - DY - Y-DIR INCREMENT
C          (10)  - PROJECTION CENTER FLAG
C          (11)  - SCANNING MODE FLAG (RIGHT ADJ COPY OF OCTET 28)
C          (12)  - LATIN 1 - FIRST LAT FROM POLE OF SECANT CONE INTER
C          (13)  - LATIN 2 - SECOND LAT FROM POLE OF SECANT CONE INTER
C     JENS         INTEGER (200) ENSEMBLE PDS PARMS FOR WHICH TO SEARCH
C                  (ONLY SEARCHED IF JPDS(23)=2)
C                  (=-1 FOR WILDCARD)
C          (1)   - APPLICATION IDENTIFIER
C          (2)   - ENSEMBLE TYPE
C          (3)   - ENSEMBLE IDENTIFIER
C          (4)   - PRODUCT IDENTIFIER
C          (5)   - SMOOTHING FLAG
C     MBUF         INTEGER LENGTH OF INDEX BUFFER IN BYTES
C     CBUF         CHARACTER*1 (MBUF) INDEX BUFFER
C                  (INITIALIZE BY SETTING J=-1)
C     NLEN         INTEGER LENGTH OF EACH INDEX RECORD IN BYTES
C                  (INITIALIZE BY SETTING J=-1)
C     NNUM         INTEGER NUMBER OF INDEX RECORDS
C                  (INITIALIZE BY SETTING J=-1)
C     MNUM         INTEGER NUMBER OF INDEX RECORDS SKIPPED
C                  (INITIALIZE BY SETTING J=-1)
C   OUTPUT ARGUMENTS:
C     CBUF         CHARACTER*1 (MBUF) INDEX BUFFER
C     NLEN         INTEGER LENGTH OF EACH INDEX RECORD IN BYTES
C     NNUM         INTEGER NUMBER OF INDEX RECORDS
C     MNUM         INTEGER NUMBER OF INDEX RECORDS SKIPPED
C     KF           INTEGER NUMBER OF DATA POINTS UNPACKED
C     K            INTEGER MESSAGE NUMBER UNPACKED
C                  (CAN BE SAME AS J IN CALLING PROGRAM
C                  IN ORDER TO FACILITATE MULTIPLE SEARCHES)
C     KPDS         INTEGER (200) UNPACKED PDS PARAMETERS
C     KGDS         INTEGER (200) UNPACKED GDS PARAMETERS
C     KENS         INTEGER (200) UNPACKED ENSEMBLE PDS PARMS
C     KPROB        INTEGER (2) PROBABILITY ENSEMBLE PARMS
C     XPROB        REAL    (2) PROBABILITY ENSEMBLE PARMS
C     KCLUST       INTEGER (16) CLUSTER ENSEMBLE PARMS
C     KMEMBR       INTEGER (8) CLUSTER ENSEMBLE PARMS
C     LB           LOGICAL*1 (KF) UNPACKED BITMAP IF PRESENT
C     F            REAL (KF) UNPACKED DATA
C     IRET         INTEGER RETURN CODE
C                    0      ALL OK
C                    96     ERROR READING INDEX FILE
C                    97     ERROR READING GRIB FILE
C                    98     NUMBER OF DATA POINTS GREATER THAN JF
C                    99     REQUEST NOT FOUND
C                    OTHER  W3FI63 GRIB UNPACKER RETURN CODE
C
C SUBPROGRAMS CALLED:
C   GETGI          READ INDEX FILE
C   GETGIR         READ INDEX BUFFER FROM GRIB FILE
C   GETGB1S        SEARCH INDEX RECORDS
C   GETGB1RE       READ AND UNPACK GRIB RECORD
C   LENGDS         RETURN THE LENGTH OF A GRID
C
C REMARKS: SPECIFY AN INDEX FILE IF FEASIBLE TO INCREASE SPEED.
C   SUBPROGRAM CAN BE CALLED FROM A MULTIPROCESSING ENVIRONMENT.
C   DO NOT ENGAGE THE SAME LOGICAL UNIT FROM MORE THAN ONE PROCESSOR.
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 77
C   MACHINE:  CRAY, WORKSTATIONS
C
C$$$
      INTEGER JPDS(200),JGDS(200),JENS(200)
      INTEGER KPDS(200),KGDS(200),KENS(200)
      INTEGER KPROB(2),KCLUST(16),KMEMBR(80)
      REAL XPROB(2)
      CHARACTER CBUF(MBUF)
      LOGICAL*1 LB(JF)
      REAL F(JF)
      PARAMETER(MSK1=32000,MSK2=4000)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  SEARCH PREVIOUS INDEX BUFFER IF POSSIBLE
      IF(J.GE.0) THEN
        IF(MNUM.GE.0) THEN
          IRGI=0
        ELSE
          MNUM=-1-MNUM
          IRGI=1
        ENDIF
        JR=J-MNUM
        IF(JR.GE.0.AND.(JR.LT.NNUM.OR.IRGI.EQ.0)) THEN
          CALL GETGB1S(CBUF,NLEN,NNUM,JR,JPDS,JGDS,JENS,
     &                 KR,KPDS,KGDS,KENS,LSKIP,LGRIB,IRGS)
          IF(IRGS.EQ.0) K=KR+MNUM
          IF(IRGI.EQ.1.AND.IRGS.EQ.0) MNUM=-1-MNUM
          IF(IRGI.EQ.1.AND.IRGS.GT.0) MNUM=MNUM+NNUM
        ELSE
          MNUM=J
          IRGI=1
          IRGS=1
        ENDIF
      ELSE
        MNUM=-1-J
        IRGI=1
        IRGS=1
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  READ AND SEARCH NEXT INDEX BUFFER
      JR=0
      DOWHILE(IRGI.EQ.1.AND.IRGS.EQ.1)
        IF(LUGI.GT.0) THEN
          CALL GETGI(LUGI,MNUM,MBUF,CBUF,NLEN,NNUM,IRGI)
        ELSE
          CALL GETGIR(LUGB,MSK1,MSK2,MNUM,MBUF,CBUF,NLEN,NNUM,IRGI)
        ENDIF
        IF(IRGI.LE.1) THEN
          CALL GETGB1S(CBUF,NLEN,NNUM,JR,JPDS,JGDS,JENS,
     &                 KR,KPDS,KGDS,KENS,LSKIP,LGRIB,IRGS)
          IF(IRGS.EQ.0) K=KR+MNUM
          IF(IRGI.EQ.1.AND.IRGS.EQ.0) MNUM=-1-MNUM
          IF(IRGI.EQ.1.AND.IRGS.GT.0) MNUM=MNUM+NNUM
        ENDIF
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  READ AND UNPACK GRIB RECORD
      IF(IRGI.GT.1) THEN
        IRET=96
      ELSEIF(IRGS.NE.0) THEN
        IRET=99
      ELSEIF(LENGDS(KGDS).GT.JF) THEN
        IRET=98
      ELSE
        CALL GETGB1RE(LUGB,LSKIP,LGRIB,KF,KPDS,KGDS,KENS,
     &                KPROB,XPROB,KCLUST,KMEMBR,LB,F,IRET)
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END
