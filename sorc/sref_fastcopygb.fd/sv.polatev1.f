C-----------------------------------------------------------------------
      SUBROUTINE POLATEV1(IPOPT,KGDSI,KGDSO,MI,MO,KM,IBI,LI,UI,VI,
     &                    NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:  POLATEV1   INTERPOLATE VECTOR FIELDS (BICUBIC)
C   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-04-10
C
C ABSTRACT: THIS SUBPROGRAM PERFORMS BICUBIC INTERPOLATION
C           FROM ANY GRID TO ANY GRID FOR SCALAR FIELDS.
C           BITMAPS ARE NOW ALLOWED EVEN WHEN INVALID POINTS ARE WITHIN
C           THE BICUBIC TEMPLATE PROVIDED THE MINIMUM WEIGHT IS REACHED. 
C           OPTIONS ALLOW CHOICES BETWEEN STRAIGHT BICUBIC (IPOPT(1)=0)
C           AND CONSTRAINED BICUBIC (IPOPT(1)=1) WHERE THE VALUE IS
C           CONFINED WITHIN THE RANGE OF THE SURROUNDING 16 POINTS.
C           ANOTHER OPTION IS THE MINIMUM PERCENTAGE FOR MASK,
C           I.E. PERCENT VALID INPUT DATA REQUIRED TO MAKE OUTPUT DATA,
C           (IPOPT(2)) WHICH DEFAULTS TO 50 (IF IPOPT(2)=-1).
C           BILINEAR USED WITHIN ONE GRID LENGTH OF BOUNDARIES.
C           ONLY HORIZONTAL INTERPOLATION IS PERFORMED.
C           THE GRIDS ARE DEFINED BY THEIR GRID DESCRIPTION SECTIONS
C           (PASSED IN INTEGER FORM AS DECODED BY SUBPROGRAM W3FI63).
C           THE CURRENT CODE RECOGNIZES THE FOLLOWING PROJECTIONS:
C             (KGDS(1)=000) EQUIDISTANT CYLINDRICAL
C             (KGDS(1)=001) MERCATOR CYLINDRICAL
C             (KGDS(1)=003) LAMBERT CONFORMAL CONICAL
C             (KGDS(1)=004) GAUSSIAN CYLINDRICAL (SPECTRAL NATIVE)
C             (KGDS(1)=005) POLAR STEREOGRAPHIC AZIMUTHAL
C             (KGDS(1)=202) ROTATED EQUIDISTANT CYLINDRICAL (ETA NATIVE)
C           WHERE KGDS COULD BE EITHER INPUT KGDSI OR OUTPUT KGDSO.
C           THE INPUT AND OUTPUT VECTORS ARE ROTATED SO THAT THEY ARE
C           EITHER RESOLVED RELATIVE TO THE DEFINED GRID
C           IN THE DIRECTION OF INCREASING X AND Y COORDINATES
C           OR RESOLVED RELATIVE TO EASTERLY AND NORTHERLY DIRECTIONS,
C           AS DESIGNATED BY THEIR RESPECTIVE GRID DESCRIPTION SECTIONS.
C           AS AN ADDED BONUS THE NUMBER OF OUTPUT GRID POINTS
C           AND THEIR LATITUDES AND LONGITUDES ARE ALSO RETURNED
C           ALONG WITH THEIR VECTOR ROTATION PARAMETERS.
C           ON THE OTHER HAND, THE OUTPUT CAN BE A SET OF STATION POINTS
C           IF KGDSO(1)<0, IN WHICH CASE THE NUMBER OF POINTS
C           AND THEIR LATITUDES AND LONGITUDES MUST BE INPUT 
C           ALONG WITH THEIR VECTOR ROTATION PARAMETERS.
C           OUTPUT BITMAPS WILL ONLY BE CREATED WHEN THE OUTPUT GRID
C           EXTENDS OUTSIDE OF THE DOMAIN OF THE INPUT GRID.
C           THE OUTPUT FIELD IS SET TO 0 WHERE THE OUTPUT BITMAP IS OFF.
C        
C PROGRAM HISTORY LOG:
C   96-04-10  IREDELL
C 1999-04-08  IREDELL  SPLIT IJKGDS INTO TWO PIECES
C 2000-02-07  GILBERT  ENSURE THAT VECTOR COMPONENTS ARE ROTATED
C                      TO THE GRID ORIENTATION.
C 2001-06-18  IREDELL  INCLUDE MINIMUM MASK PERCENTAGE OPTION
C 2002-01-17  IREDELL  SAVE DATA FROM LAST CALL FOR OPTIMIZATION
C 2007-05-22  IREDELL  EXTRAPOLATE UP TO HALF A GRID CELL
C 2007-10-30  IREDELL  CORRECT NORTH POLE INDEXING PROBLEM,
C                      UNIFY MASKED AND NON-MASKED ALGORITHMS,
C                      AND SAVE WEIGHTS FOR PERFORMANCE.
C
C USAGE:    CALL POLATEV1(IPOPT,KGDSI,KGDSO,MI,MO,KM,IBI,LI,UI,VI,
C    &                    NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
C
C   INPUT ARGUMENT LIST:
C     IPOPT    - INTEGER (20) INTERPOLATION OPTIONS
C                IPOPT(1)=0 FOR STRAIGHT BICUBIC;
C                IPOPT(1)=1 FOR CONSTRAINED BICUBIC WHERE VALUE IS
C                CONFINED WITHIN THE RANGE OF THE SURROUNDING 4 POINTS.
C                IPOPT(2) IS MINIMUM PERCENTAGE FOR MASK
C                (DEFAULTS TO 50 IF IPOPT(2)=-1)
C     KGDSI    - INTEGER (200) INPUT GDS PARAMETERS AS DECODED BY W3FI63
C     KGDSO    - INTEGER (200) OUTPUT GDS PARAMETERS
C                (KGDSO(1)<0 IMPLIES RANDOM STATION POINTS)
C     MI       - INTEGER SKIP NUMBER BETWEEN INPUT GRID FIELDS IF KM>1
C                OR DIMENSION OF INPUT GRID FIELDS IF KM=1
C     MO       - INTEGER SKIP NUMBER BETWEEN OUTPUT GRID FIELDS IF KM>1
C                OR DIMENSION OF OUTPUT GRID FIELDS IF KM=1
C     KM       - INTEGER NUMBER OF FIELDS TO INTERPOLATE
C     IBI      - INTEGER (KM) INPUT BITMAP FLAGS
C     LI       - LOGICAL*1 (MI,KM) INPUT BITMAPS (IF SOME IBI(K)=1)
C     UI       - REAL (MI,KM) INPUT U-COMPONENT FIELDS TO INTERPOLATE
C     VI       - REAL (MI,KM) INPUT V-COMPONENT FIELDS TO INTERPOLATE
C     NO       - INTEGER NUMBER OF OUTPUT POINTS (ONLY IF KGDSO(1)<0)
C     RLAT     - REAL (NO) OUTPUT LATITUDES IN DEGREES (IF KGDSO(1)<0)
C     RLON     - REAL (NO) OUTPUT LONGITUDES IN DEGREES (IF KGDSO(1)<0)
C     CROT     - REAL (NO) VECTOR ROTATION COSINES (IF KGDSO(1)<0)
C     SROT     - REAL (NO) VECTOR ROTATION SINES (IF KGDSO(1)<0)
C                (UGRID=CROT*UEARTH-SROT*VEARTH;
C                 VGRID=SROT*UEARTH+CROT*VEARTH)
C
C   OUTPUT ARGUMENT LIST:
C     NO       - INTEGER NUMBER OF OUTPUT POINTS (ONLY IF KGDSO(1)>=0)
C     RLAT     - REAL (MO) OUTPUT LATITUDES IN DEGREES (IF KGDSO(1)>=0)
C     RLON     - REAL (MO) OUTPUT LONGITUDES IN DEGREES (IF KGDSO(1)>=0)
C     CROT     - REAL (NO) VECTOR ROTATION COSINES (IF KGDSO(1)>=0)
C     SROT     - REAL (NO) VECTOR ROTATION SINES (IF KGDSO(1)>=0)
C                (UGRID=CROT*UEARTH-SROT*VEARTH;
C                 VGRID=SROT*UEARTH+CROT*VEARTH)
C     IBO      - INTEGER (KM) OUTPUT BITMAP FLAGS
C     LO       - LOGICAL*1 (MO,KM) OUTPUT BITMAPS (ALWAYS OUTPUT)
C     UO       - REAL (MO,KM) OUTPUT U-COMPONENT FIELDS INTERPOLATED
C     VO       - REAL (MO,KM) OUTPUT V-COMPONENT FIELDS INTERPOLATED
C     IRET     - INTEGER RETURN CODE
C                0    SUCCESSFUL INTERPOLATION
C                2    UNRECOGNIZED INPUT GRID OR NO GRID OVERLAP
C                3    UNRECOGNIZED OUTPUT GRID
C
C SUBPROGRAMS CALLED:
C   GDSWIZ       GRID DESCRIPTION SECTION WIZARD
C   IJKGDS0      SET UP PARAMETERS FOR IJKGDS1
C   (IJKGDS1)    RETURN FIELD POSITION FOR A GIVEN GRID POINT
C   (MOVECT)     MOVE A VECTOR ALONG A GREAT CIRCLE
C   POLFIXV      MAKE MULTIPLE POLE VECTOR VALUES CONSISTENT
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 77
C
C$$$
      IMPLICIT NONE
      INTEGER,INTENT(IN):: IPOPT(20),KGDSI(200),KGDSO(200),MI,MO,KM
      INTEGER,INTENT(IN):: IBI(KM)
      LOGICAL*1,INTENT(IN):: LI(MI,KM)
      REAL,INTENT(IN):: UI(MI,KM),VI(MI,KM)
      INTEGER,INTENT(INOUT):: NO
      REAL,INTENT(INOUT):: RLAT(MO),RLON(MO),CROT(MO),SROT(MO)
      INTEGER,INTENT(OUT):: IBO(KM)
      LOGICAL*1,INTENT(OUT):: LO(MO,KM)
      REAL,INTENT(OUT):: UO(MO,KM),VO(MO,KM)
      INTEGER,INTENT(OUT):: IRET
      REAL XPTS(MO),YPTS(MO)
      INTEGER IX(4),JY(4)
      REAL WX(4),WY(4)
      INTEGER IJKGDSA(20)
      REAL,PARAMETER:: FILL=-9999.
      INTEGER MCON,MP,N,I,J,K,NK,NV,IJKGDS1
      REAL PMP,XI,YJ,XF,YF,U,V,W,DUM,UMIN,UMAX,VMIN,VMAX
      REAL XPTI(MI),YPTI(MI),RLOI(MI),RLAI(MI),CROI(MI),SROI(MI)
      REAL CM,SM,UROT,VROT
      INTEGER,SAVE:: KGDSIX(200)=-1,KGDSOX(200)=-1,NOX=-1,IRETX=-1
      INTEGER,ALLOCATABLE,SAVE:: NXY(:,:,:),NC(:)
      REAL,ALLOCATABLE,SAVE:: RLATX(:),RLONX(:),CROTX(:),SROTX(:)
      REAL,ALLOCATABLE,SAVE:: WXY(:,:,:),CXY(:,:,:),SXY(:,:,:)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  SET PARAMETERS
      IRET=0
      NO=0
      MCON=IPOPT(1)
      MP=IPOPT(2)
      IF(MP.EQ.-1.OR.MP.EQ.0) MP=50
      IF(MP.LT.0.OR.MP.GT.100) IRET=32
      PMP=MP*0.01
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  SAVE OR SKIP WEIGHT COMPUTATION
      IF(IRET.EQ.0.AND.(KGDSO(1).LT.0.OR.
     &    ANY(KGDSI.NE.KGDSIX).OR.ANY(KGDSO.NE.KGDSOX))) THEN
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  COMPUTE NUMBER OF OUTPUT POINTS AND THEIR LATITUDES AND LONGITUDES.
        IF(KGDSO(1).GE.0) THEN
          CALL GDSWIZ(KGDSO, 0,MO,FILL,XPTS,YPTS,RLON,RLAT,NO,
     &                1,CROT,SROT)
          IF(NO.EQ.0) IRET=3
        ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  LOCATE INPUT POINTS
        CALL GDSWIZ(KGDSI,-1,NO,FILL,XPTS,YPTS,RLON,RLAT,NV,0,DUM,DUM)
        IF(IRET.EQ.0.AND.NV.EQ.0) IRET=2
        CALL GDSWIZ(KGDSI, 0,MI,FILL,XPTI,YPTI,RLOI,RLAI,NV,1,CROI,SROI)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  ALLOCATE AND SAVE GRID DATA
        KGDSIX=KGDSI
        KGDSOX=KGDSO
        IF(NOX.NE.NO) THEN
          IF(NOX.GE.0) DEALLOCATE(RLATX,RLONX,CROTX,SROTX,NC,
     &                            NXY,WXY,CXY,SXY)
          ALLOCATE(RLATX(NO),RLONX(NO),CROTX(NO),SROTX(NO),NC(NO),
     &             NXY(4,4,NO),WXY(4,4,NO),CXY(4,4,NO),SXY(4,4,NO))
          NOX=NO
        ENDIF
        IRETX=IRET
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  COMPUTE WEIGHTS
        IF(IRET.EQ.0) THEN
          CALL IJKGDS0(KGDSI,IJKGDSA)
C$OMP PARALLEL DO
C$OMP&PRIVATE(N,XI,YJ,IX,JY,XF,YF,J,I,WX,WY,CM,SM)
          DO N=1,NO
            RLONX(N)=RLON(N)
            RLATX(N)=RLAT(N)
            CROTX(N)=CROT(N)
            SROTX(N)=SROT(N)
            XI=XPTS(N)
            YJ=YPTS(N)
            IF(XI.NE.FILL.AND.YJ.NE.FILL) THEN
              IX(1:4)=FLOOR(XI-1)+(/0,1,2,3/)
              JY(1:4)=FLOOR(YJ-1)+(/0,1,2,3/)
              XF=XI-IX(2)
              YF=YJ-JY(2)
              DO J=1,4
                DO I=1,4
                  NXY(I,J,N)=IJKGDS1(IX(I),JY(J),IJKGDSA)
                ENDDO
              ENDDO
              IF(MINVAL(NXY(1:4,1:4,N)).GT.0) THEN
C  BICUBIC WHERE 16-POINT STENCIL IS AVAILABLE
                NC(N)=1
                WX(1)=XF*(1-XF)*(2-XF)/(-6.)
                WX(2)=(XF+1)*(1-XF)*(2-XF)/2.
                WX(3)=(XF+1)*XF*(2-XF)/2.
                WX(4)=(XF+1)*XF*(1-XF)/(-6.)
                WY(1)=YF*(1-YF)*(2-YF)/(-6.)
                WY(2)=(YF+1)*(1-YF)*(2-YF)/2.
                WY(3)=(YF+1)*YF*(2-YF)/2.
                WY(4)=(YF+1)*YF*(1-YF)/(-6.)
              ELSE
C  BILINEAR ELSEWHERE NEAR THE EDGE OF THE GRID
                NC(N)=2
                WX(1)=0
                WX(2)=(1-XF)
                WX(3)=XF
                WX(4)=0
                WY(1)=0
                WY(2)=(1-YF)
                WY(3)=YF
                WY(4)=0
              ENDIF
              DO J=1,4
                DO I=1,4
                  WXY(I,J,N)=WX(I)*WY(J)
                  IF(NXY(I,J,N).GT.0) THEN
                    CALL MOVECT(RLAI(NXY(I,J,N)),RLOI(NXY(I,J,N)),
     &                          RLAT(N),RLON(N),CM,SM)
                    CXY(I,J,N)=CM*CROI(NXY(I,J,N))+SM*SROI(NXY(I,J,N))
                    SXY(I,J,N)=SM*CROI(NXY(I,J,N))-CM*SROI(NXY(I,J,N))
                  ENDIF
                ENDDO
              ENDDO
            ELSE
              NC(N)=0
            ENDIF
          ENDDO
        ENDIF
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  INTERPOLATE OVER ALL FIELDS
      IF(IRET.EQ.0.AND.IRETX.EQ.0) THEN
        IF(KGDSO(1).GE.0) THEN
          NO=NOX
          DO N=1,NO
            RLON(N)=RLONX(N)
            RLAT(N)=RLATX(N)
            CROT(N)=CROTX(N)
            SROT(N)=SROTX(N)
          ENDDO
        ENDIF
C$OMP PARALLEL DO
C$OMP&PRIVATE(NK,K,N,U,V,W,UMIN,UMAX,VMIN,VMAX,UROT,VROT,J,I)
        DO NK=1,NO*KM
          K=(NK-1)/NO+1
          N=NK-NO*(K-1)
          IF(NC(N).GT.0) THEN
            U=0
            V=0
            W=0
            IF(MCON.GT.0) UMIN=HUGE(UMIN)
            IF(MCON.GT.0) UMAX=-HUGE(UMAX)
            IF(MCON.GT.0) VMIN=HUGE(VMIN)
            IF(MCON.GT.0) VMAX=-HUGE(VMAX)
            DO J=NC(N),5-NC(N)
              DO I=NC(N),5-NC(N)
                IF(NXY(I,J,N).GT.0.AND.
     &             (IBI(K).EQ.0.OR.LI(NXY(I,J,N),K))) THEN
                  UROT=CXY(I,J,N)*UI(NXY(I,J,N),K)-
     &                 SXY(I,J,N)*VI(NXY(I,J,N),K)
                  VROT=SXY(I,J,N)*UI(NXY(I,J,N),K)+
     &                 CXY(I,J,N)*VI(NXY(I,J,N),K)
                  U=U+WXY(I,J,N)*UROT
                  V=V+WXY(I,J,N)*VROT
                  W=W+WXY(I,J,N)
                  IF(MCON.GT.0) UMIN=MIN(UMIN,UROT)
                  IF(MCON.GT.0) UMAX=MAX(UMAX,UROT)
                  IF(MCON.GT.0) VMIN=MIN(VMIN,VROT)
                  IF(MCON.GT.0) VMAX=MAX(VMAX,VROT)
                ENDIF
              ENDDO
            ENDDO
            LO(N,K)=W.GE.PMP
            IF(LO(N,K)) THEN
              UROT=CROT(N)*U-SROT(N)*V
              VROT=SROT(N)*U+CROT(N)*V
              UO(N,K)=UROT/W
              VO(N,K)=VROT/W
              IF(MCON.GT.0) UO(N,K)=MIN(MAX(UO(N,K),UMIN),UMAX)
              IF(MCON.GT.0) VO(N,K)=MIN(MAX(VO(N,K),VMIN),VMAX)
            ELSE
              UO(N,K)=0.
              VO(N,K)=0.
            ENDIF
          ELSE
            LO(N,K)=.FALSE.
            UO(N,K)=0.
            VO(N,K)=0.
          ENDIF
        ENDDO
        DO K=1,KM
          IBO(K)=IBI(K)
          IF(.NOT.ALL(LO(1:NO,K))) IBO(K)=1
        ENDDO
        IF(KGDSO(1).EQ.0) CALL POLFIXV(NO,MO,KM,RLAT,RLON,IBO,LO,UO,VO)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ELSE
        IF(IRET.EQ.0) IRET=IRETX
        NO=0 
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
