C      SUBROUTINE PROF(NHB,LRSTRT,ITAG,LCLAS1)
C      SUBROUTINE PROF(ITAG,LCLAS1)
      SUBROUTINE PROF_NMM_NET(filename,startdate,ITAG,INCR)
C
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .
C   SUBROUTINE:  PROF        PROFILE SOUNDINGS
C   PRGRMMR: BLACK           ORG: W/NP22     DATE: 99-04-22
C
C ABSTRACT:  THIS ROUTINE GENERATES THE RAW PROFILE SOUNDING
C            OUTPUT FILES FROM THE FORECAST RESTRT FILE AND
C            AUXILIARY FILES
C
C PROGRAM HISTORY LOG:
C   99-04-22  T BLACK - ORIGINATOR
C   02-07-01  G MANIKIN - FIXED PROBLEM WITH DHCNVC AND DHRAIN
C                          COMPUTATIONS - SEE COMMENTS BELOW
C   03-04-01  M PYLE - BEGAN CONVERTING FOR WRF
C   04-05-26  M PYLE - MADE CHANGES FOR WRF-NMM
C
C USAGE:  CALL PROF FROM PROGRAM POST0
C
C   INPUT ARGUMENT LIST:
C     NHB    - THE UNIT NUMBER FOR READING THE NHB FILE
C     LRSTRT - THE UNIT NUMBER FOR READING THE RESTRT FILE
C     ITAG   - THE FORECAST HOUR WE ARE DEALING WITH
C     LCLAS1 - THE UNIT NUMBER FOR WRITING THE PROFILE DATA
C
C   OUTPUT ARGUMENT LIST:
C     NONE
C
C   SUBPROGRAMS CALLED:
C     UNIQUE:
C
C-----------------------------------------------------------------------
c      use vrbls3d
c      use vrbls2d
c      use soil
c      use masks
       use kinds, only             : i_llong
C
      include 'wrf_io_flags.h'
      include 'mpif.h'

!      INCLUDE "parmeta"
      INCLUDE "parmsoil"
C-----------------------------------------------------------------------
                             P A R A M E T E R
     & (NSTAT=1504,LCL1ML=15,LCL1SL=52
     &, D608=0.608)

C-----------------------------------------------------------------------
C
C    PARMS FOR HOURLY PROFILER OUTPUT
C      NSTAT - MAX NUMBER OF STATIONS
C      NWORDM - DIMENSION OF OUTPUT ARRAY, MUST BE LARGE ENOUGH
C          TO HOLD ALL VARIABLES
C          (MAX NO MULTI-LAYER VARIABLES*LM + NO OF SINGLE LAYER VARS)
C      LCL1ML - NUMBER OF MULTI-LAYER VARIABLES OUTPUT FOR CLASS 1
C      LCL1SL - NUMBER OF SINGLE LAYER VARIABLES OUTPUT FOR CLASS 1
C
C------------------------------------------------------------------------
                             P A R A M E T E R
     & (ITB=76,JTB=134)
                             P A R A M E T E R
     & (A2=17.2693882,A3=273.16,A4=35.86,PQ0=379.90516,DTR=1.74532925E-2
     &, G=9.81,GI=1./G,RD=287.04,CP=1004.6,CAPA=RD/CP,RHCRIT=0.9999)
C
      PARAMETER (GAMMA=6.5/1000.,ZSL=0.0)
      PARAMETER (TAUCR=RD*GI*290.66,CONST=0.005*G/RD)
      PARAMETER (GORD=G/RD,DP=60.E2)
C------------------------------------------------------------------------
                             R E A L
     & STNLAT(NSTAT),STNLON(NSTAT)
                             I N T E G E R
     & IDSTN(NSTAT),IHINDX(NSTAT),JHINDX(NSTAT)
     &,             IVINDX(NSTAT),JVINDX(NSTAT)
C

	REAL, ALLOCATABLE:: DETA(:),RDETA(:),AETA(:),UL(:)
     &,RES(:),FIS(:),THS(:),HBOT(:)
     &,CFRACL(:),CFRACM(:),CFRACH(:),SNO(:)
     &,SOILTB(:),SFCEXC(:),SMSTAV(:),SMSTOT(:)
     &,Z0(:),CZEN(:),CZMEAN(:),U00(:),SR(:),CPRATE(:)
     &,ACPREC(:),CUPREC(:),ACSNOW(:),ACSNOM(:)
     &,SSROFF(:),BGROFF(:),SFCSHX(:),SFCLHX(:)
     &,SUBSHX(:),SNOPCX(:),ASWIN(:),ASWOUT(:)
     &,ASWTOA(:),ALWIN(:),ALWOUT(:),ALWTOA(:)
     &,TSHLTR(:),TSHLTR_hold(:),QSHLTR(:),PSHLTR(:)
     &,TH10(:),Q10(:),U10(:),V10(:)
     &,TLMIN(:),TLMAX(:)
     &,SMC(:,:),CMC(:),STC(:,:),SH2O(:,:)
     &,VEGFRC(:),POTFLX(:),PSLP(:),PDSL1(:)
     &,EGRID2(:),SM(:),SICE(:)
     &,HBM2(:),FACTR(:)
     &,PTBL(:,:),TTBL(:,:),VEGFRA(:)
     &,T(:,:),Q(:,:),U(:,:),V(:,:),Q2(:,:)
     &,OMGALF(:,:),CWM(:,:),TRAIN(:,:),TCUCN(:,:)
     &,F_RAIN(:,:),F_ICE(:,:),CLDFRA(:,:)
     &,F_RIMEF(:,:)
     &,RSWTT(:,:),RLWTT(:,:),CCR(:,:),RTOP(:,:)
     &,HTM(:,:),OMGA(:,:)
     &,PRODAT(:),FPACK(:)
     &,STATPR(:),STACPR(:),STAEVP(:)
     &,STAPOT(:),STASHX(:),STASUB(:),STAPCX(:)
     &,STASWI(:),STASWO(:),STALWI(:),STALWO(:)
     &,STALWT(:),STASWT(:),STASNM(:),STASRF(:)
     &,STABRF(:),STASNO(:),DHCNVC(:,:),DHRAIN(:,:)
     &,STADHC(:),STADHR(:)
     &,ACPREC0(:),CUPREC0(:),SFCLHX0(:),POTFLX0(:)
     &,SFCSHX0(:),SUBSHX0(:),SNOPCX0(:),ASWIN0(:)
     &,ASWOUT0(:),ALWIN0(:),ALWOUT0(:),ALWTOA0(:)
     &,ASWTOA0(:),ACSNOW0(:),ACSNOM0(:),SSROFF0(:)
     &,BGROFF0(:)
     &,TCUCN0(:,:),TRAIN0(:,:)
C
!     & DUMSOIL(NSOIL)
!     &,DUM3D(IM+1,JM+1,LM+1),DUM3D2(IM+1,JM+1,LM+1)
!     &,DUM3D3(IM+1,JM+1,LM+1),DUMMY2(IM,JM),DUMMY(IM,JM)
!     &,PD(:),PDS(:),GDLAT(IM,JM),GDLON(IM,JM)
!       REAL:: PMID(NSTAT,LM),PINT(NSTAT,LM+1),LATSTART,LONSTART
!       REAL:: W(NSTAT,LM+1),WH(NSTAT,LM),PDTOP,PT

	real, allocatable:: DUMSOIL(:)
     &,DUM3D(:,:,:),DUM3D2(:,:,:)
     &,DUM3D3(:,:,:),DUMMY2(:,:),DUMMY(:,:)
     &,PD(:),PDS(:),GDLAT(:,:),GDLON(:,:)
     &,PMID(:,:),PINT(:,:)
     &,W(:,:),WH(:,:)
	

        real, allocatable:: CROT(:),SROT(:)
C------------------------------------------------------------------------
C

        integer, allocatable:: IW(:,:), IDUM(:,:),LMH(:,:),IDUMMY(:,:)

                             I N T E G E R
     & IDAT(3),IDAT0(3),GDS(200),SYSDEPINFO
C------------------------------------------------------------------------
                             L O G I C A L
     & RUN,RESTRT,FRST
C------------------------------------------------------------------------
                             C H A R A C T E R
     & RSTFIL*98,RESTHR*4,LABEL*32,CISTAT*8,CIDSTN(NSTAT)*8
     &,FNAME*256,ENVAR*98,BLANK*4

C	new stuff
      character(len=31) :: VarName,varin
	character(len=256) :: fileName
	character(len=256) :: fileNamehold
       character(len=256) :: newname
      integer :: Status
      character(len=19):: startdate,datestr,datestrold

	real:: rinc(5)
	integer:: IDATE(8),JDATE(8)

C------------------------------------------------------------------------
      DATA BLANK/'    '/
C------------------------------------------------------------------------
C***
C***  READ IN THE INFORMATION FILE ABOUT THE SOUNDINGS
C***

	write(6,*) 'filename= ', filename
	write(6,*) 'startedate= ', startdate
c      len=lnblnk(filename)   !total defined length 
       len=index(filename,' ')-1  !actual filename length
       write(6,*) 'actual file name LEN= ', LEN

	write(6,*) 'inside PROF_NMM_NET'

	datestr=startdate

      REWIND 19
C
      READ(19)NUMSTA,IDSTN,STNLAT,STNLON
     1,       IHINDX,JHINDX,IVINDX,JVINDX,CIDSTN
      WRITE(6,20)NUMSTA
   20 FORMAT('INIT:  NUMBER OF PROFILE STATIONS ',I5)
	if (ITAG .eq. 0) then
      WRITE(6,30)(IDSTN(N),STNLAT(N)/DTR,STNLON(N)/DTR
     1,               IHINDX(N),JHINDX(N),IVINDX(N),JVINDX(N)
     2,               CIDSTN(N),N=1,NUMSTA)
	endif
   30 FORMAT(3X,I7,2F9.3,4I9,4X,A10)

c	if (ITAG .eq. 0) then
	  FRST=.TRUE.
c	else 
c	  FRST=.FALSE.
c	endif


       if ( frst ) then
         frst = .false.
        write(0,*) 'calling ext_ncd_ioinit: '
c        CALL ext_ncd_ioinit(DataHandle, Status)
         CALL ext_ncd_ioinit(SYSDEPINFO, Status)
          write(0,*) 'CALLed ioinit', Status
c	write(6,*) 'filename early in PROF= ', trim(filename)
c         CALL ext_ncd_open_for_read(trim(fileName), 0, 0, " ",
c    &  DataHandle, Status)
	write(6,*) 'filename early in PROF= ', filename(1:len)
         CALL ext_ncd_open_for_read(fileName(1:len), 0, 0, " ",
     &  DataHandle, Status)
          print*,'CALLed open for read', Status
       else
           Status = 0
       endif
       if ( Status /= 0 ) then
         print*,'error opening ',fileName(1:len), ' Status = ', Status ; stop
       endif

C Getting start time
      CALL ext_ncd_get_dom_ti_char(DataHandle
     1 ,'START_DATE',startdate, status )
        print*,'startdate= ',startdate
      jdate=0
      idate=0
      read(startdate,15)iyear,imn,iday,ihrst
         IDATE(1)=iyear
         IDATE(2)=imn
         IDATE(3)=iday
         IDATE(5)=ihrst

 15   format(i4,1x,i2,1x,i2,1x,i2)
      print*,'start yr mo day hr =',iyear,imn,iday,ihrst

        call ext_ncd_get_dom_ti_integer(DataHandle,
     &   'WEST-EAST_GRID_DIMENSION',itmp
     + ,1,ioutcount,istatus)

        write(6,*) 'west-east dimension: ', itmp
        IM=itmp-1

        call ext_ncd_get_dom_ti_integer(DataHandle,
     &   'SOUTH-NORTH_GRID_DIMENSION',itmp
     + ,1,ioutcount,istatus)

        write(6,*) 'south-north dimension: ', itmp
        JM=itmp-1

        call ext_ncd_get_dom_ti_integer(DataHandle,
     &   'BOTTOM-TOP_GRID_DIMENSION',itmp
     + ,1,ioutcount,istatus)

        write(6,*) 'bottom-top dimension: ', itmp
        LM=itmp-1

        write(6,*) 'allocate with IM, JM, LM: ', IM, JM, LM

!  The end j row is going to be jend_2u for all variables except for V.
	JSTA_2L=1
	JEND_2U=JM
      JS=JSTA_2L
      JE=JEND_2U
      IF (JEND_2U.EQ.JM) THEN
       JEV=JEND_2U+1
      ELSE
       JEV=JEND_2U
      ENDIF

	write(6,*) 'JS, JE: ', JS, JE



! INSERT FROM INITPOST_NMM_BIN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call ext_ncd_get_dom_ti_real(DataHandle,'DX',tmp
     + ,1,ioutcount,istatus)
        dxval=nint(tmp*1000.) ! E-grid dlamda in degree
!        write(6,*) 'dxval= ', dxval

        call ext_ncd_get_dom_ti_real(DataHandle,'DY',tmp
     + ,1,ioutcount,istatus)
        dyval=nint(1000.*tmp)
!        write(6,*) 'dyval= ', dyval

        call ext_ncd_get_dom_ti_real(DataHandle,'DT',tmp
     + ,1,ioutcount,istatus)
        DT=tmp
        write(6,*) 'DT= ', DT

        call ext_ncd_get_dom_ti_real(DataHandle,'CEN_LAT',tmp
     + ,1,ioutcount,istatus)
	TPH0D=tmp
        cenlat=nint(1000.*tmp)
        write(6,*) 'cenlat= ', cenlat

        call ext_ncd_get_dom_ti_real(DataHandle,'CEN_LON',tmp
     + ,1,ioutcount,istatus)
	TLM0D=tmp
        cenlon=nint(1000.*tmp)
        write(6,*) 'cenlon= ', cenlon

        call ext_ncd_get_dom_ti_real(DataHandle,'TRUELAT1',tmp
     + ,1,ioutcount,istatus)
        truelat1=nint(1000.*tmp)
!        write(6,*) 'truelat1= ', truelat1

        call ext_ncd_get_dom_ti_real(DataHandle,'TRUELAT2',tmp
     + ,1,ioutcount,istatus)
        truelat2=nint(1000.*tmp)
!        write(6,*) 'truelat2= ', truelat2

        call ext_ncd_get_dom_ti_integer(DataHandle,'MAP_PROJ',itmp
     + ,1,ioutcount,istatus)
        maptype=itmp
!        write(6,*) 'maptype is ', maptype

! former parameter statements
        NWORDM=(LCL1ML+1)*LM+2*LCL1SL
        LRECPR=4*(8+9+LCL1ML*LM+LCL1SL)

	write(6,*) 'NWORDM, LRECPR: ', NWORDM, LRECPR
! former parameter statements

       ALLOCATE( DETA(LM),RDETA(LM),AETA(LM),UL(2*LM)
     &,RES(NUMSTA),FIS(NUMSTA),THS(NUMSTA),HBOT(NUMSTA)
     &,CFRACL(NUMSTA),CFRACM(NUMSTA),CFRACH(NUMSTA),SNO(NUMSTA)
     &,SOILTB(NUMSTA),SFCEXC(NUMSTA),SMSTAV(NUMSTA),SMSTOT(NUMSTA)
     &,Z0(NUMSTA),CZEN(NUMSTA),CZMEAN(NUMSTA),U00(NUMSTA),SR(NUMSTA)
     &,ACPREC(NUMSTA),CUPREC(NUMSTA),ACSNOW(NUMSTA),ACSNOM(NUMSTA)
     &,SSROFF(NUMSTA),BGROFF(NUMSTA),SFCSHX(NUMSTA),SFCLHX(NUMSTA)
     &,SUBSHX(NUMSTA),SNOPCX(NUMSTA),ASWIN(NUMSTA),ASWOUT(NUMSTA)
     &,ASWTOA(NUMSTA),ALWIN(NUMSTA),ALWOUT(NUMSTA),ALWTOA(NUMSTA)
     &,TSHLTR(NUMSTA),TSHLTR_hold(NUMSTA),QSHLTR(NUMSTA),PSHLTR(NUMSTA)
     &,TH10(NUMSTA),Q10(NUMSTA),U10(NUMSTA),V10(NUMSTA)
     &,TLMIN(NUMSTA),TLMAX(NUMSTA)
     &,SMC(NUMSTA,NSOIL),CMC(NUMSTA),STC(NUMSTA,NSOIL)
     &,SH2O(NUMSTA,NSOIL)
     &,VEGFRC(NUMSTA),POTFLX(NUMSTA),PSLP(NUMSTA),PDSL1(NUMSTA)
     &,EGRID2(NUMSTA),SM(NUMSTA),SICE(NUMSTA)
     &,HBM2(NUMSTA),FACTR(NUMSTA)
     &,PTBL(ITB,JTB),TTBL(JTB,ITB),VEGFRA(NUMSTA)
     &,T(NUMSTA,LM),Q(NUMSTA,LM),U(NUMSTA,LM),V(NUMSTA,LM)
     &,Q2(NUMSTA,LM)
     &,OMGALF(NUMSTA,LM),CWM(NUMSTA,LM),TRAIN(NUMSTA,LM)
     &,F_RAIN(NUMSTA,LM),F_ICE(NUMSTA,LM),CLDFRA(NUMSTA,LM)
     &,TCUCN(NUMSTA,LM),F_RIMEF(NUMSTA,LM)
     &,RSWTT(NUMSTA,LM),RLWTT(NUMSTA,LM),CCR(NUMSTA,LM)
     &,RTOP(NUMSTA,LM)
     &,HTM(NUMSTA,LM),OMGA(NUMSTA,LM)
     &,PRODAT(NWORDM),FPACK(NWORDM)
     &,STATPR(NUMSTA),STACPR(NUMSTA),STAEVP(NUMSTA)
     &,STAPOT(NUMSTA),STASHX(NUMSTA),STASUB(NUMSTA),STAPCX(NUMSTA)
     &,STASWI(NUMSTA),STASWO(NUMSTA),STALWI(NUMSTA),STALWO(NUMSTA)
     &,STALWT(NUMSTA),STASWT(NUMSTA),STASNM(NUMSTA),STASRF(NUMSTA)
     &,STABRF(NUMSTA),STASNO(NUMSTA),DHCNVC(LM,NUMSTA)
     &,DHRAIN(LM,NUMSTA)
     &,STADHC(LM),STADHR(LM),CPRATE(NUMSTA)
     &,ACPREC0(NUMSTA),CUPREC0(NUMSTA),SFCLHX0(NUMSTA),POTFLX0(NUMSTA)
     &,SFCSHX0(NUMSTA),SUBSHX0(NUMSTA),SNOPCX0(NUMSTA),ASWIN0(NUMSTA)
     &,ASWOUT0(NUMSTA),ALWIN0(NUMSTA),ALWOUT0(NUMSTA),ALWTOA0(NUMSTA)
     &,ASWTOA0(NUMSTA),ACSNOW0(NUMSTA),ACSNOM0(NUMSTA),SSROFF0(NUMSTA)
     &,BGROFF0(NUMSTA)
     &,TCUCN0(NUMSTA,LM),TRAIN0(NUMSTA,LM))

       ALLOCATE ( DUMSOIL(NSOIL)
     &,DUM3D(IM+1,JM+1,LM+1),DUM3D2(IM+1,JM+1,LM+1)
     &,DUM3D3(IM+1,JM+1,LM+1),DUMMY2(IM,JM),DUMMY(IM,JM)
     &,PD(NUMSTA),PDS(NUMSTA),GDLAT(IM,JM),GDLON(IM,JM)
     &,PMID(NUMSTA,LM),PINT(NUMSTA,LM+1)
     &,W(NUMSTA,LM+1),WH(NUMSTA,LM) )

      ALLOCATE(IW(NUMSTA,LM),IDUM(IM,JM),LMH(IM,JM),IDUMMY(IM,JM))


      VarName='HBM2'
!	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY2,
     &   IM,1,JM,1,IM,JS,JE,1)
      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
	HBM2(N)=DUMMY2(I,J)
      END DO
!	write(6,*) 'DUMMY2(20,20): ', DUMMY2(20,20)

      VarName='SM'
!	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY2,
     &   IM,1,JM,1,IM,JS,JE,1)
!	write(6,*) 'DUMMY2(20,20): ', DUMMY2(20,20)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
	SM(N)=DUMMY2(I,J)
      END DO

      VarName='SICE'
!	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY2,
     &   IM,1,JM,1,IM,JS,JE,1)
!	write(6,*) 'DUMMY2(20,20): ', DUMMY2(20,20)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
	SICE(N)=DUMMY2(I,J)
      END DO

!      VarName='HTM'
!	write(6,*) 'call getVariable for : ', VarName
!        write(6,*) 'LM: ', LM
!      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,
!     &  IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
!	write(6,*) 'DUM3D(20,20,1): ', DUM3D(20,20,1)

      DO L = 1,LM
      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            HTM(N,L) = 1.
      END DO 
      END DO

      VarName='PD'
	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY2,
     &   IM,1,JM,1,IM,JS,JE,1)
	write(6,*) 'DUMMY2(20,20): ', DUMMY2(20,20)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
	PD(N)=DUMMY2(I,J)
      END DO


      VarName='FIS'
!	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY2,
     &   IM,1,JM,1,IM,JS,JE,1)
!	write(6,*) 'DUMMY2(20,20): ', DUMMY2(20,20)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
	FIS(N)=DUMMY2(I,J)
      END DO

      VarName='RES'
!	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY2,
     &   IM,1,JM,1,IM,JS,JE,1)
!	write(6,*) 'DUMMY2(20,20): ', DUMMY2(20,20)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
	RES(N)=DUMMY2(I,J)
      END DO

      VarName='T'
!	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,
     &  IM+1,1,JM+1,LM+1,IM,JS,JE,LM)

      DO L = 1,LM
      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            T(N,L) = DUM3D(I,J,L)
      END DO 
      END DO

      VarName='Q'
!	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,
     &  IM+1,1,JM+1,LM+1,IM,JS,JE,LM)

      DO L = 1,LM
      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            Q(N,L) = DUM3D(I,J,L)
      END DO 
      END DO

      VarName='U'
!	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,
     &  IM+1,1,JM+1,LM+1,IM,JS,JE,LM)

      DO L = 1,LM
      DO N=1,NUMSTA
	I=IVINDX(N)	
	J=JVINDX(N)
            U(N,L) = DUM3D(I,J,L)
      END DO 
      END DO


      VarName='V'
!	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,
     &  IM+1,1,JM+1,LM+1,IM,JS,JE,LM)

      DO L = 1,LM
      DO N=1,NUMSTA
	I=IVINDX(N)	
	J=JVINDX(N)
            V(N,L) = DUM3D(I,J,L)
      END DO 
      END DO



        varname='PDTOP'
!	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,PDTOP,
     &  1,1,1,1,1,1,1,1)

	write(6,*) 'PDTOP: ', PDTOP

        varname='PT'
	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY2,
     &  1,1,1,1,1,1,1,1)
	write(6,*) 'DUMMY2(1,1): ', DUMMY2(1,1)
	PT=DUMMY2(1,1)


!used?

	varname='Z0'
!	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY2,
     &  IM,1,JM,1,IM,JS,JE,1)
!	write(6,*) 'DUMMY2(20,20): ', DUMMY2(20,20)
      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
	Z0(N)=DUMMY2(I,J)
      END DO

	varname='THS'
	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY2,
     &  IM,1,JM,1,IM,JS,JE,1)
	write(6,*) 'DUMMY2(20,20): ', DUMMY2(20,20)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
	THS(N)=DUMMY2(I,J)
      END DO

	varname='ACPREC' ! accum total precip
!	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY2,
     &  IM,1,JM,1,IM,JS,JE,1)
!	write(6,*) 'DUMMY2(20,20): ', DUMMY2(20,20)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
	ACPREC(N)=DUMMY2(I,J)
      END DO

	varname='CUPREC' ! accum cumulus precip
!	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY2,
     &  IM,1,JM,1,IM,JS,JE,1)
!	write(6,*) 'DUMMY2(20,20): ', DUMMY2(20,20)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
	CUPREC(N)=DUMMY2(I,J)
!        CPRATE(N)=-9999.
      END DO 

        varname='TH10'
!        write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY2,
     &  IM,1,JM,1,IM,JS,JE,1)
!        write(6,*) 'DUMMY2(20,20): ', DUMMY2(20,20)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
	TH10(N)=DUMMY2(I,J)
      END DO

       varname='Q10'
!        write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY2,
     &  IM,1,JM,1,IM,JS,JE,1)
!        write(6,*) 'DUMMY2(20,20): ', DUMMY2(20,20)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
	Q10(N)=DUMMY2(I,J)
      END DO

        

        varname='PSHLTR'
        write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY2,
     &  IM,1,JM,1,IM,JS,JE,1)
        write(6,*) 'DUMMY2(20,20): ', DUMMY2(20,20)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
	PSHLTR(N)=DUMMY2(I,J)
      END DO


	varname='TSHLTR'
	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY2,
     &  IM,1,JM,1,IM,JS,JE,1)
	write(6,*) 'DUMMY2(20,20): ', DUMMY2(20,20)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
	TSHLTR_hold(N)=DUMMY2(I,J)
      END DO

	varname='QSHLTR'
!	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY2,
     &  IM,1,JM,1,IM,JS,JE,1)
!	write(6,*) 'DUMMY2(20,20): ', DUMMY2(20,20)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
CHC CONVERT FROM MIXING RATIO TO SPECIFIC HUMIDITY
	QSHLTR(N) = DUMMY2(I,J) / (1.0 + DUMMY2(I,J))
      END DO

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Added based on the binary version: Jan. 20, 2013 by Jun Du
!
! Very confusing story ...
!
! Retrieve hbot => It is named CNVBOT in the model and
! with HBOTS (shallow conv) and HBOTD (deep conv) represent
! the 3 sets of convective cloud base/top arrays tied to the frequency
! that history files are written.
!
! IN THE *MODEL*, array HBOT is similar to CNVBOT but is
! used in radiation and is tied to the frequency of radiation updates.
!
! For historical reasons model array CNVBOT is renamed HBOT
! and manipulated throughout the sounding post.
!
!     varname='HBOT'
      VarName='CNVBOT'
!	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY2,
     &  IM,1,JM,1,IM,JS,JE,1)
!	write(6,*) 'DUMMY2(20,20): ', DUMMY2(20,20)

        DO N=1,NUMSTA
         I=IHINDX(N)
         J=JHINDX(N)
!        HBOT(N)=DUMMY2(I,J)
         HBOT(N)= float(LM)-dummy2(i,j)+1.0
        END DO

      VarName='Q2'
!	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,
     &  IM+1,1,JM+1,LM+1,IM,JS,JE,LM)

        DO L=1,LM
         DO N=1,NUMSTA
          I=IHINDX(N)
          J=JHINDX(N)
          Q2(N,L) = DUM3D (I,J,L)
         END DO
        END DO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        varname='CZEN'
!        write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY2,
     &  IM,1,JM,1,IM,JS,JE,1)
!        write(6,*) 'DUMMY2(20,20): ', DUMMY2(20,20)
      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
	CZEN(N)=DUMMY2(I,J)
      END DO

        varname='CZMEAN'
!        write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY2,
     &  IM,1,JM,1,IM,JS,JE,1)
!        write(6,*) 'DUMMY2(20,20): ', DUMMY2(20,20)
      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
	CZMEAN(N)=DUMMY2(I,J)
      END DO


	varname='GLAT'
	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY2,
     &  IM,1,JM,1,IM,JS,JE,1)
	write(6,*) 'DUMMY2(20,20): ', DUMMY2(20,20)

!       do j = jsta_2l, jend_2u
!        do i = 1, im
!	F(I,J)=.00014584*sin(DUMMY2(I,J))   ! 2*omeg*sin(phi)
!	GDLAT(I,J)=DUMMY2(I,J)*(180./acos(-1.))
!        enddo
!       enddo

	varname='GLON'
	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY2,
     &  IM,1,JM,1,IM,JS,JE,1)
	write(6,*) 'DUMMY2(20,20): ', DUMMY2(20,20)

!       do j = jsta_2l, jend_2u
!        do i = 1, im
!	GDLON(I,J)=(DUMMY2(I,J)*(180./acos(-1.)))
!        enddo
!       enddo



	varname='CWM'
	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,
     &  IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
!	write(6,*) 'DUM3D(20,20,30): ', DUM3D(20,20,30)

      DO L = 1,LM
      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            CWM(N,L) = DUM3D(I,J,L)
      END DO 
      END DO

        varname='F_ICE'
        write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,
     &  IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
!       write(6,*) 'DUM3D(20,20,30): ', DUM3D(20,20,30)

      DO L = 1,LM
      DO N=1,NUMSTA
	I=IHINDX(N)	
        J=JHINDX(N)
            F_ICE(N,L) = DUM3D(I,J,L)
      END DO 
      END DO

	varname='F_RAIN'
	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,
     &  IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
!	write(6,*) 'DUM3D(20,20,30): ', DUM3D(20,20,30)

      DO L = 1,LM
      DO N=1,NUMSTA
	I=IHINDX(N)	
        J=JHINDX(N)
            F_RAIN(N,L) = DUM3D(I,J,L)
      END DO 
      END DO

	varname='F_RIMEF'
	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,
     &  IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
!	write(6,*) 'DUM3D(20,20,30): ', DUM3D(20,20,30)

      DO L = 1,LM
      DO N=1,NUMSTA
	I=IHINDX(N)	
        J=JHINDX(N)
            F_RIMEF(N,L) = DUM3D(I,J,L)
      END DO 
      END DO

 	varname='CLDFRA'
 	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,
     &  IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
!	write(6,*) 'DUM3D(20,20,30): ', DUM3D(20,20,30)

       DO L = 1,LM
       DO N=1,NUMSTA
 	I=IHINDX(N)	
 	J=JHINDX(N)
            CLDFRA(N,L) = DUM3D(I,J,L)
       END DO 
       END DO

C Added: Jan. 20, 2013 by Jun Du
      varname='SR'
	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY2,
     &  IM,1,JM,1,IM,JS,JE,1)
         DO N=1,NUMSTA
          I=IHINDX(N)
          J=JHINDX(N)
          SR(N)=DUMMY2(I,J)
         ENDDO
	write(6,*) 'SR(1): ', SR(1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	varname='CFRACH'
	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY2,
     &  IM,1,JM,1,IM,JS,JE,1)
!	write(6,*) 'DUMMY2(20,20): ', DUMMY2(20,20)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            CFRACH(N)=DUMMY2(I,J)
      ENDDO

	varname='CFRACL'
	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY2,
     &  IM,1,JM,1,IM,JS,JE,1)
!	write(6,*) 'DUMMY2(20,20): ', DUMMY2(20,20)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            CFRACL(N)=DUMMY2(I,J)
      ENDDO

	varname='CFRACM'
!	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY2,
     &  IM,1,JM,1,IM,JS,JE,1)
!	write(6,*) 'DUMMY2(20,20): ', DUMMY2(20,20)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            CFRACM(N)=DUMMY2(I,J)
      ENDDO



! either assign SLDPTH to be the same as eta (which is original
! setup in WRF LSM) or extract thickness of soil layers from wrf
! output

! assign SLDPTH to be the same as eta

!         SLDPTH(1)=0.10
!         SLDPTH(2)=0.3
!         SLDPTH(3)=0.6
!         SLDPTH(4)=1.0

! or get SLDPTH from wrf output

      VarName='SLDPTH'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMSOIL,
     & LM,1,1,1,NSOIL,1,1,1)

! if SLDPTH in wrf output is non-zero, then use it
!      DUMCST=0.0
!      DO N=1,NSOIL
!       DUMCST=DUMCST+SLDPTH2(N)
!      END DO
!      IF(ABS(DUMCST-0.).GT.1.0E-2)THEN
!       DO N=1,NSOIL
!        SLDPTH(N)=SLDPTH2(N)
!       END DO
!      END IF
!      print*,'SLDPTH= ',(SLDPTH(N),N=1,NSOIL)

      VarName='CMC'
        write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            CMC(N) = DUMMY(i,j)
      END DO

        varname='SOILTB'
        write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY2,
     &  IM,1,JM,1,IM,JS,JE,1)
!        write(6,*) 'DUMMY2(20,20): ', DUMMY2(20,20)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
  	    SOILTB(N)=DUMMY2(I,J)
      ENDDO

	varname='VEGFRC'
	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY2,
     &  IM,1,JM,1,IM,JS,JE,1)
!	write(6,*) 'DUMMY2(20,20): ', DUMMY2(20,20)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
  	    VEGFRC(N)=DUMMY2(I,J)

	if (mod(N,20) .eq. 0) write(6,*) 'N, VEGFRC(N): ', 
     &			N, VEGFRC(N)

      ENDDO

      VarName='SH2O'
	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,
     &  IM+1,1,JM+1,LM+1,IM,JS,JE,NSOIL)

      DO L = 1,NSOIL
      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
        SH2O(N,L)=DUM3D(I,J,NSOIL-L+1)
      ENDDO
      ENDDO

      VarName='SMC'
	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,
     &  IM+1,1,JM+1,LM+1,IM,JS,JE,NSOIL)

      DO L = 1,NSOIL
      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)

! flip soil layer again because wrf soil variable vertical indexing
! is the same with eta and vertical indexing was flipped for both
! atmospheric and soil layers within getVariable

            SMC(N,L) = DUM3D(i,j,nsoil-l+1)
       END DO
       END DO
        write(6,*) 'past SMC definition, about to print SMC values'
      print*,'SMC at N=10 ',smc(10,1),smc(10,2)
     &,smc(10,3),smc(10,4)

      VarName='STC'
	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,
     &  IM+1,1,JM+1,LM+1,IM,JS,JE,NSOIL)

      DO L = 1,NSOIL
      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            STC(N,L) = DUM3D(i,j,nsoil-l+1)
       END DO
       END DO
!      print*,'STC at ', N ,' = ',stc(N,1),stc(N,2)
!     &,stc(N,3),stc(N,4)


      VarName='PINT'
	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,
     &  IM+1,1,JM+1,LM+1,IM,JS,JE,LM+1)
	write(6,*) 'PINT(20,20,25): ', DUM3D(20,20,25)

      DO L=1,LM+1
      DO N=1,NUMSTA 
	I=IHINDX(N)	
	J=JHINDX(N)
            PINT (N,L) = dum3d ( i, j, l ) 
	if (L .ge. 2) then
!tst            PMID ( N,L-1 ) = 0.5* EXP( ALOG(DUM3D(I,J,L-1))+
!tst     &               ALOG(DUM3D(I,J,L)) )
            PMID ( N,L-1 ) =  EXP( 0.5* (ALOG(DUM3D(I,J,L-1))+
     &                                   ALOG(DUM3D(I,J,L))) )

!	if (N .eq. 25) then
!	write(6,*) 'ALOG(DUM3D(I,J,L-1), ALOG(DUM3D(I,J,L): ' ,
!     &       ALOG(DUM3D(I,J,L-1)), ALOG(DUM3D(I,J,L))
!	write(6,*) 'L-1, PMID(25,L-1):  ', L-1,PMID(25,L-1)
!	endif

	endif
      END DO
      END DO 



!!!!! CONSTRAIN Q TO A PARTICULAR RH VALUE, FOLLOWING CALRH OF WRFPOST

	do N=1,NUMSTA
 	 do L=1,LM

           QC= (PQ0/PMID(N,L)) *EXP(A2*(T(N,L)-A3)/(T(N,L)-A4))
           RH=Q(N,L)/QC

           IF (RH .gt. RHCRIT) THEN
	   Q(N,L)=0.999*RHCRIT*QC
	   IF (RH .gt. 1.02) THEN
!	   write(6,*) 'reducing RH from: ', RH, ' at N,L: ', N,L
	   ENDIF
	   ENDIF

         enddo
        enddo


!!!!! END RH CONSTRAIN

!!!!! COMPUTE Z
!       do j = jsta_2l, jend_2u
!        do i = 1, im
!            ZINT(I,J,LM+1)=FIS(I,J)/G
!	if (I .eq. 1 .and. J .eq. jsta_2l) then
!                   write(6,*) 'G,ZINT: ', G,ZINT(I,J,LM+1)
!	endif
!            FI(I,J,1)=FIS(I,J)
!        end do
!       end do

! SECOND, INTEGRATE HEIGHT HYDROSTATICLY
!      DO L=LM,1,-1
!       do j = jsta_2l, jend_2u
!        do i = 1, im
!         FI(I,J,2)=HTM(I,J,L)*T(I,J,L)*(Q(I,J,L)*D608+1.0)*RD*
!     1             (ALPINT(I,J,L+1)-ALPINT(I,J,L))+FI(I,J,1)
!         ZINT(I,J,L)=FI(I,J,2)/G
!         if(i.eq.ii.and.j.eq.jj)
!     1  print*,'L,sample HTM,T,Q,ALPINT(L+1),ALPINT(l),ZINT= '
!     2  ,l,HTM(I,J,L),T(I,J,L),Q(I,J,L),ALPINT(I,J,L+1),
!     3  ALPINT(I,J,L),ZINT(I,J,L)
!         FI(I,J,1)=FI(I,J,2)
!        ENDDO
!       ENDDO
!      END DO
!      print*,'finish deriving geopotential in nmm'
!
!      DO L=1,LM
!       DO I=1,IM
!        DO J=JS,JE
!         ZMID(I,J,L)=(ZINT(I,J,L+1)+ZINT(I,J,L))*0.5  ! ave of z
!        ENDDO
!       ENDDO
!      ENDDO




      VarName='W'
!	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,
     &  IM+1,1,JM+1,LM+1,IM,JS,JE,LM+1)

      DO L = 1,LM+1
      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            W(N,L)=DUM3D(I,J,L)
      END DO
      END DO

      DO L = 1,LM
      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            WH(N,L) = (W(N,L)+W(N,L+1))*0.5
      END DO
      END DO

      VarName='SSROFF'
!        write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            SSROFF(N) = dummy(i,j)
      ENDDO

c reading UNDERGROUND RUNOFF
      VarName='BGROFF'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            BGROFF(N) = dummy ( i, j )
!	write(6,*) 'N, BGROFF(N): ', N, BGROFF(N)
      END DO




!     VarName='RLWIN'
!	write(6,*) 'call getVariable for : ', VarName
!     call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
!    &  IM,1,JM,1,IM,JS,JE,1)

!	write(6,*) 'DUMMY(20,20): ', DUMMY(20,20)

      VarName='ALWIN'
!	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            ALWIN(N)=dummy(i,j)
      END DO
!	write(6,*) 'DUMMY(20,20): ', DUMMY(20,20)

!here
      VarName='ALWOUT'
!        write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)
      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            ALWOUT(N)=dummy(i,j)
      END DO

      VarName='ALWTOA'
!        write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)
      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            ALWTOA(N)=dummy(i,j)
      END DO

      VarName='ASWIN'
!        write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            ASWIN(N)=dummy(i,j)
      END DO

      VarName='ASWOUT'
!        write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            ASWOUT(N)=dummy(i,j)
      END DO

      VarName='ASWTOA'
!        write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            ASWTOA(N)=dummy(i,j)
      END DO

      VarName='SFCSHX'
!        write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)
      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            SFCSHX(N)=dummy(i,j)
      END DO

      VarName='SFCLHX'
!        write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)
      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            SFCLHX(N)=dummy(i,j)
	if (IDSTN(N) .eq. 725840) then
	write(6,*) 'BLU value for SFCLHX: ', SFCLHX(N)
	endif

      END DO

      VarName='SUBSHX'
!        write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)
      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            SUBSHX(N)=dummy(i,j)
      END DO

      VarName='SNOPCX'
!        write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)
      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            SNOPCX(N)=dummy(i,j)
      END DO
!here2

      VarName='POTFLX'
!        write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            POTFLX(N)=dummy(i,j)
      END DO

C Added: Jan. 20, 2013 By Jun Du
      VarName='TLMIN'
!        write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)
          DO N=1,NUMSTA
           I=IHINDX(N)
           J=JHINDX(N)
           TLMIN(N)=dummy(i,j)
          END DO

      VarName='TLMAX'
!        write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)
          DO N=1,NUMSTA
           I=IHINDX(N)
           J=JHINDX(N)
           TLMAX(N)=dummy(i,j)
          END DO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      varname='RLWTT'
!	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,
     &  IM+1,1,JM+1,LM+1,IM,JS,JE,LM+1)
          DO L=1,LM
           DO N=1,NUMSTA
            I=IHINDX(N)
            J=JHINDX(N)
            RLWTT(N,L)=DUM3D(I,J,L)
           ENDDO
          ENDDO

      varname='RSWTT'
!	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,
     &  IM+1,1,JM+1,LM+1,IM,JS,JE,LM+1)
          DO L=1,LM
           DO N=1,NUMSTA
            I=IHINDX(N)
            J=JHINDX(N)
            RSWTT(N,L)=DUM3D(I,J,L)
           ENDDO
          ENDDO

      varname='TCUCN'
!	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,
     &  IM+1,1,JM+1,LM+1,IM,JS,JE,LM+1)
          DO L=1,LM
           DO N=1,NUMSTA
            I=IHINDX(N)
            J=JHINDX(N)
            TCUCN(N,L)=DUM3D(I,J,L)
           ENDDO
          ENDDO

      varname='TRAIN'
!	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,
     &  IM+1,1,JM+1,LM+1,IM,JS,JE,LM+1)
          DO L=1,LM
           DO N=1,NUMSTA
            I=IHINDX(N)
            J=JHINDX(N)
            TRAIN(N,L)=DUM3D(I,J,L)
           ENDDO
          ENDDO

      VarName='AVRAIN'
	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,AVRAIN,
     &  1,1,1,1,1,1,1,1)
c         AVRAIN=DUMMY2(1,1)

      VarName='NPHS0'
	write(6,*) 'call getVariable for : ', VarName
      call getIVariable(fileName,DateStr,DataHandle,VarName,NPHS,
     &  1,1,1,1,1,1,1,1)
        write(6,*) 'pulled NPHS from NPHS0 as: ', NPHS

      VarName='NPREC'
!        write(6,*) 'call getVariable for : ', VarName
      call getIVariable(fileName,DateStr,DataHandle,VarName,NPREC,
     &  1,1,1,1,1,1,1,1)
      write(6,*) 'NPREC= ', NPREC

      VarName='NCLOD'
!        write(6,*) 'call getVariable for : ', VarName
      call getIVariable(fileName,DateStr,DataHandle,VarName,NCLOD,
     &  1,1,1,1,1,1,1,1)
      write(6,*) 'NCLOD= ', NCLOD

      VarName='NHEAT'
!        write(6,*) 'call getIVariable for : ', VarName
      call getIVariable(fileName,DateStr,DataHandle,VarName,NHEAT,
     &  1,1,1,1,1,1,1,1)
      write(6,*) 'NHEAT= ', NHEAT

      VarName='NRDLW'
!        write(6,*) 'call getVariable for : ', VarName
      call getIVariable(fileName,DateStr,DataHandle,VarName,NRDLW,
     &  1,1,1,1,1,1,1,1)
      write(6,*) 'NRDLW= ', NRDLW

      VarName='NRDSW'
!        write(6,*) 'call getVariable for : ', VarName
      call getIVariable(fileName,DateStr,DataHandle,VarName,NRDSW,
     &  1,1,1,1,1,1,1,1)
        write(6,*) 'NRDSW= ', NRDSW

      VarName='NSRFC'
!        write(6,*) 'call getVariable for : ', VarName
      call getIVariable(fileName,DateStr,DataHandle,VarName,NSRFC,
     &  1,1,1,1,1,1,1,1)
        write(6,*) 'NSRFC= ', NSRFC

      VarName='ACUTIM'
	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,ACUTIM,
     &  1,1,1,1,1,1,1,1)
      write(6,*) 'ACUTIM= ', ACUTIM

      VarName='ARDLW'
!        write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,ARDLW,
     &  1,1,1,1,1,1,1,1)
!      write(6,*) 'ARDLW= ', ARDLW

      VarName='ARDSW'
!        write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,ARDSW,
     &  1,1,1,1,1,1,1,1)
!        write(6,*) 'ARDSW= ', ARDSW

      VarName='ASRFC'
!        write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,ASRFC,
     &  1,1,1,1,1,1,1,1)
!        write(6,*) 'ASRFC= ', ASRFC

      VarName='APHTIM'
!        write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,APHTIM,
     &  1,1,1,1,1,1,1,1)
        write(6,*) 'APHTIM= ', APHTIM

!  end nmm-core specific vars



c
c reading 10 m wind
      VarName='U10'
!	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            U10(N) = dummy ( i, j )
      END DO

      VarName='V10'
!	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            V10(N) = dummy ( i, j )
      END DO
c

c reading SMSTAV
      VarName='SMSTAV'
!	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            SMSTAV (N) = dummy ( i, j )
      ENDDO

      VarName='SMSTOT'
!        write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)
      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            SMSTOT (N) = dummy ( i, j )
      ENDDO

      VarName='VEGFRA'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)
!      print*,'VEGFRA at ',ii,jj,' = ',DUMMY(ii,jj) 

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
  	    VEGFRA(N)=DUMMY(I,J)

!	if (mod(N,20) .eq. 0) write(6,*) 'N, VEGFRA(N): ', 
!     &			N, VEGFRA(N)

      ENDDO

      VarName='SFCEXC'
!        write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            SFCEXC (N) = dummy ( i, j )
      ENDDO


      VarName='ACSNOW'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            ACSNOW (N) = dummy ( i, j )
      END DO
 
      VarName='ACSNOM'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            ACSNOM (N) = dummy ( i, j )
      END DO

! Which name is correct?
!     VarName='SNOW'
      VarName='SNO'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            SNO(N) = dummy ( i, j )
      END DO

      VarName='CPRATE'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)
          DO N=1,NUMSTA
           I=IHINDX(N)
           J=JHINDX(N)
           CPRATE(N) = dummy ( i, j )
          END DO

	write(6,*) 'past getVariable'

	write(6,*) 'past do loops'



!!!! does closing the datahandle help?
c	call ext_ncd_ioclose(SYSDEPINFO, Status)
	call ext_ncd_ioclose(DataHandale, Status)

	write(6,*) 'past ext_ncd_ioclose'


!!!! DONE GETTING

!      DO L=1,LM+1
      DO L=1,LM
      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)

            IF(ABS(T(N,L)).GT.1.0E-3)
     &        OMGA(N,L) = -W(N,L)*PINT(N,L)*G/
     &                 (RD*T(N,L)*(1.+D608*Q(N,L)))

	if (mod(N,60) .eq. 0 .and. mod(L,10) .eq. 0) then
	write(6,*) 'N, L, W, PINT, T, Q, OMGA: ', N,L,W(N,L),PINT(N,L),
     &  T(N,L),Q(N,L),OMGA(N,L)
	endif
	
      END DO
      END DO

	write(6,*) 'past OMGA def'



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DONE INSERT

C------------------------------------------------------------------------
C***
C***  READ QUANTITIES NEEDED FROM THE NHB FILE
C***



	write(6,*) 'set LMH to : ', LM
	write(6,*) 'IM,jm: ', Im,jm
	do J=1,JM
	 do I=1,IM
	   LMH(I,J)=LM
	 enddo
	enddo

C
c     DO L=1,LM
c       READ(NHB)
c     ENDDO
C
c     READ(NHB)DY,CPGFV,EN,ENT,R,PT,TDDAMP
c    1,        F4D,F4Q,EF4T,DETA,RDETA,AETA



!!!!!!!!!!!!
!!!!!!!!!!!!  Made modifications down to about here.
!!!!!!!!!!!!

	write(6,*) 'TLM0D, TPH0D: ', TLM0D, TPH0D

C
C	GET ROTATION ANGLES FOR WIND
C
	
!	write(6,*) 'hardwired DT: '
!	DT=18

      NTSPH=INT(3600./DT+0.50)

	write(6,*) 'rot angles defined, and ntsph= ', ntsph


C
C------------------------------------------------------------------------
C***

!       DO L = 1, LM
!       DO N=1,NUMSTA
!	  Q2(N,L)=0.
!	ENDDO
!	ENDDO


!	DO L=1,LM
!        DO N=1,NUMSTA
!	  TRAIN(N,L)=-9999.
!	  TCUCN(N,L)=-9999.
!        ENDDO
!      ENDDO

C

!      DO N=1,NUMSTA
!	CZEN(N)=-9999.
!      ENDDO
C

!     DO N=1,NUMSTA
!       I=IHINDX(N)	
!       J=JHINDX(N)
!       HBOT(N)=-9999.
!     ENDDO
C
!      DO N=1,NUMSTA
!	PSLP(N)=-9999.
!      ENDDO

!      DO N=1,NUMSTA
!	SOILTB(N)=-9999.
!	SMSTOT(N)=-9999.
!	SFCEXC(N)=-9999.
!      ENDDO
C
!      DO N=1,NUMSTA
!	CZMEAN(N)=-9999.
!	U00(N)=-9999.
!	SR(N)=-9999.
!      ENDDO
C

C
C
!      DO N=1,NUMSTA
!       SFCSHX(N)=-9999.
!       SFCLHX(N)=-9999.
!       SUBSHX(N)=-9999.
!       SNOPCX(N)=-9999.
!      ENDDO
C
!      DO N=1,NUMSTA
!        ASWIN(N)=-9999.
!        ASWOUT(N)=-9999.
!        ASWTOA(N)=-9999.
!        ALWOUT(N)=-9999.
!        ALWTOA(N)=-9999.
!      ENDDO

!      DO N=1,NUMSTA
!        TH10(N)=-9999.
!        Q10(N)=-9999.
!      ENDDO

C
C------------------------------------------------------------------------
!      DO N=1,NUMSTA
!	POTFLX(N)=-9999.
c	TLMIN(N)=-9999.
c	TLMAX(N)=-9999.
!      ENDDO
C------------------------------------------------------------------------
C***
C***  READ RADIATIVE TEMPERATURE TENDENCIES
C***
!     DO L=1,LM
!       DO N=1,NUMSTA
!        RSWTT(N,L)=-9999.
!        RLWTT(N,L)=-9999.
!       ENDDO
!     ENDDO
C
c     CLOSE(LRSTRT)
C------------------------------------------------------------------------
C***
C***  THE FORECAST HOUR
C***  
c     IFHR=NTSD/NTSPH


	write(6,*) 'past block of code'
	IFHR=ITAG
	write(6,*) 'IFHR: ', IFHR
C------------------------------------------------------------------------
      IF(ITAG.GT.0)THEN
C***
C***  GENERATE THE NAME OF THE PRECEDING RESTRT FILE
C***
        ITAG0=ITAG-INCR


        RINC(1)=0.
        RINC(2)=float(ITAG0)
        RINC(3)=0.
        RINC(4)=0.
        RINC(5)=0.

	write(6,*) 'RINC(2): ', rinc(2)

        call w3movdat(rinc,idate,jdate)
        write(DateStrold,301) JDATE(1),JDATE(2),JDATE(3),JDATE(5)
 301  format(i4,'-',i2.2,'-',i2.2,'_',i2.2,':00:00')


!	IOLD=ITAG0*NTSPH

	write(6,*) 'filename later in PROF: ', filename, '_END'
        len=index(filename,' ')-1
	write(6,*) 'LEN= ', LEN
	write(6,*) 'carried over part: ', filename(1:len-19)
        filename=filename(1:len-19)//DateStrold

	write(6,*) 'old filename is ', trim(filename)
	write(6,*) 'date for old file is: ', datestrold

        write(6,*) 'calling ext_ncd_ioinit'
	
         CALL ext_ncd_ioinit(SYSDEPINFO,Status)
        write(6,*) 'status from call: ', status
        write(6,*) 'trim(filename): ', trim(filename)
         CALL ext_ncd_open_for_read( trim(fileName), 0, 0, " ",
     &  DataHandle, Status)

       if ( Status /= 0 ) then
         print*,'error opening ',fileName, ' Status = ', Status 
	print*, 'quitting at open of EARLIER file'
	STOP
       endif


C***  READ THE PREVIOUS RESTRT FILE
C***
C
!          DO N=1,NUMSTA
!          DO L=1,LM
!	TRAIN0(N,L)=-9999.
!	TCUCN0(N,L)=-9999.
!          ENDDO
!          ENDDO
C


!!!! INSERT FULL READ FROM ABOVE

!!!!!
	DateStr=DateStrold
!!!!!


	write(6,*) 'FROM EARLIER FILE'
        call ext_ncd_get_dom_ti_real(DataHandle,'DX',tmp
     + ,1,ioutcount,istatus)
        dxval=nint(tmp*1000.) ! E-grid dlamda in degree

        call ext_ncd_get_dom_ti_real(DataHandle,'DY',tmp
     + ,1,ioutcount,istatus)
        dyval=nint(1000.*tmp)

        call ext_ncd_get_dom_ti_real(DataHandle,'DT',tmp
     + ,1,ioutcount,istatus)
        DT=tmp

        call ext_ncd_get_dom_ti_real(DataHandle,'CEN_LAT',tmp
     + ,1,ioutcount,istatus)
	TPH0D=tmp
        cenlat=nint(1000.*tmp)

        call ext_ncd_get_dom_ti_real(DataHandle,'CEN_LON',tmp
     + ,1,ioutcount,istatus)
	TLM0D=tmp
        cenlon=nint(1000.*tmp)

        call ext_ncd_get_dom_ti_real(DataHandle,'TRUELAT1',tmp
     + ,1,ioutcount,istatus)
        truelat1=nint(1000.*tmp)

        call ext_ncd_get_dom_ti_real(DataHandle,'TRUELAT2',tmp
     + ,1,ioutcount,istatus)
        truelat2=nint(1000.*tmp)

        call ext_ncd_get_dom_ti_integer(DataHandle,'MAP_PROJ',itmp
     + ,1,ioutcount,istatus)
        maptype=itmp

	varname='ACPREC' ! accum total precip
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY2,
     &  IM,1,JM,1,IM,JS,JE,1)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
	ACPREC0(N)=DUMMY2(I,J)
      END DO

	varname='CUPREC' ! accum cumulus precip
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY2,
     &  IM,1,JM,1,IM,JS,JE,1)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
	CUPREC0(N)=DUMMY2(I,J)
      END DO 

      VarName='SSROFF'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)
          DO N=1,NUMSTA
           I=IHINDX(N)
           J=JHINDX(N)
           SSROFF0(N) = dummy ( i, j )
          END DO

      VarName='BGROFF'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)
          DO N=1,NUMSTA
           I=IHINDX(N)
           J=JHINDX(N)
           BGROFF0(N) = dummy ( i, j )
          END DO

      VarName='ALWIN'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            ALWIN0(N)=dummy(i,j)
      END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!here
      VarName='ALWOUT'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)
      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            ALWOUT0(N)=dummy(i,j)
      END DO

      VarName='ALWTOA'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)
      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            ALWTOA0(N)=dummy(i,j)
      END DO

      VarName='RSWIN'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)

      VarName='RSWOUT'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)

      VarName='ASWIN'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            ASWIN0(N)=dummy(i,j)
      END DO

      VarName='ASWOUT'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            ASWOUT0(N)=dummy(i,j)
      END DO

      VarName='ASWTOA'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            ASWTOA0(N)=dummy(i,j)
      END DO

      VarName='SFCSHX'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)
      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            SFCSHX0(N)=dummy(i,j)
      END DO

      VarName='SFCLHX'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)
      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            SFCLHX0(N)=dummy(i,j)
      END DO

      VarName='SUBSHX'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)
      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            SUBSHX0(N)=dummy(i,j)
      END DO

      VarName='SNOPCX'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)
      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            SNOPCX0(N)=dummy(i,j)
      END DO
!here2

      VarName='SFCUVX'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)

      VarName='POTEVP'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)

      VarName='POTFLX'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            POTFLX0(N)=dummy(i,j)
      END DO

C Added: Jan. 21, 2013 by Jun Du
      varname='TCUCN'
!	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,
     &  IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
          DO L=1,LM
           DO N=1,NUMSTA
            I=IHINDX(N)
            J=JHINDX(N)
            TCUCN0(N,L)=DUM3D(I,J,L)
           ENDDO
          ENDDO

      varname='TRAIN'
!	write(6,*) 'call getVariable for : ', VarName
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM3D,
     &  IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
          DO L=1,LM
           DO N=1,NUMSTA
            I=IHINDX(N)
            J=JHINDX(N)
            TRAIN0(N,L)=DUM3D(I,J,L)
           ENDDO
          ENDDO

      VarName='AVRAIN'
      call getVariable(fileName,DateStr,DataHandle,VarName,AVRAIN0,
     &  1,1,1,1,1,1,1,1)

      VarName='AVCNVC'
      call getVariable(fileName,DateStr,DataHandle,VarName,AVCNVC,
     &  1,1,1,1,1,1,1,1)
          AVCNVC=DUMMY2(1,1)

      VarName='ACUTIM'
      call getVariable(fileName,DateStr,DataHandle,VarName,ACUTIM0,
     &  1,1,1,1,1,1,1,1)
      write(6,*) 'ACUTIM0= ', ACUTIM0

      VarName='APHTIM'
      call getVariable(fileName,DateStr,DataHandle,VarName,APHTIM0,
     &  1,1,1,1,1,1,1,1)
      write(6,*) 'APHTIM0 now : ', APHTIM0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      VarName='NCLOD'
      call getIVariable(fileName,DateStr,DataHandle,VarName,DUM0D,
     &  1,1,1,1,1,1,1,1)

      VarName='NPREC'
      call getIVariable(fileName,DateStr,DataHandle,VarName,DUM0D,
     &  1,1,1,1,1,1,1,1)

      VarName='NHEAT'
      call getIVariable(fileName,DateStr,DataHandle,VarName,DUM0D,
     &  1,1,1,1,1,1,1,1)

      VarName='NRDLW'
      call getIVariable(fileName,DateStr,DataHandle,VarName,DUM0D,
     &  1,1,1,1,1,1,1,1)

      VarName='NRDSW'
      call getIVariable(fileName,DateStr,DataHandle,VarName,DUM0D,
     &  1,1,1,1,1,1,1,1)

      VarName='NSRFC'
      call getIVariable(fileName,DateStr,DataHandle,VarName,DUM0D,
     &  1,1,1,1,1,1,1,1)

      VarName='ARDLW'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM0D,
     &  1,1,1,1,1,1,1,1)

      VarName='ARDSW'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM0D,
     &  1,1,1,1,1,1,1,1)

      VarName='ASRFC'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUM0D,
     &  1,1,1,1,1,1,1,1)

!     VarName='APHTIM'
!     call getVariable(fileName,DateStr,DataHandle,VarName,APHTIM0,
!    &  1,1,1,1,1,1,1,1)

!       write(6,*) 'APHTIM0 now : ', APHTIM0
!  end nmm-core specific vars

      VarName='ACSNOW'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            ACSNOW0(N) = dummy ( i, j )
      END DO
 
      VarName='ACSNOM'
      call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
     &  IM,1,JM,1,IM,JS,JE,1)

      DO N=1,NUMSTA
	I=IHINDX(N)	
	J=JHINDX(N)
            ACSNOM0 (N) = dummy ( i, j )
      END DO

!     VarName='SSROFF'
!     call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
!    &  IM,1,JM,1,IM,JS,JE,1)

!     DO N=1,NUMSTA
!I=IHINDX(N)	
!J=JHINDX(N)
!           SSROFF0(N) = dummy ( i, j )
!     END DO


!     VarName='BGROFF'
!     call getVariable(fileName,DateStr,DataHandle,VarName,DUMMY,
!    &  IM,1,JM,1,IM,JS,JE,1)

!     DO N=1,NUMSTA
!I=IHINDX(N)	
!J=JHINDX(N)
!           BGROFF0(N) = dummy ( i, j )
!     END DO



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! END INSERT FULL READ

!        DO N=1,NUMSTA
!          POTFLX0(N)=-9999.
!        ENDDO
C
      ENDIF

C
Cmp 	IDEALLY, WON'T NEED MANY MODS BELOW THIS POINT
C
C

C------------------------------------------------------------------------
C***
C***  ALL THE DATA IS NOW IN.
C***  CALCULATE CLOUD FRACTION AND CLOUD WATER/ICE ID NUMBER.
C***
C------------------------------------------------------------------------
      UTIM=1.
      US=1.
      CCLIMIT=1.E-3
      CLIMIT =1.E-20
C-----------------------------------------------------------------------
!$OMP parallel do 
      DO N=1,NUMSTA
!       IW(N,1)=-9999
!       CCR(N,1)=-9999.
        PDSL1(N)=PD(N)*RES(N)
      ENDDO

!	write(6,*) 'here b'
C
C------------------QW, QI AND QINT--------------------------------------
C

!!!
!!! skip section for now
!!!
 	goto 221

      DO 220 L=2,LM
C
!$OMP parallel do private(cwmkl,fiq,hh,iwkl,lml,pp,qc,qi,qint,qkl,qw,
!$OMP*                    rqkl,tkl,tmt0,tmt15,u00kl)
      DO 210 N=1,NUMSTA

!	CCR(N,L)=-9999.
!	IW(N,L)=-9999


!      LML=LM-LMH(IHINDX(N),JHINDX(N))
!	write(6,*) 'LML, IHINDX,JHINDX,LMH: ', IHINDX(N), 
!     &            JHINDX(N),LMH(IHINDX(N),JHINDX(N))
!      HH=HTM(N,L)*HBM2(N)
!      TKL=T(N,L)
!      QKL=Q(N,L)
!      CWMKL=CWM(N,L)
!      TMT0=(TKL-273.16)*HH
!      TMT15=AMIN1(TMT0,-15.)*HH
!      AI=0.008855
!      BI=1.
C
!      IF(TMT0.LT.-20.)THEN
!        AI=0.007225
!        BI=0.9674
!      ENDIF
C
Cmp      PP=PDSL1(N)*AETA(L)+PT
!      PP=PMID(N,L)
!      QW=HH*PQ0/PP*EXP(HH*A2*(TKL-A3)/(TKL-A4))
!      QI=QW*(BI+AI*AMIN1(TMT0,0.))
!      QINT=QW*(1.-0.00032*TMT15*(TMT15+15.))
!      IF(TMT0.LE.-40.)QINT=QI
C
C-------------------ICE-WATER ID NUMBER IW------------------------------
C

!	no defs for U00 or UL
!
!	so no U00KL,FIQ,IW,CWM
!
!	write(6,*) 'here c'
!	write(6,*) 'L+LML: ', L+LML
!	write(6,*) 'U00(N): ', U00(N)
!	write(6,*) 'UL(L+LML): ', UL(L+LML)


! test elim this section

!      U00KL=U00(N)+UL(L+LML)*(0.95-U00(N))*UTIM
!      IF(TMT0.LT.-15.)THEN
!        FIQ=QKL-U00KL*QI
!        IF(FIQ.GT.0..OR.CWMKL.GT.CLIMIT)THEN
!          IW(N,L)=1
!        ELSE
!          IW(N,L)=0
!        ENDIF
!      ENDIF
C
!      IF(TMT0.GE.0.)THEN
!        IW(N,L)=0
!      ENDIF
C
!      IF(TMT0.LT.0..AND.TMT0.GE.-15.)THEN
!        IW(N,L)=0
!        IF(IW(N,L-1).EQ.1.AND.CWMKL.GT.CLIMIT)IW(N,L)=1
!      ENDIF
C
!      IWKL=IW(N,L)
C----------------THE SATUATION SPECIFIC HUMIDITY------------------------
!      FIW=FLOAT(IWKL)
!      QC=(1.-FIW)*QINT+FIW*QI
C----------------THE RELATIVE HUMIDITY----------------------------------
!      IF(QC.LE.0.)THEN
!         RQKL=0.
!       ELSE
!         RQKL=QKL/QC
!      ENDIF
C----------------CLOUD COVER RATIO CCR----------------------------------
!      IF(RQKL.GE.0.9999)THEN
!        CCR(N,L)=AMIN1(US,RQKL)
!      ELSE
!        ARG=-1000.*CWMKL/(US-RQKL)
!        ARG=AMAX1(ARG,-25.)
!        CCR(N,L)= RQKL*(1.-EXP(ARG))
!      ENDIF
C----------------------------------------------------------------------
  210                 CONTINUE
  220                 CONTINUE
  221                 continue
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C***
C***  BEGIN THE PROFILE POSTING CODE.
C***
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C***
C***  USE ZERO IN ACCUMULATION ARRAYS AT APPROPRIATE TIMES
C***
       IF(ITAG .eq. 0) THEN
C
C set TLMAX,TLMIN= Lowest Layer Temp if forecast hour (ITAG)=0
C
         DO N=1,NUMSTA
           TLMIN(N)=T(N,LM)
           TLMAX(N)=T(N,LM)
         ENDDO
	write(6,*) 'here (2)...should only see once'
C
C
C what would appropriate if test be here?
C
C
        DO N=1,NUMSTA
C
C*** ZERO ACCUMLATION ARRAYS.
C
          STATPR(N)=0.
          STACPR(N)=0.
          STAEVP(N)=0.
          STAPOT(N)=0.
          STASHX(N)=0.
          STASUB(N)=0.
          STAPCX(N)=0.
          STASWI(N)=0.
          STASWO(N)=0.
          STALWI(N)=0.
          STALWO(N)=0.
          STALWT(N)=0.
          STASWT(N)=0.
          STASNM(N)=0.
          STASNO(N)=0.
          STASRF(N)=0.
          STABRF(N)=0.
          DO L=1,LM
            DHCNVC(L,N)=0.
            DHRAIN(L,N)=0.
          ENDDO
        ENDDO
C
        GO TO 300
       ENDIF
C---------------------------------------------------------------------
C***
C***  WE MUST CHECK TO SEE IF WE ARE 1 HOUR AFTER ANY OF THE 
C***  ACCUMULATION BUCKETS HAVE BEEN EMPTIED.  IF WE ARE AT SUCH A 
C***  TIME THEN WE NEED TO SET TO ZERO THE VARIABLES USED TO HOLD
C***  THE PRECEDING HOUR'S VALUES.
C***
C---------------------------------------------------------------------
C
C
C 	At this point, accumulation buckets are a foreign concept in
C	the WRF model.  
C
C

c     TIME=(NTSD-1)*DT

!	write(6,*) 'hardwired NPHS'
!
!        write(6,*) 'DT is: ', dt
!	if (DT .eq. 18) then
!	NPHS=10
!	elseif (DT .eq. 20) then
!	NPHS=20
!        else
!        NPHS=8
!	endif
 
         write(6,*) 'NPHS set to : ', NPHS

        print*,'ITAG=',ITAG
        IFCST=3600*ITAG
	TIME=IFCST
        write(6,*) 'set time to : ', TIME
	NTSD=(TIME/DT)+1

        write(6,*) 'TIME, DT, NTSD: ', TIME, DT, NTSD

      RESET0=TIME-(NTSD/NPREC)*NPREC*DT
      RESET1=(NPHS-1)*DT+3600.

        write(6,*) 'NPHS, RESET0, RESET1: ', NPHS, RESET0, RESET1

       IF(MOD(NTSD,NPREC).GE.NPHS.AND.RESET0.LE.RESET1)THEN
         DO N=1,NUMSTA
           STATPR(N)=0.
           STACPR(N)=0.
           STASNM(N)=0.
           STASNO(N)=0.
           STASRF(N)=0.
           STABRF(N)=0.
         ENDDO
       ELSE
	write(6,*) 'set STATPR'
        DO N=1,NUMSTA
          STATPR(N)=ACPREC0(N)*1.E3
          STACPR(N)=CUPREC0(N)*1.E3
!         STASNM(N)=ACSNOM0(N)*1.E3
!         STASNO(N)=ACSNOW0(N)*1.E3
!         STASRF(N)=SSROFF0(N)*1.E3
!         STABRF(N)=BGROFF0(N)*1.E3
          STASNM(N)=ACSNOM0(N)
          STASNO(N)=ACSNOW0(N)
          STASRF(N)=SSROFF0(N)
          STABRF(N)=BGROFF0(N)
        ENDDO

	write(6,*) 'past set'
       ENDIF          
C
      RESET0=TIME-(NTSD/NRDSW)*NRDSW*DT
	
	write(6,*) 'TIME, DT: ', TIME, DT
	write(6,*) 'RESET0, RESET1 : ', RESET0, RESET1
	write(6,*) 'NTSD, NRDSW, NPHS: ', NTSD, NRDSW, NPHS

!maybe      IF( MOD(NTSD,NRDSW) .GE. NPHS .AND. RESET0 .LE. RESET1)THEN
      IF( MOD(NTSD,NRDSW) .LE. NPHS .AND. RESET0 .LE. RESET1)THEN
	write(6,*) 'resetting STASWI'
       DO N=1,NUMSTA
         STASWI(N)=0.
         STASWO(N)=0.
         STASWT(N)=0.
       ENDDO
      ELSE
	write(6,*) 'setting STASWI to ASWIN0'
        DO N=1,NUMSTA
          STASWI(N)=ASWIN0(N)
          STASWO(N)=ASWOUT0(N)
          STASWT(N)=ASWTOA0(N)
        ENDDO
      ENDIF
C
      RESET0=TIME-(NTSD/NRDLW)*NRDLW*DT
	write(6,*) 'TIME, DT: ', TIME, DT
	write(6,*) 'RESET0, RESET1 : ', RESET0, RESET1
	write(6,*) 'NTSD, NRDLW, NPHS: ', NTSD, NRDLW, NPHS
!      IF(MOD(NTSD,NRDLW).GE.NPHS.AND.RESET0.LE.RESET1)THEN
      IF(MOD(NTSD,NRDLW).LE.NPHS.AND.RESET0.LE.RESET1)THEN
	write(6,*) 'resetting STALWI'
       DO N=1,NUMSTA
         STALWI(N)=0.
         STALWO(N)=0.
         STALWT(N)=0.
       ENDDO
      ELSE
	write(6,*) 'setting STALWI to ALWIN0'
        DO N=1,NUMSTA
          STALWI(N)=ALWIN0(N)
          STALWO(N)=ALWOUT0(N)
          STALWT(N)=-ALWTOA0(N)
        ENDDO
      ENDIF
C
      RESET0=TIME-(NTSD/NSRFC)*NSRFC*DT
!      IF(MOD(NTSD,NSRFC).GE.NPHS.AND.RESET0.LE.RESET1)THEN
      IF(MOD(NTSD,NSRFC).LE.NPHS.AND.RESET0.LE.RESET1)THEN
       DO N=1,NUMSTA
         STAEVP(N)=0.
         STAPOT(N)=0.
         STASHX(N)=0.
         STASUB(N)=0.
         STAPCX(N)=0.
       ENDDO
      ELSE
        DO N=1,NUMSTA
          STAEVP(N)=SFCLHX0(N)
          STAPOT(N)=POTFLX0(N)
          STASHX(N)=SFCSHX0(N)
          STASUB(N)=SUBSHX0(N)
          STAPCX(N)=SNOPCX0(N)
        ENDDO
      ENDIF
C
      RESET0=TIME-(NTSD/NHEAT)*NHEAT*DT
      IF(MOD(NTSD,NHEAT).GE.NCNVC.AND.RESET0.LE.RESET1)THEN
        DO N=1,NUMSTA
          DO L=1,LM
            DHCNVC(L,N)=0.
            DHRAIN(L,N)=0.
          ENDDO
        ENDDO
      ELSE
        DO N=1,NUMSTA
          DO L=1,LM
             DHCNVC(L,N)=TCUCN0(N,L)
             DHRAIN(L,N)=TRAIN0(N,L)
          ENDDO
        ENDDO
      ENDIF
 
C------------------------------------------------------------------
  300 CONTINUE
C------------------------------------------------------------------
C
C***  FOR ROTATION OF WINDS FROM E-GRID TO GEODETIC ORIENTATION
C***  WE NEED THE TWO QUANTITIES BELOW.
C
      SINPH0=SIN(TPH0D*DTR)
      COSPH0=COS(TPH0D*DTR)
C
C***  INITIAL CALCULATIONS/PREPARATIONS.  WE LOAD SEVERAL
C***  ARRAYS WITH PROFILE VARIABLES.
C
!$OMP parallel do
      DO N=1,NUMSTA
        IF(CZMEAN(N).GT.0.)THEN
          FACTR(N)=CZEN(N)/CZMEAN(N)
        ELSE
          FACTR(N)=0.
        ENDIF
      ENDDO
C
C***  ADJUST SHORTAVE TENDENCIES TO ACCOUNT FOR CHANGE OF SOLAR POSITION
C***  BETWEEN CALLS TO RADIATION
C
!$OMP parallel do
      DO L=1,LM
        DO N=1,NUMSTA
          RSWTT(N,L)=RSWTT(N,L)*FACTR(N)
        ENDDO
      ENDDO
C
C***  COMPUTE RTOP
C
!$OMP parallel do
      DO L=1,LM
        DO N=1,NUMSTA
!          APEL=PT+AETA(L)*PDSL1(N)
          APEL=PMID(N,L)
          RTOP(N,L)=RD*T(N,L)*(1.+0.608*Q(N,L))/APEL
        ENDDO
      ENDDO

        write(6,*) 'PD(20) (dd): ', PD(20)
C
C***  PDS IS SURFACE PRESSURE.
C
!$OMP parallel do 
	DO N=1,NUMSTA
        PDS(N)=PD(N)+PDTOP+PT
        if (PDS(N) .lt. 50000.) then
        write(6,*) 'N,PD,PDTOP,PDS: ', N,PD(N),PDTOP,PDS(N)
        endif
	ENDDO
C
C***  EGRID2 IS THE SURFACE TEMPERATURE.
C
!$OMP parallel do 
      DO N=1,NUMSTA

	TSHLTR(N)= TSHLTR_hold(N)*(PDS((N))*1.E-5)**CAPA
        EGRID2(N)= THS(N)*(PDS((N))*1.E-5)**CAPA
!	write(6,*) 'THETA, T, Psurf : ', N,THS(N),EGRID2(N),
!     &                   PDS(N)
        if (EGRID2(N) .lt. 200.) then
        write(6,*) 'THETA, PDS, CAPA, TSK: ', N,THS(N),
     &           PDS(N),CAPA,EGRID2(N)
        endif

        IF(ACPREC(N).LT.0.)ACPREC(N)=0.
        IF(CUPREC(N).LT.0.)CUPREC(N)=0.

!!! constrain surface RH

           QC=(PQ0/PSHLTR(N))*EXP(A2*(TSHLTR(N)-A3)/(TSHLTR(N)-A4))
           RH=QSHLTR(N)/QC
           IF (RH .gt. RHCRIT) THEN
           write(6,*) 'reducing surface RH from: ', RH, ' at N: ', N
           QSHLTR(N)=0.999*RHCRIT*QC
           ENDIF

      ENDDO
C
C** Compute PSLP using NMC reduction
C
      DO N=1,NUMSTA
       PSFC = PDS(N)
       ZSFC = FIS(N)*GI
       PSLP(N) = PSFC
C
C    COMPUTE LAYER TAU (VIRTUAL TEMP*RD/G).
       TVRT = T(N,LM)*(1.0+D608*Q(N,LM))
       TAU  = TVRT*RD*GI
C
C    COMPUTE TAU AT THE GROUND (Z=ZSFC) AND SEA LEVEL (Z=0)
C    ASSUMING A CONSTANT LAPSE RATE OF GAMMA=6.5DEG/KM.
       TVRSFC = TVRT + (ZSFC- ZSL)*GAMMA
       TAUSFC = TVRSFC*RD*GI
       TVRSL  = TVRT + (ZSFC- ZSL)*GAMMA
       TAUSL  = TVRSL*RD*GI
C
C    IF NEED BE APPLY SHEULL CORRECTION.
       IF ((TAUSL.GT.TAUCR).AND.(TAUSFC.LE.TAUCR)) THEN
          TAUSL=TAUCR
       ELSEIF ((TAUSL.GT.TAUCR).AND.(TAUSFC.GT.TAUCR)) THEN
          TAUSL = TAUCR-CONST*(TAUSFC-TAUCR)**2
       ENDIF
C
C    COMPUTE MEAN TAU.
       TAUAVG = 0.5*(TAUSL+TAUSFC)
C
C    COMPUTE SEA LEVEL PRESSURE.
       IF (FIS(N).GT.1.0) PSLP(N) = PSFC*EXP(ZSFC/TAUAVG)
c      print *,n,idstn(n),pslp(n),tvrt

      ENDDO

C
C***  SET CYCLE, DATE, AND FORECAST TIME.
C
c      IHR  =NTSD/NTSPH+0.5

!MAYBE      read(datestr,15)iyear,imn,iday,ihrst
         IDATE(2)=imn
         IDATE(3)=iday
         IDATE(1)=iyear
         IDATE(5)=ihrst

	IDAT(3)=IDATE(1)
	IDAT(1)=IDATE(2)
	IDAT(2)=IDATE(3)
	
      IYR  =IDAT(3)
      IMNTH=IDAT(1)
      IDAY =IDAT(2)
c      IFCST=(NTSD-1)*DT

!      IFCST=NTSPH*ITAG
      IFCST=3600*ITAG
	IHR=ITAG
	write(6,*) 'IFCST: ', IFCST
C
      WRITE(6,*)' POST PROFILE FOR ',
     1                       IYR,IMNTH,IDAY,IHR
	write(6,*) 'IHRST= ', IHRST
C
C***  SET RTSPH,RTSCU,RTSRA TO 1. OVER THE NUMBER OF TIMES THE
C***  VARIOUS PHYSICS ROUTINES HAVE BEEN
C***  CALLED SINCE LAST OUTPUT OF PROFILER DATA.  NECESSARY FOR
C***  CORRECT AVERAGING OF VARIABLES.
C
      IF(APHTIM.GT.0.)THEN
        RTSPH=1./(APHTIM-APHTIM0)
      ELSE
        RTSPH=1.
      ENDIF
	write(6,*) 'APHTIM, RTSPH: ', APHTIM-APHTIM0, RTSPH
C
C
C  Counter AVRAIN is reset every 6 or 3-h in WRF-NMM, so for stable
c  precip heating rate use counter for convection (ACUTIM)
c  since in the WRF-NMM microphysics (for stable precip physics)
c  is called at the same time as the convection.
C
      IF(ACUTIM.GT.0.)THEN
        RTSCU=1./(ACUTIM-ACUTIM0)
        RTSRA=1./(ACUTIM-ACUTIM0)
      ELSE
        RTSCU=1.
        RTSRA=1.
      ENDIF
        write(6,*) 'ACUTIM, RTSCU: ', ACUTIM, ACUTIM-ACUTIM0, RTSCU

!
!       ACUTIM=0

!     IF(ACUTIM.GT.0.)THEN
!       RTSCU=1./ACUTIM
!     ELSE
!       RTSCU=1.
!     ENDIF
C
!      IF(ARATIM.GT.0.)THEN
!        RTSRA=1./ARATIM
!      ELSE
        RTSRA=1.
!      ENDIF
C
C--------------------------------------------------------------------------
C--------------------------------------------------------------------------
C***
C***  OUTPUT PROFILE DATA.  THE FOLLOWING LOOP IS OVER ALL PROFILE SITES.
C***
C--------------------------------------------------------------------------
	LCLAS1=79

       
      OPEN(UNIT=LCLAS1,ACCESS='DIRECT',RECL=LRECPR,IOSTAT=IER)
C--------------------------------------------------------------------------
	write(6,*) 'RECORD LENGTH = ', LRECPR

      DO 1000 N=1,NUMSTA
        write(0,*) 'working station N: ', N
C
C***  ZERO OUTPUT ARRAY.
C
      DO K=1,NWORDM
        PRODAT(K)=0.
        FPACK(K) =0.
      ENDDO
C
C***  CONSTRUCT HEADER FOR CURRENT PROFILE SITE.  THE HEADER CONTAINS
C***  THE FOLLOWING INFORMATION:  PACKED CYCLE-DATE, FORECAST TIME,
C***  INTEGER STATION ID, STATION LATITUDE, STATION LONGITUDE, STATION
C***  ELEVATION, NUMBER OF VERTICAL LEVELS IN PROFILE, NUMBER OF MULTI-
C***  LEVEL PARAMETERS, NUMBER OF SINGLE LEVEL PARAMETERS, TOTAL LENGTH
C***  (IN WORDS) OF MULTI- AND SINGLE LEVEL DATA, PROFILE CLASS FLAG,
C***  AND A DUMMY WORD FOR FUTURE USE.
C
      IH=IHINDX(N)
      JH=JHINDX(N)
      LMHK     = LMH(IH,JH)
      NWORD2   = 2*LMHK
      NWORD3   = 3*LMHK
      NWORD4   = 4*LMHK
      NWORD5   = 5*LMHK
      NWORD6   = 6*LMHK
      NWORD7   = 7*LMHK
      NWORD8   = 8*LMHK
      NWORD9   = 9*LMHK
      NWORD10  = 10*LMHK
      NWORD11  = 11*LMHK
      NWORD12  = 12*LMHK
      NWORD13  = 13*LMHK
      NWORD14  = 14*LMHK
      NWORD15  = 15*LMHK
      ISTAT    = IDSTN(N)
      CISTAT   = CIDSTN(N)
C
      FPACK(1) = STNLAT(N)/DTR
      FPACK(2) = -STNLON(N)/DTR
!CHANGEDFORNMM      FPACK(2) = STNLON(N)/DTR
      IF(FPACK(2).LT.-180.)FPACK(2)=FPACK(2)+360.
      FPACK(3) = FIS(N)*GI
      FPACK(4) = FLOAT(LMHK)
      FPACK(5) = LCL1ML
      FPACK(6) = LCL1SL
      FPACK(7) = 9+FPACK(5)*FPACK(4)+FPACK(6)
      FPACK(8) = 999.
      FPACK(9) = 999.

C
C***  WIND ROTATION SINES AND COSINES
C

!	write(6,*) 'STNLON, TLM0D, COSPH0, SINPH0: ', N,STNLON(N), 
!     &                                       TLM0D, COSPH0, SINPH0

      DLM    = STNLON(N)+TLM0D*DTR
      XX     = COSPH0*COS(STNLAT(N))*COS(DLM)
     1        +SINPH0*SIN(STNLAT(N))
      YY     = -COS(STNLAT(N))*SIN(DLM)
      TLON   = ATAN(YY/XX)
      ALPHA  = ASIN(SINPH0*SIN(TLON)/COS(STNLAT(N)))
      SINALP = SIN(ALPHA)
      COSALP = COS(ALPHA)

!      SINALP = SROT(N)
!      COSALP = CROT(N)
C
C------------------------------------------------------------------
C***  EXTRACT PRESSURE AND TEMPERATURE PROFILES.
C***  EXTRACT/ROTATE U AND V WIND COMPONENT PROFILES.
C***  EXTRACT SPECIFIC HUMIDITY AND TEMPERATURE TENDENCY.
C***  EXTRACT CLOUD WATER, HEATING DUE TO CONVECTION, LARGE
C***  SCALE RAIN, SHORT WAVE RADIATION, LONG WAVE RADIATION,
C***  AND CLOUD FRACTION.
C------------------------------------------------------------------
C
      DO LV=1,LMHK
        LVL=LMHK-LV+1
!        PRODAT(LVL)      = PDSL1(N)*AETA(LV)+PT
        PRODAT(LVL)      = PMID(N,LV)
        PRODAT(LMHK+LVL) = T(N,LV)

C***  ROTATE WINDS
C
        UT     = U(N,LV)
        VT     = V(N,LV)

        PRODAT(NWORD2+LVL) = UT*COSALP+VT*SINALP
        PRODAT(NWORD3+LVL) = VT*COSALP-UT*SINALP

	if (N .eq. 1) THEN
c	WRITE(6,*) 'orig U,V: ', UT,VT	
c	write(6,*) 'COSALP,SINALP: ', COSALP,SINALP
c	WRITE(6,*) 'rotat U,V: ', PRODAT(NWORD2+LVL),PRODAT(NWORD3+LVL)
c	write(6,*) '-----------------'
	endif

C
        PRODAT(NWORD4+LVL) = Q(N,LV)
C
        IF(RTOP(N,LV).GT.1.E-12) THEN
           PRODAT(NWORD5+LVL) = OMGA(N,LV)
Cmp     1   PRODAT(NWORD5+LVL) = OMGALF(N,LV)*CP/(RTOP(N,LV)*DT)
	ENDIF

!        IF(IW(N,LV).GT.0.5)THEN
!          PRODAT(NWORD6+LVL) = -CWM(N,LV)
!        ELSE
!          PRODAT(NWORD6+LVL) = CWM(N,LV)
!        ENDIF


!--- Separate posting of cloud water (NWORD6), rain (NWORD14), and
!    ice (cloud ice + snow, NWORD13).  The ice and rain fields are
!    *APPENDED* at the end of the section dealing with vertical
!    profile fields  (Ferrier/Manikin 11/30/04)
!
        QLIQ = (1.-F_ice(N,LV))*CWM(N,LV)
        PRODAT(NWORD13+LVL) = F_rain(N,LV)*QLIQ
        PRODAT(NWORD6+LVL) = QLIQ-PRODAT(NWORD13+LVL)
        PRODAT(NWORD14+LVL) = F_ice(N,LV)*CWM(N,LV)
        PRODAT(NWORD7+LVL) = TCUCN(N,LV)
        PRODAT(NWORD8+LVL) = TRAIN(N,LV)
        PRODAT(NWORD9+LVL) = RSWTT(N,LV)
        PRODAT(NWORD10+LVL)= RLWTT(N,LV)
        PRODAT(NWORD11+LVL)= CLDFRA(N,LV)*100.
C
        IF(LV.EQ.1)THEN
          PRODAT(NWORD12+LVL)=Q2(N,LV)
        ELSE
          PRODAT(NWORD12+LVL)=(Q2(N,LV)+Q2(N,LV-1))*0.5
        ENDIF
      ENDDO

C
C***  MODIFY ACCUMLATIONS SO AS TO REPRESENT ACCUMULATED
C***  CHANGE SINCE LAST PROFILE OUTPUT TIME.
C
CGSM  MODIFIED CODE TO ACCOUNT FOR DHCNVC AND DHRAIN BEING
C       COMPUTED FROM TOP DOWN WHILE PRODAT IS FILLED FROM
C       BOTTOM UP 
C
      DO LL=1,LMHK
        LVL=LMHK-LL+1
        STADHC(LL) = PRODAT(NWORD7+LL) - DHCNVC(LVL,N)
        STADHR(LL) = PRODAT(NWORD8+LL) - DHRAIN(LVL,N)
C
        DHCNVC(LVL,N) = PRODAT(NWORD7+LL)
        DHRAIN(LVL,N) = PRODAT(NWORD8+LL)
C
Ctmp        IF(MOD(NTSD,NHEAT).LT.NCNVC)THEN
          DHCNVC(LVL,N) = 0.
          DHRAIN(LVL,N) = 0.
Ctmp        ENDIF
      ENDDO
C
C***  EXTRACT SINGLE LEVEL DATA.   EGRID2 IS SURFACE TEMPERATURE.
C
      PRODAT(NWORD15+1)  = PSLP  (N)
      PRODAT(NWORD15+2)  = PDS   (N)
      PRODAT(NWORD15+3)  = EGRID2(N)
      PRODAT(NWORD15+4)  = TLMIN (N)
      PRODAT(NWORD15+5)  = TLMAX (N)
c     print *,'tlmax, tlmin ',n,tlmax(n),tlmin(n)
      PRODAT(NWORD15+6)  = SMSTAV(N)*100.
      PRODAT(NWORD15+7)  = ACPREC(N)*1000.
      PRODAT(NWORD15+8)  = CUPREC(N)*1000.
      PRODAT(NWORD15+27) = Z0    (N)
C
      STAPRX=PRODAT(NWORD15+7)-STATPR(N)
        if (IDSTN(N) .eq. 724050) then
        write(6,*) 'STAPRX at DCA: ', STAPRX
        endif
      STACRX=PRODAT(NWORD15+8)-STACPR(N)

!	if (STAPRX .gt. 0) then
!	write(6,*) '1hr precip: ',  N,STAPRX
!	endif
C
C***  ROTATE WINDS
C
      UT     = U10(N)
      VT     = V10(N)
      PRODAT(NWORD15+28) = UT*COSALP+VT*SINALP
      PRODAT(NWORD15+29) = VT*COSALP-UT*SINALP
C
      PRODAT(NWORD15+30) = TH10  (N)
! Q10 undefined at F00; use 2-m value
      IF(ITAG .eq. 0) THEN
        PRODAT(NWORD15+31) = QSHLTR(N)
      ELSE
        PRODAT(NWORD15+31) = Q10   (N)
      ENDIF
      PRODAT(NWORD15+32) = TSHLTR(N)
      PRODAT(NWORD15+33) = QSHLTR(N)
      PRODAT(NWORD15+34) = SFCEXC(N)
      PRODAT(NWORD15+35) = VEGFRC(N)*100.
      PRODAT(NWORD15+36) = CMC   (N)*1000.
      PRODAT(NWORD15+37) = SMC   (N,1)
      PRODAT(NWORD15+38) = SMC   (N,2)
      PRODAT(NWORD15+39) = SMC   (N,3)
      PRODAT(NWORD15+40) = SMC   (N,4)
      PRODAT(NWORD15+41) = STC   (N,1)
      PRODAT(NWORD15+42) = STC   (N,2)
      PRODAT(NWORD15+43) = STC   (N,3)
      PRODAT(NWORD15+44) = STC   (N,4)
      PRODAT(NWORD15+45) = SM    (N) + SICE(N)
      PRODAT(NWORD15+46) = CFRACL(N)*100.
      PRODAT(NWORD15+47) = CFRACM(N)*100.
      PRODAT(NWORD15+48) = CFRACH(N)*100.
      PRODAT(NWORD15+49) = SR    (N)*100.
      PRODAT(NWORD15+50) = NINT(HBOT(N))
      PRODAT(NWORD15+51) = CPRATE(N)
      PRODAT(NWORD15+52) = F_RIMEF(N,LMHK)

C
      PRODAT(NWORD15+9)   = SFCLHX(N)
      PRODAT(NWORD15+10)  = POTFLX(N)
      PRODAT(NWORD15+11)  = SFCSHX(N)
      PRODAT(NWORD15+12)  = SUBSHX(N)
      PRODAT(NWORD15+13)  = SNOPCX(N)
      PRODAT(NWORD15+14)  = ASWIN (N)
      PRODAT(NWORD15+15)  = ASWOUT(N)
      PRODAT(NWORD15+16)  = ALWIN (N)
      PRODAT(NWORD15+17)  = ALWOUT(N)
      PRODAT(NWORD15+18)  =-ALWTOA(N)
      PRODAT(NWORD15+19)  = ASWTOA(N)
      PRODAT(NWORD15+20)  = ACSNOW(N)
      PRODAT(NWORD15+21)  = SMSTOT(N)
      PRODAT(NWORD15+22)  = SNO   (N)
      PRODAT(NWORD15+23)  = ACSNOM(N)
      PRODAT(NWORD15+24)  = SSROFF(N)
      PRODAT(NWORD15+25)  = BGROFF(N)
      PRODAT(NWORD15+26)  = SOILTB(N)
C
C***  ACCUMULATED CHANGE SINCE LAST PROFILE OUTPUT TIME.
C
      PSFCEVP  = PRODAT(NWORD15+9 ) - STAEVP(N)
      PPOTEVP  = PRODAT(NWORD15+10) - STAPOT(N)
      PSFCSHX  = PRODAT(NWORD15+11) - STASHX(N)
      PSFCSUB  = PRODAT(NWORD15+12) - STASUB(N)
      PSNOPCX  = PRODAT(NWORD15+13) - STAPCX(N)
      PRSWIN   = PRODAT(NWORD15+14) - STASWI(N)
	if (N .eq. 90) then
	write(6,*) 'N, PRODAT(NWORD15+14) , STASWI(N), PRSWIN: ', N,
     &     PRODAT(NWORD15+14), STASWI(N), PRSWIN
	endif
      PRSWOUT  = PRODAT(NWORD15+15) - STASWO(N)
      PRLWIN   = PRODAT(NWORD15+16) - STALWI(N)
	if (N .eq. 90) then
	write(6,*) 'N, PRODAT(NWORD53+16) , STALWI(N), PRLWIN: ', N,
     &     PRODAT(NWORD15+16), STALWI(N), PRLWIN
	endif
      PRLWOUT  = PRODAT(NWORD15+17) - STALWO(N)
      PRLWTOA  = PRODAT(NWORD15+18) - STALWT(N)
      PRSWTOA  = PRODAT(NWORD15+19) - STASWT(N)
      PACSNOW  = PRODAT(NWORD15+20) - STASNO(N)
      PACSNOM  = PRODAT(NWORD15+23) - STASNM(N)
      PSSROFF  = PRODAT(NWORD15+24) - STASRF(N)
      PBGROFF  = PRODAT(NWORD15+25) - STABRF(N)
	write(6,*) 'N, PBGROFF: ', N, PBGROFF
C***
C***  TRANSFER STATION PROFILE DATA TO "PACKED" OUTPUT ARRAY.
C***
      NN   = 0
      NLEN = FPACK(7)
        write(6,*) 'SHELT TEMP ', PRODAT(NWORD15+32) 
C
      DO NL = 10,NLEN
        NN = NL-9
        FPACK(NL) = PRODAT(NN)
      ENDDO
C
C***  REPLACE ACCUMULATED QUANTITIES WITH ACCUMULATION
C***  SINCE LAST PROFILE OUTPUT TIME.
C
      DO LL = 1,LMHK
         FPACK(9+NWORD7+LL) = STADHC(LL)*RTSCU
         FPACK(9+NWORD8+LL) = STADHR(LL)*RTSRA
!       FPACK(9+NWORD7+LL) = -9999.
!       FPACK(9+NWORD8+LL) = -9999.
      ENDDO
C
      FPACK(9+NWORD15+7)  = STAPRX
!	write(6,*) 'precip written to FPACK element: ', 9+NWORD13+7
      FPACK(9+NWORD15+8)  = STACRX
      FPACK(9+NWORD15+9)  = PSFCEVP * RTSPH
      FPACK(9+NWORD15+10) = PPOTEVP * RTSPH

!	if (IDSTN(N) .eq. 725840) then
!	write(6,*) 'PPOTEVP, RTSPH: ', PPOTEVP, RTSPH
!	write(6,*) 'packing PPOTEVP*RTSPH= ', FPACK(9+NWORD13+10)
!	endif

      FPACK(9+NWORD15+11) = PSFCSHX * RTSPH
      FPACK(9+NWORD15+12) =-PSFCSUB * RTSPH
      FPACK(9+NWORD15+13) = PSNOPCX * RTSPH
      FPACK(9+NWORD15+14) = PRSWIN  * RTSPH
!	if (mod(N,NUMSTA/5) .eq. 0) then
!	write(6,*) 'N, RTSPH, SWRD: ', N, RTSPH, FPACK(9+NWORD13+14)
!	endif
      FPACK(9+NWORD15+15) = PRSWOUT * RTSPH
      FPACK(9+NWORD15+16) = PRLWIN  * RTSPH
      FPACK(9+NWORD15+17) = PRLWOUT * RTSPH
      FPACK(9+NWORD15+18) = PRLWTOA * RTSPH
      FPACK(9+NWORD15+19) = PRSWTOA * RTSPH
      FPACK(9+NWORD15+20) = PACSNOW
      FPACK(9+NWORD15+23) = PACSNOM
      FPACK(9+NWORD15+24) = PSSROFF
      FPACK(9+NWORD15+25) = PBGROFF
C
!      IF(RESTRT)THEN
      IF(ITAG .eq. 0)THEN
        DO LL = 1,LMHK
          FPACK(9+NWORD7+LL) = 0.
          FPACK(9+NWORD8+LL) = 0.
        ENDDO
C
        FPACK(9+NWORD15+7)  = 0.
        FPACK(9+NWORD15+8)  = 0.
        FPACK(9+NWORD15+9)  = 0.
        FPACK(9+NWORD15+10) = 0.
        FPACK(9+NWORD15+11) = 0.
        FPACK(9+NWORD15+12) = 0.
        FPACK(9+NWORD15+13) = 0.
        FPACK(9+NWORD15+14) = 0.
        FPACK(9+NWORD15+15) = 0.
        FPACK(9+NWORD15+16) = 0.
        FPACK(9+NWORD15+17) = 0.
        FPACK(9+NWORD15+18) = 0.
        FPACK(9+NWORD15+19) = 0.
        FPACK(9+NWORD15+20) = 0.
        FPACK(9+NWORD15+23) = 0.
        FPACK(9+NWORD15+24) = 0.
        FPACK(9+NWORD15+25) = 0.
      ENDIF
C---------------------------------------------------------------------
C***
C***  WRITE PROFILE DATA
C***
      
!      NREC=(IFHR/INCR)*NUMSTA+N
      NREC=N
      write(6,*) 'NREC,IHRST,IDAT: ', NREC,IHRST,IDAT,IFCST

        if (mod(NREC,5) .eq. 0) then
        write(6,*) 'IDAT, IFCST, ISTAT, CISTAT: ',
     &  IDAT, IFCST, ISTAT, CISTAT
        endif

      WRITE(LCLAS1,REC=NREC)IHRST,IDAT,IFCST,ISTAT,CISTAT
     1,                    (FPACK(NL),NL=1,NLEN)
	write(6,*) 'FPACK(458) written: ', FPACK(458)

C---------------------------------------------------------------------
 1000 CONTINUE
        write(6,*) 'past 1000 continue'
      CLOSE(LCLAS1)
C
        write(6,*) 'past file close'
C***  END OF PROFILE SITE LOOP
C
C***  END PROFILE POSTING CODE.
C---------------------------------------------------------------------
        write(6,*) 'to deallocate'
        DEALLOCATE( DETA,RDETA,AETA,UL
     &,RES,FIS,THS,HBOT
     &,CFRACL,CFRACM,CFRACH,SNO
     &,SOILTB,SFCEXC,SMSTAV,SMSTOT)
        write(6,*) 'past 1st dealloc'
        DEALLOCATE(Z0,CZEN,CZMEAN,U00,SR,CPRATE
     &,ACPREC,CUPREC,ACSNOW,ACSNOM
     &,SSROFF,BGROFF,SFCSHX,SFCLHX
     &,SUBSHX,SNOPCX,ASWIN,ASWOUT
     &,ASWTOA,ALWIN,ALWOUT,ALWTOA
     &,TSHLTR,TSHLTR_hold,QSHLTR,PSHLTR)
        write(6,*) 'past 2nd dealloc'

      DEALLOCATE(TH10,Q10,U10,V10
     &,TLMIN,TLMAX
     &,SMC,CMC,STC,SH2O
     &,VEGFRC,POTFLX,PSLP,PDSL1
     &,EGRID2,SM,SICE
     &,HBM2,FACTR)
        write(6,*) 'past 3rd dealloc'
      DEALLOCATE(PTBL,TTBL,VEGFRA
     &,T,Q,U,V,Q2
     &,OMGALF,CWM,TRAIN,TCUCN)

        write(6,*) 'past 4th dealloc'
      DEALLOCATE(RSWTT,RLWTT,CCR,RTOP
     &,HTM,OMGA
     &,PRODAT,FPACK)
        write(6,*) 'past 5th dealloc'

 
      DEALLOCATE(STATPR,STACPR,STAEVP
     &,STAPOT,STASHX,STASUB,STAPCX
     &,STASWI,STASWO,STALWI,STALWO
     &,STALWT,STASWT,STASNM,STASRF
     &,STABRF,STASNO,DHCNVC,DHRAIN
     &,STADHC,STADHR)
        write(6,*) 'past 6th dealloc'

       DEALLOCATE(ACPREC0,CUPREC0,SFCLHX0,POTFLX0
     &,SFCSHX0,SUBSHX0,SNOPCX0,ASWIN0)

        write(6,*) 'past 7th dealloc'
       DEALLOCATE(ASWOUT0,ALWIN0,ALWOUT0,ALWTOA0)
        write(6,*) 'past 8th dealloc'
       DEALLOCATE(ASWTOA0,ACSNOW0,ACSNOM0,SSROFF0)
        write(6,*) 'past 9th dealloc'
        write(6,*) 'to dealloc of BGROFF0'
        write(6,*) 'min/max of BGROFF0: ', minval(BGROFF0),
     &                                     maxval(BGROFF0)
        DEALLOCATE(BGROFF0)
         write(6,*) 'past 10th dealloc'
        write(6,*) 'past deallocate, to return'

        write(6,*) 'End of PROF_NMM_NET'

      RETURN
      END
