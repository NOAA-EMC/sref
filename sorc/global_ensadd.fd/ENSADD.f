C$$$  MAIN PROGRAM DOCUMENTATION BLOCK
C
C MAIN PROGRAM: ENSADD
C   PRGMMR: ZHU              ORG: NP23        DATE: 1999-08-31
C
C ABSTRACT: THIS PROGRAM WILL EXTEND PDS MESSAGE WHICH WILL  
C           INCLUDE ENS(5) MESSAGE
C
C PROGRAM HISTORY LOG:
C   96-10-??   MARK IREDELL     - Originator
C   97-03-17   YUEJIAN ZHU      - Added DOCBLOACK
C   99-07-26   YUEJIAN ZHU      - Modified to IBM-SP
C   15-02-24   JUN DU           - Expanded and adopted for SREF
C                                 as well as added a precision control
C
C USAGE:
C
C   INPUT FILES:
C     UNIT  11  GRIB FILE ( maximum 512*256 )
C     UNIT  21  GRIB INDEX FILE
C     UNIT   5  READ *, IENST,IENSI (For ensemble message)
C
C   OUTPUT FILES:
C     UNIT  51  GRIB FILE ( maximum 512*256 )
C
C   SUBPROGRAMS CALLED:
C     GETGB  -- W3LIB ROUTINE
C     PUTGBE -- W3LIB ROUTINE
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN
C
C$$$
      program      ensadd
C     parameter    (jf=512*256)
      parameter    (jf=697*553)
c     integer      jpds(200),jgds(200),kpds(200),kgds(200),kens(200)
      integer      jpds(25),jgds(22),kpds(25),kgds(22),kens(20)
      integer      ienst,iensi                                        
      real         f(jf)
      character*80 cpgb,cpgi,cpge
      logical*1      lb(jf)
      namelist     /namin/ ienst,iensi,cpgb,cpgi,cpge      
c
c
c     CALL W3LOG('$S97118.73','ENSADD')
      CALL W3TAGB('ENSADD',1999,0243,0068,'NP23')

      read (5,namin)
      lpgb=len_trim(cpgb)
      lpgi=len_trim(cpgi)
      lpge=len_trim(cpge)
      print *, cpgb(1:lpgb),'  ',cpgi(1:lpgi),'  ',cpge(1:lpge)
      call baopenr(11,cpgb(1:lpgb),iretb)
      call baopenr(21,cpgi(1:lpgi),ireti)
      call baopen (51,cpge(1:lpge),irete)

      j=0
      jpds=-1
      call getgb(11,21,jf,j,jpds,jgds,kf,j,kpds,kgds,lb,f,iret)
      dowhile(iret.eq.0)
        kpds(23)=2   !ensemble
C Precision control (J.Du)
        kpds(22)=2    !default

c H, mslp, slp, pressure
        if(kpds(5).eq.7.or.kpds(5).eq.2) kpds(22)=0
        if(kpds(5).eq.130.or.kpds(5).eq.1) kpds(22)=0
c cetegorical precip types
        if(kpds(5).eq.140.or.kpds(5).eq.141) kpds(22)=0
        if(kpds(5).eq.142.or.kpds(5).eq.143) kpds(22)=0
c visibility
        if(kpds(5).eq.20) kpds(22)=0
c categorical icing, turbulence
        if(kpds(5).eq.186.or.kpds(5).eq.185) kpds(22)=0
c land-sea mask
        if(kpds(5).eq.81) kpds(22)=0
c sea-ice mask
        if(kpds(5).eq.91) kpds(22)=0

c radar reflectivity and echo top
        if(kpds(5).eq.212.or.kpds(5).eq.240) kpds(22)=1
c u, v
        if(kpds(5).eq.33.or.kpds(5).eq.34) kpds(22)=1
c cape, cin
        if(kpds(5).eq.157.or.kpds(5).eq.156) kpds(22)=1
c storm relative helicity
        if(kpds(5).eq.190) kpds(22)=1
c total cloud cover
        if(kpds(5).eq.71) kpds(22)=1
c PBL height
        if(kpds(5).eq.221) kpds(22)=1
c relative humidity and temperature
        if(kpds(5).eq.52.or.kpds(5).eq.11) kpds(22)=1
c dew point temperature
        if(kpds(5).eq.17) kpds(22)=1

c TKE
        if(kpds(5).eq.158) kpds(22)=3
c soil moisture fraction
        if(kpds(5).eq.144) kpds(22)=4
c precip rate
        if(kpds(5).eq.59) kpds(22)=5
c specific humidity
        if(kpds(5).eq.51) kpds(22)=5
c vorticity
        if(kpds(5).eq.41) kpds(22)=6
c moisture divergence
        if(kpds(5).eq.135) kpds(22)=6
c cloud ice and water
        if(kpds(5).eq.58.or.kpds(5).eq.153) kpds(22)=6

c to avoid zero in relative humidity and specific humidity
        if(kpds(5).eq.52) f=f+0.1
        if(kpds(5).eq.51) f=f+0.00001

        kens(1)=1
        kens(2)=ienst
        kens(3)=iensi
        kens(4)=1
        kens(5)=255
        if(kpds(6).eq.100.and.kpds(5).ne.222.and.kpds(5).ne.52)
     &  kens(5)=kgds(3)/2-1
        call putgbe(51,kf,kpds,kgds,kens,lb,f,iret)
        call getgb(11,21,jf,j,jpds,jgds,kf,j,kpds,kgds,lb,f,iret)
      enddo
 
      call baclose(11,iretb)
      call baclose(21,ireti)
      call baclose(51,irete)

c     CALL W3LOG('$E')
      CALL W3TAGE('ENSADD')

      stop
      end
