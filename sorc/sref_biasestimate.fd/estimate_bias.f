      program biasestimate

C$$  MAIN PROGRAM DOCUMENTATION BLOCK
C      
C MAIN PROGRAM: bias
C PRGMMR: Jun Du          ORG: NP21       DATE: Nov. 15, 2004
C
C ABSTRACT: 
c  1. decode grib data 
c  2. calculate ensemble mean and spread
c  3. estimate bias with various methods for all SREF members
c     method=1: simple average
c     method=2: decaying average
c     method=3: regime-selective average
c  4. write out ensemble mean, spread and bias for next day use
c  5. produce new debiased ensemble forecasts
c
c PROGRAM HISTORY LOG:
c  11/05/2004, Jun Du, initial program
c  11/23/2004, Jun Du, added multi-method functionality
c  03/03/2006, Jun Du, add normalized spread (by mean spread in past) 
c                      as an indicator of forecast confidence
c  05/20/2007, Jun Du, modified for operational implementation
c  07/28/2008, Jun Du, increased wrf membership from 6 to 10 and reduced
c                      eta membership from 10 to 6
c  11/03/2011, Jun Du, modified to the new 2012 NEMS-based version of SREF
c                      (v6.0.0) by reducing grouping from 4 (Eta, RSM, NMM 
c                      and ARW) to 3 (NMMB, NMM and ARW)
c  03/10/2015, Jun Du, modified to the new 2015 2-model 26-member SREFv7.0.0
c  03/25/2015, Jun Du, added ensemble section
c
c INPUT:
c
c OUTPUT:
c
c SUBPROGRAMS CALLED: 
c  getname
c  grange
c  mean_spread
c  correlation
c  bias_averaging
c
C ATTRIBUTES:
c   language: Fortran 77
c   machine: IBM-SP
C
C$$$

C OUTPUT:
C KPDS(25) INTEGER*4
C ARRAY TO CONTAIN THE ELEMENTS OF THE PRODUCT DEFINITION SEC
C (VERSION 1)
C KPDS(1)  - ID OF CENTER
C KPDS(2)  - MODEL IDENTIFICATION (SEE "GRIB" TABLE 1)
C KPDS(3)  - GRID IDENTIFICATION (SEE "GRIB" TABLE 2)
C KPDS(4)  - GDS/BMS FLAG
C            BIT       DEFINITION
C            25        0 - GDS OMITTED
C                      1 - GDS INCLUDED
C            26        0 - BMS OMITTED
C                      1 - BMS INCLUDED
C            NOTE:- LEFTMOST BIT = 1,
C                   RIGHTMOST BIT = 32
C KPDS(5)  - INDICATOR OF PARAMETER (SEE "GRIB" TABLE 5)
C KPDS(6)  - TYPE OF LEVEL (SEE "GRIB" TABLES 6 & 7)
C KPDS(7)  - HEIGHT,PRESSURE,ETC  OF LEVEL
C KPDS(8)  - YEAR INCLUDING CENTURY
C KPDS(9)  - MONTH OF YEAR
C KPDS(10) - DAY OF MONTH
C KPDS(11) - HOUR OF DAY
C KPDS(12) - MINUTE OF HOUR
C KPDS(13) - INDICATOR OF FORECAST TIME UNIT (SEE "GRIB"
C            TABLE 8)
C KPDS(14) - TIME 1               (SEE "GRIB" TABLE 8A)
C KPDS(15) - TIME 2               (SEE "GRIB" TABLE 8A)
C KPDS(16) - TIME RANGE INDICATOR (SEE "GRIB" TABLE 8A)
C KPDS(17) - NUMBER INCLUDED IN AVERAGE
C KPDS(18) - EDITION NR OF GRIB SPECIFICATION
C KPDS(19) - VERSION NR OF PARAMETER TABLE
C
C KPDS(20)  - NR MISSING FROM AVERAGE/ACCUMULATION
C KPDS(21)  - CENTURY OF REFERENCE TIME OF DATA
C KPDS(22)  - UNITS DECIMAL SCALE FACTOR
C KPDS(23)  - SUBCENTER NUMBER
C KPDS(24)  - PDS BYTE 29, FOR NMC ENSEMBLE PRODUCTS
C         128 IF FORECAST FIELD ERROR
C          64 IF BIAS CORRECTED FCST FIELD
C          32 IF SMOOTHED FIELD
C         WARNING: CAN BE COMBINATION OF MORE THAN 1
C KPDS(25)  - PDS BYTE 30, NOT USED
C KPDS(26-35)  - RESERVED
C KPDS(36-N)   - CONSECUTIVE BYTES EXTRACTED FROM PROGRAM
C                DEFINITION SECTION (PDS) OF GRIB MESSAGE
c     kpds(1)=7
c     kpds(2)=130
c     kpds(3)=212
c     kpds(4)=0
c     kpds(5)=variable dependent
c     kpds(6)=variable dependent
c     kpds(7)=variable dependent
c     kpds(8)=iyr-(iyr/100)*100
c     kpds(9)=imon
c     kpds(10)=idy
c     kpds(11)=ihr
c     kpds(12)=0
c     kpds(13)=1
c     kpds(14)=variable dependent
c     kpds(15)=variable dependent
c     kpds(16)=variable dependent
c     kpds(17)=-1
c     kpds(18)=1
c     kpds(19)=2
c     kpds(20)=-1
c     kpds(21)=21
c     kpds(22)=variable dependent (units decimal scale factor/precision)
c     kpds(23)=2
c     kpds(24)=128
c     kpds(25)=-1
       
c     kgds(20)=255
CCCCCCCCCCCCCCCCCCCC

 
C nfile - 1-21: ens mems
      parameter(mem=26,nfile=mem,ntime=29,nvmax=150)
      parameter(nSREF=2,npastday=20) !NMB, NMM, ARW
      parameter(jf=185*129,lon=185,lat=129,ngrid=lon*lat)  !for Grid 212
c     parameter(jf=512*256,lon=185,lat=129,ngrid=lon*lat)

      integer mem_sta(nSREF),mem_end(nSREF)
      dimension f(jf),data(ngrid,nfile,nvmax),x(ngrid),y(ngrid),
     &spread(ngrid),avesprd(ngrid),sprd_norm(ngrid),fcst(nfile),
     &xmean(ngrid),corr(npastday),fmean(ngrid,nSREF,nvmax),
     &pasterro(ngrid,npastday),pastfcst(ngrid,npastday),
     &pastsprd(ngrid,npastday),pastbias(ngrid),bias(ngrid),w(npastday),
     &fcst_BC(ngrid),weight(npastday)
      dimension jpds5(nvmax),jpds6(nvmax),jpds7(nvmax),kkpds(25,nvmax)
      dimension kftime(nfile,ntime),jpds(25),jgds(22),kpds(25),kgds(22)
      dimension kens(20)
      integer yy1,mm1,dd1,hh1,method,type
      integer pyy(npastday),pmm(npastday),pdd(npastday)
      integer yesyy,yesmm,yesdd
      logical lb(jf)
      character*2  date
      character*2  sys
      character*1  meth
      character*13  fname1
      character*15  fname2
      character*15 file1,file3,file5,file7
      character*17 file2,file4,file6,file8
      character*16 bias_today
      character*22 fcst_new1,fcst_new2,fcst_new3,fcst_new4,fcst_new5,
     &fcst_new6,fcst_new7,fcst_new8,fcst_new9,fcst_new10,fcst_new11,
     &fcst_new12,fcst_new13,fcst_new14,fcst_new15,fcst_new16,fcst_new17,
     &fcst_new18,fcst_new19,fcst_new20,fcst_new21,fcst_new22,fcst_new23,
     &fcst_new24,fcst_new25,fcst_new26
      namelist/namin/yy1,mm1,dd1,hh1,method,type,nvar0,
     &nvar
      data (mem_sta(i),i=1,nSREF) /01,14/
      data (mem_end(i),i=1,nSREF) /13,26/

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C decoding current-day forecast grib data
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c Passing over date information
      read(5,namin,end=1000)
1000  continue
      print*,yy1,mm1,dd1,hh1
      print*,method,type,nvar0,nvar
      read(7,*) ((kftime(i,j),j=1,ntime),i=1,nfile)
      do i=1,nfile
       print *,(kftime(i,j),j=1,ntime)
      enddo

      read(8,*) (pyy(i),i=1,npastday),(pmm(i),i=1,npastday),
     &(pdd(i),i=1,npastday)
      print*,(pyy(i),i=1,npastday)
      print*,(pmm(i),i=1,npastday)
      print*,(pdd(i),i=1,npastday)

      read(9,*) yesyy,yesmm,yesdd
      print*,'yesyy,yesmm,yesdd=',yesyy,yesmm,yesdd

      write(meth,'(i1.1)') method

       nunit=101
       open (nunit, file='variable.tbl', status='old')
       read (nunit,*)
       read (nunit,*)
102    format(9x,3i8)
      do nv=nvar0,nvar
        read(nunit,102) jpds5(nv),jpds6(nv),jpds7(nv)
        print*,jpds5(nv), ' ',jpds6(nv),' ',jpds7(nv)
      enddo

      do 10 nt=1,ntime

        print *,' '
        print *,'----------------------------------------------- '
        print *,'NEW TIME: nt = ',nt
        print *,' '

        j=0 
        data=0.0

        do 5 nf=1,nfile
        lugb=8       
        lugi=9     
        jpds=-1
        jgds=-1
C specify date information here:
       if(nf.le.mem) then
c       jpds(8)=yy1-(yy1/100)*100
c       jpds(8)=yy1
c       jpds(9)=mm1
        jpds(10)=dd1
        jpds(11)=hh1
       endif

        print *,' '
        print *,'----------'
        print *,'New File: nt = ',nt,' nf= ',nf
        print *,' '

c define variables to be decoded
c     do 111 nv=nvar0,nvar
      do 111 nv=nvar0+1,nvar   !skip a duable variable T2m in the tbl
       jpds(5)=jpds5(nv)
       jpds(6)=jpds6(nv)
       jpds(7)=jpds7(nv)
       print*,jpds(5), ' ',jpds(6),' ',jpds(7)

       jpds(14)=kftime(nf,nt)
       jpds(15)=-1
c      print*, jpds(10)
c      print*, jpds(15)

      lugb=lugb+1
      lugi=lugi+1
      print *,'lugb= ',lugb,' lugi= ',lugi
      call getname(nf,fname1,fname2)
c     print *,' fname1= ',fname1
c     print *,' fname2= ',fname2
      call baopenr(lugb,fname1,ierr)
      call baopenr(lugi,fname2,ierr)
      call getgb(lugb,lugi,jf,j,jpds,jgds,
     &             kf,k,kpds,kgds,lb,f,iret)
      call baclose(lugb,ierr)
      call baclose(lugi,ierr)
c     print *,'lugb= ',lugb,' lugi= ',lugi
 
      call grange(kf,lb,f,dmin,dmax) 
c     print *,'immediately after call to getgb, iret= ',iret
c    &,' j= ',j,' k= ',k,' nf= ',nf,' nt= ',nt
c    &,' j5= ',jpds(5),' j6= ',jpds(6),' j7= ',jpds(7)
c    &,' j8= ',jpds(8),' j9= ',jpds(9),' j10= ',jpds(10)
c    &,' j11= ',jpds(11),' j14= ',jpds(14),' j15= ',jpds(15)
c     print '(i4,2x,9i5,i8,2g12.4)',
c    &k,(kpds(i),i=5,11),kpds(14),kpds(15),kf,dmin,dmax

      if (iret.eq.0) then
       do ipt=1,ngrid
        data(ipt,nf,nv) = f(ipt)
        if(nt.eq.1) then
         do i=1,25
          kkpds(i,nv)=kpds(i)
         enddo
        endif
       enddo
      else
       print *, "Something wrong in decoding ensemble forecasts!!!!"
      endif
111   continue    !diff variables (nv)
5     continue    !files (nf: mems)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C bias estimation
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        write(date,'(i2.2)') (nt-1)*3+3

         fcst_new1='BC' // meth // '_nmb.pgrb212.c1.f' // date
         fcst_new2='BC' // meth // '_nmb.pgrb212.n1.f' // date
         fcst_new3='BC' // meth // '_nmb.pgrb212.p1.f' // date
         fcst_new4='BC' // meth // '_nmb.pgrb212.n2.f' // date
         fcst_new5='BC' // meth // '_nmb.pgrb212.p2.f' // date
         fcst_new6='BC' // meth // '_nmb.pgrb212.n3.f' // date
         fcst_new7='BC' // meth // '_nmb.pgrb212.p3.f' // date
         fcst_new8='BC' // meth // '_nmb.pgrb212.n4.f' // date
         fcst_new9='BC' // meth // '_nmb.pgrb212.p4.f' // date
        fcst_new10='BC' // meth // '_nmb.pgrb212.n5.f' // date
        fcst_new11='BC' // meth // '_nmb.pgrb212.p5.f' // date
        fcst_new12='BC' // meth // '_nmb.pgrb212.n6.f' // date
        fcst_new13='BC' // meth // '_nmb.pgrb212.p6.f' // date
        fcst_new14='BC' // meth // '_arw.pgrb212.c1.f' // date
        fcst_new15='BC' // meth // '_arw.pgrb212.n1.f' // date
        fcst_new16='BC' // meth // '_arw.pgrb212.p1.f' // date
        fcst_new17='BC' // meth // '_arw.pgrb212.n2.f' // date
        fcst_new18='BC' // meth // '_arw.pgrb212.p2.f' // date
        fcst_new19='BC' // meth // '_arw.pgrb212.n3.f' // date
        fcst_new20='BC' // meth // '_arw.pgrb212.p3.f' // date
        fcst_new21='BC' // meth // '_arw.pgrb212.n4.f' // date
        fcst_new22='BC' // meth // '_arw.pgrb212.p4.f' // date
        fcst_new23='BC' // meth // '_arw.pgrb212.n5.f' // date
        fcst_new24='BC' // meth // '_arw.pgrb212.p5.f' // date
        fcst_new25='BC' // meth // '_arw.pgrb212.n6.f' // date
        fcst_new26='BC' // meth // '_arw.pgrb212.p6.f' // date
        call baopen(10,fcst_new1,ierr)
        call baopen(11,fcst_new2,ierr)
        call baopen(12,fcst_new3,ierr)
        call baopen(13,fcst_new4,ierr)
        call baopen(14,fcst_new5,ierr)
        call baopen(15,fcst_new6,ierr)
        call baopen(16,fcst_new7,ierr)
        call baopen(17,fcst_new8,ierr)
        call baopen(18,fcst_new9,ierr)
        call baopen(19,fcst_new10,ierr)
        call baopen(20,fcst_new11,ierr)
        call baopen(21,fcst_new12,ierr)
        call baopen(22,fcst_new13,ierr)
        call baopen(23,fcst_new14,ierr)
        call baopen(24,fcst_new15,ierr)
        call baopen(25,fcst_new16,ierr)
        call baopen(26,fcst_new17,ierr)
        call baopen(27,fcst_new18,ierr)
        call baopen(28,fcst_new19,ierr)
        call baopen(29,fcst_new20,ierr)
        call baopen(30,fcst_new21,ierr)
        call baopen(31,fcst_new22,ierr)
        call baopen(32,fcst_new23,ierr)
        call baopen(33,fcst_new24,ierr)
        call baopen(34,fcst_new25,ierr)
        call baopen(35,fcst_new26,ierr)

       do nens=1,nSREF  !how many sub-SREFs
         write(sys,'(i2.2)') nens
         file1='grb_pasterro.' // sys
         file2='grb_pasterro.' // sys // '.i'
         file3='grb_pastfcst.' // sys
         file4='grb_pastfcst.' // sys // '.i'
         file5='grb_pastsprd.' // sys
         file6='grb_pastsprd.' // sys // '.i'
         file7='grb_pastbias.' // sys
         file8='grb_pastbias.' // sys // '.i'

         bias_today='grb_bias' // meth // '.' // sys // '.f' // date
         lugb_newbc=70+nens
         call baopen(lugb_newbc,bias_today,ierr)

c     open(97,file='asc_bias' // sys // '.f' // date,
c    &status='unknown')
      if (method.eq.3) then
       open(98,file='asc_wgts' // sys // '.f' // date,
     & status='unknown')
      endif

C if the number of subgroup increases, these unit values need to be modifed.
          lugb_err=40+nens
          lugi_err=41+nens
          call baopenr(lugb_err,file1,ierr)
          call baopenr(lugi_err,file2,ierr)
          lugb_fcst=45+nens
          lugi_fcst=46+nens
          call baopenr(lugb_fcst,file3,ierr)
          call baopenr(lugi_fcst,file4,ierr)
          lugb_sprd=50+nens
          lugi_sprd=51+nens
          call baopenr(lugb_sprd,file5,ierr)
          call baopenr(lugi_sprd,file6,ierr)
          lugb_bias=55+nens
          lugi_bias=56+nens
          call baopenr(lugb_bias,file7,ierr)
          call baopenr(lugi_bias,file8,ierr)

        do nv=nvar0+1,nvar  !for selected variables
         jpds=-1
         jgds=-1
         jpds(5)=jpds5(nv)
         jpds(6)=jpds6(nv)
         jpds(7)=jpds7(nv)
         jpds(14)=kftime(1,nt)   !same as ensemble member 1
         jpds(15)=-1
         jpds(21)=21

         do npast=1,npastday
          jpds(8)=pyy(npast)-(pyy(npast)/100)*100
          jpds(9)=pmm(npast)
          jpds(10)=pdd(npast)
          jpds(11)=hh1
C------------error--------------------------
          if(method.eq.2.and.npast.gt.1) goto 888
          call getgb(lugb_err,lugi_err,jf,j,jpds,jgds,kf,k,kpds,
     &               kgds,lb,f,iret)
          call grange(kf,lb,f,dmin,dmax) 
          if (nens.eq.1.and.npast.eq.1) then
          print *,' lugb_err= ',lugb_err
          print *,' lugi_err= ',lugi_err
c         print *,' file1= ',file1
c         print *,' file2= ',file2
          print *,'read pasterr, iret=, nt= ',iret,nt,dmin,dmax
          endif
c         print *,'read pasterr, iret=, nt= ',iret,nt,dmin,dmax
c         print *,'reading pastday error, iret= ',iret
c    &,' j= ',j,' k= ',k,' nsub= ',nens,' nt= ',nt
c    &,' j5= ',jpds(5),' j6= ',jpds(6),' j7= ',jpds(7)
c    &,' j8= ',jpds(8),' j9= ',jpds(9),' j10= ',jpds(10)
c    &,' j11= ',jpds(11),' j14= ',jpds(14),' j15= ',jpds(15)
c         print '(i4,2x,9i5,i8,2g12.4)',
c    &k,(kpds(i),i=5,11),kpds(14),kpds(15),kf,dmin,dmax
          if (iret.eq.0) then
           do ipt=1,ngrid
            pasterro(ipt,npast) = f(ipt)
           enddo
          else
           print *,"Something wrong in reading pastday error!"
          endif
888       continue
c----------------fcst---------------------
          if(method.eq.3) then
          call getgb(lugb_fcst,lugi_fcst,jf,j,jpds,jgds,kf,k,kpds,
     &               kgds,lb,f,iret)
          call grange(kf,lb,f,dmin,dmax) 
          if (nens.eq.1.and.npast.eq.1) then
c         print *,' file3= ',file3
c         print *,' file4= ',file4
          print *,'read pastfcst, iret=, nt= ',iret,nt,dmin,dmax
          endif
c         print *,'read pastfcst, iret=, nt= ',iret,nt,dmin,dmax
          if (iret.eq.0) then
           do ipt=1,ngrid
            pastfcst(ipt,npast) = f(ipt)
           enddo
          else
           print *,"Something wrong in reading pastday fcst!"
          endif
          endif
c----------------sprd---------------------
          if(method.eq.4) then
          call getgb(lugb_sprd,lugi_sprd,jf,j,jpds,jgds,kf,k,kpds,
     &               kgds,lb,f,iret)
          call grange(kf,lb,f,dmin,dmax) 
          if (nens.eq.1.and.npast.eq.1) then
c         print *,' file5= ',file5
c         print *,' file6= ',file6
          print *,'read pastsprd, iret=, nt= ',iret,nt,dmin,dmax
          endif
c         print *,'read pastsprd, iret=, nt= ',iret,nt,dmin,dmax
          if (iret.eq.0) then
           do ipt=1,ngrid
            pastsprd(ipt,npast) = f(ipt)
           enddo
          else
           print *,"Something wrong in reading pastday spread!"
          endif
          endif
c----------------bias---------------------
          if(type.eq.1.and.method.eq.2.and.npast.eq.1) then
           jpds(8)=yesyy-yesyy/100*100
           jpds(9)=yesmm
           jpds(10)=yesdd
           call getgb(lugb_bias,lugi_bias,jf,j,jpds,jgds,kf,k,kpds,
     &                kgds,lb,f,iret)
           call grange(kf,lb,f,dmin,dmax) 
c          print *,' file7= ',file7
c          print *,' file8= ',file8
           print *,'read pastbias, iret=, nt= ',iret,nt,dmin,dmax
           if (iret.eq.0) then
            do ipt=1,ngrid
             pastbias(ipt) = f(ipt)
            enddo
           else
            print *,"Something wrong in reading pastday bias!"
           endif
          endif

         enddo  !npast

C to calculate and write out
         do i=1,25
          kpds(i)=kkpds(i,nv)
         enddo
c        kpds(2)=111  !different models, deal with this later
c        kpds(4)=0  !bitmap
         kpds(8)=yy1-(yy1/100)*100
         kpds(9)=mm1
         kpds(10)=dd1
         kpds(11)=hh1
         kpds(12)=0
         kpds(14)=(nt-1)*3+3
         kpds(15)=0
         kpds(21)=21
         kpds(22)=4  !precision
         kpds(5)=jpds5(nv)
         kpds(6)=jpds6(nv)
         kpds(7)=jpds7(nv)

c        kpds(23)=0  !not ensemble
         kpds(23)=2  !ensemble

c find ensemble mean and spread
         if(method.eq.3.or.method.eq.4) then
         do ipt=1,ngrid
          do nf=mem_sta(nens),mem_end(nens) 
           fcst(nf)=data(ipt,nf,nv)
          enddo
      call mean_spread(nfile,mem_sta(nens),mem_end(nens),fcst,ave,std)
c         fmean(ipt,nens,nv)=ave
          xmean(ipt)=ave
          spread(ipt)=std
         enddo
         endif

         if(method.eq.4) then
         avesprd=0.0
         do ipt=1,ngrid
          do nday=1,npastday
           avesprd(ipt)=avesprd(ipt)+pastsprd(ipt,nday)/float(npastday)
          enddo
          sprd_norm(ipt)=spread(ipt)/avesprd(ipt)   !????
         enddo
         endif

c corr coeff between today's fcst and past fcsts (wrt subgroup ens mean)
c (currently over the CONUS and plan to divide it into subregions in future)
         if (method.eq.3) then
         sum_r=0.0000000001
         do npast=1,npastday
          do ipt=1,ngrid
           x(ipt)=xmean(ipt)
           y(ipt)=pastfcst(ipt,npast)
          enddo
          call correlation(x,y,r,lat,lon)
          corr(npast)=r
c --------------------------
c         if(nv.eq.1) then    !use SLP pattern to identify flow regime
c or use SLP+500H+the variable itself three fields (?)
           if(corr(npast).gt.0.) weight(npast)=corr(npast)
           if(corr(npast).le.0.) weight(npast)=0.0
c         endif
          sum_r=sum_r+weight(npast)
c --------------------------
         enddo 
         endif

c assign weights
         do npast=1,npastday
          if(method.eq.1) w(npast)=1.0/npastday
          if(method.eq.2) w(npast)=0.05
c         if(method.eq.2) w(npast)=0.02
c         if(method.eq.2) w(npast)=0.1
          if(method.eq.3) w(npast)=weight(npast)/sum_r
         enddo
         if(method.eq.3) write(98,90) (corr(npast),npast=1,npastday)
         if(method.eq.3) write(98,91) (w(npast),npast=1,npastday)
90       format(1x,20f5.1)
91       format(1x,20f5.2)
c estimate bias with averaging methods 
         call bias_averaging(method,type,lon,lat,npastday,pasterro,w,
     &                       pastbias,bias)
c        write(97,*) nens
c        write(97,*) (bias(ipt),ipt=1,ngrid)
         call putgb(lugb_newbc,jf,kpds,kgds,lb,bias,iret)
c        call putgbe(lugb_newbc,jf,kpds,kgds,kens,lb,bias,iret)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C write out debiased new ensemble fcsts:
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c         if (nens.eq.1) kpds(2)=111  !nmmb
c         if (nens.eq.2) kpds(2)=112  !nmm
c         if (nens.eq.3) kpds(2)=116  !arw
          if (nens.eq.1) kpds(2)=111  !nmmb
          if (nens.eq.2) kpds(2)=116  !arw
          do nf=mem_sta(nens),mem_end(nens)
           do ipt=1,ngrid
            fcst_BC(ipt)=data(ipt,nf,nv)-bias(ipt)
            if(jpds5(nv).eq.51.or.jpds5(nv).eq.52.or.
     &jpds5(nv).eq.61) then
             if(fcst_BC(ipt).lt.0.) fcst_BC(ipt)=0.  
            endif
            if(jpds5(nv).eq.52) then
             if(fcst_BC(ipt).gt.100.) fcst_BC(ipt)=100.  
            endif
           enddo

            if(kpds(23).eq.0) then
           if(nf.eq.1) call putgb(10,jf,kpds,kgds,lb,fcst_BC,iret)
           if(nf.eq.2) call putgb(11,jf,kpds,kgds,lb,fcst_BC,iret)
           if(nf.eq.3) call putgb(12,jf,kpds,kgds,lb,fcst_BC,iret)
           if(nf.eq.4) call putgb(13,jf,kpds,kgds,lb,fcst_BC,iret)
           if(nf.eq.5) call putgb(14,jf,kpds,kgds,lb,fcst_BC,iret)
           if(nf.eq.6) call putgb(15,jf,kpds,kgds,lb,fcst_BC,iret)
           if(nf.eq.7) call putgb(16,jf,kpds,kgds,lb,fcst_BC,iret)
           if(nf.eq.8) call putgb(17,jf,kpds,kgds,lb,fcst_BC,iret)
           if(nf.eq.9) call putgb(18,jf,kpds,kgds,lb,fcst_BC,iret)
           if(nf.eq.10) call putgb(19,jf,kpds,kgds,lb,fcst_BC,iret)
           if(nf.eq.11) call putgb(20,jf,kpds,kgds,lb,fcst_BC,iret)
           if(nf.eq.12) call putgb(21,jf,kpds,kgds,lb,fcst_BC,iret)
           if(nf.eq.13) call putgb(22,jf,kpds,kgds,lb,fcst_BC,iret)
           if(nf.eq.14) call putgb(23,jf,kpds,kgds,lb,fcst_BC,iret)
           if(nf.eq.15) call putgb(24,jf,kpds,kgds,lb,fcst_BC,iret)
           if(nf.eq.16) call putgb(25,jf,kpds,kgds,lb,fcst_BC,iret)
           if(nf.eq.17) call putgb(26,jf,kpds,kgds,lb,fcst_BC,iret)
           if(nf.eq.18) call putgb(27,jf,kpds,kgds,lb,fcst_BC,iret)
           if(nf.eq.19) call putgb(28,jf,kpds,kgds,lb,fcst_BC,iret)
           if(nf.eq.20) call putgb(29,jf,kpds,kgds,lb,fcst_BC,iret)
           if(nf.eq.21) call putgb(30,jf,kpds,kgds,lb,fcst_BC,iret)
           if(nf.eq.22) call putgb(31,jf,kpds,kgds,lb,fcst_BC,iret)
           if(nf.eq.23) call putgb(32,jf,kpds,kgds,lb,fcst_BC,iret)
           if(nf.eq.24) call putgb(33,jf,kpds,kgds,lb,fcst_BC,iret)
           if(nf.eq.25) call putgb(34,jf,kpds,kgds,lb,fcst_BC,iret)
           if(nf.eq.26) call putgb(35,jf,kpds,kgds,lb,fcst_BC,iret)

            else

c make an ensemble section
           kens=-1
           kens(1)=1
           kens(4)=1
           kens(5)=255

c nmb
           if(nf.eq.1) then   
            kens(2)=1         !low res control
            kens(3)=2         !ctl
            call putgbe(10,jf,kpds,kgds,kens,lb,fcst_BC,iret)
           endif
           if(nf.eq.2) then   
            kens(2)=2         !negatively perturbed member
            kens(3)=1         !n1
            call putgbe(11,jf,kpds,kgds,kens,lb,fcst_BC,iret)
           endif
           if(nf.eq.3) then   
            kens(2)=3         !positively perturbed member
            kens(3)=1         !p1
            call putgbe(12,jf,kpds,kgds,kens,lb,fcst_BC,iret)
           endif
           if(nf.eq.4) then   
            kens(2)=2
            kens(3)=2         !n2
            call putgbe(13,jf,kpds,kgds,kens,lb,fcst_BC,iret)
           endif
           if(nf.eq.5) then   
            kens(2)=3
            kens(3)=2         !p2
            call putgbe(14,jf,kpds,kgds,kens,lb,fcst_BC,iret)
           endif
           if(nf.eq.6) then 
            kens(2)=2
            kens(3)=3         !n3
            call putgbe(15,jf,kpds,kgds,kens,lb,fcst_BC,iret)
           endif
           if(nf.eq.7) then 
            kens(2)=3
            kens(3)=3         !p3
            call putgbe(16,jf,kpds,kgds,kens,lb,fcst_BC,iret)
           endif
           if(nf.eq.8) then 
            kens(2)=2
            kens(3)=4         !n4
            call putgbe(17,jf,kpds,kgds,kens,lb,fcst_BC,iret)
           endif
           if(nf.eq.9) then 
            kens(2)=3
            kens(3)=4         !p4
            call putgbe(18,jf,kpds,kgds,kens,lb,fcst_BC,iret)
           endif
           if(nf.eq.10) then 
            kens(2)=2
            kens(3)=5         !n5
            call putgbe(19,jf,kpds,kgds,kens,lb,fcst_BC,iret)
           endif
           if(nf.eq.11) then 
            kens(2)=3
            kens(3)=5         !p5
            call putgbe(20,jf,kpds,kgds,kens,lb,fcst_BC,iret)
           endif
           if(nf.eq.12) then 
            kens(2)=2
            kens(3)=6         !n6
            call putgbe(21,jf,kpds,kgds,kens,lb,fcst_BC,iret)
           endif
           if(nf.eq.13) then 
            kens(2)=3
            kens(3)=6         !p6
            call putgbe(22,jf,kpds,kgds,kens,lb,fcst_BC,iret)
           endif

c arw
           if(nf.eq.14) then   
            kens(2)=1         !low res control
            kens(3)=2         !ctl
            call putgbe(23,jf,kpds,kgds,kens,lb,fcst_BC,iret)
           endif
           if(nf.eq.15) then   
            kens(2)=2         !negatively perturbed member
            kens(3)=1         !n1
            call putgbe(24,jf,kpds,kgds,kens,lb,fcst_BC,iret)
           endif
           if(nf.eq.16) then   
            kens(2)=3         !positively perturbed member
            kens(3)=1         !p1
            call putgbe(25,jf,kpds,kgds,kens,lb,fcst_BC,iret)
           endif
           if(nf.eq.17) then   
            kens(2)=2
            kens(3)=2         !n2
            call putgbe(26,jf,kpds,kgds,kens,lb,fcst_BC,iret)
           endif
           if(nf.eq.18) then   
            kens(2)=3
            kens(3)=2         !p2
            call putgbe(27,jf,kpds,kgds,kens,lb,fcst_BC,iret)
           endif
           if(nf.eq.19) then 
            kens(2)=2
            kens(3)=3         !n3
            call putgbe(28,jf,kpds,kgds,kens,lb,fcst_BC,iret)
           endif
           if(nf.eq.20) then 
            kens(2)=3
            kens(3)=3         !p3
            call putgbe(29,jf,kpds,kgds,kens,lb,fcst_BC,iret)
           endif
           if(nf.eq.21) then 
            kens(2)=2
            kens(3)=4         !n4
            call putgbe(30,jf,kpds,kgds,kens,lb,fcst_BC,iret)
           endif
           if(nf.eq.22) then 
            kens(2)=3
            kens(3)=4         !p4
            call putgbe(31,jf,kpds,kgds,kens,lb,fcst_BC,iret)
           endif
           if(nf.eq.23) then 
            kens(2)=2
            kens(3)=5         !n5
            call putgbe(32,jf,kpds,kgds,kens,lb,fcst_BC,iret)
           endif
           if(nf.eq.24) then 
            kens(2)=3
            kens(3)=5         !p5
            call putgbe(33,jf,kpds,kgds,kens,lb,fcst_BC,iret)
           endif
           if(nf.eq.25) then 
            kens(2)=2
            kens(3)=6         !n6
            call putgbe(34,jf,kpds,kgds,kens,lb,fcst_BC,iret)
           endif
           if(nf.eq.26) then 
            kens(2)=3
            kens(3)=6         !p6
            call putgbe(35,jf,kpds,kgds,kens,lb,fcst_BC,iret)
           endif
            endif   !kpds(23)
          enddo

        enddo    !variables
          call baclose(lugb_err,ierr)
          call baclose(lugi_err,ierr)
          call baclose(lugb_fcst,ierr)
          call baclose(lugi_fcst,ierr)
          call baclose(lugb_sprd,ierr)
          call baclose(lugi_sprd,ierr)
          call baclose(lugb_bias,ierr)
          call baclose(lugi_bias,ierr)
          call baclose(lugb_newbc,ierr)
       enddo  !nSREF
        call baclose(10,ierr)
        call baclose(11,ierr)
        call baclose(12,ierr)
        call baclose(13,ierr)
        call baclose(14,ierr)
        call baclose(15,ierr)
        call baclose(16,ierr)
        call baclose(17,ierr)
        call baclose(18,ierr)
        call baclose(19,ierr)
        call baclose(20,ierr)
        call baclose(21,ierr)
        call baclose(22,ierr)
        call baclose(23,ierr)
        call baclose(24,ierr)
        call baclose(25,ierr)
        call baclose(26,ierr)
        call baclose(27,ierr)
        call baclose(28,ierr)
        call baclose(29,ierr)
        call baclose(30,ierr)
        call baclose(31,ierr)
        call baclose(32,ierr)
        call baclose(33,ierr)
        call baclose(34,ierr)
        call baclose(35,ierr)

10    continue  !time

9999	END

C************************************************
      subroutine grange(n,ld,d,dmin,dmax)
      logical ld
      dimension ld(n),d(n)
      dmin=1.e38
      dmax=-1.e38
      do i=1,n
       if(ld(i)) then
        dmin=min(dmin,d(i))
        dmax=max(dmax,d(i))
       endif
      enddo
      return
      end
 
C***************************************************
      subroutine getname (fnum,fname1,fname2)
      character*13  fname1
      character*15  fname2
      integer    fnum
      if (fnum.eq.1)  fname1='r_gribawips01'
      if (fnum.eq.2)  fname1='r_gribawips02'
      if (fnum.eq.3)  fname1='r_gribawips03'
      if (fnum.eq.4)  fname1='r_gribawips04'
      if (fnum.eq.5)  fname1='r_gribawips05'
      if (fnum.eq.6)  fname1='r_gribawips06'
      if (fnum.eq.7)  fname1='r_gribawips07'
      if (fnum.eq.8)  fname1='r_gribawips08'
      if (fnum.eq.9)  fname1='r_gribawips09'
      if (fnum.eq.10) fname1='r_gribawips10'
      if (fnum.eq.11) fname1='r_gribawips11'
      if (fnum.eq.12) fname1='r_gribawips12'
      if (fnum.eq.13) fname1='r_gribawips13'
      if (fnum.eq.14) fname1='r_gribawips14'
      if (fnum.eq.15) fname1='r_gribawips15'
      if (fnum.eq.16) fname1='r_gribawips16'
      if (fnum.eq.17) fname1='r_gribawips17'
      if (fnum.eq.18) fname1='r_gribawips18'
      if (fnum.eq.19) fname1='r_gribawips19'
      if (fnum.eq.20) fname1='r_gribawips20'
      if (fnum.eq.21) fname1='r_gribawips21'
      if (fnum.eq.22) fname1='r_gribawips22'
      if (fnum.eq.23) fname1='r_gribawips23'
      if (fnum.eq.24) fname1='r_gribawips24'
      if (fnum.eq.25) fname1='r_gribawips25'
      if (fnum.eq.26) fname1='r_gribawips26'
      if (fnum.eq.27) fname1='r_gribawips27'
      if (fnum.eq.28) fname1='r_gribawips28'
      if (fnum.eq.29) fname1='r_gribawips29'
      if (fnum.eq.30) fname1='r_gribawips30'

      if (fnum.eq.1)  fname2='r_gribawips01.i'
      if (fnum.eq.2)  fname2='r_gribawips02.i'
      if (fnum.eq.3)  fname2='r_gribawips03.i'
      if (fnum.eq.4)  fname2='r_gribawips04.i'
      if (fnum.eq.5)  fname2='r_gribawips05.i'
      if (fnum.eq.6)  fname2='r_gribawips06.i'
      if (fnum.eq.7)  fname2='r_gribawips07.i'
      if (fnum.eq.8)  fname2='r_gribawips08.i'
      if (fnum.eq.9)  fname2='r_gribawips09.i'
      if (fnum.eq.10) fname2='r_gribawips10.i'
      if (fnum.eq.11) fname2='r_gribawips11.i'
      if (fnum.eq.12) fname2='r_gribawips12.i'
      if (fnum.eq.13) fname2='r_gribawips13.i'
      if (fnum.eq.14) fname2='r_gribawips14.i'
      if (fnum.eq.15) fname2='r_gribawips15.i'
      if (fnum.eq.16) fname2='r_gribawips16.i'
      if (fnum.eq.17) fname2='r_gribawips17.i'
      if (fnum.eq.18) fname2='r_gribawips18.i'
      if (fnum.eq.19) fname2='r_gribawips19.i'
      if (fnum.eq.20) fname2='r_gribawips20.i'
      if (fnum.eq.21) fname2='r_gribawips21.i'
      if (fnum.eq.22) fname2='r_gribawips22.i'
      if (fnum.eq.23) fname2='r_gribawips23.i'
      if (fnum.eq.24) fname2='r_gribawips24.i'
      if (fnum.eq.25) fname2='r_gribawips25.i'
      if (fnum.eq.26) fname2='r_gribawips26.i'
      if (fnum.eq.27) fname2='r_gribawips27.i'
      if (fnum.eq.28) fname2='r_gribawips28.i'
      if (fnum.eq.29) fname2='r_gribawips29.i'
      if (fnum.eq.30) fname2='r_gribawips30.i'
      return
      end

C******************************************************
        subroutine correlation(x,y,r,lat,lon)
        dimension x(lon*lat),y(lon*lat)
        r=0.0
        r1=0.0
        r2=0.0
        xm=0.0
        ym=0.0
C calculate mean values over the selected region
        do ipt=1,lat*lon
         xm=xm+x(ipt)/(lat*lon)
         ym=ym+y(ipt)/(lat*lon)
        enddo
c calculate anamoally correlation coefficient
         do ipt=1,lat*lon
        r1=r1+(x(ipt)-xm)**2
        r2=r2+(y(ipt)-ym)**2
        r=r+(x(ipt)-xm)*(y(ipt)-ym)
         enddo
        if(r1.eq.0.0.and.r2.eq.0.0) r=100.0
        if(r1.eq.0.0.and.r2.ne.0.0) r1=0.01
        if(r2.eq.0.0.and.r1.ne.0.0) r2=0.01
        if(r1.ne.0.0.and.r2.ne.0.0)
     &  r=r/sqrt(r1)/sqrt(r2)*100.
        return
        end

C******************************************************
      subroutine mean_spread(nmax,mem_sta,mem_end,f,ave,std)

c  Purpose: computing ensemble mean and spread among x number of members
c           at a given grid point (could be a full set or a subset)
c  PRGMMR: Jun Du          ORG: NP21       DATE: Apr. 05, 2004
c
c  Input:
c	 nmax=maximum number of ensemble members allowed
c	 mem_sta=the number of the first member in a sub-ensemble set
c	 mem_end=the number of the last member in a sub-ensemble set
c	 f=ensemble forecasts data
c  Output:
c	 ave=ensemble mean
c	 std=ensemble spread

      dimension f(nmax)
      ave=0.
      std=0.
      do i=mem_sta,mem_end
       ave=ave+f(i)/float(mem_end-mem_sta+1)
      enddo
      do i=mem_sta,mem_end
       std=std+(f(i)-ave)**2.0
      enddo
      std=sqrt(std/float(mem_end-mem_sta+1))
      return
      end

      subroutine bias_averaging(nopt,type,ix,jy,kd,err,w,b_old,b_new)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C This subroutine estimates model forecast mean Bias using popular
C Averaging Method based on past forecast error information in
C three different ways (option to chose from):
C  (1) simple average over a period of some past days (equally weighting)
C  (2) decaying average (more weight toward most recent forecast errors)
C  (3) regime-selective average (use corr coeff between today's forecast
C      and past forecasts as weights, Du 2004)
C
C Change log:
C  inital program: 2004-09-15, Jun Du
C
C INPUT:
C  nopt:  1 for simple, 2 decaying, 3 regime-selective averaging
C  ix:    x dimention
C  jy:    y dimension
C  kd:    number of past days used
C  err:   forecast error(s) of past days
C  w:     weights (w(1)+w(2)+w(3)+....=1.)
C  b_old: previous accumulated bias (for decaying average method only)
C
C OUTPUT:
C  b_new: new bias estimation
C
c SUBPROGRAMS CALLED:
C  none
C ATTRIBUTES:
c  language: Fortran 77
c  machine: IBM-SP
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
      dimension err(ix*jy,kd),w(kd),b_old(ix*jy),b_new(ix*jy)
      integer type
 
      print*,"inside subroutine Bias_averaging"
      print*,"method = ",nopt
      print*,"type = ",type
      print*,"npastday = ",kd
      print*,"weights = ", (w(i),i=1,kd)

c initialization
      b_new=0.0
 
c choose a method
      if (nopt.eq.1.or.nopt.eq.3) go to 111
      if (nopt.eq.2) go to 222
 
111   continue
C simple average or regime-weighted average method (kd>=1)
      if(nopt.eq.1) print*,'simple-average method applied'
      if(nopt.eq.3) print*,'regime-selective approach applied'
      do k=1,kd
       do i=1,ix*jy
        b_new(i)=b_new(i)+w(k)*err(i,k)
       enddo
      enddo
      goto 9999
 
222   continue
C decaying average method (kd=1)
      print*,'decaying-average method applied'
      do k=1,1
       do i=1,ix*jy
c       if(b_old(i).lt.-9999.0) then
        if(type.eq.2) then
         if(i.eq.1) print*,"no revious bias existed"
         b_new(i)=0.25*err(i,k)        !if previous accumulated bias not available
        else
         if(i.eq.1) print*,"previous bias found"
         b_new(i)=(1-w(k))*b_old(i)+w(k)*err(i,k)
        endif
       enddo
      enddo
 
9999  return
      end

